#include "spacewind/covariant_em.h"
#include "spacewind/em_pic3d.h"
#include "spacewind/four_momentum_delta.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);
static double worst_mapping_relative = 0.0;
static double worst_boost_covariance = 0.0;
static double worst_dynamic_residual = 0.0;
static double worst_stress_momentum_relative = 0.0;

static double max2(double a, double b) { return a > b ? a : b; }
static double max3(double a, double b, double c) {
    return max2(max2(a, b), c);
}

static void check_true(int condition, const char *name) {
    ++checks;
    transcript_hash ^= swa_fnv1a64(name, strlen(name));
    transcript_hash *= UINT64_C(1099511628211);
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}

static size_t index3(
    const swa_em_pic3d *state,
    size_t i,
    size_t j,
    size_t k
) {
    return (k * state->ny + j) * state->nx + i;
}

static swa_em_pic3d_limits limits(void) {
    swa_em_pic3d_limits gate;
    memset(&gate, 0, sizeof(gate));
    gate.continuity_relative_tolerance = 2.0e-10;
    gate.gauss_relative_tolerance = 2.0e-11;
    gate.magnetic_divergence_relative_tolerance = 2.0e-11;
    gate.local_current_difference_relative_tolerance = 0.85;
    gate.spectral_oracle_relative_tolerance = 2.0e-10;
    gate.spectral_oracle_curl_relative_tolerance = 2.0e-10;
    gate.energy_relative_tolerance = 2.0e-10;
    gate.particle_work_relative_tolerance = 5.0e-6;
    gate.field_work_relative_tolerance = 2.0e-10;
    gate.momentum_relative_tolerance = 2.0e-4;
    gate.nonlinear_relative_tolerance = 2.0e-10;
    gate.maximum_courant = 0.8;
    gate.maximum_particle_cells_per_step = 0.8;
    gate.maximum_nonlinear_iterations = 16U;
    return gate;
}

static double rest_energy_J(const swa_em_pic3d *state) {
    swa_kahan energy = {0.0, 0.0};
    size_t p;
    for (p = 0U; p < state->particle_count; ++p) {
        swa_kahan_add(
            &energy,
            state->particles[p].macro_weight *
            state->particles[p].mass_kg * SWA_C * SWA_C
        );
    }
    return energy.sum;
}

static swa_four zero_four(void) {
    return swa_four_make(0.0, 0.0, 0.0, 0.0);
}

static swa_four_momentum_ledger ledger_from_result(
    const swa_em_pic3d_result *result,
    double particle_rest_energy_J
) {
    swa_four_momentum_ledger ledger;
    ledger.field_initial_Ns = swa_four_make(
        result->initial_field_energy_J / SWA_C,
        result->initial_field_momentum_Ns.x,
        result->initial_field_momentum_Ns.y,
        result->initial_field_momentum_Ns.z
    );
    ledger.field_final_Ns = swa_four_make(
        result->final_field_energy_J / SWA_C,
        result->final_field_momentum_Ns.x,
        result->final_field_momentum_Ns.y,
        result->final_field_momentum_Ns.z
    );
    ledger.matter_initial_Ns = swa_four_make(
        (particle_rest_energy_J + result->initial_particle_energy_J) /
            SWA_C,
        result->initial_particle_momentum_Ns.x,
        result->initial_particle_momentum_Ns.y,
        result->initial_particle_momentum_Ns.z
    );
    ledger.matter_final_Ns = swa_four_make(
        (particle_rest_energy_J + result->final_particle_energy_J) /
            SWA_C,
        result->final_particle_momentum_Ns.x,
        result->final_particle_momentum_Ns.y,
        result->final_particle_momentum_Ns.z
    );
    ledger.external_impulse_Ns = zero_four();
    ledger.momentum_in_Ns = zero_four();
    ledger.momentum_out_Ns = zero_four();
    return ledger;
}

static swa_four_momentum_delta_ledger delta_ledger_from_result(
    const swa_em_pic3d_result *result
) {
    swa_four_momentum_delta_ledger ledger;
    ledger.field_change_Ns = swa_four_make(
        result->field_energy_change_J / SWA_C,
        result->final_field_momentum_Ns.x -
            result->initial_field_momentum_Ns.x,
        result->final_field_momentum_Ns.y -
            result->initial_field_momentum_Ns.y,
        result->final_field_momentum_Ns.z -
            result->initial_field_momentum_Ns.z
    );
    ledger.matter_change_Ns = swa_four_make(
        result->particle_energy_change_J / SWA_C,
        result->final_particle_momentum_Ns.x -
            result->initial_particle_momentum_Ns.x,
        result->final_particle_momentum_Ns.y -
            result->initial_particle_momentum_Ns.y,
        result->final_particle_momentum_Ns.z -
            result->initial_particle_momentum_Ns.z
    );
    ledger.external_impulse_Ns = zero_four();
    ledger.momentum_in_Ns = zero_four();
    ledger.momentum_out_Ns = zero_four();
    return ledger;
}

static double dynamic_scale_Ns(const swa_em_pic3d_result *result) {
    const double initial_energy = fabs(
        result->initial_particle_energy_J +
        result->initial_field_energy_J
    ) / SWA_C;
    const double final_energy = fabs(
        result->final_particle_energy_J +
        result->final_field_energy_J
    ) / SWA_C;
    const double initial_momentum = swa_vnorm(swa_vadd(
        result->initial_particle_momentum_Ns,
        result->initial_field_momentum_Ns
    ));
    const double final_momentum = swa_vnorm(swa_vadd(
        result->final_particle_momentum_Ns,
        result->final_field_momentum_Ns
    ));
    return max3(
        max2(initial_energy, final_energy),
        max2(initial_momentum, final_momentum),
        DBL_MIN
    );
}

static double four_relative(
    swa_four observed,
    swa_four expected,
    double physical_scale
) {
    return swa_four_euclidean_norm(swa_four_sub(observed, expected)) /
        max3(
            swa_four_euclidean_norm(observed),
            swa_four_euclidean_norm(expected),
            physical_scale
        );
}

static int transform_ledger(
    swa_vec3 frame_velocity,
    const swa_four_momentum_ledger *input,
    swa_four_momentum_ledger *output
) {
    return swa_lorentz_transform_four(
               frame_velocity, input->field_initial_Ns,
               &output->field_initial_Ns) &&
           swa_lorentz_transform_four(
               frame_velocity, input->field_final_Ns,
               &output->field_final_Ns) &&
           swa_lorentz_transform_four(
               frame_velocity, input->matter_initial_Ns,
               &output->matter_initial_Ns) &&
           swa_lorentz_transform_four(
               frame_velocity, input->matter_final_Ns,
               &output->matter_final_Ns) &&
           swa_lorentz_transform_four(
               frame_velocity, input->external_impulse_Ns,
               &output->external_impulse_Ns) &&
           swa_lorentz_transform_four(
               frame_velocity, input->momentum_in_Ns,
               &output->momentum_in_Ns) &&
           swa_lorentz_transform_four(
               frame_velocity, input->momentum_out_Ns,
               &output->momentum_out_Ns);
}

static swa_four integrate_cell_stress_energy(
    const swa_em_pic3d *state,
    double *maximum_dominant_ratio,
    double *trace_integral_J
) {
    const double dx = state->length_x_m / (double)state->nx;
    const double dy = state->length_y_m / (double)state->ny;
    const double dz = state->length_z_m / (double)state->nz;
    const double volume = dx * dy * dz;
    swa_kahan p0 = {0.0, 0.0};
    swa_kahan px = {0.0, 0.0};
    swa_kahan py = {0.0, 0.0};
    swa_kahan pz = {0.0, 0.0};
    swa_kahan trace = {0.0, 0.0};
    size_t i;
    size_t j;
    size_t k;
    *maximum_dominant_ratio = 0.0;
    for (k = 0U; k < state->nz; ++k) {
        const size_t kp = (k + 1U) % state->nz;
        for (j = 0U; j < state->ny; ++j) {
            const size_t jp = (j + 1U) % state->ny;
            for (i = 0U; i < state->nx; ++i) {
                const size_t ip = (i + 1U) % state->nx;
                const size_t center = index3(state, i, j, k);
                swa_em_fields fields;
                swa_tensor4 tensor;
                swa_covariant_em_diagnostics diagnostics;
                fields.electric_field_V_m = swa_v3(
                    0.5 * (
                        state->ex_V_m[center] +
                        state->ex_V_m[index3(state, ip, j, k)]
                    ),
                    0.5 * (
                        state->ey_V_m[center] +
                        state->ey_V_m[index3(state, i, jp, k)]
                    ),
                    0.5 * (
                        state->ez_V_m[center] +
                        state->ez_V_m[index3(state, i, j, kp)]
                    )
                );
                fields.magnetic_field_T = swa_v3(
                    0.25 * (
                        state->bx_T[center] +
                        state->bx_T[index3(state, i, jp, k)] +
                        state->bx_T[index3(state, i, j, kp)] +
                        state->bx_T[index3(state, i, jp, kp)]
                    ),
                    0.25 * (
                        state->by_T[center] +
                        state->by_T[index3(state, ip, j, k)] +
                        state->by_T[index3(state, i, j, kp)] +
                        state->by_T[index3(state, ip, j, kp)]
                    ),
                    0.25 * (
                        state->bz_T[center] +
                        state->bz_T[index3(state, ip, j, k)] +
                        state->bz_T[index3(state, i, jp, k)] +
                        state->bz_T[index3(state, ip, jp, k)]
                    )
                );
                check_true(swa_em_stress_energy_tensor(
                               &fields, &tensor, &diagnostics),
                           "3D3V cell stress-energy tensor executes");
                check_true(diagnostics.dominant_energy_condition,
                           "3D3V cell stress-energy obeys dominant-energy bound");
                swa_kahan_add(&p0,
                    tensor.component[0][0] * volume / SWA_C);
                swa_kahan_add(&px,
                    tensor.component[0][1] * volume / SWA_C);
                swa_kahan_add(&py,
                    tensor.component[0][2] * volume / SWA_C);
                swa_kahan_add(&pz,
                    tensor.component[0][3] * volume / SWA_C);
                swa_kahan_add(&trace,
                    diagnostics.tensor_trace_J_m3 * volume);
                *maximum_dominant_ratio = max2(
                    *maximum_dominant_ratio,
                    diagnostics.dominant_energy_ratio
                );
            }
        }
    }
    *trace_integral_J = trace.sum;
    return swa_four_make(p0.sum, px.sum, py.sum, pz.sum);
}

static void test_covariant_3d3v_step(void) {
    swa_em_pic3d state;
    swa_em_pic3d_result result;
    swa_em_pic3d_limits gate = limits();
    const double dx = 1.0 / 4.0;
    const double dy = 1.2 / 5.0;
    const double dz = 1.4 / 6.0;
    const double dt = 0.08 / (
        SWA_C * sqrt(
            1.0 / (dx * dx) +
            1.0 / (dy * dy) +
            1.0 / (dz * dz)
        )
    );
    const swa_vec3 boost = swa_v3(
        0.07 * SWA_C,
        -0.04 * SWA_C,
        0.03 * SWA_C
    );
    double rest_energy;
    size_t p;
    size_t index;
    size_t step;
    FILE *fp;
    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, 4U, 5U, 6U, 8U,
                   1.0, 1.2, 1.4, dt),
               "covariant 3D3V grid initializes");
    for (p = 0U; p < state.particle_count; ++p) {
        const double d = (double)p;
        const double sign = p % 2U == 0U ? 1.0 : -1.0;
        check_true(swa_em_pic3d_set_particle(
                       &state, p,
                       swa_v3(
                           0.10 + 0.11 * d,
                           0.08 + 0.13 * d,
                           0.06 + 0.17 * d
                       ),
                       swa_v3(
                           sign * 0.005 * SWA_C,
                           (0.001 + 0.0001 * d) * SWA_C,
                           -sign * 0.0007 * SWA_C
                       ),
                       sign * 1.0e-16,
                       2.0e-27,
                       1.0),
                   "covariant 3D3V particle configures");
    }
    check_true(swa_em_pic3d_initialize_gauss_field(
                   &state, swa_v3(2.0e5, -1.0e5, 0.5e5),
                   1.0e-25),
               "covariant 3D3V Gauss field initializes");
    for (index = 0U; index < state.cell_count; ++index) {
        state.bx_T[index] = 2.0e-3;
        state.by_T[index] = -1.0e-3;
        state.bz_T[index] = 0.5e-3;
    }
    rest_energy = rest_energy_J(&state);
    check_true(isfinite(rest_energy) && rest_energy > 0.0,
               "covariant 3D3V rest energy finite");

    for (step = 0U; step < 16U; ++step) {
        swa_four_momentum_ledger ledger;
        swa_four_momentum_ledger boosted_ledger;
        swa_four_momentum_ledger_result audit;
        swa_four_momentum_ledger_result boosted_audit;
        swa_four_momentum_delta_ledger delta_ledger;
        swa_four_momentum_delta_ledger boosted_delta_ledger;
        swa_four_momentum_delta_result delta_audit;
        swa_four_momentum_delta_result boosted_delta_audit;
        swa_four expected_residual;
        swa_four expected_boosted_residual;
        swa_four stress_four;
        const double dynamic_floor = 1.0e-30;
        double dynamic_scale;
        double mapping_error;
        double boost_error;
        double dominant_ratio;
        double trace_integral;
        double stress_momentum_error;

        check_true(swa_em_pic3d_step(&state, &gate, &result),
                   "covariant 3D3V step executes");
        check_true(result.passes,
                   "covariant 3D3V solver gates pass before promotion");
        check_true(result.local_current_passes &&
                   result.spectral_oracle_passes,
                   "covariant 3D3V step uses charge-conserving current independently accepted by spectral oracle");
        ledger = ledger_from_result(&result, rest_energy);
        delta_ledger = delta_ledger_from_result(&result);
        dynamic_scale = max2(dynamic_scale_Ns(&result), dynamic_floor);
        check_true(swa_audit_four_momentum_ledger(
                       &ledger, 0.0, 3.0e-4 * dynamic_scale, &audit),
                   "same-step total 3D3V four-momentum audit executes");
        check_true(audit.finite && audit.passes,
                   "same-step total 3D3V four-momentum closes on dynamic scale");
        check_true(swa_audit_four_momentum_delta_ledger(
                       &delta_ledger, dynamic_scale,
                       3.0e-4, 0.0, &delta_audit),
                   "incremental 3D3V four-momentum audit executes");
        check_true(delta_audit.finite && delta_audit.passes,
                   "incremental 3D3V four-momentum closes without total cancellation");
        expected_residual = swa_four_make(
            result.total_energy_residual_J / SWA_C,
            result.total_momentum_residual_Ns.x,
            result.total_momentum_residual_Ns.y,
            result.total_momentum_residual_Ns.z
        );
        mapping_error = four_relative(
            delta_audit.residual_Ns, expected_residual, dynamic_scale
        );
        check_true(mapping_error < 8.0e-15,
                   "incremental four-ledger exactly reproduces solver ledgers");
        check_true(transform_ledger(boost, &ledger, &boosted_ledger),
                   "boost complete total 3D3V four-momentum ledger");
        check_true(swa_audit_four_momentum_ledger(
                       &boosted_ledger, 1.0, DBL_MIN, &boosted_audit),
                   "boosted total 3D3V four-momentum audit executes");
        check_true(swa_lorentz_transform_four_momentum_delta_ledger(
                       boost, &delta_ledger, &boosted_delta_ledger),
                   "boost incremental 3D3V four-momentum ledger");
        check_true(swa_audit_four_momentum_delta_ledger(
                       &boosted_delta_ledger, dynamic_scale,
                       1.0, DBL_MIN, &boosted_delta_audit),
                   "boosted incremental 3D3V four-momentum audit executes");
        check_true(swa_lorentz_transform_four(
                       boost, delta_audit.residual_Ns,
                       &expected_boosted_residual),
                   "boost incremental 3D3V numerical four-residual");
        boost_error = four_relative(
            boosted_delta_audit.residual_Ns,
            expected_boosted_residual,
            max2(dynamic_scale,
                 swa_four_euclidean_norm(expected_boosted_residual))
        );
        check_true(boost_error < 3.0e-14,
                   "incremental 3D3V numerical residual transforms covariantly");

        stress_four = integrate_cell_stress_energy(
            &state, &dominant_ratio, &trace_integral
        );
        stress_momentum_error = swa_vnorm(swa_vsub(
            swa_v3(
                stress_four.component[1],
                stress_four.component[2],
                stress_four.component[3]
            ),
            result.final_field_momentum_Ns
        )) / max3(
            swa_vnorm(result.final_field_momentum_Ns),
            swa_vnorm(swa_v3(
                stress_four.component[1],
                stress_four.component[2],
                stress_four.component[3]
            )),
            result.final_field_energy_J / SWA_C
        );
        check_true(stress_momentum_error < 2.0e-12,
                   "cell stress-energy momentum matches 3D Yee field momentum");
        check_true(dominant_ratio <= 1.0 + 1.0e-13,
                   "integrated 3D field obeys dominant-energy bound");
        check_true(fabs(trace_integral) <
                   2.0e-12 * max2(result.final_field_energy_J, 1.0),
                   "integrated 3D electromagnetic stress-energy is traceless");

        worst_mapping_relative = max2(
            worst_mapping_relative, mapping_error
        );
        worst_boost_covariance = max2(
            worst_boost_covariance, boost_error
        );
        worst_dynamic_residual = max2(
            worst_dynamic_residual,
            delta_audit.residual_norm_Ns / dynamic_scale
        );
        worst_stress_momentum_relative = max2(
            worst_stress_momentum_relative,
            stress_momentum_error
        );
    }

    {
        swa_four_momentum_delta_ledger corrupted =
            delta_ledger_from_result(&result);
        swa_four_momentum_delta_result corrupted_audit;
        const double dynamic_scale = max2(
            dynamic_scale_Ns(&result), 1.0e-30
        );
        corrupted.field_change_Ns.component[2] +=
            0.01 * dynamic_scale;
        check_true(swa_audit_four_momentum_delta_ledger(
                       &corrupted, dynamic_scale,
                       3.0e-4, 0.0, &corrupted_audit),
                   "corrupted incremental 3D3V four-ledger audit executes");
        check_true(!corrupted_audit.passes,
                   "untracked incremental 3D3V four-impulse is rejected");
    }

    fp = fopen("output/em_pic3d_covariant_receipt.json", "w");
    check_true(fp != NULL, "open 3D3V covariant receipt");
    if (fp != NULL) {
        check_true(fprintf(
            fp,
            "{\n"
            "  \"schema\": \"spacewind.em-pic3d-covariant/v2\",\n"
            "  \"steps\": 16,\n"
            "  \"worst_residual_mapping_relative\": %.17g,\n"
            "  \"worst_boost_covariance_relative\": %.17g,\n"
            "  \"worst_dynamic_four_residual\": %.17g,\n"
            "  \"worst_stress_momentum_relative\": %.17g,\n"
            "  \"dynamic_scale_excludes_constant_rest_energy\": true,\n"
            "  \"passes\": %s,\n"
            "  \"nonclaim\": \"this bounded test verifies same-step 3D3V incremental four-momentum mapping, stress-energy consistency, Lorentz covariance, and the accepted coordinate-split charge-conserving current; it does not establish equivalence to a final production Esirkepov deposition, open-boundary plasma-wing force, chamber validity, or propulsion\"\n"
            "}\n",
            worst_mapping_relative,
            worst_boost_covariance,
            worst_dynamic_residual,
            worst_stress_momentum_relative,
            failures == 0U ? "true" : "false"
        ) > 0, "write 3D3V covariant receipt");
        (void)fclose(fp);
    }
    swa_em_pic3d_destroy(&state);
}

int main(void) {
    test_covariant_3d3v_step();
    printf(
        "spacewind 3D3V/covariant integration: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
