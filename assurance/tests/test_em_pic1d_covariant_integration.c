#include "spacewind/covariant_em.h"
#include "spacewind/em_pic1d.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);

static double max2(double a, double b) { return a > b ? a : b; }
static double max3(double a, double b, double c) { return max2(max2(a, b), c); }

static void check_true(int condition, const char *name) {
    ++checks;
    transcript_hash ^= swa_fnv1a64(name, strlen(name));
    transcript_hash *= UINT64_C(1099511628211);
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}

static swa_em_pic1d_limits limits(void) {
    swa_em_pic1d_limits result;
    result.continuity_relative_tolerance = 2.0e-11;
    result.gauss_relative_tolerance = 2.0e-11;
    result.magnetic_divergence_relative_tolerance = 1.0e-15;
    result.energy_relative_tolerance = 6.0e-2;
    result.momentum_relative_tolerance = 2.5e-1;
    result.maximum_courant = 0.8;
    result.maximum_particle_cells_per_step = 0.8;
    return result;
}

static void initialize_right_wave(
    swa_em_pic1d *state,
    unsigned mode,
    double amplitude_V_m
) {
    const double dx = state->domain_length_m / (double)state->cells;
    const double wave_number = 2.0 * SWA_PI * (double)mode /
                               state->domain_length_m;
    size_t i;
    for (i = 0U; i < state->cells; ++i) {
        const double x_face = (double)i * dx;
        const double x_cell = ((double)i + 0.5) * dx;
        state->ey_face_V_m[i] += amplitude_V_m * sin(wave_number * x_face);
        state->bz_cell_T[i] += amplitude_V_m / SWA_C *
                               sin(wave_number * x_cell);
    }
}

static double weighted_rest_energy_J(const swa_em_pic1d *state) {
    swa_kahan sum = {0.0, 0.0};
    size_t i;
    for (i = 0U; i < state->particle_count; ++i) {
        swa_kahan_add(
            &sum,
            state->particles[i].macro_weight *
            state->particles[i].mass_kg * SWA_C * SWA_C
        );
    }
    return sum.sum;
}

static swa_four field_initial_four(const swa_em_pic1d_result *r) {
    return swa_four_make(
        r->initial_field_energy_J / SWA_C,
        r->initial_field_momentum_Ns.x,
        r->initial_field_momentum_Ns.y,
        r->initial_field_momentum_Ns.z
    );
}

static swa_four field_final_four(const swa_em_pic1d_result *r) {
    return swa_four_make(
        r->final_field_energy_J / SWA_C,
        r->final_field_momentum_Ns.x,
        r->final_field_momentum_Ns.y,
        r->final_field_momentum_Ns.z
    );
}

static swa_four matter_initial_four(
    const swa_em_pic1d_result *r,
    double rest_energy_J
) {
    return swa_four_make(
        (rest_energy_J + r->initial_particle_energy_J) / SWA_C,
        r->initial_particle_momentum_Ns.x,
        r->initial_particle_momentum_Ns.y,
        r->initial_particle_momentum_Ns.z
    );
}

static swa_four matter_final_four(
    const swa_em_pic1d_result *r,
    double rest_energy_J
) {
    return swa_four_make(
        (rest_energy_J + r->final_particle_energy_J) / SWA_C,
        r->final_particle_momentum_Ns.x,
        r->final_particle_momentum_Ns.y,
        r->final_particle_momentum_Ns.z
    );
}

static swa_four zero_four(void) {
    return swa_four_make(0.0, 0.0, 0.0, 0.0);
}

static swa_four_momentum_ledger ledger_from_step(
    const swa_em_pic1d_result *r,
    double rest_energy_J
) {
    swa_four_momentum_ledger ledger;
    ledger.field_initial_Ns = field_initial_four(r);
    ledger.field_final_Ns = field_final_four(r);
    ledger.matter_initial_Ns = matter_initial_four(r, rest_energy_J);
    ledger.matter_final_Ns = matter_final_four(r, rest_energy_J);
    ledger.external_impulse_Ns = zero_four();
    ledger.momentum_in_Ns = zero_four();
    ledger.momentum_out_Ns = zero_four();
    return ledger;
}

static double dynamic_four_scale_Ns(const swa_em_pic1d_result *r) {
    const double initial_dynamic_energy =
        fabs(r->initial_particle_energy_J + r->initial_field_energy_J) / SWA_C;
    const double final_dynamic_energy =
        fabs(r->final_particle_energy_J + r->final_field_energy_J) / SWA_C;
    const double initial_momentum = swa_vnorm(swa_vadd(
        r->initial_particle_momentum_Ns,
        r->initial_field_momentum_Ns
    ));
    const double final_momentum = swa_vnorm(swa_vadd(
        r->final_particle_momentum_Ns,
        r->final_field_momentum_Ns
    ));
    return max3(
        max2(initial_dynamic_energy, final_dynamic_energy),
        max2(initial_momentum, final_momentum),
        DBL_MIN
    );
}

static double four_difference_relative(
    swa_four observed,
    swa_four expected,
    double physical_scale
) {
    return swa_four_euclidean_norm(swa_four_sub(observed, expected)) /
           max3(swa_four_euclidean_norm(observed),
                swa_four_euclidean_norm(expected),
                physical_scale);
}

static int transform_ledger(
    swa_vec3 frame_velocity_mps,
    const swa_four_momentum_ledger *input,
    swa_four_momentum_ledger *output
) {
    return swa_lorentz_transform_four(frame_velocity_mps,
                                      input->field_initial_Ns,
                                      &output->field_initial_Ns) &&
           swa_lorentz_transform_four(frame_velocity_mps,
                                      input->field_final_Ns,
                                      &output->field_final_Ns) &&
           swa_lorentz_transform_four(frame_velocity_mps,
                                      input->matter_initial_Ns,
                                      &output->matter_initial_Ns) &&
           swa_lorentz_transform_four(frame_velocity_mps,
                                      input->matter_final_Ns,
                                      &output->matter_final_Ns) &&
           swa_lorentz_transform_four(frame_velocity_mps,
                                      input->external_impulse_Ns,
                                      &output->external_impulse_Ns) &&
           swa_lorentz_transform_four(frame_velocity_mps,
                                      input->momentum_in_Ns,
                                      &output->momentum_in_Ns) &&
           swa_lorentz_transform_four(frame_velocity_mps,
                                      input->momentum_out_Ns,
                                      &output->momentum_out_Ns);
}

static void test_same_step_four_momentum(void) {
    enum { PARTICLES = 8, STEPS = 64 };
    swa_em_pic1d state;
    swa_em_pic1d_result step_result;
    swa_em_pic1d_limits gate = limits();
    const size_t cells = 64U;
    const double length = 2.0;
    const double area = 1.5;
    const double dx = length / (double)cells;
    const double dt = 0.10 * dx / SWA_C;
    const double charge = 2.0e-18;
    const double mass = 5.0e-22;
    const swa_vec3 boost_velocity = swa_v3(
        0.08 * SWA_C,
        -0.03 * SWA_C,
        0.02 * SWA_C
    );
    double rest_energy_J;
    double worst_mapping_relative = 0.0;
    double worst_boost_residual_covariance = 0.0;
    double worst_dynamic_four_residual = 0.0;
    size_t i;
    size_t step;
    FILE *fp;

    check_true(
        swa_em_pic1d_init(&state, cells, PARTICLES,
                          length, area, dt),
        "initialize covariant 1D3V PIC integration state"
    );
    state.bx_uniform_T = 2.0e-4;
    for (i = 0U; i < PARTICLES; ++i) {
        const double sign = (i % 2U) == 0U ? 1.0 : -1.0;
        const double x = length * ((double)i + 0.37) /
                         (double)PARTICLES;
        const swa_vec3 velocity = swa_v3(
            sign * (0.025 + 0.001 * (double)i) * SWA_C,
            (0.010 - 0.0005 * (double)i) * SWA_C,
            sign * 0.006 * SWA_C
        );
        check_true(
            swa_em_pic1d_set_particle(
                &state, i, x, velocity,
                sign * charge, mass, 1.0
            ),
            "configure covariant integration particle"
        );
    }
    check_true(
        swa_em_pic1d_initialize_gauss_field(&state, 0.0, 1.0e-28),
        "initialize covariant integration Gauss field"
    );
    initialize_right_wave(&state, 2U, 40.0);
    rest_energy_J = weighted_rest_energy_J(&state);
    check_true(isfinite(rest_energy_J) && rest_energy_J > 0.0,
               "weighted rest energy finite and positive");

    for (step = 0U; step < STEPS; ++step) {
        swa_four_momentum_ledger ledger;
        swa_four_momentum_ledger boosted_ledger;
        swa_four_momentum_ledger_result ledger_result;
        swa_four_momentum_ledger_result boosted_result;
        swa_four expected_residual;
        swa_four boosted_expected_residual;
        const double dynamic_scale_floor = 1.0e-30;
        double dynamic_scale;
        double dynamic_relative;
        double mapping_relative;
        double boost_covariance;
        double absolute_gate;

        check_true(
            swa_em_pic1d_step(&state, &gate, &step_result),
            "covariant integration EM PIC step executes"
        );
        check_true(step_result.finite,
                   "covariant integration EM PIC result finite");
        check_true(step_result.continuity_passes &&
                   step_result.initial_gauss_passes &&
                   step_result.final_gauss_passes,
                   "covariant integration charge and Gauss gates pass");
        check_true(step_result.subluminal_passes,
                   "covariant integration particles remain subluminal");
        check_true(step_result.energy_passes && step_result.momentum_passes,
                   "EM PIC energy and momentum gates pass before covariant promotion");

        ledger = ledger_from_step(&step_result, rest_energy_J);
        dynamic_scale = max2(dynamic_four_scale_Ns(&step_result),
                             dynamic_scale_floor);
        absolute_gate = 0.30 * dynamic_scale;
        check_true(
            swa_audit_four_momentum_ledger(
                &ledger, 0.0, absolute_gate, &ledger_result
            ),
            "same-step four-momentum ledger executes"
        );
        dynamic_relative = ledger_result.residual_norm_Ns / dynamic_scale;
        check_true(ledger_result.finite &&
                   dynamic_relative <= 0.30,
                   "same-step four-momentum closes on dynamic rather than rest-energy scale");

        expected_residual = swa_four_make(
            step_result.total_energy_residual_J / SWA_C,
            step_result.total_momentum_residual_Ns.x,
            step_result.total_momentum_residual_Ns.y,
            step_result.total_momentum_residual_Ns.z
        );
        mapping_relative = four_difference_relative(
            ledger_result.residual_Ns,
            expected_residual,
            dynamic_scale
        );
        check_true(mapping_relative < 5.0e-13,
                   "covariant ledger residual matches EM PIC energy and momentum residuals");

        check_true(transform_ledger(boost_velocity, &ledger,
                                    &boosted_ledger),
                   "boost complete same-step four-momentum ledger");
        check_true(
            swa_audit_four_momentum_ledger(
                &boosted_ledger, 1.0, DBL_MIN, &boosted_result
            ),
            "boosted same-step four-momentum ledger executes"
        );
        check_true(
            swa_lorentz_transform_four(
                boost_velocity,
                ledger_result.residual_Ns,
                &boosted_expected_residual
            ),
            "boost laboratory numerical residual"
        );
        boost_covariance = four_difference_relative(
            boosted_result.residual_Ns,
            boosted_expected_residual,
            max2(dynamic_scale, swa_four_euclidean_norm(
                boosted_expected_residual))
        );
        check_true(boost_covariance < 2.0e-12,
                   "four-momentum numerical residual transforms covariantly");

        worst_mapping_relative = max2(worst_mapping_relative,
                                      mapping_relative);
        worst_boost_residual_covariance = max2(
            worst_boost_residual_covariance,
            boost_covariance
        );
        worst_dynamic_four_residual = max2(
            worst_dynamic_four_residual,
            dynamic_relative
        );
    }

    {
        swa_four_momentum_ledger corrupted =
            ledger_from_step(&step_result, rest_energy_J);
        swa_four_momentum_ledger_result corrupted_result;
        const double dynamic_scale = max2(
            dynamic_four_scale_Ns(&step_result), 1.0e-30
        );
        corrupted.field_final_Ns.component[1] += 0.5 * dynamic_scale;
        check_true(
            swa_audit_four_momentum_ledger(
                &corrupted, 0.0, 0.30 * dynamic_scale,
                &corrupted_result
            ),
            "corrupted same-step four-momentum ledger executes"
        );
        check_true(!corrupted_result.passes,
                   "untracked spatial four-impulse is rejected");
    }

    fp = fopen("output/em_pic1d_covariant_receipt.json", "w");
    check_true(fp != NULL, "open EM PIC covariant integration receipt");
    if (fp != NULL) {
        check_true(
            fprintf(
                fp,
                "{\n"
                "  \"schema\": \"spacewind.em-pic1d-covariant/v1\",\n"
                "  \"steps\": %u,\n"
                "  \"worst_residual_mapping_relative\": %.17g,\n"
                "  \"worst_boost_residual_covariance_relative\": %.17g,\n"
                "  \"worst_dynamic_four_momentum_relative\": %.17g,\n"
                "  \"dynamic_scale_excludes_constant_rest_energy_from_tolerance\": true,\n"
                "  \"passes\": %s,\n"
                "  \"nonclaim\": \"this verifies same-step four-momentum bookkeeping and covariance for the bounded periodic 1D3V reference; it does not establish exact momentum conservation, multidimensional plasma-wing force, chamber similarity, or propulsion\"\n"
                "}\n",
                (unsigned)STEPS,
                worst_mapping_relative,
                worst_boost_residual_covariance,
                worst_dynamic_four_residual,
                failures == 0U ? "true" : "false"
            ) > 0,
            "write EM PIC covariant integration receipt"
        );
        (void)fclose(fp);
    }

    swa_em_pic1d_destroy(&state);
}

int main(void) {
    test_same_step_four_momentum();
    printf(
        "spacewind EM PIC/covariant integration: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
