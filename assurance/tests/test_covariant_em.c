#include "spacewind/covariant_em.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);
static double worst_field_invariant_relative = 0.0;
static double worst_tensor_covariance_relative = 0.0;
static double worst_four_force_covariance_relative = 0.0;
static double worst_mass_shell_relative = 0.0;

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

static double relative_scalar(double a, double b) {
    return fabs(a - b) / max3(fabs(a), fabs(b), DBL_MIN);
}

static double relative_four(swa_four a, swa_four b) {
    return swa_four_euclidean_norm(swa_four_sub(a, b)) /
           max3(swa_four_euclidean_norm(a),
                swa_four_euclidean_norm(b), DBL_MIN);
}

static double relative_tensor(
    const swa_tensor4 *a,
    const swa_tensor4 *b
) {
    double difference2 = 0.0;
    double a2 = 0.0;
    double b2 = 0.0;
    size_t mu;
    size_t nu;
    for (mu = 0U; mu < 4U; ++mu) {
        for (nu = 0U; nu < 4U; ++nu) {
            const double difference =
                a->component[mu][nu] - b->component[mu][nu];
            difference2 += difference * difference;
            a2 += a->component[mu][nu] * a->component[mu][nu];
            b2 += b->component[mu][nu] * b->component[mu][nu];
        }
    }
    return sqrt(difference2) /
           max3(sqrt(a2), sqrt(b2), DBL_MIN);
}

static swa_vec3 deterministic_direction(size_t i,
                                        unsigned a,
                                        unsigned b,
                                        unsigned c) {
    swa_vec3 direction = swa_v3(
        2.0 * swa_halton(i, a) - 1.0,
        2.0 * swa_halton(i, b) - 1.0,
        2.0 * swa_halton(i, c) - 1.0
    );
    const double norm = swa_vnorm(direction);
    if (!(norm > 1.0e-12)) {
        direction = swa_v3(1.0, 0.0, 0.0);
    } else {
        direction = swa_vscale(direction, 1.0 / norm);
    }
    return direction;
}

static void test_invalid_contracts(void) {
    swa_tensor4 matrix;
    swa_four four;
    swa_em_fields fields;
    swa_covariant_em_diagnostics diagnostics;
    double residual;
    double relative;
    int passes;
    memset(&fields, 0, sizeof(fields));
    check_true(!swa_lorentz_boost_matrix(
                   swa_v3(SWA_C, 0.0, 0.0), &matrix),
               "luminal boost rejected");
    check_true(!swa_lorentz_boost_matrix(
                   swa_v3(NAN, 0.0, 0.0), &matrix),
               "nonfinite boost rejected");
    check_true(!swa_particle_four_momentum(
                   -1.0, swa_v3(0.0, 0.0, 0.0), &four),
               "negative rest mass rejected");
    check_true(!swa_particle_four_momentum(
                   1.0, swa_v3(1.1 * SWA_C, 0.0, 0.0), &four),
               "superluminal particle velocity rejected");
    fields.electric_field_V_m.x = NAN;
    check_true(!swa_em_stress_energy_tensor(
                   &fields, &matrix, &diagnostics),
               "nonfinite field rejected");
    check_true(!swa_audit_particle_mass_shell(
                   1.0, swa_four_make(NAN, 0.0, 0.0, 0.0),
                   1.0e-12, 0.0,
                   &residual, &relative, &passes),
               "nonfinite four-momentum rejected");
}

static void test_field_covariance(void) {
    size_t i;
    for (i = 1U; i <= 2500U; ++i) {
        const swa_vec3 e_direction =
            deterministic_direction(i, 2U, 3U, 5U);
        const swa_vec3 b_direction =
            deterministic_direction(i, 7U, 11U, 13U);
        const swa_vec3 boost_direction =
            deterministic_direction(i, 17U, 19U, 23U);
        const double e_scale = 1.0e2 + 9.0e6 * swa_halton(i, 29U);
        const double b_scale = 1.0e-10 + 2.0e-2 * swa_halton(i, 31U);
        const double boost_speed =
            0.90 * SWA_C * swa_halton(i, 37U);
        const swa_vec3 boost =
            swa_vscale(boost_direction, boost_speed);
        const swa_em_fields fields = {
            swa_vscale(e_direction, e_scale),
            swa_vscale(b_direction, b_scale)
        };
        swa_em_fields transformed_fields;
        swa_tensor4 field_tensor;
        swa_tensor4 transformed_field_tensor;
        swa_tensor4 field_tensor_from_fields;
        swa_tensor4 stress;
        swa_tensor4 transformed_stress;
        swa_tensor4 stress_from_fields;
        swa_covariant_em_diagnostics diagnostics;
        swa_covariant_em_diagnostics transformed_diagnostics;
        double invariant_error;
        double tensor_error;

        check_true(swa_lorentz_transform_em_fields(
                       boost, &fields, &transformed_fields),
                   "field Lorentz transform executes");
        check_true(swa_em_field_tensor(&fields, &field_tensor),
                   "field tensor construction executes");
        check_true(swa_lorentz_transform_tensor(
                       boost, &field_tensor,
                       &transformed_field_tensor),
                   "field tensor boost executes");
        check_true(swa_em_field_tensor(
                       &transformed_fields,
                       &field_tensor_from_fields),
                   "transformed field tensor reconstructs");
        tensor_error = relative_tensor(
            &transformed_field_tensor,
            &field_tensor_from_fields
        );
        worst_tensor_covariance_relative = max2(
            worst_tensor_covariance_relative, tensor_error
        );
        check_true(tensor_error < 2.0e-11,
                   "field tensor transforms covariantly");

        check_true(swa_em_stress_energy_tensor(
                       &fields, &stress, &diagnostics),
                   "stress-energy construction executes");
        check_true(swa_em_stress_energy_tensor(
                       &transformed_fields,
                       &stress_from_fields,
                       &transformed_diagnostics),
                   "transformed stress-energy reconstructs");
        check_true(swa_lorentz_transform_tensor(
                       boost, &stress, &transformed_stress),
                   "stress-energy tensor boost executes");
        tensor_error = relative_tensor(
            &transformed_stress, &stress_from_fields
        );
        worst_tensor_covariance_relative = max2(
            worst_tensor_covariance_relative, tensor_error
        );
        check_true(tensor_error < 4.0e-11,
                   "stress-energy transforms covariantly");

        invariant_error = max2(
            relative_scalar(
                diagnostics.invariant_B2_minus_E2_over_c2_T2,
                transformed_diagnostics.
                    invariant_B2_minus_E2_over_c2_T2
            ),
            relative_scalar(
                diagnostics.invariant_E_dot_B_over_c,
                transformed_diagnostics.invariant_E_dot_B_over_c
            )
        );
        worst_field_invariant_relative = max2(
            worst_field_invariant_relative, invariant_error
        );
        check_true(invariant_error < 8.0e-9,
                   "electromagnetic scalar invariants are boost invariant");
        check_true(diagnostics.dominant_energy_condition &&
                   transformed_diagnostics.dominant_energy_condition,
                   "electromagnetic dominant energy condition holds");
        check_true(fabs(diagnostics.tensor_trace_J_m3) <=
                   5.0e-12 * max2(diagnostics.energy_density_J_m3, 1.0),
                   "electromagnetic stress-energy trace vanishes");
    }
}

static void test_four_force_covariance(void) {
    size_t i;
    for (i = 1U; i <= 1500U; ++i) {
        const swa_vec3 e = swa_vscale(
            deterministic_direction(i, 2U, 5U, 7U),
            1.0e3 + 2.0e6 * swa_halton(i, 11U)
        );
        const swa_vec3 b = swa_vscale(
            deterministic_direction(i, 13U, 17U, 19U),
            1.0e-9 + 5.0e-3 * swa_halton(i, 23U)
        );
        const swa_vec3 current = swa_vscale(
            deterministic_direction(i, 29U, 31U, 37U),
            1.0e-3 + 4.0e2 * swa_halton(i, 41U)
        );
        const double rho =
            -2.0e-6 + 4.0e-6 * swa_halton(i, 43U);
        const swa_vec3 boost = swa_vscale(
            deterministic_direction(i, 47U, 53U, 59U),
            0.82 * SWA_C * swa_halton(i, 61U)
        );
        const swa_em_fields fields = {e, b};
        swa_em_fields transformed_fields;
        swa_four current_four;
        swa_four transformed_current_four;
        swa_four force;
        swa_four transformed_force;
        swa_four expected_force;
        double transformed_rho;
        swa_vec3 transformed_current;
        double error;

        check_true(swa_four_current(rho, current, &current_four),
                   "four-current construction executes");
        check_true(swa_lorentz_transform_four(
                       boost, current_four,
                       &transformed_current_four),
                   "four-current boost executes");
        transformed_rho =
            transformed_current_four.component[0] / SWA_C;
        transformed_current = swa_v3(
            transformed_current_four.component[1],
            transformed_current_four.component[2],
            transformed_current_four.component[3]
        );
        check_true(swa_lorentz_transform_em_fields(
                       boost, &fields, &transformed_fields),
                   "four-force field transform executes");
        check_true(swa_lorentz_force_density(
                       &fields, rho, current, &force),
                   "four-force density construction executes");
        check_true(swa_lorentz_transform_four(
                       boost, force, &expected_force),
                   "four-force density boost executes");
        check_true(swa_lorentz_force_density(
                       &transformed_fields, transformed_rho,
                       transformed_current, &transformed_force),
                   "transformed four-force density reconstructs");
        error = relative_four(expected_force, transformed_force);
        worst_four_force_covariance_relative = max2(
            worst_four_force_covariance_relative, error
        );
        check_true(error < 5.0e-11,
                   "Lorentz four-force density is covariant");
        check_true(relative_scalar(
                       force.component[0],
                       swa_vdot(e, current) / SWA_C) < 2.0e-14,
                   "four-force time component equals J dot E over c");
        check_true(relative_scalar(
                       force.component[1],
                       rho * e.x + swa_vcross(current, b).x) < 2.0e-14,
                   "four-force spatial x equals Lorentz force density");
    }
}

static void test_particle_mass_shell(void) {
    size_t i;
    for (i = 1U; i <= 2000U; ++i) {
        const double mass = 1.0e-30 +
                            1.0e4 * swa_halton(i, 2U);
        const swa_vec3 velocity = swa_vscale(
            deterministic_direction(i, 3U, 5U, 7U),
            0.97 * SWA_C * swa_halton(i, 11U)
        );
        const swa_vec3 boost = swa_vscale(
            deterministic_direction(i, 13U, 17U, 19U),
            0.88 * SWA_C * swa_halton(i, 23U)
        );
        swa_four momentum;
        swa_four transformed;
        swa_four recovered;
        double residual;
        double relative;
        int passes;
        check_true(swa_particle_four_momentum(
                       mass, velocity, &momentum),
                   "particle four-momentum constructs");
        check_true(swa_audit_particle_mass_shell(
                       mass, momentum, 5.0e-13, 1.0e-50,
                       &residual, &relative, &passes),
                   "particle mass-shell audit executes");
        worst_mass_shell_relative = max2(
            worst_mass_shell_relative, relative
        );
        check_true(passes, "particle mass shell closes");
        check_true(swa_lorentz_transform_four(
                       boost, momentum, &transformed),
                   "particle four-momentum boost executes");
        check_true(swa_audit_particle_mass_shell(
                       mass, transformed, 3.0e-12, 1.0e-50,
                       &residual, &relative, &passes),
                   "boosted particle mass-shell audit executes");
        worst_mass_shell_relative = max2(
            worst_mass_shell_relative, relative
        );
        check_true(passes, "boosted particle mass shell closes");
        check_true(swa_lorentz_transform_four(
                       swa_vscale(boost, -1.0),
                       transformed, &recovered),
                   "inverse particle boost executes");
        check_true(relative_four(recovered, momentum) < 2.0e-12,
                   "inverse boost recovers particle four-momentum");
    }
}

static void test_plane_wave_and_ledger(void) {
    const double electric = 4.0e5;
    const swa_em_fields wave = {
        {0.0, electric, 0.0},
        {0.0, 0.0, electric / SWA_C}
    };
    swa_tensor4 stress;
    swa_covariant_em_diagnostics diagnostics;
    swa_four_momentum_ledger ledger;
    swa_four_momentum_ledger_result result;
    memset(&ledger, 0, sizeof(ledger));
    check_true(swa_em_stress_energy_tensor(
                   &wave, &stress, &diagnostics),
               "plane-wave stress-energy constructs");
    check_true(relative_scalar(stress.component[0][0],
                               stress.component[0][1]) < 2.0e-14,
               "plane-wave energy density equals c momentum density");
    check_true(diagnostics.dominant_energy_ratio > 1.0 - 1.0e-14 &&
               diagnostics.dominant_energy_ratio <= 1.0 + 1.0e-14,
               "plane wave saturates dominant energy bound");
    check_true(fabs(diagnostics.
                       invariant_B2_minus_E2_over_c2_T2) < 1.0e-24,
               "plane-wave first field invariant vanishes");
    check_true(fabs(diagnostics.invariant_E_dot_B_over_c) < 1.0e-24,
               "plane-wave second field invariant vanishes");

    ledger.field_initial_Ns = swa_four_make(10.0, 3.0, 2.0, 1.0);
    ledger.matter_initial_Ns = swa_four_make(5.0, -1.0, 0.5, 0.0);
    ledger.external_impulse_Ns = swa_four_make(1.0, 0.5, -0.25, 0.75);
    ledger.momentum_in_Ns = swa_four_make(0.4, 0.2, 0.1, -0.1);
    ledger.momentum_out_Ns = swa_four_make(0.1, 0.05, 0.0, 0.0);
    ledger.field_final_Ns = swa_four_make(8.0, 1.5, 1.0, 0.5);
    ledger.matter_final_Ns = swa_four_add(
        swa_four_add(
            swa_four_add(ledger.field_initial_Ns,
                         ledger.matter_initial_Ns),
            swa_four_add(
                ledger.external_impulse_Ns,
                swa_four_sub(ledger.momentum_in_Ns,
                             ledger.momentum_out_Ns)
            )
        ),
        swa_four_scale(ledger.field_final_Ns, -1.0)
    );
    check_true(swa_audit_four_momentum_ledger(
                   &ledger, 1.0e-14, 1.0e-14, &result) &&
               result.passes,
               "closed four-momentum ledger passes");
    ledger.matter_final_Ns.component[2] += 1.0e-3;
    check_true(swa_audit_four_momentum_ledger(
                   &ledger, 1.0e-14, 1.0e-14, &result) &&
               !result.passes,
               "corrupted four-momentum ledger fails");
}

int main(void) {
    FILE *fp;
    test_invalid_contracts();
    test_field_covariance();
    test_four_force_covariance();
    test_particle_mass_shell();
    test_plane_wave_and_ledger();
    fp = fopen("output/covariant_em_receipt.json", "w");
    check_true(fp != NULL, "open covariant EM receipt");
    if (fp != NULL) {
        check_true(swa_write_covariant_em_receipt(
                       fp, checks, failures, transcript_hash,
                       worst_field_invariant_relative,
                       worst_tensor_covariance_relative,
                       worst_four_force_covariance_relative,
                       worst_mass_shell_relative),
                   "write covariant EM receipt");
        (void)fclose(fp);
    }
    printf(
        "spacewind covariant EM assurance: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
