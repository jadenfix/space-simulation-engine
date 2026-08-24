#include "spacewind/assurance.h"
#include "spacewind/evidence.h"
#include "spacewind/experiment.h"
#include "spacewind/invariants.h"
#include "spacewind/pde_checks.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t receipt_hash = UINT64_C(14695981039346656037);

static void check_true(int condition, const char *name) {
    ++checks;
    receipt_hash ^= swa_fnv1a64(name, strlen(name));
    receipt_hash *= UINT64_C(1099511628211);
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}

static int near(double a, double b, double relative, double absolute) {
    return fabs(a - b) <= fmax(absolute, relative * fmax(fabs(a), fabs(b)));
}

static swa_plasma_state plasma_from_index(size_t i) {
    swa_plasma_state s;
    memset(&s, 0, sizeof(s));
    s.ion_density_m3 = pow(10.0, 4.0 + 6.0 * swa_halton(i, 2U));
    s.electron_density_m3 = s.ion_density_m3;
    s.ion_temperature_K = pow(10.0, 3.0 + 4.0 * swa_halton(i, 3U));
    s.electron_temperature_K = pow(10.0, 3.0 + 4.0 * swa_halton(i, 5U));
    s.ion_mass_kg = SWA_MP;
    s.ion_charge_C = SWA_QE;
    s.neutral_density_m3 = pow(10.0, -2.0 + 4.0 * swa_halton(i, 7U));
    s.collision_cross_section_m2 = 1e-19;
    s.electrode_radius_m = pow(10.0, -6.0 + 3.0 * swa_halton(i, 11U));
    s.interaction_length_m = pow(10.0, 0.0 + 5.0 * swa_halton(i, 13U));
    s.electrode_voltage_V = 100.0 + 100000.0 * swa_halton(i, 17U);
    s.circuit_capacitance_F = pow(10.0, -9.0 + 5.0 * swa_halton(i, 19U));
    s.available_power_W = 1.0 + 10000.0 * swa_halton(i, 23U);
    s.collection_current_A = 1e-6 + 0.1 * swa_halton(i, 29U);
    s.bulk_velocity_mps = swa_v3(1e4 + 9e5 * swa_halton(i, 31U),
                                 1e4 * (swa_halton(i, 37U) - 0.5),
                                 1e4 * (swa_halton(i, 41U) - 0.5));
    s.magnetic_field_T = swa_v3(1e-10 + 2e-8 * swa_halton(i, 43U),
                                1e-8 * (swa_halton(i, 47U) - 0.5),
                                1e-8 * (swa_halton(i, 53U) - 0.5));
    return s;
}

static void test_plasma_scaling(void) {
    size_t i;
    for (i = 1U; i <= 5000U; ++i) {
        swa_plasma_state base = plasma_from_index(i);
        swa_plasma_state density_scaled = base;
        swa_plasma_state temperature_scaled = base;
        swa_plasma_regime a, b, c;
        density_scaled.ion_density_m3 *= 4.0;
        density_scaled.electron_density_m3 *= 4.0;
        temperature_scaled.electron_temperature_K *= 4.0;
        check_true(swa_compute_plasma_regime(&base, &a) && a.finite,
                   "random plasma regime finite");
        check_true(swa_compute_plasma_regime(&density_scaled, &b) && b.finite,
                   "density-scaled plasma regime finite");
        check_true(swa_compute_plasma_regime(&temperature_scaled, &c) && c.finite,
                   "temperature-scaled plasma regime finite");
        check_true(near(b.debye_length_m, 0.5 * a.debye_length_m, 3e-13, 1e-30),
                   "Debye length follows inverse square-root density law");
        check_true(near(b.ion_plasma_frequency_rad_s,
                        2.0 * a.ion_plasma_frequency_rad_s, 3e-13, 1e-30),
                   "ion plasma frequency follows square-root density law");
        check_true(near(b.ion_inertial_length_m,
                        0.5 * a.ion_inertial_length_m, 3e-13, 1e-30),
                   "ion inertial length follows inverse square-root density law");
        check_true(near(b.alfven_speed_mps, 0.5 * a.alfven_speed_mps, 4e-13, 1e-30),
                   "Alfven speed follows inverse square-root density law");
        check_true(near(c.debye_length_m, 2.0 * a.debye_length_m, 3e-13, 1e-30),
                   "Debye length follows square-root electron-temperature law");
        check_true(a.debye_length_m > 0.0 && a.fast_magnetosonic_speed_mps > 0.0 &&
                   a.power_margin > 0.0,
                   "random plasma positive physical scales");
    }
}

static void test_force_bound_monotonicity(void) {
    size_t i;
    for (i = 1U; i <= 5000U; ++i) {
        swa_force_bound a, b, c;
        const double irradiance = 1.0 + 2000.0 * swa_halton(i, 2U);
        const double area = 0.1 + 1e5 * swa_halton(i, 3U);
        const double reflectivity = swa_halton(i, 5U);
        check_true(swa_audit_photon_force(irradiance, area, reflectivity,
                                          0.0, 0.0, &a),
                   "photon bound base evaluates");
        check_true(swa_audit_photon_force(2.0 * irradiance, area, reflectivity,
                                          0.0, 0.0, &b),
                   "photon bound irradiance scaling evaluates");
        check_true(swa_audit_photon_force(irradiance, 3.0 * area, reflectivity,
                                          0.0, 0.0, &c),
                   "photon bound area scaling evaluates");
        check_true(near(b.upper_bound_force_N, 2.0 * a.upper_bound_force_N,
                        2e-15, 1e-30),
                   "photon force bound linear in irradiance");
        check_true(near(c.upper_bound_force_N, 3.0 * a.upper_bound_force_N,
                        2e-15, 1e-30),
                   "photon force bound linear in area");
        {
            const double rho = 1e-22 + 1e-18 * swa_halton(i, 7U);
            const double speed = 1e4 + 1e6 * swa_halton(i, 11U);
            check_true(swa_audit_plasma_momentum_force(rho, speed, area, 2.0,
                                                        0.0, 0.0, &a),
                       "plasma force bound base evaluates");
            check_true(swa_audit_plasma_momentum_force(rho, 2.0 * speed, area, 2.0,
                                                        0.0, 0.0, &b),
                       "plasma force speed scaling evaluates");
            check_true(near(b.upper_bound_force_N, 4.0 * a.upper_bound_force_N,
                            3e-15, 1e-30),
                       "plasma momentum bound quadratic in relative speed");
        }
    }
}

static swa_claim_gate_input progressive_input(unsigned level,
                                              swa_reversal_result *reversal,
                                              swa_replication_result *replication,
                                              swa_cycle_result *cycle) {
    swa_claim_gate_input i;
    memset(&i, 0, sizeof(i));
    i.mission_net_gain_interval_J = swa_interval_make(-1.0, 1.0);
    i.required_deployment_survival_lower_bound = 0.99;
    if (level >= 1U) {
        i.dimensional_contract_passes = 1;
        i.conservation_ledgers_pass = 1;
        i.deterministic_replay_passes = 1;
        i.memory_and_undefined_behavior_checks_pass = 1;
    }
    if (level >= 2U) {
        i.convergence_study_passes = 1;
        i.independent_implementation_agrees = 1;
    }
    if (level >= 3U) {
        i.similarity_contract_passes = 1;
        i.available_to_required_power_ratio = 1.2;
        i.thermal_margin_fraction = 0.2;
        i.structural_envelope_passes = 1;
        i.fault_containment_passes = 1;
        reversal->finite = 1;
        reversal->passes = 1;
        reversal->reversed_effect_mean_N = 1.0;
        reversal->effect_interval_N = swa_interval_make(0.8, 1.2);
        replication->passing_replicates = 2U;
        replication->same_sign = 1;
        replication->interval_overlap = 1;
        replication->passes = 1;
        i.reversal_experiment = reversal;
        i.independent_replication = replication;
    }
    if (level >= 4U) {
        cycle->numerical_cycle_passes = 1;
        cycle->robust_positive_gain = 1;
        cycle->physical_promotion_passes = 1;
        i.closed_cycle = cycle;
        i.mission_net_gain_interval_J = swa_interval_make(1.0, 2.0);
    }
    if (level >= 5U) {
        i.deployment_survival_lower_bound = 0.995;
        i.flight_demonstration_passes = 1;
    }
    return i;
}

static void test_claim_gate_monotonicity(void) {
    size_t sample;
    for (sample = 1U; sample <= 5000U; ++sample) {
        swa_reversal_result reversal;
        swa_replication_result replication;
        swa_cycle_result cycle;
        swa_claim_gate_result previous;
        unsigned level;
        memset(&reversal, 0, sizeof(reversal));
        memset(&replication, 0, sizeof(replication));
        memset(&cycle, 0, sizeof(cycle));
        memset(&previous, 0, sizeof(previous));
        for (level = 0U; level <= 5U; ++level) {
            swa_claim_gate_input input = progressive_input(level, &reversal,
                                                            &replication, &cycle);
            swa_claim_gate_result current;
            check_true(swa_evaluate_claim_gate(&input, &current),
                       "progressive claim gate evaluates");
            check_true(current.highest_tier >= previous.highest_tier,
                       "adding evidence never lowers claim tier");
            check_true(current.highest_tier <= (swa_claim_tier)level,
                       "claim tier never outruns supplied stage");
            previous = current;
        }
    }
}

static void test_invalid_inputs_fail_closed(void) {
    swa_plasma_state plasma = plasma_from_index(1U);
    swa_plasma_regime regime;
    swa_force_bound bound;
    swa_four_vector u;
    swa_discretization_error error;
    swa_claim_gate_result claim;
    check_true(!swa_compute_plasma_regime(NULL, &regime),
               "null plasma pointer rejected");
    plasma.ion_density_m3 = NAN;
    check_true(!swa_compute_plasma_regime(&plasma, &regime),
               "NaN plasma density rejected");
    check_true(!swa_audit_photon_force(-1.0, 1.0, 1.0, 0.0, 0.0, &bound),
               "negative irradiance rejected");
    check_true(!swa_audit_photon_force(1.0, 1.0, 1.1, 0.0, 0.0, &bound),
               "reflectivity above one rejected");
    check_true(!swa_four_velocity_from_three_velocity(swa_v3(INFINITY, 0.0, 0.0), &u),
               "infinite velocity rejected");
    check_true(!swa_poisson_manufactured_error(2U, &error),
               "undersized Poisson grid rejected");
    check_true(!swa_continuity_manufactured_error(32U, 0.0, &error),
               "zero continuity CFL rejected");
    check_true(!swa_evaluate_claim_gate(NULL, &claim),
               "null claim gate rejected");
}

int main(void) {
    FILE *fp;
    test_plasma_scaling();
    test_force_bound_monotonicity();
    test_claim_gate_monotonicity();
    test_invalid_inputs_fail_closed();
    fp = fopen("output/robustness_assurance_receipt.json", "w");
    if (fp != NULL) {
        (void)swa_write_assurance_receipt(fp, checks, failures, receipt_hash);
        (void)fclose(fp);
    } else {
        ++failures;
    }
    printf("spacewind robustness assurance: %zu checks, %zu failures, hash=%016llx\n",
           checks, failures, (unsigned long long)receipt_hash);
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
