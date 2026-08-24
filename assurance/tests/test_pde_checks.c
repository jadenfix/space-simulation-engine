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

static void test_poisson_mms(void) {
    swa_discretization_error coarse, medium, fine;
    swa_convergence_result convergence;
    check_true(swa_poisson_manufactured_error(16U, &coarse),
               "Poisson manufactured coarse grid");
    check_true(swa_poisson_manufactured_error(32U, &medium),
               "Poisson manufactured medium grid");
    check_true(swa_poisson_manufactured_error(64U, &fine),
               "Poisson manufactured fine grid");
    check_true(coarse.finite && medium.finite && fine.finite,
               "Poisson manufactured errors finite");
    check_true(coarse.l2_error > medium.l2_error &&
               medium.l2_error > fine.l2_error,
               "Poisson manufactured error decreases monotonically");
    check_true(swa_richardson_three_level(coarse.l2_error,
                                         medium.l2_error,
                                         fine.l2_error,
                                         2.0, 2.0, 0.12,
                                         &convergence),
               "Poisson Richardson audit executes");
    check_true(convergence.passes,
               "Poisson five-point stencil recovers second order");
}

static void test_continuity_mms(void) {
    swa_discretization_error coarse, medium, fine;
    swa_convergence_result convergence;
    check_true(swa_continuity_manufactured_error(32U, 0.4, &coarse),
               "continuity manufactured coarse grid");
    check_true(swa_continuity_manufactured_error(64U, 0.4, &medium),
               "continuity manufactured medium grid");
    check_true(swa_continuity_manufactured_error(128U, 0.4, &fine),
               "continuity manufactured fine grid");
    check_true(coarse.l2_error > medium.l2_error &&
               medium.l2_error > fine.l2_error,
               "continuity residual decreases monotonically");
    check_true(swa_richardson_three_level(coarse.l2_error,
                                         medium.l2_error,
                                         fine.l2_error,
                                         2.0, 2.0, 0.15,
                                         &convergence),
               "continuity Richardson audit executes");
    check_true(convergence.passes,
               "centered continuity residual recovers second order");
}

static void test_maxwell_wave(void) {
    swa_maxwell_wave_result coarse, medium, fine;
    swa_convergence_result convergence;
    check_true(swa_vacuum_maxwell_period(32U, 0.5, &coarse),
               "vacuum Maxwell coarse wave");
    check_true(swa_vacuum_maxwell_period(64U, 0.5, &medium),
               "vacuum Maxwell medium wave");
    check_true(swa_vacuum_maxwell_period(128U, 0.5, &fine),
               "vacuum Maxwell fine wave");
    check_true(coarse.stable && medium.stable && fine.stable,
               "Yee waves remain stable below CFL limit");
    check_true(coarse.l2_electric_error_V_m > medium.l2_electric_error_V_m &&
               medium.l2_electric_error_V_m > fine.l2_electric_error_V_m,
               "Yee wave error decreases monotonically");
    check_true(swa_richardson_three_level(coarse.l2_electric_error_V_m,
                                         medium.l2_electric_error_V_m,
                                         fine.l2_electric_error_V_m,
                                         2.0, 2.0, 0.25,
                                         &convergence),
               "Yee-wave Richardson audit executes");
    check_true(convergence.passes,
               "Yee vacuum wave recovers second-order convergence");
    check_true(fine.rms_amplitude_relative_error < 5e-3,
               "fine Yee wave retains RMS amplitude");
    check_true(!swa_vacuum_maxwell_period(64U, 1.01, &fine),
               "Yee audit rejects super-CFL request");
}

static void test_boris_orbit(void) {
    const double charge_to_mass = SWA_QE / SWA_MP;
    const double magnetic_field = 5e-6;
    const double speed = 2e5;
    swa_boris_orbit_result coarse, medium, fine;
    swa_convergence_result velocity_convergence;
    swa_convergence_result position_convergence;
    check_true(swa_uniform_boris_gyroperiod(32U, charge_to_mass,
                                            magnetic_field, speed, &coarse),
               "Boris coarse gyroperiod");
    check_true(swa_uniform_boris_gyroperiod(64U, charge_to_mass,
                                            magnetic_field, speed, &medium),
               "Boris medium gyroperiod");
    check_true(swa_uniform_boris_gyroperiod(128U, charge_to_mass,
                                            magnetic_field, speed, &fine),
               "Boris fine gyroperiod");
    check_true(coarse.speed_relative_error < 1e-14 &&
               medium.speed_relative_error < 1e-14 &&
               fine.speed_relative_error < 1e-14,
               "Boris magnetic update preserves speed");
    check_true(coarse.velocity_closure_error_mps >
               medium.velocity_closure_error_mps &&
               medium.velocity_closure_error_mps >
               fine.velocity_closure_error_mps,
               "Boris velocity closure improves monotonically");
    check_true(coarse.position_closure_error_m >
               medium.position_closure_error_m &&
               medium.position_closure_error_m >
               fine.position_closure_error_m,
               "Boris position closure improves monotonically");
    check_true(swa_richardson_three_level(coarse.velocity_closure_error_mps,
                                         medium.velocity_closure_error_mps,
                                         fine.velocity_closure_error_mps,
                                         2.0, 2.0, 0.20,
                                         &velocity_convergence),
               "Boris velocity Richardson audit executes");
    check_true(velocity_convergence.passes,
               "Boris velocity phase error recovers second order");
    check_true(swa_richardson_three_level(coarse.position_closure_error_m,
                                         medium.position_closure_error_m,
                                         fine.position_closure_error_m,
                                         2.0, 2.0, 0.25,
                                         &position_convergence),
               "Boris position Richardson audit executes");
    check_true(position_convergence.passes,
               "Boris orbit closure recovers second order");
}

static void test_parker_implicit_relation(void) {
    double residual;
    check_true(swa_parker_implicit_residual(1.0, 1.0, &residual),
               "Parker critical residual evaluates");
    check_true(fabs(residual) < 1e-15,
               "Parker critical point satisfies implicit relation");
    check_true(!swa_parker_implicit_residual(0.0, 1.0, &residual),
               "Parker residual rejects nonpositive radius");
    check_true(!swa_parker_implicit_residual(1.0, 0.0, &residual),
               "Parker residual rejects nonpositive speed");
}

int main(void) {
    FILE *fp;
    test_poisson_mms();
    test_continuity_mms();
    test_maxwell_wave();
    test_boris_orbit();
    test_parker_implicit_relation();
    fp = fopen("output/pde_assurance_receipt.json", "w");
    if (fp != NULL) {
        (void)swa_write_assurance_receipt(fp, checks, failures, receipt_hash);
        (void)fclose(fp);
    } else {
        ++failures;
    }
    printf("spacewind PDE assurance: %zu checks, %zu failures, hash=%016llx\n",
           checks, failures, (unsigned long long)receipt_hash);
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
