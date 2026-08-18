#include "spacewind/charge_deposition.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);

static void check_true(int condition, const char *name) {
    ++checks;
    transcript_hash ^= swa_fnv1a64(name, strlen(name));
    transcript_hash *= UINT64_C(1099511628211);
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}

static int near(double a, double b, double relative, double absolute) {
    return fabs(a - b) <=
           fmax(absolute, relative * fmax(fabs(a), fabs(b)));
}

static swa_periodic_cic_transport_input base_input(
    size_t cells,
    size_t particles,
    const double *initial,
    const double *final,
    const double *charge
) {
    swa_periodic_cic_transport_input input;
    memset(&input, 0, sizeof(input));
    input.cells = cells;
    input.particles = particles;
    input.domain_length_m = 8.0;
    input.cross_section_area_m2 = 0.5;
    input.dt_s = 0.2;
    input.initial_position_unwrapped_m = initial;
    input.final_position_unwrapped_m = final;
    input.particle_charge_C = charge;
    input.continuity_relative_tolerance = 2.0e-12;
    input.continuity_absolute_tolerance_A_m3 = 1.0e-18;
    input.continuity_global_absolute_tolerance_C_m = 1.0e-18;
    input.charge_absolute_tolerance_C = 1.0e-24;
    input.mean_current_relative_tolerance = 2.0e-12;
    input.mean_current_absolute_tolerance_A_m2 = 1.0e-18;
    return input;
}

static void test_zero_motion(void) {
    enum { CELLS = 32, PARTICLES = 128 };
    double initial[PARTICLES];
    double final[PARTICLES];
    double charge[PARTICLES];
    double rho_initial[CELLS];
    double rho_final[CELLS];
    double current[CELLS + 1];
    swa_periodic_cic_transport_input input;
    swa_periodic_cic_transport_buffers buffers;
    swa_periodic_cic_transport_result result;
    size_t i;

    for (i = 0U; i < PARTICLES; ++i) {
        initial[i] = -8.0 + 24.0 * swa_halton(i + 1U, 2U);
        final[i] = initial[i];
        charge[i] = 1.0e-15 *
                    (0.5 + swa_halton(i + 1U, 3U));
    }
    input = base_input(CELLS, PARTICLES,
                       initial, final, charge);
    buffers.rho_initial_C_m3 = rho_initial;
    buffers.rho_final_C_m3 = rho_final;
    buffers.rho_capacity = CELLS;
    buffers.current_face_A_m2 = current;
    buffers.current_face_capacity = CELLS + 1U;

    check_true(
        swa_deposit_periodic_cic_transport(
            &input, &buffers, &result
        ),
        "zero-motion deposition executes"
    );
    check_true(result.passes,
               "zero-motion deposition passes");
    check_true(near(result.target_mean_current_A_m2,
                    0.0, 0.0, 1.0e-30),
               "zero motion has zero target current");
    for (i = 0U; i <= CELLS; ++i) {
        check_true(near(current[i], 0.0, 0.0, 1.0e-25),
                   "zero motion has zero reconstructed current");
    }
    for (i = 0U; i < CELLS; ++i) {
        check_true(near(rho_initial[i], rho_final[i],
                        1.0e-15, 1.0e-25),
                   "zero motion preserves deposited density");
    }
}

static void test_multiwrap_transport(void) {
    enum { CELLS = 16, PARTICLES = 1 };
    double initial[PARTICLES] = {1.7};
    double final[PARTICLES] = {17.7};
    double charge[PARTICLES] = {1.0e-9};
    double rho_initial[CELLS];
    double rho_final[CELLS];
    double current[CELLS + 1];
    swa_periodic_cic_transport_input input =
        base_input(CELLS, PARTICLES,
                   initial, final, charge);
    swa_periodic_cic_transport_buffers buffers;
    swa_periodic_cic_transport_result result;
    const double expected = 2.0e-8;
    size_t i;

    buffers.rho_initial_C_m3 = rho_initial;
    buffers.rho_final_C_m3 = rho_final;
    buffers.rho_capacity = CELLS;
    buffers.current_face_A_m2 = current;
    buffers.current_face_capacity = CELLS + 1U;
    input.charge_absolute_tolerance_C = 1.0e-20;

    check_true(
        swa_deposit_periodic_cic_transport(
            &input, &buffers, &result
        ),
        "multiwrap deposition executes"
    );
    check_true(result.passes,
               "multiwrap deposition passes");
    check_true(near(result.transported_charge_turns_C,
                    2.0e-9, 1.0e-14, 1.0e-24),
               "multiwrap transported charge-turns are retained");
    check_true(near(result.target_mean_current_A_m2,
                    expected, 1.0e-14, 1.0e-24),
               "multiwrap target mean current is physical");
    for (i = 0U; i <= CELLS; ++i) {
        check_true(near(current[i], expected,
                        1.0e-13, 1.0e-22),
                   "multiwrap harmonic current is uniform");
    }
}

static void test_random_transport_and_reversal(void) {
    enum { CELLS = 64, PARTICLES = 1024 };
    double initial[PARTICLES];
    double final[PARTICLES];
    double charge[PARTICLES];
    double rho_initial[CELLS];
    double rho_final[CELLS];
    double current[CELLS + 1];
    double reverse_rho_initial[CELLS];
    double reverse_rho_final[CELLS];
    double reverse_current[CELLS + 1];
    swa_periodic_cic_transport_input input;
    swa_periodic_cic_transport_input reverse_input;
    swa_periodic_cic_transport_buffers buffers;
    swa_periodic_cic_transport_buffers reverse_buffers;
    swa_periodic_cic_transport_result result;
    swa_periodic_cic_transport_result reverse_result;
    size_t i;

    for (i = 0U; i < PARTICLES; ++i) {
        const double base =
            -16.0 + 40.0 * swa_halton(i + 1U, 2U);
        const double displacement =
            40.0 * (swa_halton(i + 1U, 3U) - 0.5);
        initial[i] = base;
        final[i] = base + displacement;
        charge[i] = ((i & 1U) == 0U ? -1.0 : 1.0) *
                    1.0e-15 *
                    (0.5 + swa_halton(i + 1U, 5U));
    }
    input = base_input(CELLS, PARTICLES,
                       initial, final, charge);
    input.charge_absolute_tolerance_C = 5.0e-25;
    buffers.rho_initial_C_m3 = rho_initial;
    buffers.rho_final_C_m3 = rho_final;
    buffers.rho_capacity = CELLS;
    buffers.current_face_A_m2 = current;
    buffers.current_face_capacity = CELLS + 1U;

    check_true(
        swa_deposit_periodic_cic_transport(
            &input, &buffers, &result
        ),
        "random transport deposition executes"
    );
    check_true(result.passes,
               "random transport closes charge and current");
    check_true(result.maximum_particle_displacement_cells > 100.0,
               "random transport includes long multi-cell motion");
    check_true(result.continuity.residual_relative_rms < 2.0e-12,
               "random transport local continuity closes");

    reverse_input = base_input(CELLS, PARTICLES,
                               final, initial, charge);
    reverse_input.charge_absolute_tolerance_C =
        input.charge_absolute_tolerance_C;
    reverse_buffers.rho_initial_C_m3 = reverse_rho_initial;
    reverse_buffers.rho_final_C_m3 = reverse_rho_final;
    reverse_buffers.rho_capacity = CELLS;
    reverse_buffers.current_face_A_m2 = reverse_current;
    reverse_buffers.current_face_capacity = CELLS + 1U;
    check_true(
        swa_deposit_periodic_cic_transport(
            &reverse_input, &reverse_buffers, &reverse_result
        ),
        "time-reversed deposition executes"
    );
    check_true(reverse_result.passes,
               "time-reversed deposition passes");
    for (i = 0U; i < CELLS; ++i) {
        check_true(
            near(rho_initial[i], reverse_rho_final[i],
                 5.0e-14, 1.0e-24),
            "time reversal swaps initial and final charge"
        );
        check_true(
            near(rho_final[i], reverse_rho_initial[i],
                 5.0e-14, 1.0e-24),
            "time reversal swaps final and initial charge"
        );
    }
    for (i = 0U; i <= CELLS; ++i) {
        check_true(
            near(current[i], -reverse_current[i],
                 5.0e-12, 1.0e-24),
            "time reversal negates reconstructed current"
        );
    }
}

static void test_permutation_stability(void) {
    enum { CELLS = 48, PARTICLES = 257 };
    double initial[PARTICLES];
    double final[PARTICLES];
    double charge[PARTICLES];
    double initial_reversed[PARTICLES];
    double final_reversed[PARTICLES];
    double charge_reversed[PARTICLES];
    double rho_a[CELLS];
    double rho_b[CELLS];
    double current_a[CELLS + 1];
    double rho_a_rev[CELLS];
    double rho_b_rev[CELLS];
    double current_rev[CELLS + 1];
    swa_periodic_cic_transport_input input;
    swa_periodic_cic_transport_input reversed;
    swa_periodic_cic_transport_buffers buffers;
    swa_periodic_cic_transport_buffers reverse_buffers;
    swa_periodic_cic_transport_result result;
    swa_periodic_cic_transport_result reverse_result;
    size_t i;

    for (i = 0U; i < PARTICLES; ++i) {
        initial[i] = 20.0 * swa_halton(i + 1U, 2U) - 6.0;
        final[i] = initial[i] +
                   12.0 * (swa_halton(i + 1U, 3U) - 0.5);
        charge[i] = 1.0e-14 *
                    (0.25 + swa_halton(i + 1U, 7U));
    }
    for (i = 0U; i < PARTICLES; ++i) {
        const size_t source = PARTICLES - 1U - i;
        initial_reversed[i] = initial[source];
        final_reversed[i] = final[source];
        charge_reversed[i] = charge[source];
    }

    input = base_input(CELLS, PARTICLES,
                       initial, final, charge);
    input.charge_absolute_tolerance_C = 1.0e-23;
    buffers.rho_initial_C_m3 = rho_a;
    buffers.rho_final_C_m3 = rho_b;
    buffers.rho_capacity = CELLS;
    buffers.current_face_A_m2 = current_a;
    buffers.current_face_capacity = CELLS + 1U;
    check_true(
        swa_deposit_periodic_cic_transport(
            &input, &buffers, &result
        ),
        "forward-order deposition executes"
    );

    reversed = base_input(CELLS, PARTICLES,
                          initial_reversed, final_reversed,
                          charge_reversed);
    reversed.charge_absolute_tolerance_C =
        input.charge_absolute_tolerance_C;
    reverse_buffers.rho_initial_C_m3 = rho_a_rev;
    reverse_buffers.rho_final_C_m3 = rho_b_rev;
    reverse_buffers.rho_capacity = CELLS;
    reverse_buffers.current_face_A_m2 = current_rev;
    reverse_buffers.current_face_capacity = CELLS + 1U;
    check_true(
        swa_deposit_periodic_cic_transport(
            &reversed, &reverse_buffers, &reverse_result
        ),
        "reverse-order deposition executes"
    );
    check_true(result.passes && reverse_result.passes,
               "both particle orders pass");
    for (i = 0U; i < CELLS; ++i) {
        check_true(near(rho_a[i], rho_a_rev[i],
                        2.0e-14, 1.0e-23),
                   "initial density is permutation stable");
        check_true(near(rho_b[i], rho_b_rev[i],
                        2.0e-14, 1.0e-23),
                   "final density is permutation stable");
    }
    for (i = 0U; i <= CELLS; ++i) {
        check_true(near(current_a[i], current_rev[i],
                        2.0e-12, 1.0e-22),
                   "current is permutation stable");
    }
}

static void test_invalid_inputs(void) {
    enum { CELLS = 8, PARTICLES = 2 };
    double initial[PARTICLES] = {0.0, 1.0};
    double final[PARTICLES] = {0.1, 1.1};
    double charge[PARTICLES] = {1.0, -1.0};
    double rho_initial[CELLS];
    double rho_final[CELLS];
    double current[CELLS + 1];
    swa_periodic_cic_transport_input input =
        base_input(CELLS, PARTICLES,
                   initial, final, charge);
    swa_periodic_cic_transport_buffers buffers;
    swa_periodic_cic_transport_result result;

    buffers.rho_initial_C_m3 = rho_initial;
    buffers.rho_final_C_m3 = rho_final;
    buffers.rho_capacity = CELLS;
    buffers.current_face_A_m2 = current;
    buffers.current_face_capacity = CELLS + 1U;

    input.dt_s = 0.0;
    check_true(
        !swa_deposit_periodic_cic_transport(
            &input, &buffers, &result
        ),
        "zero transport timestep rejected"
    );
    input = base_input(CELLS, PARTICLES,
                       initial, final, charge);
    buffers.current_face_capacity = CELLS;
    check_true(
        !swa_deposit_periodic_cic_transport(
            &input, &buffers, &result
        ),
        "undersized current buffer rejected"
    );
    buffers.current_face_capacity = CELLS + 1U;
    final[0] = NAN;
    check_true(
        !swa_deposit_periodic_cic_transport(
            &input, &buffers, &result
        ),
        "nonfinite particle position rejected"
    );
}

int main(void) {
    enum { CELLS = 32, PARTICLES = 64 };
    double initial[PARTICLES];
    double final[PARTICLES];
    double charge[PARTICLES];
    double rho_initial[CELLS];
    double rho_final[CELLS];
    double current[CELLS + 1];
    swa_periodic_cic_transport_input input;
    swa_periodic_cic_transport_buffers buffers;
    swa_periodic_cic_transport_result result;
    FILE *fp;
    size_t i;

    test_zero_motion();
    test_multiwrap_transport();
    test_random_transport_and_reversal();
    test_permutation_stability();
    test_invalid_inputs();

    for (i = 0U; i < PARTICLES; ++i) {
        initial[i] = 8.0 * swa_halton(i + 1U, 2U);
        final[i] = initial[i] +
                   6.0 * (swa_halton(i + 1U, 3U) - 0.5);
        charge[i] = 1.0e-14 *
                    (0.2 + swa_halton(i + 1U, 5U));
    }
    input = base_input(CELLS, PARTICLES,
                       initial, final, charge);
    input.charge_absolute_tolerance_C = 1.0e-23;
    buffers.rho_initial_C_m3 = rho_initial;
    buffers.rho_final_C_m3 = rho_final;
    buffers.rho_capacity = CELLS;
    buffers.current_face_A_m2 = current;
    buffers.current_face_capacity = CELLS + 1U;
    check_true(
        swa_deposit_periodic_cic_transport(
            &input, &buffers, &result
        ),
        "receipt deposition executes"
    );
    check_true(result.passes,
               "receipt deposition passes");
    fp = fopen("output/charge_deposition_receipt.json", "w");
    check_true(fp != NULL, "open charge-deposition receipt");
    if (fp != NULL) {
        check_true(
            swa_write_periodic_cic_transport_receipt(
                fp, &input, &result
            ),
            "write charge-deposition receipt"
        );
        (void)fclose(fp);
    }

    printf(
        "spacewind charge deposition assurance: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
