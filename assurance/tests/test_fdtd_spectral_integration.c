#include "spacewind/constants.h"
#include "spacewind/maxwell1d_reference.h"
#include "spacewind/plasma.h"

#include <float.h>
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

static double right_wave(double x_m, double time_s,
                         double length_m, unsigned mode) {
    const double phase = 2.0 * SWA_PI * (double)mode *
                         (x_m - SWA_C * time_s) / length_m;
    return sin(phase);
}

static int run_level(
    size_t cells,
    unsigned mode,
    double *relative_l2_error,
    double *spectral_relative_l2_error,
    double *relative_energy_change,
    double *courant
) {
    const double length_m = 1.0;
    const double dx_m = length_m / (double)cells;
    const double sigma = 0.5;
    const double dt_s = sigma * dx_m / SWA_C;
    const size_t steps = cells;
    const double final_time_s = (double)steps * dt_s;
    sw_fdtd1d grid;
    swa_maxwell1d_state reference;
    swa_maxwell1d_step_result reference_result;
    double initial_energy;
    double final_energy;
    double error_sum = 0.0;
    double spectral_error_sum = 0.0;
    double exact_sum = 0.0;
    size_t i;
    size_t step;

    memset(&grid, 0, sizeof(grid));
    memset(&reference, 0, sizeof(reference));
    if (relative_l2_error == NULL || spectral_relative_l2_error == NULL ||
        relative_energy_change == NULL || courant == NULL ||
        cells < 5U || (cells & 1U) == 0U ||
        !sw_fdtd1d_init(&grid, cells, length_m, dt_s) ||
        !swa_maxwell1d_init(&reference, cells, length_m, 1.0)) {
        sw_fdtd1d_destroy(&grid);
        swa_maxwell1d_destroy(&reference);
        return 0;
    }

    for (i = 0U; i < cells; ++i) {
        const double electric_x_m = ((double)i + 0.5) * dx_m;
        const double magnetic_x_m = ((double)i + 1.0) * dx_m;
        const double electric = right_wave(
            electric_x_m, 0.0, length_m, mode
        );
        const double magnetic_staggered = right_wave(
            magnetic_x_m, -0.5 * dt_s, length_m, mode
        ) / SWA_C;
        grid.electric_y_v_m[i] = electric;
        grid.magnetic_z_t[i] = magnetic_staggered;
        reference.electric_y_V_m[i] = electric;
        reference.magnetic_z_T[i] = electric / SWA_C;
    }

    initial_energy = sw_fdtd1d_energy_j_m2(&grid);
    for (step = 0U; step < steps; ++step) {
        sw_fdtd1d_step(&grid);
    }
    final_energy = sw_fdtd1d_energy_j_m2(&grid);
    if (!swa_maxwell1d_advance(
            &reference, NULL, NULL, final_time_s,
            4.0e-12, 1.0e-20, 1.0e-30, &reference_result) ||
        !reference_result.passes) {
        sw_fdtd1d_destroy(&grid);
        swa_maxwell1d_destroy(&reference);
        return 0;
    }

    for (i = 0U; i < cells; ++i) {
        const double x_m = ((double)i + 0.5) * dx_m;
        const double exact = right_wave(
            x_m, final_time_s, length_m, mode
        );
        const double fdtd_error = grid.electric_y_v_m[i] - exact;
        const double spectral_error =
            reference.electric_y_V_m[i] - exact;
        error_sum += fdtd_error * fdtd_error;
        spectral_error_sum += spectral_error * spectral_error;
        exact_sum += exact * exact;
    }

    *relative_l2_error = sqrt(error_sum / fmax(exact_sum, DBL_MIN));
    *spectral_relative_l2_error =
        sqrt(spectral_error_sum / fmax(exact_sum, DBL_MIN));
    *relative_energy_change =
        (final_energy - initial_energy) / fmax(fabs(initial_energy), DBL_MIN);
    *courant = SWA_C * dt_s / dx_m;

    sw_fdtd1d_destroy(&grid);
    swa_maxwell1d_destroy(&reference);
    return isfinite(*relative_l2_error) &&
           isfinite(*spectral_relative_l2_error) &&
           isfinite(*relative_energy_change) && isfinite(*courant);
}

static double observed_order(double coarse_error, double fine_error,
                             size_t coarse_cells, size_t fine_cells) {
    return log(coarse_error / fine_error) /
           log((double)fine_cells / (double)coarse_cells);
}

static void test_wrong_staggering_is_detected(void) {
    enum { CELLS = 33, ADVERSARY_STEPS = 7 };
    const double length_m = 1.0;
    const double dx_m = length_m / (double)CELLS;
    const double dt_s = 0.5 * dx_m / SWA_C;
    sw_fdtd1d correct;
    sw_fdtd1d wrong;
    double correct_error = 0.0;
    double wrong_error = 0.0;
    double exact_sum = 0.0;
    size_t i;
    size_t step;

    memset(&correct, 0, sizeof(correct));
    memset(&wrong, 0, sizeof(wrong));
    check_true(sw_fdtd1d_init(&correct, CELLS, length_m, dt_s),
               "correct-staggering grid initializes");
    check_true(sw_fdtd1d_init(&wrong, CELLS, length_m, dt_s),
               "wrong-staggering grid initializes");
    for (i = 0U; i < CELLS; ++i) {
        const double electric_x_m = ((double)i + 0.5) * dx_m;
        const double magnetic_x_m = ((double)i + 1.0) * dx_m;
        const double electric = right_wave(electric_x_m, 0.0, length_m, 3U);
        correct.electric_y_v_m[i] = electric;
        correct.magnetic_z_t[i] = right_wave(
            magnetic_x_m, -0.5 * dt_s, length_m, 3U
        ) / SWA_C;
        wrong.electric_y_v_m[i] = electric;
        wrong.magnetic_z_t[i] = -correct.magnetic_z_t[i];
    }
    for (step = 0U; step < ADVERSARY_STEPS; ++step) {
        sw_fdtd1d_step(&correct);
        sw_fdtd1d_step(&wrong);
    }
    for (i = 0U; i < CELLS; ++i) {
        const double x_m = ((double)i + 0.5) * dx_m;
        const double exact = right_wave(
            x_m, (double)ADVERSARY_STEPS * dt_s, length_m, 3U
        );
        const double dc = correct.electric_y_v_m[i] - exact;
        const double dw = wrong.electric_y_v_m[i] - exact;
        correct_error += dc * dc;
        wrong_error += dw * dw;
        exact_sum += exact * exact;
    }
    correct_error = sqrt(correct_error / exact_sum);
    wrong_error = sqrt(wrong_error / exact_sum);
    check_true(wrong_error > 20.0 * correct_error,
               "wrong magnetic sign is rejected by spectral comparison");
    sw_fdtd1d_destroy(&correct);
    sw_fdtd1d_destroy(&wrong);
}

int main(void) {
    const size_t levels[3] = {33U, 65U, 129U};
    double error[3];
    double spectral_error[3];
    double energy_change[3];
    double courant[3];
    double order_coarse_medium;
    double order_medium_fine;
    FILE *fp;
    size_t level;

    for (level = 0U; level < 3U; ++level) {
        check_true(run_level(
                       levels[level], 3U,
                       &error[level],
                       &spectral_error[level],
                       &energy_change[level],
                       &courant[level]),
                   "FDTD and spectral comparison level executes");
        check_true(fabs(courant[level] - 0.5) < 1.0e-15,
                   "declared Courant number is preserved");
        check_true(spectral_error[level] < 5.0e-12,
                   "spectral oracle matches represented analytic wave");
    }
    check_true(error[2] < error[1] && error[1] < error[0],
               "Yee electric-field error decreases monotonically");
    order_coarse_medium = observed_order(
        error[0], error[1], levels[0], levels[1]
    );
    order_medium_fine = observed_order(
        error[1], error[2], levels[1], levels[2]
    );
    check_true(order_coarse_medium > 1.75 && order_coarse_medium < 2.25,
               "coarse-to-medium Yee convergence is second order");
    check_true(order_medium_fine > 1.75 && order_medium_fine < 2.25,
               "medium-to-fine Yee convergence is second order");
    check_true(fabs(energy_change[2]) < fabs(energy_change[0]),
               "Yee staggered energy diagnostic improves with refinement");
    test_wrong_staggering_is_detected();

    fp = fopen("output/fdtd_spectral_integration_receipt.json", "w");
    check_true(fp != NULL, "open FDTD spectral receipt");
    if (fp != NULL) {
        check_true(fprintf(
            fp,
            "{\n"
            "  \"schema\": \"spacewind.fdtd-spectral-integration/v1\",\n"
            "  \"cells\": [%zu, %zu, %zu],\n"
            "  \"courant\": [%.17g, %.17g, %.17g],\n"
            "  \"relative_l2_error\": [%.17g, %.17g, %.17g],\n"
            "  \"spectral_relative_l2_error\": [%.17g, %.17g, %.17g],\n"
            "  \"relative_energy_change\": [%.17g, %.17g, %.17g],\n"
            "  \"observed_order\": [%.17g, %.17g],\n"
            "  \"passes\": %s,\n"
            "  \"nonclaim\": \"this bounded cross-method test verifies the one-polarization periodic Yee update against an independent exact spectral Maxwell oracle; it does not establish electromagnetic PIC conservation or plasma-wing thrust\"\n"
            "}\n",
            levels[0], levels[1], levels[2],
            courant[0], courant[1], courant[2],
            error[0], error[1], error[2],
            spectral_error[0], spectral_error[1], spectral_error[2],
            energy_change[0], energy_change[1], energy_change[2],
            order_coarse_medium, order_medium_fine,
            failures == 0U ? "true" : "false"
        ) > 0, "write FDTD spectral receipt");
        (void)fclose(fp);
    }

    printf(
        "spacewind FDTD/spectral integration: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
