#include "spacewind/maxwell1d_reference.h"

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

static double max2(double a, double b) { return a > b ? a : b; }

static double manufactured(size_t i, size_t n) {
    const double x = 2.0 * SWA_PI * (double)i / (double)n;
    return 0.7 * cos(3.0 * x) +
           0.2 * sin(5.0 * x) -
           0.11 * cos(9.0 * x);
}

static double manufactured_two(size_t i, size_t n) {
    const double x = 2.0 * SWA_PI * (double)i / (double)n;
    return -0.31 * sin(2.0 * x) +
            0.17 * cos(7.0 * x);
}

static double maximum_state_error(
    const swa_maxwell1d_state *a,
    const swa_maxwell1d_state *b
) {
    double maximum = 0.0;
    size_t i;
    for (i = 0U; i < a->samples; ++i) {
        maximum = max2(maximum,
                       fabs(a->electric_y_V_m[i] - b->electric_y_V_m[i]));
        maximum = max2(maximum,
                       fabs(a->electric_z_V_m[i] - b->electric_z_V_m[i]));
        maximum = max2(maximum,
                       SWA_C * fabs(a->magnetic_y_T[i] - b->magnetic_y_T[i]));
        maximum = max2(maximum,
                       SWA_C * fabs(a->magnetic_z_T[i] - b->magnetic_z_T[i]));
    }
    return maximum;
}

static void test_invalid_contracts(void) {
    swa_maxwell1d_state state;
    swa_maxwell1d_step_result result;
    double nan_current[5] = {0.0, 0.0, NAN, 0.0, 0.0};
    memset(&state, 0, sizeof(state));
    check_true(!swa_maxwell1d_init(&state, 8U, 1.0, 1.0),
               "even collocated grid rejected to avoid Nyquist ambiguity");
    check_true(!swa_maxwell1d_init(&state, 3U, 1.0, 1.0),
               "undersized spectral grid rejected");
    check_true(!swa_maxwell1d_init(&state, 5U, -1.0, 1.0),
               "negative domain length rejected");
    check_true(swa_maxwell1d_init(&state, 5U, 1.0, 1.0),
               "small valid odd grid initializes");
    check_true(!swa_maxwell1d_advance(
                   &state, nan_current, NULL, 1.0e-12,
                   1.0e-12, 1.0e-20, 1.0e-30, &result),
               "nonfinite source current rejected");
    check_true(!swa_maxwell1d_advance(
                   &state, state.electric_y_V_m, NULL, 1.0e-12,
                   1.0e-12, 1.0e-20, 1.0e-30, &result),
               "source aliasing mutable field rejected");
    check_true(!swa_maxwell1d_advance(
                   &state, NULL, NULL, NAN,
                   1.0e-12, 1.0e-20, 1.0e-30, &result),
               "nonfinite timestep rejected");
    swa_maxwell1d_destroy(&state);
}

static void test_exact_traveling_modes(
    swa_maxwell1d_step_result *receipt_result,
    swa_maxwell1d_state *receipt_state
) {
    enum { N = 63 };
    const double length = 63.0;
    const double dx = length / (double)N;
    const double dt = dx / SWA_C;
    swa_maxwell1d_state right;
    swa_maxwell1d_state left;
    swa_maxwell1d_step_result right_result;
    swa_maxwell1d_step_result left_result;
    double maximum_right_error = 0.0;
    double maximum_left_error = 0.0;
    size_t i;
    check_true(swa_maxwell1d_init(&right, N, length, 2.0),
               "right-wave state initializes");
    check_true(swa_maxwell1d_init(&left, N, length, 2.0),
               "left-wave state initializes");
    for (i = 0U; i < N; ++i) {
        const double f = manufactured(i, N);
        const double g = manufactured_two(i, N);
        right.electric_y_V_m[i] = f;
        right.magnetic_z_T[i] = f / SWA_C;
        right.electric_z_V_m[i] = g;
        right.magnetic_y_T[i] = -g / SWA_C;
        left.electric_y_V_m[i] = f;
        left.magnetic_z_T[i] = -f / SWA_C;
        left.electric_z_V_m[i] = g;
        left.magnetic_y_T[i] = g / SWA_C;
    }
    check_true(swa_maxwell1d_advance(
                   &right, NULL, NULL, dt,
                   2.0e-12, 1.0e-22, 1.0e-30, &right_result),
               "right-wave exact spectral step executes");
    check_true(right_result.passes,
               "right-wave energy and momentum ledgers close");
    check_true(swa_maxwell1d_advance(
                   &left, NULL, NULL, dt,
                   2.0e-12, 1.0e-22, 1.0e-30, &left_result),
               "left-wave exact spectral step executes");
    check_true(left_result.passes,
               "left-wave energy and momentum ledgers close");
    for (i = 0U; i < N; ++i) {
        const size_t right_source = (i + N - 1U) % N;
        const size_t left_source = (i + 1U) % N;
        maximum_right_error = max2(
            maximum_right_error,
            fabs(right.electric_y_V_m[i] -
                 manufactured(right_source, N))
        );
        maximum_right_error = max2(
            maximum_right_error,
            fabs(right.electric_z_V_m[i] -
                 manufactured_two(right_source, N))
        );
        maximum_left_error = max2(
            maximum_left_error,
            fabs(left.electric_y_V_m[i] -
                 manufactured(left_source, N))
        );
        maximum_left_error = max2(
            maximum_left_error,
            fabs(left.electric_z_V_m[i] -
                 manufactured_two(left_source, N))
        );
    }
    check_true(maximum_right_error < 2.0e-12,
               "right-moving finite modes translate one exact cell");
    check_true(maximum_left_error < 2.0e-12,
               "left-moving finite modes translate one exact cell");
    check_true(fabs(right_result.field_momentum_initial_Ns -
                    right_result.field_energy_initial_J / SWA_C) < 1.0e-26,
               "pure right wave has momentum plus energy over c");
    check_true(fabs(left_result.field_momentum_initial_Ns +
                    left_result.field_energy_initial_J / SWA_C) < 1.0e-26,
               "pure left wave has momentum minus energy over c");
    *receipt_result = right_result;
    check_true(swa_maxwell1d_copy(&right, receipt_state),
               "receipt state copied");
    swa_maxwell1d_destroy(&right);
    swa_maxwell1d_destroy(&left);
}

static void test_time_reversal_and_replay(void) {
    enum { N = 41 };
    swa_maxwell1d_state initial;
    swa_maxwell1d_state evolved;
    swa_maxwell1d_state replay;
    swa_maxwell1d_step_result result_a;
    swa_maxwell1d_step_result result_b;
    const double dt = 0.173 * 41.0 / SWA_C;
    size_t i;
    check_true(swa_maxwell1d_init(&initial, N, 41.0, 1.0),
               "reversal initial state initializes");
    for (i = 0U; i < N; ++i) {
        initial.electric_y_V_m[i] = manufactured(i, N);
        initial.electric_z_V_m[i] = manufactured_two(i, N);
        initial.magnetic_y_T[i] =
            0.37 * manufactured((i + 3U) % N, N) / SWA_C;
        initial.magnetic_z_T[i] =
            -0.22 * manufactured_two((i + 5U) % N, N) / SWA_C;
    }
    check_true(swa_maxwell1d_copy(&initial, &evolved),
               "reversal evolved state copied");
    check_true(swa_maxwell1d_copy(&initial, &replay),
               "deterministic replay state copied");
    check_true(swa_maxwell1d_advance(
                   &evolved, NULL, NULL, dt,
                   3.0e-12, 1.0e-22, 1.0e-30, &result_a) &&
               result_a.passes,
               "forward arbitrary vacuum state passes");
    check_true(swa_maxwell1d_advance(
                   &replay, NULL, NULL, dt,
                   3.0e-12, 1.0e-22, 1.0e-30, &result_b) &&
               result_b.passes,
               "deterministic replay executes");
    check_true(maximum_state_error(&evolved, &replay) == 0.0,
               "vacuum spectral replay is bitwise deterministic");
    check_true(swa_maxwell1d_advance(
                   &evolved, NULL, NULL, -dt,
                   3.0e-12, 1.0e-22, 1.0e-30, &result_a) &&
               result_a.passes,
               "negative-time inverse step passes");
    check_true(maximum_state_error(&evolved, &initial) < 5.0e-12,
               "vacuum spectral propagation reverses to initial state");
    swa_maxwell1d_destroy(&initial);
    swa_maxwell1d_destroy(&evolved);
    swa_maxwell1d_destroy(&replay);
}

static void test_polarization_covariance(void) {
    enum { N = 37 };
    const double angle = 0.731;
    const double c = cos(angle);
    const double s = sin(angle);
    swa_maxwell1d_state base;
    swa_maxwell1d_state rotated;
    swa_maxwell1d_step_result base_result;
    swa_maxwell1d_step_result rotated_result;
    double maximum_error = 0.0;
    size_t i;
    check_true(swa_maxwell1d_init(&base, N, 23.0, 0.7),
               "polarization base initializes");
    check_true(swa_maxwell1d_init(&rotated, N, 23.0, 0.7),
               "polarization rotated initializes");
    for (i = 0U; i < N; ++i) {
        const double ey = manufactured(i, N);
        const double ez = manufactured_two(i, N);
        const double by = -ez / SWA_C;
        const double bz = ey / SWA_C;
        base.electric_y_V_m[i] = ey;
        base.electric_z_V_m[i] = ez;
        base.magnetic_y_T[i] = by;
        base.magnetic_z_T[i] = bz;
        rotated.electric_y_V_m[i] = c * ey - s * ez;
        rotated.electric_z_V_m[i] = s * ey + c * ez;
        rotated.magnetic_y_T[i] = c * by - s * bz;
        rotated.magnetic_z_T[i] = s * by + c * bz;
    }
    check_true(swa_maxwell1d_advance(
                   &base, NULL, NULL, 0.41 * 23.0 / SWA_C,
                   3.0e-12, 1.0e-22, 1.0e-30, &base_result) &&
               base_result.passes,
               "base polarization propagates");
    check_true(swa_maxwell1d_advance(
                   &rotated, NULL, NULL, 0.41 * 23.0 / SWA_C,
                   3.0e-12, 1.0e-22, 1.0e-30, &rotated_result) &&
               rotated_result.passes,
               "rotated polarization propagates");
    for (i = 0U; i < N; ++i) {
        maximum_error = max2(
            maximum_error,
            fabs(rotated.electric_y_V_m[i] -
                 (c * base.electric_y_V_m[i] -
                  s * base.electric_z_V_m[i]))
        );
        maximum_error = max2(
            maximum_error,
            fabs(rotated.electric_z_V_m[i] -
                 (s * base.electric_y_V_m[i] +
                  c * base.electric_z_V_m[i]))
        );
        maximum_error = max2(
            maximum_error,
            SWA_C * fabs(rotated.magnetic_y_T[i] -
                 (c * base.magnetic_y_T[i] -
                  s * base.magnetic_z_T[i]))
        );
        maximum_error = max2(
            maximum_error,
            SWA_C * fabs(rotated.magnetic_z_T[i] -
                 (s * base.magnetic_y_T[i] +
                  c * base.magnetic_z_T[i]))
        );
    }
    check_true(maximum_error < 5.0e-12,
               "transverse polarization rotation commutes with propagation");
    swa_maxwell1d_destroy(&base);
    swa_maxwell1d_destroy(&rotated);
}

static void test_prescribed_current_energy_ledger(void) {
    enum { N = 35 };
    swa_maxwell1d_state state;
    swa_maxwell1d_step_result result;
    double jy[N];
    double jz[N];
    double corrupted_residual;
    double corrupted_relative;
    int corrupted_passes;
    size_t i;
    check_true(swa_maxwell1d_init(&state, N, 17.0, 1.3),
               "current-driven state initializes");
    for (i = 0U; i < N; ++i) {
        const double x = 2.0 * SWA_PI * (double)i / (double)N;
        state.electric_y_V_m[i] = 2.0e3 * manufactured(i, N);
        state.electric_z_V_m[i] = 1.5e3 * manufactured_two(i, N);
        state.magnetic_y_T[i] = 2.0e-6 * cos(4.0 * x);
        state.magnetic_z_T[i] = -1.0e-6 * sin(6.0 * x);
        jy[i] = 2.0e-7 * sin(3.0 * x);
        jz[i] = -1.5e-7 * cos(5.0 * x);
    }
    check_true(swa_maxwell1d_advance(
                   &state, jy, jz, 2.0e-10,
                   8.0e-12, 1.0e-20, 1.0e-30, &result),
               "prescribed-current Strang step executes");
    check_true(result.energy_closes && result.passes,
               "field energy plus midpoint source work closes");
    check_true(!result.vacuum_step,
               "current-driven result is marked non-vacuum");
    check_true(swa_maxwell1d_audit_energy_balance(
                   result.field_energy_initial_J,
                   result.field_energy_final_J,
                   result.current_work_J + 1.0e-3,
                   8.0e-12, 1.0e-20,
                   &corrupted_residual,
                   &corrupted_relative,
                   &corrupted_passes),
               "corrupted source-work audit executes");
    check_true(!corrupted_passes && fabs(corrupted_residual) > 1.0e-4,
               "untracked source energy is rejected");
    swa_maxwell1d_destroy(&state);
}

int main(void) {
    swa_maxwell1d_step_result receipt_result;
    swa_maxwell1d_state receipt_state;
    FILE *fp;
    memset(&receipt_result, 0, sizeof(receipt_result));
    memset(&receipt_state, 0, sizeof(receipt_state));
    test_invalid_contracts();
    test_exact_traveling_modes(&receipt_result, &receipt_state);
    test_time_reversal_and_replay();
    test_polarization_covariance();
    test_prescribed_current_energy_ledger();
    fp = fopen("output/maxwell1d_reference_receipt.json", "w");
    check_true(fp != NULL, "open Maxwell reference receipt");
    if (fp != NULL) {
        check_true(swa_write_maxwell1d_reference_receipt(
                       fp, &receipt_state, &receipt_result),
                   "write Maxwell reference receipt");
        (void)fclose(fp);
    }
    swa_maxwell1d_destroy(&receipt_state);
    printf(
        "spacewind Maxwell 1D reference: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
