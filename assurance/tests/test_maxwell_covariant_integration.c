#include "spacewind/covariant_em.h"
#include "spacewind/maxwell1d_reference.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);
static double worst_energy_relative = 0.0;
static double worst_momentum_relative = 0.0;
static double worst_characteristic_relative = 0.0;
static double worst_vacuum_ledger_relative = 0.0;

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

static int integrate_stress_energy(
    const swa_maxwell1d_state *state,
    swa_four *out_field_four_momentum_Ns,
    double *out_trace_integral_J,
    double *out_maximum_dominant_energy_ratio
) {
    const double dx = state != NULL ?
        state->length_m / (double)state->samples : 0.0;
    const double volume = state != NULL ?
        state->cross_section_area_m2 * dx : 0.0;
    swa_kahan p0 = {0.0, 0.0};
    swa_kahan px = {0.0, 0.0};
    swa_kahan py = {0.0, 0.0};
    swa_kahan pz = {0.0, 0.0};
    swa_kahan trace = {0.0, 0.0};
    double maximum_dominant = 0.0;
    size_t i;
    if (!swa_maxwell1d_state_is_finite(state) ||
        out_field_four_momentum_Ns == NULL ||
        out_trace_integral_J == NULL ||
        out_maximum_dominant_energy_ratio == NULL ||
        !isfinite(dx) || !(dx > 0.0) ||
        !isfinite(volume) || !(volume > 0.0)) {
        return 0;
    }
    for (i = 0U; i < state->samples; ++i) {
        const swa_em_fields fields = {
            {0.0, state->electric_y_V_m[i],
             state->electric_z_V_m[i]},
            {0.0, state->magnetic_y_T[i],
             state->magnetic_z_T[i]}
        };
        swa_tensor4 tensor;
        swa_covariant_em_diagnostics diagnostics;
        if (!swa_em_stress_energy_tensor(
                &fields, &tensor, &diagnostics) ||
            !diagnostics.dominant_energy_condition) {
            return 0;
        }
        swa_kahan_add(&p0, tensor.component[0][0] * volume / SWA_C);
        swa_kahan_add(&px, tensor.component[0][1] * volume / SWA_C);
        swa_kahan_add(&py, tensor.component[0][2] * volume / SWA_C);
        swa_kahan_add(&pz, tensor.component[0][3] * volume / SWA_C);
        swa_kahan_add(&trace, diagnostics.tensor_trace_J_m3 * volume);
        maximum_dominant = max2(
            maximum_dominant,
            diagnostics.dominant_energy_ratio
        );
    }
    *out_field_four_momentum_Ns =
        swa_four_make(p0.sum, px.sum, py.sum, pz.sum);
    *out_trace_integral_J = trace.sum;
    *out_maximum_dominant_energy_ratio = maximum_dominant;
    return swa_four_is_finite(*out_field_four_momentum_Ns) &&
           isfinite(*out_trace_integral_J) &&
           isfinite(*out_maximum_dominant_energy_ratio);
}

static void fill_arbitrary_fields(swa_maxwell1d_state *state) {
    size_t i;
    for (i = 0U; i < state->samples; ++i) {
        const double x = 2.0 * SWA_PI *
                         (double)i / (double)state->samples;
        state->electric_y_V_m[i] =
            3.0e5 * (0.7 * sin(2.0 * x) +
                     0.2 * cos(5.0 * x));
        state->electric_z_V_m[i] =
            2.0e5 * (-0.4 * cos(3.0 * x) +
                     0.15 * sin(7.0 * x));
        state->magnetic_y_T[i] =
            6.0e-4 * (0.3 * sin(4.0 * x) -
                      0.1 * cos(8.0 * x));
        state->magnetic_z_T[i] =
            5.0e-4 * (-0.5 * cos(2.0 * x) +
                      0.12 * sin(6.0 * x));
    }
}

static void test_global_diagnostic_equivalence(
    swa_four *receipt_momentum
) {
    enum { N = 47 };
    swa_maxwell1d_state state;
    swa_four integrated;
    double trace_integral;
    double maximum_dominant;
    double energy;
    double momentum;
    double right_energy;
    double left_energy;
    double energy_error;
    double momentum_error;
    double characteristic_energy_error;
    double characteristic_momentum_error;
    memset(&state, 0, sizeof(state));
    check_true(swa_maxwell1d_init(&state, N, 19.0, 2.3),
               "global diagnostic state initializes");
    fill_arbitrary_fields(&state);
    check_true(integrate_stress_energy(
                   &state, &integrated, &trace_integral,
                   &maximum_dominant),
               "stress-energy integration executes");
    energy = swa_maxwell1d_field_energy_J(&state);
    momentum = swa_maxwell1d_field_momentum_Ns(&state);
    right_energy = 0.5 * (energy + SWA_C * momentum);
    left_energy = 0.5 * (energy - SWA_C * momentum);
    energy_error = relative_scalar(
        integrated.component[0], energy / SWA_C
    );
    momentum_error = relative_scalar(
        integrated.component[1], momentum
    );
    characteristic_energy_error = relative_scalar(
        energy, right_energy + left_energy
    );
    characteristic_momentum_error = relative_scalar(
        momentum, (right_energy - left_energy) / SWA_C
    );
    worst_energy_relative = max2(worst_energy_relative, energy_error);
    worst_momentum_relative = max2(
        worst_momentum_relative, momentum_error
    );
    worst_characteristic_relative = max2(
        worst_characteristic_relative,
        max2(characteristic_energy_error,
             characteristic_momentum_error)
    );
    check_true(energy_error < 5.0e-14,
               "integrated T00 equals Maxwell reference energy over c");
    check_true(momentum_error < 5.0e-14,
               "integrated T0x equals Maxwell reference field momentum");
    check_true(fabs(integrated.component[2]) <
               1.0e-14 * max2(fabs(integrated.component[0]), 1.0),
               "transverse y field momentum vanishes in one-dimensional geometry");
    check_true(fabs(integrated.component[3]) <
               1.0e-14 * max2(fabs(integrated.component[0]), 1.0),
               "transverse z field momentum vanishes in one-dimensional geometry");
    check_true(characteristic_energy_error < 5.0e-15 &&
               characteristic_momentum_error < 5.0e-15,
               "right-left energy decomposition reconstructs four-momentum");
    check_true(fabs(trace_integral) <
               2.0e-12 * max2(energy, 1.0),
               "integrated electromagnetic stress-energy remains traceless");
    check_true(maximum_dominant <= 1.0 + 1.0e-13,
               "cellwise dominant-energy bound holds");
    *receipt_momentum = integrated;
    swa_maxwell1d_destroy(&state);
}

static void test_vacuum_four_momentum_ledger(void) {
    enum { N = 53, STEPS = 125 };
    swa_maxwell1d_state state;
    swa_maxwell1d_step_result step_result;
    swa_four initial;
    swa_four final;
    swa_four_momentum_ledger ledger;
    swa_four_momentum_ledger_result result;
    double trace_integral;
    double maximum_dominant;
    const double dt = 0.137 * 31.0 / SWA_C;
    size_t step;
    memset(&state, 0, sizeof(state));
    memset(&ledger, 0, sizeof(ledger));
    check_true(swa_maxwell1d_init(&state, N, 31.0, 0.8),
               "vacuum ledger state initializes");
    fill_arbitrary_fields(&state);
    check_true(integrate_stress_energy(
                   &state, &initial, &trace_integral,
                   &maximum_dominant),
               "initial vacuum four-momentum integrates");
    for (step = 0U; step < STEPS; ++step) {
        check_true(swa_maxwell1d_advance(
                       &state, NULL, NULL, dt,
                       5.0e-12, 1.0e-20, 1.0e-29,
                       &step_result) && step_result.passes,
                   "vacuum Maxwell step passes local ledgers");
    }
    check_true(integrate_stress_energy(
                   &state, &final, &trace_integral,
                   &maximum_dominant),
               "final vacuum four-momentum integrates");
    ledger.field_initial_Ns = initial;
    ledger.field_final_Ns = final;
    ledger.matter_initial_Ns = swa_four_make(0.0, 0.0, 0.0, 0.0);
    ledger.matter_final_Ns = swa_four_make(0.0, 0.0, 0.0, 0.0);
    ledger.external_impulse_Ns = swa_four_make(0.0, 0.0, 0.0, 0.0);
    ledger.momentum_in_Ns = swa_four_make(0.0, 0.0, 0.0, 0.0);
    ledger.momentum_out_Ns = swa_four_make(0.0, 0.0, 0.0, 0.0);
    check_true(swa_audit_four_momentum_ledger(
                   &ledger, 8.0e-12, 1.0e-28, &result),
               "vacuum field four-momentum ledger executes");
    worst_vacuum_ledger_relative = max2(
        worst_vacuum_ledger_relative, result.relative_residual
    );
    check_true(result.passes,
               "vacuum spectral propagation conserves field four-momentum");
    ledger.external_impulse_Ns.component[1] = 1.0e-6;
    check_true(swa_audit_four_momentum_ledger(
                   &ledger, 8.0e-12, 1.0e-28, &result) &&
               !result.passes,
               "invented external four-impulse is rejected");
    swa_maxwell1d_destroy(&state);
}

static void test_prescribed_current_time_component(void) {
    enum { N = 39 };
    swa_maxwell1d_state state;
    swa_maxwell1d_step_result result;
    swa_four initial;
    swa_four final;
    double jy[N];
    double jz[N];
    double trace_integral;
    double maximum_dominant;
    double p0_residual;
    size_t i;
    memset(&state, 0, sizeof(state));
    check_true(swa_maxwell1d_init(&state, N, 13.0, 1.4),
               "source-driven stress state initializes");
    fill_arbitrary_fields(&state);
    for (i = 0U; i < N; ++i) {
        const double x = 2.0 * SWA_PI * (double)i / (double)N;
        jy[i] = 2.0e-7 * sin(3.0 * x);
        jz[i] = -1.0e-7 * cos(5.0 * x);
    }
    check_true(integrate_stress_energy(
                   &state, &initial, &trace_integral,
                   &maximum_dominant),
               "initial source-driven four-momentum integrates");
    check_true(swa_maxwell1d_advance(
                   &state, jy, jz, 1.0e-10,
                   8.0e-12, 1.0e-20, 1.0e-29, &result) &&
               result.energy_closes,
               "source-driven Maxwell step closes energy ledger");
    check_true(integrate_stress_energy(
                   &state, &final, &trace_integral,
                   &maximum_dominant),
               "final source-driven four-momentum integrates");
    p0_residual = final.component[0] - initial.component[0] +
                  result.current_work_J / SWA_C;
    check_true(fabs(p0_residual) <
               2.0e-11 * max3(fabs(initial.component[0]),
                               fabs(final.component[0]), DBL_MIN),
               "stress-energy time component closes prescribed-current work");
    check_true(fabs(p0_residual + 1.0e-3 / SWA_C) >
               1.0e-12 / SWA_C,
               "untracked source energy would break four-momentum time ledger");
    swa_maxwell1d_destroy(&state);
}

int main(void) {
    swa_four receipt_momentum;
    FILE *fp;
    memset(&receipt_momentum, 0, sizeof(receipt_momentum));
    test_global_diagnostic_equivalence(&receipt_momentum);
    test_vacuum_four_momentum_ledger();
    test_prescribed_current_time_component();
    fp = fopen("output/maxwell_covariant_integration_receipt.json", "w");
    check_true(fp != NULL, "open Maxwell covariant receipt");
    if (fp != NULL) {
        check_true(fprintf(
            fp,
            "{\n"
            "  \"schema\": \"spacewind.maxwell-covariant-integration/v1\",\n"
            "  \"checks\": %zu,\n"
            "  \"failures\": %zu,\n"
            "  \"deterministic_hash\": \"%016llx\",\n"
            "  \"receipt_field_four_momentum_Ns\": [%.17g, %.17g, %.17g, %.17g],\n"
            "  \"worst_energy_relative\": %.17g,\n"
            "  \"worst_momentum_relative\": %.17g,\n"
            "  \"worst_characteristic_relative\": %.17g,\n"
            "  \"worst_vacuum_ledger_relative\": %.17g,\n"
            "  \"passes\": %s,\n"
            "  \"nonclaim\": \"this bounded integration proves consistency between two implemented electromagnetic diagnostics; it is not a self-consistent particle-plus-field propulsion result\"\n"
            "}\n",
            checks,
            failures,
            (unsigned long long)transcript_hash,
            receipt_momentum.component[0],
            receipt_momentum.component[1],
            receipt_momentum.component[2],
            receipt_momentum.component[3],
            worst_energy_relative,
            worst_momentum_relative,
            worst_characteristic_relative,
            worst_vacuum_ledger_relative,
            failures == 0U ? "true" : "false"
        ) > 0, "write Maxwell covariant receipt");
        (void)fclose(fp);
    }
    printf(
        "spacewind Maxwell/covariant integration: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
