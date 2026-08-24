#include "spacewind/four_momentum_delta.h"

#include <float.h>
#include <math.h>
#include <string.h>

static double max2(double a, double b) { return a > b ? a : b; }
static double max3(double a, double b, double c) {
    return max2(max2(a, b), c);
}

int swa_audit_four_momentum_delta_ledger(
    const swa_four_momentum_delta_ledger *ledger,
    double physical_scale_floor_Ns,
    double relative_tolerance,
    double absolute_tolerance_Ns,
    swa_four_momentum_delta_result *out
) {
    swa_four supplied;
    swa_four expected;
    if (ledger == NULL || out == NULL ||
        !swa_four_is_finite(ledger->field_change_Ns) ||
        !swa_four_is_finite(ledger->matter_change_Ns) ||
        !swa_four_is_finite(ledger->external_impulse_Ns) ||
        !swa_four_is_finite(ledger->momentum_in_Ns) ||
        !swa_four_is_finite(ledger->momentum_out_Ns) ||
        !isfinite(physical_scale_floor_Ns) ||
        physical_scale_floor_Ns < 0.0 ||
        !isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
        !isfinite(absolute_tolerance_Ns) ||
        absolute_tolerance_Ns < 0.0) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    supplied = swa_four_add(
        ledger->field_change_Ns,
        ledger->matter_change_Ns
    );
    expected = swa_four_add(
        ledger->external_impulse_Ns,
        swa_four_sub(
            ledger->momentum_in_Ns,
            ledger->momentum_out_Ns
        )
    );
    out->residual_Ns = swa_four_sub(supplied, expected);
    out->residual_norm_Ns = swa_four_euclidean_norm(out->residual_Ns);
    out->physical_scale_Ns = max3(
        swa_four_euclidean_norm(ledger->field_change_Ns),
        swa_four_euclidean_norm(ledger->matter_change_Ns),
        max3(
            swa_four_euclidean_norm(expected),
            physical_scale_floor_Ns,
            DBL_MIN
        )
    );
    out->relative_residual =
        out->residual_norm_Ns / out->physical_scale_Ns;
    out->finite = swa_four_is_finite(out->residual_Ns) &&
        isfinite(out->residual_norm_Ns) &&
        isfinite(out->physical_scale_Ns) &&
        isfinite(out->relative_residual);
    out->passes = out->finite &&
        (out->residual_norm_Ns <= absolute_tolerance_Ns ||
         out->relative_residual <= relative_tolerance);
    return 1;
}

int swa_lorentz_transform_four_momentum_delta_ledger(
    swa_vec3 frame_velocity_mps,
    const swa_four_momentum_delta_ledger *input,
    swa_four_momentum_delta_ledger *output
) {
    if (input == NULL || output == NULL) {
        return 0;
    }
    return swa_lorentz_transform_four(
               frame_velocity_mps, input->field_change_Ns,
               &output->field_change_Ns) &&
           swa_lorentz_transform_four(
               frame_velocity_mps, input->matter_change_Ns,
               &output->matter_change_Ns) &&
           swa_lorentz_transform_four(
               frame_velocity_mps, input->external_impulse_Ns,
               &output->external_impulse_Ns) &&
           swa_lorentz_transform_four(
               frame_velocity_mps, input->momentum_in_Ns,
               &output->momentum_in_Ns) &&
           swa_lorentz_transform_four(
               frame_velocity_mps, input->momentum_out_Ns,
               &output->momentum_out_Ns);
}
