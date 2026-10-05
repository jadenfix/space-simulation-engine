#include "spacewind/sail_validation.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

bool sw_sail_screen_cold_passive_stream(
    double cd, double cl, sw_sail_stream_screen *out
) {
    double speed_ratio;
    double power_ratio;
    if (out == NULL) {
        return false;
    }
    memset(out, 0, sizeof(*out));
    if (!isfinite(cd) || !isfinite(cl)) {
        return false;
    }
    speed_ratio = hypot(1.0 - cd, cl);
    power_ratio = speed_ratio * speed_ratio;
    if (!isfinite(power_ratio)) {
        return false;
    }
    out->valid = true;
    out->outgoing_kinetic_power_lower_bound_ratio = power_ratio;
    out->additional_power_lower_bound_over_incoming = fmax(0.0, power_ratio - 1.0);
    /* Only roundoff allowance, not an adjustable physical acceptance budget. */
    out->passive_necessary_bound_satisfied = power_ratio <= 1.0 + 64.0 * DBL_EPSILON;
    return true;
}
