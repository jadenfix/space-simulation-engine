#ifndef SPACEWIND_SAIL_VALIDATION_H
#define SPACEWIND_SAIL_VALIDATION_H

#include <stdbool.h>

/* Necessary condition ONLY: steady passive device, cold uniform stream,
 * unchanged stored energy, no extra power or field-boundary momentum, and
 * A_eff equal to the actual intercepted mass-flux area. Coefficients use
 * F = rho*u^2*A_eff*(cd*drag_hat + cl*lift_hat), NOT the half-rho convention.
 * Passing is not evidence of a realizable device. Failure requests a complete
 * power/particle/field ledger; it does not falsify an active plasma sail. */
typedef struct {
    bool valid;
    bool passive_necessary_bound_satisfied;
    double outgoing_kinetic_power_lower_bound_ratio;
    double additional_power_lower_bound_over_incoming;
} sw_sail_stream_screen;

/* Signed cl is supported. Invalid/overflowing inputs return false and clear
 * the output, so a previous successful result cannot be reused accidentally. */
bool sw_sail_screen_cold_passive_stream(
    double cd, double cl, sw_sail_stream_screen *out
);

#endif
