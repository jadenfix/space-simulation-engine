#ifndef SPACEWIND_METRIC_H
#define SPACEWIND_METRIC_H

#include <stdbool.h>

#include "spacewind/math3.h"

typedef void (*sw_metric_eval_fn)(const double x[4], double g_cov[4][4], const void *ctx);

typedef struct {
    sw_metric_eval_fn eval;
    const void *ctx;
    double derivative_step[4];
} sw_metric;

typedef struct {
    double mass_kg;
} sw_schwarzschild_metric;

typedef struct {
    double mass_kg;
    double dimensionless_spin;
} sw_kerr_metric;

typedef struct {
    double scale_factor_at_epoch;
    double hubble_rate_per_s;
} sw_flrw_metric;

typedef struct {
    double x[4];
    double u[4];
} sw_geodesic_state;

void sw_metric_minkowski_cartesian(const double x[4], double g_cov[4][4], const void *ctx);
void sw_metric_schwarzschild(const double x[4], double g_cov[4][4], const void *ctx);
void sw_metric_kerr_boyer_lindquist(const double x[4], double g_cov[4][4], const void *ctx);
void sw_metric_flrw_de_sitter(const double x[4], double g_cov[4][4], const void *ctx);
bool sw_kerr_equatorial_circular_state(
    const sw_kerr_metric *metric,
    double radius_m,
    int orbital_direction,
    sw_geodesic_state *state
);

bool sw_metric_christoffel(
    const sw_metric *metric,
    const double x[4],
    double gamma[4][4][4]
);

bool sw_geodesic_rhs(
    const sw_metric *metric,
    const sw_geodesic_state *state,
    sw_geodesic_state *derivative
);

bool sw_geodesic_rk4_step(
    const sw_metric *metric,
    sw_geodesic_state *state,
    double affine_step
);

double sw_metric_interval(
    const sw_metric *metric,
    const double x[4],
    const double u[4]
);

double sw_weak_field_proper_time_rate(double gravitational_mu, sw_vec3 position, sw_vec3 velocity);
sw_vec3 sw_schwarzschild_1pn_acceleration(double gravitational_mu, sw_vec3 position, sw_vec3 velocity);

#endif
