#include "spacewind/metric.h"
#include "spacewind/constants.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

void sw_metric_minkowski_cartesian(const double x[4], double g_cov[4][4], const void *ctx) {
    (void)x;
    (void)ctx;
    memset(g_cov, 0, sizeof(double) * 16U);
    g_cov[0][0] = -1.0;
    g_cov[1][1] = 1.0;
    g_cov[2][2] = 1.0;
    g_cov[3][3] = 1.0;
}

void sw_metric_schwarzschild(const double x[4], double g_cov[4][4], const void *ctx) {
    const sw_schwarzschild_metric *params = (const sw_schwarzschild_metric *)ctx;
    const double r = x[1];
    const double theta = x[2];
    const double rs = 2.0 * SW_G * params->mass_kg / (SW_C * SW_C);
    const double f = 1.0 - rs / r;

    memset(g_cov, 0, sizeof(double) * 16U);
    g_cov[0][0] = -f;
    g_cov[1][1] = 1.0 / f;
    g_cov[2][2] = r * r;
    g_cov[3][3] = r * r * sin(theta) * sin(theta);
}

void sw_metric_kerr_boyer_lindquist(const double x[4], double g_cov[4][4], const void *ctx) {
    const sw_kerr_metric *params = (const sw_kerr_metric *)ctx;
    const double mass_length = SW_G * params->mass_kg / (SW_C * SW_C);
    const double spin = fmax(-1.0, fmin(1.0, params->dimensionless_spin)) * mass_length;
    const double r = x[1];
    const double theta = x[2];
    const double sin_theta = sin(theta);
    const double cos_theta = cos(theta);
    const double sigma = r * r + spin * spin * cos_theta * cos_theta;
    const double delta = r * r - 2.0 * mass_length * r + spin * spin;
    const double sin2 = sin_theta * sin_theta;

    memset(g_cov, 0, sizeof(double) * 16U);
    g_cov[0][0] = -(1.0 - 2.0 * mass_length * r / sigma);
    g_cov[0][3] = -2.0 * mass_length * spin * r * sin2 / sigma;
    g_cov[3][0] = g_cov[0][3];
    g_cov[1][1] = sigma / delta;
    g_cov[2][2] = sigma;
    g_cov[3][3] = (
        r * r + spin * spin
        + 2.0 * mass_length * spin * spin * r * sin2 / sigma
    ) * sin2;
}


bool sw_kerr_equatorial_circular_state(
    const sw_kerr_metric *metric,
    double radius_m,
    int orbital_direction,
    sw_geodesic_state *state
) {
    double x[4];
    double g[4][4];
    double mass_length;
    double spin_length;
    double root_mass;
    double direction;
    double denominator;
    double omega_per_m;
    double normalization_coefficient;
    double u0;

    if (metric == NULL || state == NULL || !(metric->mass_kg > 0.0)
        || !(radius_m > 0.0) || (orbital_direction != 1 && orbital_direction != -1)) {
        return false;
    }
    mass_length = SW_G * metric->mass_kg / (SW_C * SW_C);
    spin_length = fmax(-1.0, fmin(1.0, metric->dimensionless_spin)) * mass_length;
    direction = (double)orbital_direction;
    root_mass = sqrt(mass_length);
    denominator = pow(radius_m, 1.5) + direction * spin_length * root_mass;
    if (!(denominator > 0.0)) {
        return false;
    }
    omega_per_m = direction * root_mass / denominator;
    x[0] = 0.0;
    x[1] = radius_m;
    x[2] = 0.5 * SW_PI;
    x[3] = 0.0;
    sw_metric_kerr_boyer_lindquist(x, g, metric);
    normalization_coefficient = g[0][0]
        + 2.0 * g[0][3] * omega_per_m
        + g[3][3] * omega_per_m * omega_per_m;
    if (!(normalization_coefficient < 0.0)) {
        return false;
    }
    u0 = 1.0 / sqrt(-normalization_coefficient);
    memset(state, 0, sizeof(*state));
    memcpy(state->x, x, sizeof(x));
    state->u[0] = u0;
    state->u[3] = omega_per_m * u0;
    return true;
}

void sw_metric_flrw_de_sitter(const double x[4], double g_cov[4][4], const void *ctx) {
    const sw_flrw_metric *params = (const sw_flrw_metric *)ctx;
    const double time_s = x[0] / SW_C;
    const double scale = params->scale_factor_at_epoch * exp(params->hubble_rate_per_s * time_s);
    (void)x[1];
    (void)x[2];
    (void)x[3];
    memset(g_cov, 0, sizeof(double) * 16U);
    g_cov[0][0] = -1.0;
    g_cov[1][1] = scale * scale;
    g_cov[2][2] = scale * scale;
    g_cov[3][3] = scale * scale;
}

static double sw_metric_step_for(const sw_metric *metric, const double x[4], size_t axis) {
    const double configured = metric->derivative_step[axis];
    if (configured > 0.0) {
        return configured;
    }
    {
        const double scale = fmax(fabs(x[axis]), 1.0);
        return cbrt(DBL_EPSILON) * scale;
    }
}

bool sw_metric_christoffel(
    const sw_metric *metric,
    const double x[4],
    double gamma[4][4][4]
) {
    double g[4][4];
    double g_inv[4][4];
    double dg[4][4][4];
    size_t derivative_axis;
    size_t i;
    size_t j;
    size_t mu;
    size_t alpha;
    size_t beta;
    size_t nu;

    if (metric == NULL || metric->eval == NULL) {
        return false;
    }

    metric->eval(x, g, metric->ctx);
    if (!sw_mat4_inverse(&g[0][0], &g_inv[0][0])) {
        return false;
    }

    memset(dg, 0, sizeof(dg));
    for (derivative_axis = 0U; derivative_axis < 4U; ++derivative_axis) {
        double xp[4];
        double xm[4];
        double gp[4][4];
        double gm[4][4];
        const double h = sw_metric_step_for(metric, x, derivative_axis);
        memcpy(xp, x, sizeof(xp));
        memcpy(xm, x, sizeof(xm));
        xp[derivative_axis] += h;
        xm[derivative_axis] -= h;
        metric->eval(xp, gp, metric->ctx);
        metric->eval(xm, gm, metric->ctx);
        for (i = 0U; i < 4U; ++i) {
            for (j = 0U; j < 4U; ++j) {
                dg[derivative_axis][i][j] = (gp[i][j] - gm[i][j]) / (2.0 * h);
            }
        }
    }

    memset(gamma, 0, sizeof(double) * 64U);
    for (mu = 0U; mu < 4U; ++mu) {
        for (alpha = 0U; alpha < 4U; ++alpha) {
            for (beta = 0U; beta < 4U; ++beta) {
                double value = 0.0;
                for (nu = 0U; nu < 4U; ++nu) {
                    value += 0.5 * g_inv[mu][nu] * (
                        dg[alpha][nu][beta]
                        + dg[beta][nu][alpha]
                        - dg[nu][alpha][beta]
                    );
                }
                gamma[mu][alpha][beta] = value;
            }
        }
    }
    return true;
}

bool sw_geodesic_rhs(
    const sw_metric *metric,
    const sw_geodesic_state *state,
    sw_geodesic_state *derivative
) {
    double gamma[4][4][4];
    size_t mu;
    size_t alpha;
    size_t beta;

    if (!sw_metric_christoffel(metric, state->x, gamma)) {
        return false;
    }

    for (mu = 0U; mu < 4U; ++mu) {
        derivative->x[mu] = state->u[mu];
        derivative->u[mu] = 0.0;
        for (alpha = 0U; alpha < 4U; ++alpha) {
            for (beta = 0U; beta < 4U; ++beta) {
                derivative->u[mu] -= gamma[mu][alpha][beta]
                    * state->u[alpha] * state->u[beta];
            }
        }
    }
    return true;
}

static sw_geodesic_state sw_geodesic_add_scaled(
    const sw_geodesic_state *a,
    const sw_geodesic_state *b,
    double scale
) {
    sw_geodesic_state out;
    size_t i;
    for (i = 0U; i < 4U; ++i) {
        out.x[i] = a->x[i] + scale * b->x[i];
        out.u[i] = a->u[i] + scale * b->u[i];
    }
    return out;
}

bool sw_geodesic_rk4_step(
    const sw_metric *metric,
    sw_geodesic_state *state,
    double affine_step
) {
    sw_geodesic_state k1;
    sw_geodesic_state k2;
    sw_geodesic_state k3;
    sw_geodesic_state k4;
    sw_geodesic_state tmp;
    size_t i;

    if (!sw_geodesic_rhs(metric, state, &k1)) {
        return false;
    }
    tmp = sw_geodesic_add_scaled(state, &k1, 0.5 * affine_step);
    if (!sw_geodesic_rhs(metric, &tmp, &k2)) {
        return false;
    }
    tmp = sw_geodesic_add_scaled(state, &k2, 0.5 * affine_step);
    if (!sw_geodesic_rhs(metric, &tmp, &k3)) {
        return false;
    }
    tmp = sw_geodesic_add_scaled(state, &k3, affine_step);
    if (!sw_geodesic_rhs(metric, &tmp, &k4)) {
        return false;
    }

    for (i = 0U; i < 4U; ++i) {
        state->x[i] += affine_step * (
            k1.x[i] + 2.0 * k2.x[i] + 2.0 * k3.x[i] + k4.x[i]
        ) / 6.0;
        state->u[i] += affine_step * (
            k1.u[i] + 2.0 * k2.u[i] + 2.0 * k3.u[i] + k4.u[i]
        ) / 6.0;
    }
    return true;
}

double sw_metric_interval(
    const sw_metric *metric,
    const double x[4],
    const double u[4]
) {
    double g[4][4];
    double interval = 0.0;
    size_t i;
    size_t j;
    metric->eval(x, g, metric->ctx);
    for (i = 0U; i < 4U; ++i) {
        for (j = 0U; j < 4U; ++j) {
            interval += g[i][j] * u[i] * u[j];
        }
    }
    return interval;
}

double sw_weak_field_proper_time_rate(double gravitational_mu, sw_vec3 position, sw_vec3 velocity) {
    const double r = sw_v3_norm(position);
    const double v2 = sw_v3_norm2(velocity);
    if (r <= 0.0) {
        return 0.0;
    }
    {
        const double inside = 1.0 - 2.0 * gravitational_mu / (r * SW_C * SW_C)
            - v2 / (SW_C * SW_C);
        return sqrt(fmax(inside, 0.0));
    }
}

sw_vec3 sw_schwarzschild_1pn_acceleration(double gravitational_mu, sw_vec3 position, sw_vec3 velocity) {
    const double r = sw_v3_norm(position);
    if (r <= 0.0) {
        return sw_v3(0.0, 0.0, 0.0);
    }
    {
        const double r2 = r * r;
        const double r3 = r2 * r;
        const double v2 = sw_v3_norm2(velocity);
        const double rv = sw_v3_dot(position, velocity);
        const double factor = gravitational_mu / (SW_C * SW_C * r3);
        const sw_vec3 radial = sw_v3_scale(
            position,
            4.0 * gravitational_mu / r - v2
        );
        const sw_vec3 velocity_term = sw_v3_scale(velocity, 4.0 * rv);
        return sw_v3_scale(sw_v3_add(radial, velocity_term), factor);
    }
}
