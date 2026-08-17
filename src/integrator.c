#include "spacewind/integrator.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

bool sw_ode_workspace_init(sw_ode_workspace *workspace, size_t dimension, size_t vectors) {
    if (workspace == NULL || dimension == 0U || vectors == 0U
        || dimension > SIZE_MAX / vectors
        || dimension * vectors > SIZE_MAX / sizeof(double)) {
        return false;
    }
    workspace->dimension = dimension;
    workspace->vectors = vectors;
    workspace->storage = (double *)calloc(dimension * vectors, sizeof(double));
    return workspace->storage != NULL;
}

void sw_ode_workspace_destroy(sw_ode_workspace *workspace) {
    if (workspace != NULL) {
        free(workspace->storage);
        workspace->storage = NULL;
        workspace->dimension = 0U;
        workspace->vectors = 0U;
    }
}

static bool sw_workspace_has(const sw_ode_workspace *workspace, size_t dimension, size_t vectors) {
    return workspace != NULL
        && workspace->storage != NULL
        && workspace->dimension >= dimension
        && workspace->vectors >= vectors;
}

bool sw_rk4_step(
    sw_ode_rhs_fn rhs,
    void *context,
    double time,
    double step,
    double *state,
    size_t dimension,
    sw_ode_workspace *workspace
) {
    double *k1;
    double *k2;
    double *k3;
    double *k4;
    double *tmp;
    size_t i;

    if (rhs == NULL || state == NULL || !sw_workspace_has(workspace, dimension, 5U)) {
        return false;
    }
    k1 = workspace->storage;
    k2 = k1 + dimension;
    k3 = k2 + dimension;
    k4 = k3 + dimension;
    tmp = k4 + dimension;

    if (!rhs(time, state, k1, dimension, context)) {
        return false;
    }
    for (i = 0U; i < dimension; ++i) {
        tmp[i] = state[i] + 0.5 * step * k1[i];
    }
    if (!rhs(time + 0.5 * step, tmp, k2, dimension, context)) {
        return false;
    }
    for (i = 0U; i < dimension; ++i) {
        tmp[i] = state[i] + 0.5 * step * k2[i];
    }
    if (!rhs(time + 0.5 * step, tmp, k3, dimension, context)) {
        return false;
    }
    for (i = 0U; i < dimension; ++i) {
        tmp[i] = state[i] + step * k3[i];
    }
    if (!rhs(time + step, tmp, k4, dimension, context)) {
        return false;
    }
    for (i = 0U; i < dimension; ++i) {
        state[i] += step * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    return true;
}

bool sw_dopri54_step(
    sw_ode_rhs_fn rhs,
    void *context,
    double time,
    double requested_step,
    double absolute_tolerance,
    double relative_tolerance,
    double *state,
    size_t dimension,
    sw_ode_workspace *workspace,
    double *accepted_step,
    double *suggested_next_step,
    unsigned *rejected_attempts
) {
    double *k1;
    double *k2;
    double *k3;
    double *k4;
    double *k5;
    double *k6;
    double *k7;
    double *tmp;
    double *y5;
    double step = requested_step;
    unsigned attempt;
    size_t i;

    if (rhs == NULL || state == NULL || accepted_step == NULL || suggested_next_step == NULL
        || !(requested_step > 0.0) || !(absolute_tolerance > 0.0)
        || !(relative_tolerance >= 0.0)
        || !sw_workspace_has(workspace, dimension, 9U)) {
        return false;
    }
    if (rejected_attempts != NULL) {
        *rejected_attempts = 0U;
    }

    k1 = workspace->storage;
    k2 = k1 + dimension;
    k3 = k2 + dimension;
    k4 = k3 + dimension;
    k5 = k4 + dimension;
    k6 = k5 + dimension;
    k7 = k6 + dimension;
    tmp = k7 + dimension;
    y5 = tmp + dimension;

    for (attempt = 0U; attempt < 16U; ++attempt) {
        double err2 = 0.0;
        if (!rhs(time, state, k1, dimension, context)) {
            return false;
        }
        for (i = 0U; i < dimension; ++i) {
            tmp[i] = state[i] + step * (1.0 / 5.0) * k1[i];
        }
        if (!rhs(time + step * (1.0 / 5.0), tmp, k2, dimension, context)) {
            return false;
        }
        for (i = 0U; i < dimension; ++i) {
            tmp[i] = state[i] + step * ((3.0 / 40.0) * k1[i] + (9.0 / 40.0) * k2[i]);
        }
        if (!rhs(time + step * (3.0 / 10.0), tmp, k3, dimension, context)) {
            return false;
        }
        for (i = 0U; i < dimension; ++i) {
            tmp[i] = state[i] + step * (
                (44.0 / 45.0) * k1[i]
                + (-56.0 / 15.0) * k2[i]
                + (32.0 / 9.0) * k3[i]
            );
        }
        if (!rhs(time + step * (4.0 / 5.0), tmp, k4, dimension, context)) {
            return false;
        }
        for (i = 0U; i < dimension; ++i) {
            tmp[i] = state[i] + step * (
                (19372.0 / 6561.0) * k1[i]
                + (-25360.0 / 2187.0) * k2[i]
                + (64448.0 / 6561.0) * k3[i]
                + (-212.0 / 729.0) * k4[i]
            );
        }
        if (!rhs(time + step * (8.0 / 9.0), tmp, k5, dimension, context)) {
            return false;
        }
        for (i = 0U; i < dimension; ++i) {
            tmp[i] = state[i] + step * (
                (9017.0 / 3168.0) * k1[i]
                + (-355.0 / 33.0) * k2[i]
                + (46732.0 / 5247.0) * k3[i]
                + (49.0 / 176.0) * k4[i]
                + (-5103.0 / 18656.0) * k5[i]
            );
        }
        if (!rhs(time + step, tmp, k6, dimension, context)) {
            return false;
        }
        for (i = 0U; i < dimension; ++i) {
            y5[i] = state[i] + step * (
                (35.0 / 384.0) * k1[i]
                + (500.0 / 1113.0) * k3[i]
                + (125.0 / 192.0) * k4[i]
                + (-2187.0 / 6784.0) * k5[i]
                + (11.0 / 84.0) * k6[i]
            );
        }
        if (!rhs(time + step, y5, k7, dimension, context)) {
            return false;
        }

        for (i = 0U; i < dimension; ++i) {
            const double y4 = state[i] + step * (
                (5179.0 / 57600.0) * k1[i]
                + (7571.0 / 16695.0) * k3[i]
                + (393.0 / 640.0) * k4[i]
                + (-92097.0 / 339200.0) * k5[i]
                + (187.0 / 2100.0) * k6[i]
                + (1.0 / 40.0) * k7[i]
            );
            const double scale = absolute_tolerance
                + relative_tolerance * fmax(fabs(state[i]), fabs(y5[i]));
            const double e = (y5[i] - y4) / fmax(scale, DBL_MIN);
            err2 += e * e;
        }

        {
            const double error_norm = sqrt(err2 / (double)dimension);
            const double safe_error = fmax(error_norm, 1.0e-16);
            const double factor = 0.9 * pow(safe_error, -0.2);
            const double bounded_factor = fmax(0.1, fmin(5.0, factor));
            if (error_norm <= 1.0) {
                memcpy(state, y5, dimension * sizeof(double));
                *accepted_step = step;
                *suggested_next_step = step * bounded_factor;
                return true;
            }
            if (rejected_attempts != NULL) {
                *rejected_attempts += 1U;
            }
            step *= fmax(0.1, fmin(0.5, bounded_factor));
        }
    }
    return false;
}

void sw_velocity_verlet_step(
    sw_acceleration_fn acceleration,
    void *context,
    double time,
    double step,
    sw_vec3 *position,
    sw_vec3 *velocity
) {
    const sw_vec3 a0 = acceleration(time, *position, context);
    const sw_vec3 half_velocity = sw_v3_add(*velocity, sw_v3_scale(a0, 0.5 * step));
    const sw_vec3 next_position = sw_v3_add(*position, sw_v3_scale(half_velocity, step));
    const sw_vec3 a1 = acceleration(time + step, next_position, context);
    *position = next_position;
    *velocity = sw_v3_add(half_velocity, sw_v3_scale(a1, 0.5 * step));
}
