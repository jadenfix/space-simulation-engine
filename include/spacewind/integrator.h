#ifndef SPACEWIND_INTEGRATOR_H
#define SPACEWIND_INTEGRATOR_H

#include <stdbool.h>
#include <stddef.h>

#include "spacewind/math3.h"

typedef bool (*sw_ode_rhs_fn)(
    double time,
    const double *state,
    double *derivative,
    size_t dimension,
    void *context
);

typedef sw_vec3 (*sw_acceleration_fn)(double time, sw_vec3 position, void *context);

typedef struct {
    size_t dimension;
    size_t vectors;
    double *storage;
} sw_ode_workspace;

bool sw_ode_workspace_init(sw_ode_workspace *workspace, size_t dimension, size_t vectors);
void sw_ode_workspace_destroy(sw_ode_workspace *workspace);

bool sw_rk4_step(
    sw_ode_rhs_fn rhs,
    void *context,
    double time,
    double step,
    double *state,
    size_t dimension,
    sw_ode_workspace *workspace
);

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
);

void sw_velocity_verlet_step(
    sw_acceleration_fn acceleration,
    void *context,
    double time,
    double step,
    sw_vec3 *position,
    sw_vec3 *velocity
);

#endif
