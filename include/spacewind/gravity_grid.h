#ifndef SPACEWIND_GRAVITY_GRID_H
#define SPACEWIND_GRAVITY_GRID_H

#include <stdbool.h>
#include <stddef.h>

#include "spacewind/math3.h"

typedef struct {
    size_t nx;
    size_t ny;
    size_t nz;
    double spacing_m;
    sw_vec3 origin_m;
    double *density_kg_m3;
    double *potential_m2_s2;
} sw_gravity_grid;

bool sw_gravity_grid_init(
    sw_gravity_grid *grid,
    size_t nx,
    size_t ny,
    size_t nz,
    double spacing_m,
    sw_vec3 origin_m
);
void sw_gravity_grid_destroy(sw_gravity_grid *grid);
void sw_gravity_grid_clear(sw_gravity_grid *grid);
size_t sw_gravity_grid_index(const sw_gravity_grid *grid, size_t i, size_t j, size_t k);
sw_vec3 sw_gravity_grid_position(const sw_gravity_grid *grid, size_t i, size_t j, size_t k);

bool sw_gravity_grid_deposit_cic(
    sw_gravity_grid *grid,
    double mass_kg,
    sw_vec3 position_m
);
void sw_gravity_grid_set_monopole_boundary(
    sw_gravity_grid *grid,
    double total_mass_kg,
    sw_vec3 center_m,
    double minimum_radius_m
);

bool sw_gravity_grid_solve_sor(
    sw_gravity_grid *grid,
    double relaxation,
    double relative_tolerance,
    size_t max_iterations,
    size_t *iterations_used,
    double *final_relative_residual
);

double sw_gravity_grid_relative_residual(const sw_gravity_grid *grid);
bool sw_gravity_grid_sample_potential(
    const sw_gravity_grid *grid,
    sw_vec3 position_m,
    double *potential_m2_s2
);
bool sw_gravity_grid_sample_acceleration(
    const sw_gravity_grid *grid,
    sw_vec3 position_m,
    sw_vec3 *acceleration_m_s2
);
bool sw_gravity_grid_weak_field_metric(
    const sw_gravity_grid *grid,
    sw_vec3 position_m,
    double g_cov[4][4]
);

#endif
