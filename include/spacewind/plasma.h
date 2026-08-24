#ifndef SPACEWIND_PLASMA_H
#define SPACEWIND_PLASMA_H

#include <stdbool.h>
#include <stddef.h>

#include "spacewind/math3.h"

typedef struct {
    sw_vec3 position_m;
    sw_vec3 velocity_m_s;
    double charge_c;
    double mass_kg;
    double macro_weight;
} sw_particle;

void sw_boris_push(
    sw_particle *particle,
    sw_vec3 electric_field_v_m,
    sw_vec3 magnetic_field_t,
    double step_s
);

void sw_relativistic_boris_push(
    sw_particle *particle,
    sw_vec3 electric_field_v_m,
    sw_vec3 magnetic_field_t,
    double step_s
);

typedef struct {
    size_t grid_count;
    size_t particle_count;
    double length_m;
    double dx_m;
    double step_s;
    sw_particle *particles;
    double *charge_density_c_m3;
    double *electric_field_v_m;
    double *potential_v;
    double neutralizing_background_c_m3;
    double *previous_unwrapped_position_x_m;
    double *unwrapped_position_x_m;
    bool unwrapped_positions_initialized;
} sw_pic1d;

bool sw_pic1d_init(
    sw_pic1d *pic,
    size_t grid_count,
    size_t particle_count,
    double length_m,
    double step_s
);

void sw_pic1d_destroy(sw_pic1d *pic);
void sw_pic1d_quiet_start(
    sw_pic1d *pic,
    double charge_c,
    double mass_kg,
    double number_density_m3,
    double drift_velocity_m_s,
    double perturbation_fraction,
    unsigned mode_number
);
bool sw_pic1d_sync_unwrapped_positions(sw_pic1d *pic);
bool sw_pic1d_unwrapped_positions_consistent(
    const sw_pic1d *pic,
    double absolute_tolerance_m
);
void sw_pic1d_deposit_charge(sw_pic1d *pic);
void sw_pic1d_solve_poisson_spectral(sw_pic1d *pic);
double sw_pic1d_interpolate_electric(const sw_pic1d *pic, double x_m);
void sw_pic1d_step(sw_pic1d *pic);
double sw_pic1d_kinetic_energy_j(const sw_pic1d *pic);
double sw_pic1d_field_energy_j(const sw_pic1d *pic);
double sw_pic1d_plasma_frequency_rad_s(const sw_pic1d *pic);

typedef struct {
    size_t cell_count;
    double length_m;
    double dx_m;
    double step_s;
    double *electric_y_v_m;
    double *magnetic_z_t;
} sw_fdtd1d;

bool sw_fdtd1d_init(sw_fdtd1d *grid, size_t cell_count, double length_m, double step_s);
void sw_fdtd1d_destroy(sw_fdtd1d *grid);
void sw_fdtd1d_gaussian_right_wave(
    sw_fdtd1d *grid,
    double center_m,
    double width_m,
    double amplitude_v_m
);
void sw_fdtd1d_step(sw_fdtd1d *grid);
double sw_fdtd1d_energy_j_m2(const sw_fdtd1d *grid);

#endif
