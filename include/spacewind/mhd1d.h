#ifndef SPACEWIND_MHD1D_H
#define SPACEWIND_MHD1D_H

#include <stdbool.h>
#include <stddef.h>

#define SW_MHD1D_VARIABLES 7U

typedef struct {
    double density;
    double velocity_x;
    double velocity_y;
    double velocity_z;
    double pressure;
    double magnetic_y;
    double magnetic_z;
} sw_mhd_primitive;

typedef struct {
    size_t cell_count;
    double length;
    double dx;
    double gamma;
    double magnetic_permeability;
    double magnetic_x;
    double density_floor;
    double pressure_floor;
    double *conserved;
    double *scratch;
    double *interface_flux;
} sw_mhd1d;

bool sw_mhd1d_init(
    sw_mhd1d *grid,
    size_t cell_count,
    double length,
    double gamma,
    double magnetic_permeability,
    double magnetic_x
);
void sw_mhd1d_destroy(sw_mhd1d *grid);
void sw_mhd1d_set_cell(sw_mhd1d *grid, size_t index, sw_mhd_primitive primitive);
sw_mhd_primitive sw_mhd1d_get_cell(const sw_mhd1d *grid, size_t index);
void sw_mhd1d_init_brio_wu(sw_mhd1d *grid);
double sw_mhd1d_max_signal_speed(const sw_mhd1d *grid);
double sw_mhd1d_cfl_step(const sw_mhd1d *grid, double cfl);
bool sw_mhd1d_step(sw_mhd1d *grid, double step);
double sw_mhd1d_total_mass(const sw_mhd1d *grid);
double sw_mhd1d_total_energy(const sw_mhd1d *grid);

#endif
