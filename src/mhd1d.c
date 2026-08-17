#include "spacewind/mhd1d.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

enum {
    SW_MHD_RHO = 0,
    SW_MHD_MX = 1,
    SW_MHD_MY = 2,
    SW_MHD_MZ = 3,
    SW_MHD_BY = 4,
    SW_MHD_BZ = 5,
    SW_MHD_ENERGY = 6
};

static double *sw_mhd_cell(double *array, size_t index) {
    return array + index * SW_MHD1D_VARIABLES;
}

static const double *sw_mhd_const_cell(const double *array, size_t index) {
    return array + index * SW_MHD1D_VARIABLES;
}

bool sw_mhd1d_init(
    sw_mhd1d *grid,
    size_t cell_count,
    double length,
    double gamma,
    double magnetic_permeability,
    double magnetic_x
) {
    if (grid == NULL || cell_count < 8U || length <= 0.0 || gamma <= 1.0
        || magnetic_permeability <= 0.0) {
        return false;
    }
    memset(grid, 0, sizeof(*grid));
    grid->cell_count = cell_count;
    grid->length = length;
    grid->dx = length / (double)cell_count;
    grid->gamma = gamma;
    grid->magnetic_permeability = magnetic_permeability;
    grid->magnetic_x = magnetic_x;
    grid->density_floor = 1.0e-12;
    grid->pressure_floor = 1.0e-12;
    grid->conserved = (double *)calloc(cell_count * SW_MHD1D_VARIABLES, sizeof(double));
    grid->scratch = (double *)calloc(cell_count * SW_MHD1D_VARIABLES, sizeof(double));
    grid->interface_flux = (double *)calloc((cell_count + 1U) * SW_MHD1D_VARIABLES, sizeof(double));
    if (grid->conserved == NULL || grid->scratch == NULL || grid->interface_flux == NULL) {
        sw_mhd1d_destroy(grid);
        return false;
    }
    return true;
}

void sw_mhd1d_destroy(sw_mhd1d *grid) {
    if (grid == NULL) {
        return;
    }
    free(grid->conserved);
    free(grid->scratch);
    free(grid->interface_flux);
    memset(grid, 0, sizeof(*grid));
}

static void sw_mhd_primitive_to_conserved(
    const sw_mhd1d *grid,
    sw_mhd_primitive p,
    double u[SW_MHD1D_VARIABLES]
) {
    const double v2 = p.velocity_x * p.velocity_x
        + p.velocity_y * p.velocity_y
        + p.velocity_z * p.velocity_z;
    const double b2 = grid->magnetic_x * grid->magnetic_x
        + p.magnetic_y * p.magnetic_y
        + p.magnetic_z * p.magnetic_z;
    p.density = fmax(p.density, grid->density_floor);
    p.pressure = fmax(p.pressure, grid->pressure_floor);
    u[SW_MHD_RHO] = p.density;
    u[SW_MHD_MX] = p.density * p.velocity_x;
    u[SW_MHD_MY] = p.density * p.velocity_y;
    u[SW_MHD_MZ] = p.density * p.velocity_z;
    u[SW_MHD_BY] = p.magnetic_y;
    u[SW_MHD_BZ] = p.magnetic_z;
    u[SW_MHD_ENERGY] = p.pressure / (grid->gamma - 1.0)
        + 0.5 * p.density * v2
        + 0.5 * b2 / grid->magnetic_permeability;
}

static sw_mhd_primitive sw_mhd_conserved_to_primitive(
    const sw_mhd1d *grid,
    const double u[SW_MHD1D_VARIABLES]
) {
    sw_mhd_primitive p;
    const double density = fmax(u[SW_MHD_RHO], grid->density_floor);
    const double vx = u[SW_MHD_MX] / density;
    const double vy = u[SW_MHD_MY] / density;
    const double vz = u[SW_MHD_MZ] / density;
    const double v2 = vx * vx + vy * vy + vz * vz;
    const double b2 = grid->magnetic_x * grid->magnetic_x
        + u[SW_MHD_BY] * u[SW_MHD_BY]
        + u[SW_MHD_BZ] * u[SW_MHD_BZ];
    const double internal = u[SW_MHD_ENERGY]
        - 0.5 * density * v2
        - 0.5 * b2 / grid->magnetic_permeability;
    p.density = density;
    p.velocity_x = vx;
    p.velocity_y = vy;
    p.velocity_z = vz;
    p.pressure = fmax((grid->gamma - 1.0) * internal, grid->pressure_floor);
    p.magnetic_y = u[SW_MHD_BY];
    p.magnetic_z = u[SW_MHD_BZ];
    return p;
}

void sw_mhd1d_set_cell(sw_mhd1d *grid, size_t index, sw_mhd_primitive primitive) {
    if (index >= grid->cell_count) {
        return;
    }
    sw_mhd_primitive_to_conserved(grid, primitive, sw_mhd_cell(grid->conserved, index));
}

sw_mhd_primitive sw_mhd1d_get_cell(const sw_mhd1d *grid, size_t index) {
    sw_mhd_primitive zero;
    memset(&zero, 0, sizeof(zero));
    if (index >= grid->cell_count) {
        return zero;
    }
    return sw_mhd_conserved_to_primitive(grid, sw_mhd_const_cell(grid->conserved, index));
}

void sw_mhd1d_init_brio_wu(sw_mhd1d *grid) {
    size_t i;
    for (i = 0U; i < grid->cell_count; ++i) {
        sw_mhd_primitive p;
        memset(&p, 0, sizeof(p));
        if (i < grid->cell_count / 2U) {
            p.density = 1.0;
            p.pressure = 1.0;
            p.magnetic_y = 1.0;
        } else {
            p.density = 0.125;
            p.pressure = 0.1;
            p.magnetic_y = -1.0;
        }
        sw_mhd1d_set_cell(grid, i, p);
    }
}

static double sw_mhd_fast_speed(const sw_mhd1d *grid, sw_mhd_primitive p) {
    const double rho = fmax(p.density, grid->density_floor);
    const double a2 = grid->gamma * p.pressure / rho;
    const double b2 = (
        grid->magnetic_x * grid->magnetic_x
        + p.magnetic_y * p.magnetic_y
        + p.magnetic_z * p.magnetic_z
    ) / (grid->magnetic_permeability * rho);
    const double bx2 = grid->magnetic_x * grid->magnetic_x
        / (grid->magnetic_permeability * rho);
    const double discriminant = fmax((a2 + b2) * (a2 + b2) - 4.0 * a2 * bx2, 0.0);
    return sqrt(fmax(0.5 * (a2 + b2 + sqrt(discriminant)), 0.0));
}

static void sw_mhd_flux(
    const sw_mhd1d *grid,
    const double u[SW_MHD1D_VARIABLES],
    double flux[SW_MHD1D_VARIABLES]
) {
    const sw_mhd_primitive p = sw_mhd_conserved_to_primitive(grid, u);
    const double bdotv = grid->magnetic_x * p.velocity_x
        + p.magnetic_y * p.velocity_y
        + p.magnetic_z * p.velocity_z;
    const double b2 = grid->magnetic_x * grid->magnetic_x
        + p.magnetic_y * p.magnetic_y
        + p.magnetic_z * p.magnetic_z;
    const double total_pressure = p.pressure + 0.5 * b2 / grid->magnetic_permeability;

    flux[SW_MHD_RHO] = u[SW_MHD_MX];
    flux[SW_MHD_MX] = u[SW_MHD_MX] * p.velocity_x
        + total_pressure
        - grid->magnetic_x * grid->magnetic_x / grid->magnetic_permeability;
    flux[SW_MHD_MY] = u[SW_MHD_MY] * p.velocity_x
        - grid->magnetic_x * p.magnetic_y / grid->magnetic_permeability;
    flux[SW_MHD_MZ] = u[SW_MHD_MZ] * p.velocity_x
        - grid->magnetic_x * p.magnetic_z / grid->magnetic_permeability;
    flux[SW_MHD_BY] = p.magnetic_y * p.velocity_x - grid->magnetic_x * p.velocity_y;
    flux[SW_MHD_BZ] = p.magnetic_z * p.velocity_x - grid->magnetic_x * p.velocity_z;
    flux[SW_MHD_ENERGY] = (u[SW_MHD_ENERGY] + total_pressure) * p.velocity_x
        - grid->magnetic_x * bdotv / grid->magnetic_permeability;
}

static void sw_mhd_hll_flux(
    const sw_mhd1d *grid,
    const double left[SW_MHD1D_VARIABLES],
    const double right[SW_MHD1D_VARIABLES],
    double out[SW_MHD1D_VARIABLES]
) {
    const sw_mhd_primitive pl = sw_mhd_conserved_to_primitive(grid, left);
    const sw_mhd_primitive pr = sw_mhd_conserved_to_primitive(grid, right);
    const double cfl = sw_mhd_fast_speed(grid, pl);
    const double cfr = sw_mhd_fast_speed(grid, pr);
    const double sl = fmin(pl.velocity_x - cfl, pr.velocity_x - cfr);
    const double sr = fmax(pl.velocity_x + cfl, pr.velocity_x + cfr);
    double fl[SW_MHD1D_VARIABLES];
    double fr[SW_MHD1D_VARIABLES];
    size_t v;
    sw_mhd_flux(grid, left, fl);
    sw_mhd_flux(grid, right, fr);
    if (sl >= 0.0) {
        memcpy(out, fl, sizeof(fl));
    } else if (sr <= 0.0) {
        memcpy(out, fr, sizeof(fr));
    } else {
        const double denom = fmax(sr - sl, DBL_MIN);
        for (v = 0U; v < SW_MHD1D_VARIABLES; ++v) {
            out[v] = (sr * fl[v] - sl * fr[v] + sl * sr * (right[v] - left[v])) / denom;
        }
    }
}

double sw_mhd1d_max_signal_speed(const sw_mhd1d *grid) {
    double maximum = 0.0;
    size_t i;
    for (i = 0U; i < grid->cell_count; ++i) {
        const sw_mhd_primitive p = sw_mhd1d_get_cell(grid, i);
        maximum = fmax(maximum, fabs(p.velocity_x) + sw_mhd_fast_speed(grid, p));
    }
    return maximum;
}

double sw_mhd1d_cfl_step(const sw_mhd1d *grid, double cfl) {
    const double maximum = sw_mhd1d_max_signal_speed(grid);
    if (maximum <= DBL_MIN) {
        return INFINITY;
    }
    return cfl * grid->dx / maximum;
}

bool sw_mhd1d_step(sw_mhd1d *grid, double step) {
    size_t interface_index;
    size_t cell_index;
    size_t v;
    if (step <= 0.0 || !isfinite(step)) {
        return false;
    }

    for (interface_index = 0U; interface_index <= grid->cell_count; ++interface_index) {
        const size_t left_index = interface_index == 0U ? 0U : interface_index - 1U;
        const size_t right_index = interface_index == grid->cell_count
            ? grid->cell_count - 1U
            : interface_index;
        sw_mhd_hll_flux(
            grid,
            sw_mhd_const_cell(grid->conserved, left_index),
            sw_mhd_const_cell(grid->conserved, right_index),
            sw_mhd_cell(grid->interface_flux, interface_index)
        );
    }

    for (cell_index = 0U; cell_index < grid->cell_count; ++cell_index) {
        const double *u = sw_mhd_const_cell(grid->conserved, cell_index);
        const double *flux_left = sw_mhd_const_cell(grid->interface_flux, cell_index);
        const double *flux_right = sw_mhd_const_cell(grid->interface_flux, cell_index + 1U);
        double *next = sw_mhd_cell(grid->scratch, cell_index);
        for (v = 0U; v < SW_MHD1D_VARIABLES; ++v) {
            next[v] = u[v] - step / grid->dx * (flux_right[v] - flux_left[v]);
            if (!isfinite(next[v])) {
                return false;
            }
        }
        {
            sw_mhd_primitive p = sw_mhd_conserved_to_primitive(grid, next);
            sw_mhd_primitive_to_conserved(grid, p, next);
        }
    }

    {
        double *tmp = grid->conserved;
        grid->conserved = grid->scratch;
        grid->scratch = tmp;
    }
    return true;
}

double sw_mhd1d_total_mass(const sw_mhd1d *grid) {
    double total = 0.0;
    size_t i;
    for (i = 0U; i < grid->cell_count; ++i) {
        total += sw_mhd_const_cell(grid->conserved, i)[SW_MHD_RHO] * grid->dx;
    }
    return total;
}

double sw_mhd1d_total_energy(const sw_mhd1d *grid) {
    double total = 0.0;
    size_t i;
    for (i = 0U; i < grid->cell_count; ++i) {
        total += sw_mhd_const_cell(grid->conserved, i)[SW_MHD_ENERGY] * grid->dx;
    }
    return total;
}
