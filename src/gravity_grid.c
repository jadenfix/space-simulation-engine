#include "spacewind/gravity_grid.h"
#include "spacewind/constants.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static bool sw_gravity_grid_dimensions_valid(size_t nx, size_t ny, size_t nz) {
    if (nx < 3U || ny < 3U || nz < 3U) {
        return false;
    }
    if (nx > SIZE_MAX / ny) {
        return false;
    }
    if (nx * ny > SIZE_MAX / nz) {
        return false;
    }
    return true;
}

bool sw_gravity_grid_init(
    sw_gravity_grid *grid,
    size_t nx,
    size_t ny,
    size_t nz,
    double spacing_m,
    sw_vec3 origin_m
) {
    size_t count;
    if (grid == NULL || !sw_gravity_grid_dimensions_valid(nx, ny, nz)
        || !(spacing_m > 0.0) || !isfinite(spacing_m) || !sw_v3_isfinite(origin_m)) {
        return false;
    }
    memset(grid, 0, sizeof(*grid));
    count = nx * ny * nz;
    grid->density_kg_m3 = (double *)calloc(count, sizeof(double));
    grid->potential_m2_s2 = (double *)calloc(count, sizeof(double));
    if (grid->density_kg_m3 == NULL || grid->potential_m2_s2 == NULL) {
        sw_gravity_grid_destroy(grid);
        return false;
    }
    grid->nx = nx;
    grid->ny = ny;
    grid->nz = nz;
    grid->spacing_m = spacing_m;
    grid->origin_m = origin_m;
    return true;
}

void sw_gravity_grid_destroy(sw_gravity_grid *grid) {
    if (grid == NULL) {
        return;
    }
    free(grid->density_kg_m3);
    free(grid->potential_m2_s2);
    memset(grid, 0, sizeof(*grid));
}

void sw_gravity_grid_clear(sw_gravity_grid *grid) {
    size_t count;
    if (grid == NULL || grid->density_kg_m3 == NULL || grid->potential_m2_s2 == NULL) {
        return;
    }
    count = grid->nx * grid->ny * grid->nz;
    memset(grid->density_kg_m3, 0, count * sizeof(double));
    memset(grid->potential_m2_s2, 0, count * sizeof(double));
}

size_t sw_gravity_grid_index(const sw_gravity_grid *grid, size_t i, size_t j, size_t k) {
    return (k * grid->ny + j) * grid->nx + i;
}

sw_vec3 sw_gravity_grid_position(const sw_gravity_grid *grid, size_t i, size_t j, size_t k) {
    return sw_v3(
        grid->origin_m.x + (double)i * grid->spacing_m,
        grid->origin_m.y + (double)j * grid->spacing_m,
        grid->origin_m.z + (double)k * grid->spacing_m
    );
}

bool sw_gravity_grid_deposit_cic(
    sw_gravity_grid *grid,
    double mass_kg,
    sw_vec3 position_m
) {
    double gx;
    double gy;
    double gz;
    size_t i0;
    size_t j0;
    size_t k0;
    double fx;
    double fy;
    double fz;
    size_t di;
    size_t dj;
    size_t dk;
    double inverse_volume;

    if (grid == NULL || grid->density_kg_m3 == NULL || !(mass_kg >= 0.0)
        || !isfinite(mass_kg) || !sw_v3_isfinite(position_m)) {
        return false;
    }
    gx = (position_m.x - grid->origin_m.x) / grid->spacing_m;
    gy = (position_m.y - grid->origin_m.y) / grid->spacing_m;
    gz = (position_m.z - grid->origin_m.z) / grid->spacing_m;
    inverse_volume = 1.0 / (
        grid->spacing_m * grid->spacing_m * grid->spacing_m
    );
    if (gx < 0.0 || gy < 0.0 || gz < 0.0
        || gx >= (double)(grid->nx - 1U)
        || gy >= (double)(grid->ny - 1U)
        || gz >= (double)(grid->nz - 1U)) {
        return false;
    }
    i0 = (size_t)floor(gx);
    j0 = (size_t)floor(gy);
    k0 = (size_t)floor(gz);
    fx = gx - (double)i0;
    fy = gy - (double)j0;
    fz = gz - (double)k0;
    for (dk = 0U; dk < 2U; ++dk) {
        const double wz = dk == 0U ? 1.0 - fz : fz;
        for (dj = 0U; dj < 2U; ++dj) {
            const double wy = dj == 0U ? 1.0 - fy : fy;
            for (di = 0U; di < 2U; ++di) {
                const double wx = di == 0U ? 1.0 - fx : fx;
                const size_t index = sw_gravity_grid_index(grid, i0 + di, j0 + dj, k0 + dk);
                grid->density_kg_m3[index] += mass_kg * wx * wy * wz * inverse_volume;
            }
        }
    }
    return true;
}

void sw_gravity_grid_set_monopole_boundary(
    sw_gravity_grid *grid,
    double total_mass_kg,
    sw_vec3 center_m,
    double minimum_radius_m
) {
    size_t i;
    size_t j;
    size_t k;
    if (grid == NULL || grid->potential_m2_s2 == NULL || !(total_mass_kg >= 0.0)
        || !isfinite(total_mass_kg)) {
        return;
    }
    for (k = 0U; k < grid->nz; ++k) {
        for (j = 0U; j < grid->ny; ++j) {
            for (i = 0U; i < grid->nx; ++i) {
                if (i == 0U || j == 0U || k == 0U
                    || i + 1U == grid->nx || j + 1U == grid->ny || k + 1U == grid->nz) {
                    const sw_vec3 position = sw_gravity_grid_position(grid, i, j, k);
                    const double radius = fmax(
                        sw_v3_norm(sw_v3_sub(position, center_m)),
                        minimum_radius_m
                    );
                    grid->potential_m2_s2[sw_gravity_grid_index(grid, i, j, k)] =
                        radius > 0.0 ? -SW_G * total_mass_kg / radius : 0.0;
                }
            }
        }
    }
}

static double sw_gravity_grid_source(const sw_gravity_grid *grid, size_t index) {
    return 4.0 * SW_PI * SW_G * grid->density_kg_m3[index];
}

double sw_gravity_grid_relative_residual(const sw_gravity_grid *grid) {
    double residual2 = 0.0;
    double scale2 = 0.0;
    double inverse_h2;
    size_t i;
    size_t j;
    size_t k;
    size_t count = 0U;
    if (grid == NULL || grid->density_kg_m3 == NULL || grid->potential_m2_s2 == NULL) {
        return INFINITY;
    }
    inverse_h2 = 1.0 / (grid->spacing_m * grid->spacing_m);
    for (k = 1U; k + 1U < grid->nz; ++k) {
        for (j = 1U; j + 1U < grid->ny; ++j) {
            for (i = 1U; i + 1U < grid->nx; ++i) {
                const size_t index = sw_gravity_grid_index(grid, i, j, k);
                const double laplacian = (
                    grid->potential_m2_s2[sw_gravity_grid_index(grid, i - 1U, j, k)]
                    + grid->potential_m2_s2[sw_gravity_grid_index(grid, i + 1U, j, k)]
                    + grid->potential_m2_s2[sw_gravity_grid_index(grid, i, j - 1U, k)]
                    + grid->potential_m2_s2[sw_gravity_grid_index(grid, i, j + 1U, k)]
                    + grid->potential_m2_s2[sw_gravity_grid_index(grid, i, j, k - 1U)]
                    + grid->potential_m2_s2[sw_gravity_grid_index(grid, i, j, k + 1U)]
                    - 6.0 * grid->potential_m2_s2[index]
                ) * inverse_h2;
                const double source = sw_gravity_grid_source(grid, index);
                const double residual = laplacian - source;
                residual2 += residual * residual;
                scale2 += source * source;
                ++count;
            }
        }
    }
    if (count == 0U) {
        return INFINITY;
    }
    if (scale2 <= 0.0) {
        return sqrt(residual2 / (double)count);
    }
    return sqrt(residual2 / scale2);
}

bool sw_gravity_grid_solve_sor(
    sw_gravity_grid *grid,
    double relaxation,
    double relative_tolerance,
    size_t max_iterations,
    size_t *iterations_used,
    double *final_relative_residual
) {
    double h2;
    size_t iteration;
    if (grid == NULL || grid->density_kg_m3 == NULL || grid->potential_m2_s2 == NULL
        || !(relaxation > 0.0 && relaxation < 2.0)
        || !(relative_tolerance > 0.0) || max_iterations == 0U) {
        return false;
    }
    h2 = grid->spacing_m * grid->spacing_m;
    for (iteration = 0U; iteration < max_iterations; ++iteration) {
        unsigned parity;
        for (parity = 0U; parity < 2U; ++parity) {
            size_t i;
            size_t j;
            size_t k;
            for (k = 1U; k + 1U < grid->nz; ++k) {
                for (j = 1U; j + 1U < grid->ny; ++j) {
                    for (i = 1U; i + 1U < grid->nx; ++i) {
                        const size_t index = sw_gravity_grid_index(grid, i, j, k);
                        double target;
                        if (((i + j + k) & 1U) != parity) {
                            continue;
                        }
                        target = (
                            grid->potential_m2_s2[sw_gravity_grid_index(grid, i - 1U, j, k)]
                            + grid->potential_m2_s2[sw_gravity_grid_index(grid, i + 1U, j, k)]
                            + grid->potential_m2_s2[sw_gravity_grid_index(grid, i, j - 1U, k)]
                            + grid->potential_m2_s2[sw_gravity_grid_index(grid, i, j + 1U, k)]
                            + grid->potential_m2_s2[sw_gravity_grid_index(grid, i, j, k - 1U)]
                            + grid->potential_m2_s2[sw_gravity_grid_index(grid, i, j, k + 1U)]
                            - h2 * sw_gravity_grid_source(grid, index)
                        ) / 6.0;
                        grid->potential_m2_s2[index] += relaxation * (
                            target - grid->potential_m2_s2[index]
                        );
                    }
                }
            }
        }
        if ((iteration + 1U) % 4U == 0U || iteration + 1U == max_iterations) {
            const double residual = sw_gravity_grid_relative_residual(grid);
            if (!isfinite(residual)) {
                return false;
            }
            if (residual <= relative_tolerance) {
                if (iterations_used != NULL) {
                    *iterations_used = iteration + 1U;
                }
                if (final_relative_residual != NULL) {
                    *final_relative_residual = residual;
                }
                return true;
            }
        }
    }
    if (iterations_used != NULL) {
        *iterations_used = max_iterations;
    }
    if (final_relative_residual != NULL) {
        *final_relative_residual = sw_gravity_grid_relative_residual(grid);
    }
    return false;
}

static bool sw_gravity_grid_cell_coordinates(
    const sw_gravity_grid *grid,
    sw_vec3 position_m,
    size_t *i0,
    size_t *j0,
    size_t *k0,
    double *fx,
    double *fy,
    double *fz
) {
    const double gx = (position_m.x - grid->origin_m.x) / grid->spacing_m;
    const double gy = (position_m.y - grid->origin_m.y) / grid->spacing_m;
    const double gz = (position_m.z - grid->origin_m.z) / grid->spacing_m;
    if (gx < 0.0 || gy < 0.0 || gz < 0.0
        || gx >= (double)(grid->nx - 1U)
        || gy >= (double)(grid->ny - 1U)
        || gz >= (double)(grid->nz - 1U)) {
        return false;
    }
    *i0 = (size_t)floor(gx);
    *j0 = (size_t)floor(gy);
    *k0 = (size_t)floor(gz);
    *fx = gx - (double)(*i0);
    *fy = gy - (double)(*j0);
    *fz = gz - (double)(*k0);
    return true;
}

bool sw_gravity_grid_sample_potential(
    const sw_gravity_grid *grid,
    sw_vec3 position_m,
    double *potential_m2_s2
) {
    size_t i0;
    size_t j0;
    size_t k0;
    double fx;
    double fy;
    double fz;
    double value = 0.0;
    size_t di;
    size_t dj;
    size_t dk;
    if (grid == NULL || grid->potential_m2_s2 == NULL || potential_m2_s2 == NULL
        || !sw_v3_isfinite(position_m)) {
        return false;
    }
    if (!sw_gravity_grid_cell_coordinates(
        grid, position_m, &i0, &j0, &k0, &fx, &fy, &fz
    )) {
        return false;
    }
    for (dk = 0U; dk < 2U; ++dk) {
        const double wz = dk == 0U ? 1.0 - fz : fz;
        for (dj = 0U; dj < 2U; ++dj) {
            const double wy = dj == 0U ? 1.0 - fy : fy;
            for (di = 0U; di < 2U; ++di) {
                const double wx = di == 0U ? 1.0 - fx : fx;
                value += wx * wy * wz * grid->potential_m2_s2[
                    sw_gravity_grid_index(grid, i0 + di, j0 + dj, k0 + dk)
                ];
            }
        }
    }
    *potential_m2_s2 = value;
    return isfinite(value);
}

bool sw_gravity_grid_sample_acceleration(
    const sw_gravity_grid *grid,
    sw_vec3 position_m,
    sw_vec3 *acceleration_m_s2
) {
    double h;
    double xp;
    double xm;
    double yp;
    double ym;
    double zp;
    double zm;
    if (grid == NULL || acceleration_m_s2 == NULL) {
        return false;
    }
    h = 0.5 * grid->spacing_m;
    if (!sw_gravity_grid_sample_potential(grid, sw_v3_add(position_m, sw_v3(h, 0.0, 0.0)), &xp)
        || !sw_gravity_grid_sample_potential(grid, sw_v3_sub(position_m, sw_v3(h, 0.0, 0.0)), &xm)
        || !sw_gravity_grid_sample_potential(grid, sw_v3_add(position_m, sw_v3(0.0, h, 0.0)), &yp)
        || !sw_gravity_grid_sample_potential(grid, sw_v3_sub(position_m, sw_v3(0.0, h, 0.0)), &ym)
        || !sw_gravity_grid_sample_potential(grid, sw_v3_add(position_m, sw_v3(0.0, 0.0, h)), &zp)
        || !sw_gravity_grid_sample_potential(grid, sw_v3_sub(position_m, sw_v3(0.0, 0.0, h)), &zm)) {
        return false;
    }
    *acceleration_m_s2 = sw_v3(
        -(xp - xm) / (2.0 * h),
        -(yp - ym) / (2.0 * h),
        -(zp - zm) / (2.0 * h)
    );
    return sw_v3_isfinite(*acceleration_m_s2);
}

bool sw_gravity_grid_weak_field_metric(
    const sw_gravity_grid *grid,
    sw_vec3 position_m,
    double g_cov[4][4]
) {
    double potential;
    double temporal;
    double spatial;
    if (g_cov == NULL || !sw_gravity_grid_sample_potential(grid, position_m, &potential)) {
        return false;
    }
    temporal = -(1.0 + 2.0 * potential / (SW_C * SW_C));
    spatial = 1.0 - 2.0 * potential / (SW_C * SW_C);
    memset(g_cov, 0, sizeof(double) * 16U);
    g_cov[0][0] = temporal;
    g_cov[1][1] = spatial;
    g_cov[2][2] = spatial;
    g_cov[3][3] = spatial;
    return true;
}
