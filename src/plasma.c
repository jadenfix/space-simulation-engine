#include "spacewind/plasma.h"
#include "spacewind/constants.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double sw_wrap_periodic(double x, double length) {
    double wrapped = fmod(x, length);
    if (wrapped < 0.0) {
        wrapped += length;
    }
    return wrapped;
}

void sw_boris_push(
    sw_particle *particle,
    sw_vec3 electric_field_v_m,
    sw_vec3 magnetic_field_t,
    double step_s
) {
    const double qmdt2 = particle->charge_c * step_s / (2.0 * particle->mass_kg);
    const sw_vec3 v_minus = sw_v3_add(
        particle->velocity_m_s,
        sw_v3_scale(electric_field_v_m, qmdt2)
    );
    const sw_vec3 t = sw_v3_scale(magnetic_field_t, qmdt2);
    const sw_vec3 s = sw_v3_scale(t, 2.0 / (1.0 + sw_v3_norm2(t)));
    const sw_vec3 v_prime = sw_v3_add(v_minus, sw_v3_cross(v_minus, t));
    const sw_vec3 v_plus = sw_v3_add(v_minus, sw_v3_cross(v_prime, s));
    particle->velocity_m_s = sw_v3_add(v_plus, sw_v3_scale(electric_field_v_m, qmdt2));
    particle->position_m = sw_v3_add(
        particle->position_m,
        sw_v3_scale(particle->velocity_m_s, step_s)
    );
}

void sw_relativistic_boris_push(
    sw_particle *particle,
    sw_vec3 electric_field_v_m,
    sw_vec3 magnetic_field_t,
    double step_s
) {
    const double c2 = SW_C * SW_C;
    const double v2 = sw_v3_norm2(particle->velocity_m_s);
    const double gamma0 = 1.0 / sqrt(fmax(1.0 - fmin(v2 / c2, 1.0 - 1.0e-15), 1.0e-15));
    sw_vec3 u = sw_v3_scale(particle->velocity_m_s, gamma0);
    const double qmdt2 = particle->charge_c * step_s / (2.0 * particle->mass_kg);
    sw_vec3 u_minus = sw_v3_add(u, sw_v3_scale(electric_field_v_m, qmdt2));
    const double gamma_minus = sqrt(1.0 + sw_v3_norm2(u_minus) / c2);
    const sw_vec3 t = sw_v3_scale(magnetic_field_t, qmdt2 / gamma_minus);
    const sw_vec3 s = sw_v3_scale(t, 2.0 / (1.0 + sw_v3_norm2(t)));
    const sw_vec3 u_prime = sw_v3_add(u_minus, sw_v3_cross(u_minus, t));
    const sw_vec3 u_plus = sw_v3_add(u_minus, sw_v3_cross(u_prime, s));
    u = sw_v3_add(u_plus, sw_v3_scale(electric_field_v_m, qmdt2));
    {
        const double gamma_new = sqrt(1.0 + sw_v3_norm2(u) / c2);
        particle->velocity_m_s = sw_v3_scale(u, 1.0 / gamma_new);
    }
    particle->position_m = sw_v3_add(
        particle->position_m,
        sw_v3_scale(particle->velocity_m_s, step_s)
    );
}

bool sw_pic1d_init(
    sw_pic1d *pic,
    size_t grid_count,
    size_t particle_count,
    double length_m,
    double step_s
) {
    if (pic == NULL || grid_count < 4U || particle_count == 0U || length_m <= 0.0 || step_s <= 0.0) {
        return false;
    }
    memset(pic, 0, sizeof(*pic));
    pic->grid_count = grid_count;
    pic->particle_count = particle_count;
    pic->length_m = length_m;
    pic->dx_m = length_m / (double)grid_count;
    pic->step_s = step_s;
    pic->particles = (sw_particle *)calloc(particle_count, sizeof(sw_particle));
    pic->charge_density_c_m3 = (double *)calloc(grid_count, sizeof(double));
    pic->electric_field_v_m = (double *)calloc(grid_count, sizeof(double));
    pic->potential_v = (double *)calloc(grid_count, sizeof(double));
    if (pic->particles == NULL || pic->charge_density_c_m3 == NULL
        || pic->electric_field_v_m == NULL || pic->potential_v == NULL) {
        sw_pic1d_destroy(pic);
        return false;
    }
    return true;
}

void sw_pic1d_destroy(sw_pic1d *pic) {
    if (pic == NULL) {
        return;
    }
    free(pic->particles);
    free(pic->charge_density_c_m3);
    free(pic->electric_field_v_m);
    free(pic->potential_v);
    memset(pic, 0, sizeof(*pic));
}

void sw_pic1d_quiet_start(
    sw_pic1d *pic,
    double charge_c,
    double mass_kg,
    double number_density_m3,
    double drift_velocity_m_s,
    double perturbation_fraction,
    unsigned mode_number
) {
    size_t i;
    const double macro_weight = number_density_m3 * pic->length_m / (double)pic->particle_count;
    const double perturbation_velocity = perturbation_fraction * drift_velocity_m_s;
    for (i = 0U; i < pic->particle_count; ++i) {
        const double phase = SW_TWO_PI * (double)i / (double)pic->particle_count;
        pic->particles[i].position_m = sw_v3(
            pic->length_m * ((double)i + 0.5) / (double)pic->particle_count,
            0.0,
            0.0
        );
        pic->particles[i].velocity_m_s = sw_v3(
            drift_velocity_m_s + perturbation_velocity * sin((double)mode_number * phase),
            0.0,
            0.0
        );
        pic->particles[i].charge_c = charge_c;
        pic->particles[i].mass_kg = mass_kg;
        pic->particles[i].macro_weight = macro_weight;
    }
    pic->neutralizing_background_c_m3 = -number_density_m3 * charge_c;
}

void sw_pic1d_deposit_charge(sw_pic1d *pic) {
    size_t i;
    memset(pic->charge_density_c_m3, 0, pic->grid_count * sizeof(double));
    for (i = 0U; i < pic->particle_count; ++i) {
        const sw_particle *p = &pic->particles[i];
        const double x = sw_wrap_periodic(p->position_m.x, pic->length_m);
        const double coordinate = x / pic->dx_m;
        const size_t left = ((size_t)floor(coordinate)) % pic->grid_count;
        const size_t right = (left + 1U) % pic->grid_count;
        const double fraction = coordinate - floor(coordinate);
        const double charge_density = p->charge_c * p->macro_weight / pic->dx_m;
        pic->charge_density_c_m3[left] += charge_density * (1.0 - fraction);
        pic->charge_density_c_m3[right] += charge_density * fraction;
    }
    for (i = 0U; i < pic->grid_count; ++i) {
        pic->charge_density_c_m3[i] += pic->neutralizing_background_c_m3;
    }
}

void sw_pic1d_solve_poisson_spectral(sw_pic1d *pic) {
    const size_t n = pic->grid_count;
    size_t j;
    size_t k;
    memset(pic->electric_field_v_m, 0, n * sizeof(double));
    memset(pic->potential_v, 0, n * sizeof(double));

    for (k = 1U; k < n; ++k) {
        const long signed_mode = (k <= n / 2U) ? (long)k : (long)k - (long)n;
        const double wave_number = SW_TWO_PI * (double)signed_mode / pic->length_m;
        double rho_real = 0.0;
        double rho_imag = 0.0;
        double phi_real;
        double phi_imag;
        double e_real;
        double e_imag;
        for (j = 0U; j < n; ++j) {
            const double angle = -SW_TWO_PI * (double)(k * j) / (double)n;
            rho_real += pic->charge_density_c_m3[j] * cos(angle);
            rho_imag += pic->charge_density_c_m3[j] * sin(angle);
        }
        phi_real = rho_real / (SW_EPS0 * wave_number * wave_number);
        phi_imag = rho_imag / (SW_EPS0 * wave_number * wave_number);
        e_real = wave_number * phi_imag;
        e_imag = -wave_number * phi_real;
        for (j = 0U; j < n; ++j) {
            const double angle = SW_TWO_PI * (double)(k * j) / (double)n;
            pic->potential_v[j] += (phi_real * cos(angle) - phi_imag * sin(angle)) / (double)n;
            pic->electric_field_v_m[j] += (e_real * cos(angle) - e_imag * sin(angle)) / (double)n;
        }
    }
}

double sw_pic1d_interpolate_electric(const sw_pic1d *pic, double x_m) {
    const double x = sw_wrap_periodic(x_m, pic->length_m);
    const double coordinate = x / pic->dx_m;
    const size_t left = ((size_t)floor(coordinate)) % pic->grid_count;
    const size_t right = (left + 1U) % pic->grid_count;
    const double fraction = coordinate - floor(coordinate);
    return pic->electric_field_v_m[left] * (1.0 - fraction)
        + pic->electric_field_v_m[right] * fraction;
}

void sw_pic1d_step(sw_pic1d *pic) {
    size_t i;
    sw_pic1d_deposit_charge(pic);
    sw_pic1d_solve_poisson_spectral(pic);
    for (i = 0U; i < pic->particle_count; ++i) {
        sw_particle *p = &pic->particles[i];
        const double e = sw_pic1d_interpolate_electric(pic, p->position_m.x);
        p->velocity_m_s.x += 0.5 * pic->step_s * p->charge_c * e / p->mass_kg;
        p->position_m.x = sw_wrap_periodic(
            p->position_m.x + pic->step_s * p->velocity_m_s.x,
            pic->length_m
        );
    }
    sw_pic1d_deposit_charge(pic);
    sw_pic1d_solve_poisson_spectral(pic);
    for (i = 0U; i < pic->particle_count; ++i) {
        sw_particle *p = &pic->particles[i];
        const double e = sw_pic1d_interpolate_electric(pic, p->position_m.x);
        p->velocity_m_s.x += 0.5 * pic->step_s * p->charge_c * e / p->mass_kg;
    }
}

double sw_pic1d_kinetic_energy_j(const sw_pic1d *pic) {
    double energy = 0.0;
    size_t i;
    for (i = 0U; i < pic->particle_count; ++i) {
        const sw_particle *p = &pic->particles[i];
        energy += 0.5 * p->mass_kg * p->macro_weight * sw_v3_norm2(p->velocity_m_s);
    }
    return energy;
}

double sw_pic1d_field_energy_j(const sw_pic1d *pic) {
    double energy = 0.0;
    size_t i;
    for (i = 0U; i < pic->grid_count; ++i) {
        energy += 0.5 * SW_EPS0 * pic->electric_field_v_m[i]
            * pic->electric_field_v_m[i] * pic->dx_m;
    }
    return energy;
}

double sw_pic1d_plasma_frequency_rad_s(const sw_pic1d *pic) {
    if (pic->particle_count == 0U || pic->particles[0].mass_kg <= 0.0) {
        return 0.0;
    }
    {
        const double number_density = fabs(pic->neutralizing_background_c_m3 / pic->particles[0].charge_c);
        return sqrt(
            number_density * pic->particles[0].charge_c * pic->particles[0].charge_c
            / (SW_EPS0 * pic->particles[0].mass_kg)
        );
    }
}

bool sw_fdtd1d_init(sw_fdtd1d *grid, size_t cell_count, double length_m, double step_s) {
    if (grid == NULL || cell_count < 4U || length_m <= 0.0 || step_s <= 0.0) {
        return false;
    }
    memset(grid, 0, sizeof(*grid));
    grid->cell_count = cell_count;
    grid->length_m = length_m;
    grid->dx_m = length_m / (double)cell_count;
    grid->step_s = step_s;
    if (SW_C * step_s / grid->dx_m > 1.0) {
        return false;
    }
    grid->electric_y_v_m = (double *)calloc(cell_count, sizeof(double));
    grid->magnetic_z_t = (double *)calloc(cell_count, sizeof(double));
    if (grid->electric_y_v_m == NULL || grid->magnetic_z_t == NULL) {
        sw_fdtd1d_destroy(grid);
        return false;
    }
    return true;
}

void sw_fdtd1d_destroy(sw_fdtd1d *grid) {
    if (grid == NULL) {
        return;
    }
    free(grid->electric_y_v_m);
    free(grid->magnetic_z_t);
    memset(grid, 0, sizeof(*grid));
}

void sw_fdtd1d_gaussian_right_wave(
    sw_fdtd1d *grid,
    double center_m,
    double width_m,
    double amplitude_v_m
) {
    size_t i;
    for (i = 0U; i < grid->cell_count; ++i) {
        const double x = ((double)i + 0.5) * grid->dx_m;
        double delta = x - center_m;
        if (delta > 0.5 * grid->length_m) {
            delta -= grid->length_m;
        } else if (delta < -0.5 * grid->length_m) {
            delta += grid->length_m;
        }
        grid->electric_y_v_m[i] = amplitude_v_m * exp(-0.5 * delta * delta / (width_m * width_m));
        grid->magnetic_z_t[i] = grid->electric_y_v_m[i] / SW_C;
    }
}

void sw_fdtd1d_step(sw_fdtd1d *grid) {
    size_t i;
    const double b_factor = grid->step_s / grid->dx_m;
    const double e_factor = SW_C * SW_C * grid->step_s / grid->dx_m;
    for (i = 0U; i < grid->cell_count; ++i) {
        const size_t next = (i + 1U) % grid->cell_count;
        grid->magnetic_z_t[i] -= b_factor * (
            grid->electric_y_v_m[next] - grid->electric_y_v_m[i]
        );
    }
    for (i = 0U; i < grid->cell_count; ++i) {
        const size_t previous = (i + grid->cell_count - 1U) % grid->cell_count;
        grid->electric_y_v_m[i] -= e_factor * (
            grid->magnetic_z_t[i] - grid->magnetic_z_t[previous]
        );
    }
}

double sw_fdtd1d_energy_j_m2(const sw_fdtd1d *grid) {
    double energy = 0.0;
    size_t i;
    for (i = 0U; i < grid->cell_count; ++i) {
        energy += (
            0.5 * SW_EPS0 * grid->electric_y_v_m[i] * grid->electric_y_v_m[i]
            + 0.5 / SW_MU0 * grid->magnetic_z_t[i] * grid->magnetic_z_t[i]
        ) * grid->dx_m;
    }
    return energy;
}
