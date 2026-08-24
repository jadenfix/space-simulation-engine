#include "spacewind/pde_checks.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double max2(double a, double b) { return a > b ? a : b; }

int swa_poisson_manufactured_error(size_t intervals,
                                   swa_discretization_error *out) {
    const double kx = SWA_PI;
    const double ky = 2.0 * SWA_PI;
    double h;
    swa_kahan squared = {0.0, 0.0};
    double maximum = 0.0;
    size_t count = 0U;
    size_t i, j;
    if (out == NULL || intervals < 4U) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    h = 1.0 / (double)intervals;
    for (i = 1U; i < intervals; ++i) {
        const double x = (double)i * h;
        for (j = 1U; j < intervals; ++j) {
            const double y = (double)j * h;
            const double phi = sin(kx * x) * sin(ky * y);
            const double xp = sin(kx * (x + h)) * sin(ky * y);
            const double xm = sin(kx * (x - h)) * sin(ky * y);
            const double yp = sin(kx * x) * sin(ky * (y + h));
            const double ym = sin(kx * x) * sin(ky * (y - h));
            const double laplacian = (xp + xm + yp + ym - 4.0 * phi) / (h * h);
            const double source = -(kx * kx + ky * ky) * phi;
            const double error = laplacian - source;
            swa_kahan_add(&squared, error * error);
            if (fabs(error) > maximum) {
                maximum = fabs(error);
            }
            ++count;
        }
    }
    out->resolution = intervals;
    out->spacing = h;
    out->l2_error = sqrt(squared.sum / (double)count);
    out->linf_error = maximum;
    out->finite = isfinite(out->l2_error) && isfinite(out->linf_error);
    return 1;
}

int swa_continuity_manufactured_error(size_t intervals,
                                      double cfl,
                                      swa_discretization_error *out) {
    const double velocity = 0.3;
    const double amplitude = 0.2;
    const double background = 1.0;
    const double wave_number = 2.0 * SWA_PI;
    const double time = 0.371;
    double dx, dt;
    swa_kahan squared = {0.0, 0.0};
    double maximum = 0.0;
    size_t i;
    if (out == NULL || intervals < 8U || !(cfl > 0.0 && cfl <= 1.0)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    dx = 1.0 / (double)intervals;
    dt = cfl * dx / velocity;
    for (i = 0U; i < intervals; ++i) {
        const double x = (double)i * dx;
        const double rho_tp = background + amplitude * cos(wave_number * (x - velocity * (time + dt)));
        const double rho_tm = background + amplitude * cos(wave_number * (x - velocity * (time - dt)));
        const double rho_xp = background + amplitude * cos(wave_number * ((x + dx) - velocity * time));
        const double rho_xm = background + amplitude * cos(wave_number * ((x - dx) - velocity * time));
        const double time_derivative = (rho_tp - rho_tm) / (2.0 * dt);
        const double current_derivative = velocity * (rho_xp - rho_xm) / (2.0 * dx);
        const double residual = time_derivative + current_derivative;
        swa_kahan_add(&squared, residual * residual);
        if (fabs(residual) > maximum) {
            maximum = fabs(residual);
        }
    }
    out->resolution = intervals;
    out->spacing = dx;
    out->l2_error = sqrt(squared.sum / (double)intervals);
    out->linf_error = maximum;
    out->finite = isfinite(out->l2_error) && isfinite(out->linf_error);
    return 1;
}

int swa_vacuum_maxwell_period(size_t cells,
                              double cfl,
                              swa_maxwell_wave_result *out) {
    const double length = 1.0;
    const double amplitude = 1.0;
    const double wave_number = 2.0 * SWA_PI / length;
    double *electric;
    double *magnetic;
    double *initial;
    double dx, dt, target_steps, time;
    size_t steps, step, i;
    swa_kahan initial_squared = {0.0, 0.0};
    swa_kahan final_squared = {0.0, 0.0};
    swa_kahan error_squared = {0.0, 0.0};
    double maximum_error = 0.0;
    if (out == NULL || cells < 16U || !(cfl > 0.0 && cfl <= 1.0)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    target_steps = (double)cells / cfl;
    steps = (size_t)llround(target_steps);
    if (fabs((double)steps - target_steps) > 1e-12 * target_steps) {
        return 0;
    }
    electric = (double *)calloc(cells, sizeof(*electric));
    magnetic = (double *)calloc(cells, sizeof(*magnetic));
    initial = (double *)calloc(cells, sizeof(*initial));
    if (electric == NULL || magnetic == NULL || initial == NULL) {
        free(electric);
        free(magnetic);
        free(initial);
        return 0;
    }
    dx = length / (double)cells;
    dt = cfl * dx / SWA_C;
    for (i = 0U; i < cells; ++i) {
        const double x_e = (double)i * dx;
        const double x_b = ((double)i + 0.5) * dx;
        electric[i] = amplitude * sin(wave_number * x_e);
        initial[i] = electric[i];
        magnetic[i] = -amplitude / SWA_C *
                      sin(wave_number * (x_b + 0.5 * SWA_C * dt));
        swa_kahan_add(&initial_squared, electric[i] * electric[i]);
    }
    for (step = 0U; step < steps; ++step) {
        for (i = 0U; i < cells; ++i) {
            const size_t ip = (i + 1U) % cells;
            magnetic[i] += (dt / dx) * (electric[ip] - electric[i]);
        }
        for (i = 0U; i < cells; ++i) {
            const size_t im = (i + cells - 1U) % cells;
            electric[i] += SWA_C * SWA_C * (dt / dx) *
                           (magnetic[i] - magnetic[im]);
        }
    }
    time = (double)steps * dt;
    for (i = 0U; i < cells; ++i) {
        const double x = (double)i * dx;
        const double exact = amplitude * sin(wave_number * (x - SWA_C * time));
        const double error = electric[i] - exact;
        swa_kahan_add(&final_squared, electric[i] * electric[i]);
        swa_kahan_add(&error_squared, error * error);
        if (fabs(error) > maximum_error) {
            maximum_error = fabs(error);
        }
    }
    out->resolution = cells;
    out->steps = steps;
    out->cfl = cfl;
    out->l2_electric_error_V_m = sqrt(error_squared.sum / (double)cells);
    out->linf_electric_error_V_m = maximum_error;
    out->initial_electric_rms_V_m = sqrt(initial_squared.sum / (double)cells);
    out->final_electric_rms_V_m = sqrt(final_squared.sum / (double)cells);
    out->rms_amplitude_relative_error =
        fabs(out->final_electric_rms_V_m - out->initial_electric_rms_V_m) /
        max2(out->initial_electric_rms_V_m, DBL_MIN);
    out->finite = isfinite(out->l2_electric_error_V_m) &&
                  isfinite(out->linf_electric_error_V_m) &&
                  isfinite(out->rms_amplitude_relative_error);
    out->stable = out->finite && out->final_electric_rms_V_m <
                  2.0 * out->initial_electric_rms_V_m;
    free(electric);
    free(magnetic);
    free(initial);
    return 1;
}

int swa_uniform_boris_gyroperiod(size_t steps_per_period,
                                 double charge_to_mass_C_kg,
                                 double magnetic_field_T,
                                 double initial_speed_mps,
                                 swa_boris_orbit_result *out) {
    const double omega = fabs(charge_to_mass_C_kg * magnetic_field_T);
    double period, dt, t_scalar, s_scalar;
    swa_vec3 velocity;
    swa_vec3 position = {0.0, 0.0, 0.0};
    const swa_vec3 initial_velocity = {1.0, 0.0, 0.0};
    size_t step;
    if (out == NULL || steps_per_period < 8U || !(omega > 0.0) ||
        !(initial_speed_mps > 0.0)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    period = 2.0 * SWA_PI / omega;
    dt = period / (double)steps_per_period;
    t_scalar = charge_to_mass_C_kg * magnetic_field_T * dt * 0.5;
    s_scalar = 2.0 * t_scalar / (1.0 + t_scalar * t_scalar);
    velocity = swa_vscale(initial_velocity, initial_speed_mps);
    for (step = 0U; step < steps_per_period; ++step) {
        const swa_vec3 old_velocity = velocity;
        const swa_vec3 t = swa_v3(0.0, 0.0, t_scalar);
        const swa_vec3 s = swa_v3(0.0, 0.0, s_scalar);
        const swa_vec3 prime = swa_vadd(velocity, swa_vcross(velocity, t));
        velocity = swa_vadd(velocity, swa_vcross(prime, s));
        position = swa_vadd(position,
                            swa_vscale(swa_vadd(old_velocity, velocity), 0.5 * dt));
    }
    out->steps_per_gyroperiod = steps_per_period;
    out->speed_initial_mps = initial_speed_mps;
    out->speed_final_mps = swa_vnorm(velocity);
    out->speed_relative_error = fabs(out->speed_final_mps - initial_speed_mps) /
                                initial_speed_mps;
    out->velocity_closure_error_mps =
        swa_vnorm(swa_vsub(velocity, swa_vscale(initial_velocity, initial_speed_mps)));
    out->position_closure_error_m = swa_vnorm(position);
    out->finite = isfinite(out->speed_relative_error) &&
                  isfinite(out->velocity_closure_error_mps) &&
                  isfinite(out->position_closure_error_m);
    return 1;
}

int swa_parker_implicit_residual(double radius_over_critical,
                                 double speed_over_sound,
                                 double *residual) {
    double y;
    if (residual == NULL || !(radius_over_critical > 0.0) ||
        !(speed_over_sound > 0.0)) {
        return 0;
    }
    y = speed_over_sound * speed_over_sound;
    *residual = y - log(y) -
                (4.0 * log(radius_over_critical) +
                 4.0 / radius_over_critical - 3.0);
    return isfinite(*residual);
}
