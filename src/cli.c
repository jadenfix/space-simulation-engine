#include "spacewind/constants.h"
#include "spacewind/fields.h"
#include "spacewind/gravity_grid.h"
#include "spacewind/metric.h"
#include "spacewind/mhd1d.h"
#include "spacewind/nbody.h"
#include "spacewind/plasma.h"
#include "spacewind/simulation.h"
#include "spacewind/version.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void sw_usage(FILE *stream) {
    fprintf(stream,
        "spacewind %s: first-principles native space/plasma simulator\n\n"
        "usage:\n"
        "  spacewind simulate <config.cfg> <trajectory.csv> [receipt.json]\n"
        "  spacewind wind <profile.csv>\n"
        "  spacewind geodesic <orbit.csv>\n"
        "  spacewind pic <plasma.csv>\n"
        "  spacewind fdtd <wave.csv>\n"
        "  spacewind mhd <shock.csv>\n"
        "  spacewind boris <particle.csv>\n"
        "  spacewind nbody <orbits.csv>\n"
        "  spacewind poisson <gravity.csv>\n"
        "  spacewind kerr <orbit.csv>\n",
        SPACEWIND_VERSION_STRING
    );
}

static int sw_command_simulate(int argc, char **argv) {
    sw_simulation_config config;
    sw_simulation simulation;
    char error[512];
    FILE *csv;
    FILE *receipt = NULL;
    bool hash_ok = false;
    uint64_t config_hash;

    if (argc < 4 || argc > 5) {
        sw_usage(stderr);
        return 2;
    }
    if (!sw_simulation_config_load_file(argv[2], &config, error, sizeof(error))) {
        fprintf(stderr, "configuration error: %s\n", error);
        return 2;
    }
    config_hash = sw_fnv1a64_file(argv[2], &hash_ok);
    if (!sw_simulation_init(&simulation, &config)) {
        fprintf(stderr, "failed to initialize simulation\n");
        return 1;
    }
    csv = fopen(argv[3], "w");
    if (csv == NULL) {
        fprintf(stderr, "cannot open %s: %s\n", argv[3], strerror(errno));
        sw_simulation_destroy(&simulation);
        return 1;
    }
    if (argc == 5) {
        receipt = fopen(argv[4], "w");
        if (receipt == NULL) {
            fprintf(stderr, "cannot open %s: %s\n", argv[4], strerror(errno));
            fclose(csv);
            sw_simulation_destroy(&simulation);
            return 1;
        }
    }

    if (!sw_simulation_run_csv(&simulation, csv, receipt)) {
        fprintf(stderr, "simulation failed at t=%.17g after %" PRIu64 " accepted steps\n",
            simulation.time_s,
            simulation.accepted_steps
        );
        if (receipt != NULL) {
            fclose(receipt);
        }
        fclose(csv);
        sw_simulation_destroy(&simulation);
        return 1;
    }
    if (receipt != NULL) {
        fprintf(stderr, "config_fnv1a64=%s%016" PRIx64 "\n", hash_ok ? "" : "unavailable:", config_hash);
        fclose(receipt);
    }
    fclose(csv);
    fprintf(stdout,
        "completed: steps=%" PRIu64 " final_r=%.6f AU final_speed=%.6f km/s\n",
        simulation.accepted_steps,
        sw_v3_norm(simulation.state.position_m) / SW_AU,
        sw_v3_norm(simulation.state.velocity_m_s) / 1000.0
    );
    sw_simulation_destroy(&simulation);
    return 0;
}

static int sw_command_wind(int argc, char **argv) {
    sw_solar_wind_model model;
    FILE *out;
    size_t i;
    const size_t samples = 400U;
    if (argc != 3) {
        sw_usage(stderr);
        return 2;
    }
    sw_solar_wind_default(&model);
    sw_solar_wind_init_turbulence(&model, 7U, 0U, 0.0, 1.0, 1.0);
    out = fopen(argv[2], "w");
    if (out == NULL) {
        return 1;
    }
    fprintf(out,
        "radius_m,radius_au,speed_m_s,density_m3,br_t,bphi_t,electric_v_m,"
        "dynamic_pressure_pa,thermal_pressure_pa,magnetic_pressure_pa,debye_m,alfven_m_s\n"
    );
    for (i = 0U; i < samples; ++i) {
        const double u = (double)i / (double)(samples - 1U);
        const double radius = 0.03 * SW_AU * pow(1000.0, u);
        sw_environment_sample sample;
        if (!sw_environment_evaluate(&model, 0.0, sw_v3(radius, 0.0, 0.0), &sample)) {
            continue;
        }
        fprintf(out,
            "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            radius,
            radius / SW_AU,
            sw_v3_norm(sample.wind_velocity_m_s),
            sample.proton_number_density_m3,
            sample.magnetic_field_t.x,
            sample.magnetic_field_t.y,
            sw_v3_norm(sample.electric_field_v_m),
            sample.dynamic_pressure_pa,
            sample.thermal_pressure_pa,
            sample.magnetic_pressure_pa,
            sample.debye_length_m,
            sample.alfven_speed_m_s
        );
    }
    fclose(out);
    return 0;
}

static int sw_command_geodesic(int argc, char **argv) {
    const double mass = 10.0 * SW_SOLAR_MASS;
    const double geometric_mass = SW_G * mass / (SW_C * SW_C);
    const double radius = 20.0 * geometric_mass;
    const double orbital_factor = sqrt(1.0 - 3.0 * geometric_mass / radius);
    const double u_phi = sqrt(geometric_mass / (radius * radius * radius)) / orbital_factor;
    const double affine_period = SW_TWO_PI / u_phi;
    const double step = affine_period / 20000.0;
    sw_schwarzschild_metric parameters = {mass};
    sw_metric metric;
    sw_geodesic_state state;
    FILE *out;
    size_t i;
    const size_t steps = 20000U;

    if (argc != 3) {
        sw_usage(stderr);
        return 2;
    }
    memset(&metric, 0, sizeof(metric));
    metric.eval = sw_metric_schwarzschild;
    metric.ctx = &parameters;
    metric.derivative_step[0] = 100.0;
    metric.derivative_step[1] = fmax(0.01, radius * 1.0e-6);
    metric.derivative_step[2] = 1.0e-6;
    metric.derivative_step[3] = 1.0e-6;
    memset(&state, 0, sizeof(state));
    state.x[0] = 0.0;
    state.x[1] = radius;
    state.x[2] = 0.5 * SW_PI;
    state.x[3] = 0.0;
    state.u[0] = 1.0 / orbital_factor;
    state.u[3] = u_phi;

    out = fopen(argv[2], "w");
    if (out == NULL) {
        return 1;
    }
    fprintf(out, "step,affine_m,ct_m,r_m,theta_rad,phi_rad,x_m,y_m,interval\n");
    for (i = 0U; i <= steps; ++i) {
        if (i % 20U == 0U) {
            const double x = state.x[1] * sin(state.x[2]) * cos(state.x[3]);
            const double y = state.x[1] * sin(state.x[2]) * sin(state.x[3]);
            fprintf(out,
                "%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                i,
                (double)i * step,
                state.x[0],
                state.x[1],
                state.x[2],
                state.x[3],
                x,
                y,
                sw_metric_interval(&metric, state.x, state.u)
            );
        }
        if (i < steps && !sw_geodesic_rk4_step(&metric, &state, step)) {
            fprintf(stderr, "geodesic integration failed at step %zu\n", i);
            fclose(out);
            return 1;
        }
    }
    fclose(out);
    return 0;
}


static int sw_command_kerr(int argc, char **argv) {
    sw_kerr_metric parameters = {10.0 * SW_SOLAR_MASS, 0.8};
    sw_metric metric;
    sw_geodesic_state state;
    FILE *out;
    const double mass_length = SW_G * parameters.mass_kg / (SW_C * SW_C);
    const double radius = 10.0 * mass_length;
    double affine_period;
    double step;
    const size_t steps = 20000U;
    size_t i;

    if (argc != 3) {
        sw_usage(stderr);
        return 2;
    }
    if (!sw_kerr_equatorial_circular_state(&parameters, radius, 1, &state)) {
        fprintf(stderr, "failed to initialize Kerr circular orbit\n");
        return 1;
    }
    memset(&metric, 0, sizeof(metric));
    metric.eval = sw_metric_kerr_boyer_lindquist;
    metric.ctx = &parameters;
    metric.derivative_step[0] = fmax(1.0, radius * 1.0e-7);
    metric.derivative_step[1] = fmax(0.01, radius * 1.0e-6);
    metric.derivative_step[2] = 1.0e-6;
    metric.derivative_step[3] = 1.0e-6;
    affine_period = SW_TWO_PI / fabs(state.u[3]);
    step = affine_period / (double)steps;
    out = fopen(argv[2], "w");
    if (out == NULL) {
        return 1;
    }
    fprintf(out,
        "step,affine_m,ct_m,r_m,theta_rad,phi_rad,x_m,y_m,interval,gtt,gtphi,gphiphi\n"
    );
    for (i = 0U; i <= steps; ++i) {
        if (i % 20U == 0U) {
            double g[4][4];
            const double x = state.x[1] * sin(state.x[2]) * cos(state.x[3]);
            const double y = state.x[1] * sin(state.x[2]) * sin(state.x[3]);
            sw_metric_kerr_boyer_lindquist(state.x, g, &parameters);
            fprintf(out,
                "%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                i,
                (double)i * step,
                state.x[0],
                state.x[1],
                state.x[2],
                state.x[3],
                x,
                y,
                sw_metric_interval(&metric, state.x, state.u),
                g[0][0],
                g[0][3],
                g[3][3]
            );
        }
        if (i < steps && !sw_geodesic_rk4_step(&metric, &state, step)) {
            fprintf(stderr, "Kerr geodesic integration failed at step %zu\n", i);
            fclose(out);
            return 1;
        }
    }
    fclose(out);
    fprintf(stdout,
        "kerr completed: spin=%.3f radius/M=%.3f final_radius/M=%.9f interval=%.12g\n",
        parameters.dimensionless_spin,
        radius / mass_length,
        state.x[1] / mass_length,
        sw_metric_interval(&metric, state.x, state.u)
    );
    return 0;
}

static int sw_command_pic(int argc, char **argv) {
    sw_pic1d pic;
    FILE *out;
    double omega;
    size_t step;
    const size_t steps = 1200U;
    if (argc != 3) {
        sw_usage(stderr);
        return 2;
    }
    if (!sw_pic1d_init(&pic, 64U, 4096U, 1000.0, 1.0e-7)) {
        return 1;
    }
    sw_pic1d_quiet_start(
        &pic,
        -SW_QE,
        SW_ME,
        1.0e6,
        1.0e5,
        0.02,
        1U
    );
    omega = sw_pic1d_plasma_frequency_rad_s(&pic);
    if (omega * pic.step_s > 0.2) {
        fprintf(stderr, "PIC time step is too large for plasma frequency\n");
        sw_pic1d_destroy(&pic);
        return 1;
    }
    out = fopen(argv[2], "w");
    if (out == NULL) {
        sw_pic1d_destroy(&pic);
        return 1;
    }
    fprintf(out, "step,time_s,kinetic_j,field_j,total_j,max_abs_e_v_m,plasma_frequency_rad_s\n");
    for (step = 0U; step <= steps; ++step) {
        if (step % 5U == 0U) {
            size_t i;
            double max_e = 0.0;
            sw_pic1d_deposit_charge(&pic);
            sw_pic1d_solve_poisson_spectral(&pic);
            for (i = 0U; i < pic.grid_count; ++i) {
                max_e = fmax(max_e, fabs(pic.electric_field_v_m[i]));
            }
            {
                const double kinetic = sw_pic1d_kinetic_energy_j(&pic);
                const double field = sw_pic1d_field_energy_j(&pic);
                fprintf(out,
                    "%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    step,
                    (double)step * pic.step_s,
                    kinetic,
                    field,
                    kinetic + field,
                    max_e,
                    omega
                );
            }
        }
        if (step < steps) {
            sw_pic1d_step(&pic);
        }
    }
    fclose(out);
    sw_pic1d_destroy(&pic);
    return 0;
}

static int sw_command_fdtd(int argc, char **argv) {
    sw_fdtd1d grid;
    FILE *out;
    size_t step;
    const size_t steps = 1200U;
    const size_t cells = 512U;
    const double length = 1000.0;
    const double dx = length / (double)cells;
    const double dt = 0.95 * dx / SW_C;
    if (argc != 3) {
        sw_usage(stderr);
        return 2;
    }
    if (!sw_fdtd1d_init(&grid, cells, length, dt)) {
        return 1;
    }
    sw_fdtd1d_gaussian_right_wave(&grid, 250.0, 30.0, 1.0);
    out = fopen(argv[2], "w");
    if (out == NULL) {
        sw_fdtd1d_destroy(&grid);
        return 1;
    }
    fprintf(out, "step,time_s,index,x_m,ey_v_m,bz_t,total_energy_j_m2\n");
    for (step = 0U; step <= steps; ++step) {
        if (step % 20U == 0U) {
            size_t i;
            const double energy = sw_fdtd1d_energy_j_m2(&grid);
            for (i = 0U; i < cells; i += 4U) {
                fprintf(out,
                    "%zu,%.17g,%zu,%.17g,%.17g,%.17g,%.17g\n",
                    step,
                    (double)step * dt,
                    i,
                    ((double)i + 0.5) * dx,
                    grid.electric_y_v_m[i],
                    grid.magnetic_z_t[i],
                    energy
                );
            }
        }
        if (step < steps) {
            sw_fdtd1d_step(&grid);
        }
    }
    fclose(out);
    sw_fdtd1d_destroy(&grid);
    return 0;
}

static int sw_command_mhd(int argc, char **argv) {
    sw_mhd1d grid;
    FILE *out;
    double time = 0.0;
    size_t step = 0U;
    const double end_time = 0.2;
    if (argc != 3) {
        sw_usage(stderr);
        return 2;
    }
    if (!sw_mhd1d_init(&grid, 512U, 1.0, 2.0, 1.0, 0.75)) {
        return 1;
    }
    sw_mhd1d_init_brio_wu(&grid);
    out = fopen(argv[2], "w");
    if (out == NULL) {
        sw_mhd1d_destroy(&grid);
        return 1;
    }
    fprintf(out, "time,step,index,x,density,vx,vy,vz,pressure,by,bz,total_mass,total_energy\n");
    while (time < end_time) {
        double dt = sw_mhd1d_cfl_step(&grid, 0.35);
        if (time + dt > end_time) {
            dt = end_time - time;
        }
        if (!sw_mhd1d_step(&grid, dt)) {
            fprintf(stderr, "MHD step failed at t=%g\n", time);
            fclose(out);
            sw_mhd1d_destroy(&grid);
            return 1;
        }
        time += dt;
        ++step;
        if (step % 20U == 0U || time >= end_time) {
            size_t i;
            const double mass = sw_mhd1d_total_mass(&grid);
            const double energy = sw_mhd1d_total_energy(&grid);
            for (i = 0U; i < grid.cell_count; i += 2U) {
                const sw_mhd_primitive p = sw_mhd1d_get_cell(&grid, i);
                fprintf(out,
                    "%.17g,%zu,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    time,
                    step,
                    i,
                    ((double)i + 0.5) * grid.dx,
                    p.density,
                    p.velocity_x,
                    p.velocity_y,
                    p.velocity_z,
                    p.pressure,
                    p.magnetic_y,
                    p.magnetic_z,
                    mass,
                    energy
                );
            }
        }
    }
    fclose(out);
    sw_mhd1d_destroy(&grid);
    return 0;
}

static int sw_command_boris(int argc, char **argv) {
    sw_particle particle;
    FILE *out;
    size_t step;
    const double magnetic = 5.0e-9;
    const double gyrofrequency = SW_QE * magnetic / SW_MP;
    const double dt = 0.02 / gyrofrequency;
    const size_t steps = 2000U;
    if (argc != 3) {
        sw_usage(stderr);
        return 2;
    }
    memset(&particle, 0, sizeof(particle));
    particle.velocity_m_s = sw_v3(400000.0, 0.0, 0.0);
    particle.charge_c = SW_QE;
    particle.mass_kg = SW_MP;
    particle.macro_weight = 1.0;
    out = fopen(argv[2], "w");
    if (out == NULL) {
        return 1;
    }
    fprintf(out, "step,time_s,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,speed_m_s\n");
    for (step = 0U; step <= steps; ++step) {
        if (step % 5U == 0U) {
            fprintf(out,
                "%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                step,
                (double)step * dt,
                particle.position_m.x,
                particle.position_m.y,
                particle.position_m.z,
                particle.velocity_m_s.x,
                particle.velocity_m_s.y,
                particle.velocity_m_s.z,
                sw_v3_norm(particle.velocity_m_s)
            );
        }
        if (step < steps) {
            sw_boris_push(&particle, sw_v3(0.0, 0.0, 0.0), sw_v3(0.0, 0.0, magnetic), dt);
        }
    }
    fclose(out);
    return 0;
}



static double sw_uniform_sphere_potential(
    double mass_kg,
    double radius_m,
    double distance_m
) {
    if (distance_m >= radius_m) {
        return -SW_G * mass_kg / distance_m;
    }
    return -SW_G * mass_kg * (
        3.0 * radius_m * radius_m - distance_m * distance_m
    ) / (2.0 * radius_m * radius_m * radius_m);
}

static int sw_command_poisson(int argc, char **argv) {
    const size_t cells = 41U;
    const double sphere_radius_m = 2.0e7;
    const double domain_half_width_m = 4.0 * sphere_radius_m;
    const double spacing_m = 2.0 * domain_half_width_m / (double)(cells - 1U);
    const double mass_kg = 1.0e25;
    const double density = mass_kg / (
        (4.0 / 3.0) * SW_PI * sphere_radius_m * sphere_radius_m * sphere_radius_m
    );
    sw_gravity_grid grid;
    FILE *out;
    size_t i;
    size_t j;
    size_t k;
    size_t iterations = 0U;
    double residual = INFINITY;
    double represented_mass = 0.0;
    const sw_vec3 origin = sw_v3(
        -domain_half_width_m,
        -domain_half_width_m,
        -domain_half_width_m
    );

    if (argc != 3) {
        sw_usage(stderr);
        return 2;
    }
    if (!sw_gravity_grid_init(&grid, cells, cells, cells, spacing_m, origin)) {
        return 1;
    }
    for (k = 0U; k < cells; ++k) {
        for (j = 0U; j < cells; ++j) {
            for (i = 0U; i < cells; ++i) {
                const sw_vec3 position = sw_gravity_grid_position(&grid, i, j, k);
                const size_t index = sw_gravity_grid_index(&grid, i, j, k);
                if (sw_v3_norm(position) <= sphere_radius_m) {
                    grid.density_kg_m3[index] = density;
                    represented_mass += density * spacing_m * spacing_m * spacing_m;
                }
            }
        }
    }
    sw_gravity_grid_set_monopole_boundary(
        &grid, represented_mass, sw_v3(0.0, 0.0, 0.0), sphere_radius_m
    );
    if (!sw_gravity_grid_solve_sor(&grid, 1.85, 1.0e-8, 20000U, &iterations, &residual)) {
        fprintf(stderr, "Poisson solve did not converge: iterations=%zu residual=%g\n", iterations, residual);
        sw_gravity_grid_destroy(&grid);
        return 1;
    }
    out = fopen(argv[2], "w");
    if (out == NULL) {
        sw_gravity_grid_destroy(&grid);
        return 1;
    }
    fprintf(out,
        "i,j,x_m,y_m,r_m,density_kg_m3,potential_numeric_m2_s2,potential_exact_m2_s2,"
        "relative_error,ax_m_s2,ay_m_s2,g00,gxx,iterations,residual,represented_mass_kg\n"
    );
    k = cells / 2U;
    for (j = 1U; j + 1U < cells; ++j) {
        for (i = 1U; i + 1U < cells; ++i) {
            const sw_vec3 position = sw_gravity_grid_position(&grid, i, j, k);
            const size_t index = sw_gravity_grid_index(&grid, i, j, k);
            const double radius = sw_v3_norm(position);
            const double exact = sw_uniform_sphere_potential(
                represented_mass, sphere_radius_m, radius
            );
            const double numeric = grid.potential_m2_s2[index];
            sw_vec3 acceleration = sw_v3(0.0, 0.0, 0.0);
            double metric[4][4];
            (void)sw_gravity_grid_sample_acceleration(&grid, position, &acceleration);
            (void)sw_gravity_grid_weak_field_metric(&grid, position, metric);
            fprintf(out,
                "%zu,%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                "%.17g,%.17g,%zu,%.17g,%.17g\n",
                i,
                j,
                position.x,
                position.y,
                radius,
                grid.density_kg_m3[index],
                numeric,
                exact,
                (numeric - exact) / fabs(exact),
                acceleration.x,
                acceleration.y,
                metric[0][0],
                metric[1][1],
                iterations,
                residual,
                represented_mass
            );
        }
    }
    fclose(out);
    fprintf(stdout,
        "poisson completed: grid=%zux%zux%zu iterations=%zu residual=%.3e represented_mass=%.6e kg\n",
        cells, cells, cells, iterations, residual, represented_mass
    );
    sw_gravity_grid_destroy(&grid);
    return 0;
}

static const char *sw_nbody_name(unsigned id) {
    switch (id) {
        case 0U:
            return "sun";
        case 1U:
            return "earth";
        case 2U:
            return "jupiter";
        default:
            return "body";
    }
}

static void sw_nbody_initialize_demo(sw_nbody_system *system) {
    const double earth_radius = SW_AU;
    const double jupiter_radius = 5.2044 * SW_AU;
    const double earth_speed = sqrt(SW_G * SW_SOLAR_MASS / earth_radius);
    const double jupiter_speed = sqrt(SW_G * SW_SOLAR_MASS / jupiter_radius);
    sw_vec3 center;
    sw_vec3 momentum;
    size_t i;

    memset(system, 0, sizeof(*system));
    system->count = 3U;
    system->softening_m = 1000.0;

    system->objects[0].id = 0U;
    system->objects[0].mass_kg = SW_SOLAR_MASS;
    system->objects[0].radius_m = SW_SOLAR_RADIUS;

    system->objects[1].id = 1U;
    system->objects[1].mass_kg = SW_EARTH_MASS;
    system->objects[1].radius_m = SW_EARTH_RADIUS;
    system->objects[1].position_m = sw_v3(earth_radius, 0.0, 0.0);
    system->objects[1].velocity_m_s = sw_v3(0.0, earth_speed, 0.0);

    system->objects[2].id = 2U;
    system->objects[2].mass_kg = SW_JUPITER_MASS;
    system->objects[2].radius_m = SW_JUPITER_RADIUS;
    system->objects[2].position_m = sw_v3(0.0, jupiter_radius, 0.0);
    system->objects[2].velocity_m_s = sw_v3(-jupiter_speed, 0.0, 0.0);

    momentum = sw_nbody_total_momentum_kg_m_s(system);
    system->objects[0].velocity_m_s = sw_v3_scale(momentum, -1.0 / SW_SOLAR_MASS);
    center = sw_nbody_center_of_mass_m(system);
    for (i = 0U; i < system->count; ++i) {
        system->objects[i].position_m = sw_v3_sub(system->objects[i].position_m, center);
    }
}

static int sw_command_nbody(int argc, char **argv) {
    sw_nbody_system system;
    FILE *out;
    const double step_s = 6.0 * 3600.0;
    const double duration_s = 12.0 * SW_YEAR;
    const size_t steps = (size_t)llround(duration_s / step_s);
    const size_t output_stride = 20U;
    double reference_energy;
    size_t step;

    if (argc != 3) {
        sw_usage(stderr);
        return 2;
    }
    sw_nbody_initialize_demo(&system);
    reference_energy = sw_nbody_total_energy_j(&system);
    out = fopen(argv[2], "w");
    if (out == NULL) {
        fprintf(stderr, "cannot open %s: %s\n", argv[2], strerror(errno));
        return 1;
    }
    fprintf(out,
        "step,time_s,time_years,id,name,mass_kg,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,"
        "system_energy_j,relative_energy_error,momentum_x,momentum_y,momentum_z,com_x,com_y,com_z\n"
    );
    for (step = 0U; step <= steps; ++step) {
        if (step % output_stride == 0U || step == steps) {
            const double energy = sw_nbody_total_energy_j(&system);
            const double relative_error = (energy - reference_energy) / fabs(reference_energy);
            const sw_vec3 momentum = sw_nbody_total_momentum_kg_m_s(&system);
            const sw_vec3 center = sw_nbody_center_of_mass_m(&system);
            size_t i;
            for (i = 0U; i < system.count; ++i) {
                const sw_nbody_object *object = &system.objects[i];
                fprintf(out,
                    "%zu,%.17g,%.17g,%u,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                    "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                    step,
                    (double)step * step_s,
                    (double)step * step_s / SW_YEAR,
                    object->id,
                    sw_nbody_name(object->id),
                    object->mass_kg,
                    object->position_m.x,
                    object->position_m.y,
                    object->position_m.z,
                    object->velocity_m_s.x,
                    object->velocity_m_s.y,
                    object->velocity_m_s.z,
                    energy,
                    relative_error,
                    momentum.x,
                    momentum.y,
                    momentum.z,
                    center.x,
                    center.y,
                    center.z
                );
            }
        }
        if (step < steps && !sw_nbody_leapfrog_step(&system, step_s)) {
            fprintf(stderr, "N-body integration failed at step %zu\n", step);
            fclose(out);
            return 1;
        }
    }
    fclose(out);
    fprintf(stdout,
        "nbody completed: bodies=%zu years=%.3f relative_energy_error=%.6e\n",
        system.count,
        duration_s / SW_YEAR,
        (sw_nbody_total_energy_j(&system) - reference_energy) / fabs(reference_energy)
    );
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        sw_usage(stderr);
        return 2;
    }
    if (strcmp(argv[1], "simulate") == 0) {
        return sw_command_simulate(argc, argv);
    }
    if (strcmp(argv[1], "wind") == 0) {
        return sw_command_wind(argc, argv);
    }
    if (strcmp(argv[1], "geodesic") == 0) {
        return sw_command_geodesic(argc, argv);
    }
    if (strcmp(argv[1], "kerr") == 0) {
        return sw_command_kerr(argc, argv);
    }
    if (strcmp(argv[1], "pic") == 0) {
        return sw_command_pic(argc, argv);
    }
    if (strcmp(argv[1], "fdtd") == 0) {
        return sw_command_fdtd(argc, argv);
    }
    if (strcmp(argv[1], "mhd") == 0) {
        return sw_command_mhd(argc, argv);
    }
    if (strcmp(argv[1], "boris") == 0) {
        return sw_command_boris(argc, argv);
    }
    if (strcmp(argv[1], "nbody") == 0) {
        return sw_command_nbody(argc, argv);
    }
    if (strcmp(argv[1], "poisson") == 0) {
        return sw_command_poisson(argc, argv);
    }
    sw_usage(stderr);
    return 2;
}
