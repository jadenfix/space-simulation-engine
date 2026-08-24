#include "spacewind/constants.h"
#include "spacewind/fields.h"
#include "spacewind/gravity_grid.h"
#include "spacewind/integrator.h"
#include "spacewind/math3.h"
#include "spacewind/metric.h"
#include "spacewind/mhd1d.h"
#include "spacewind/nbody.h"
#include "spacewind/plasma.h"
#include "spacewind/simulation.h"
#include "spacewind/spacecraft.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static unsigned sw_tests_run = 0U;
static unsigned sw_tests_failed = 0U;

#define SW_CHECK(condition) do { \
    ++sw_tests_run; \
    if (!(condition)) { \
        ++sw_tests_failed; \
        fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #condition); \
    } \
} while (0)

static bool sw_near(double a, double b, double relative, double absolute) {
    return fabs(a - b) <= absolute + relative * fmax(fabs(a), fabs(b));
}

static void test_vectors_and_matrix(void) {
    const sw_vec3 a = sw_v3(1.0, 2.0, 3.0);
    const sw_vec3 b = sw_v3(-4.0, 5.0, -6.0);
    const sw_vec3 c = sw_v3_cross(a, b);
    const double m[4][4] = {
        {4.0, 1.0, 0.0, 0.0},
        {1.0, 3.0, 0.0, 0.0},
        {0.0, 0.0, 2.0, 0.5},
        {0.0, 0.0, 0.5, 1.0}
    };
    double inv[4][4];
    double product[4][4];
    size_t i;
    size_t j;
    SW_CHECK(sw_near(sw_v3_dot(c, a), 0.0, 0.0, 1.0e-14));
    SW_CHECK(sw_near(sw_v3_dot(c, b), 0.0, 0.0, 1.0e-14));
    SW_CHECK(sw_near(sw_v3_norm(sw_v3_normalize(a)), 1.0, 1.0e-14, 1.0e-14));
    SW_CHECK(sw_mat4_inverse(&m[0][0], &inv[0][0]));
    sw_mat4_mul(&m[0][0], &inv[0][0], &product[0][0]);
    for (i = 0U; i < 4U; ++i) {
        for (j = 0U; j < 4U; ++j) {
            SW_CHECK(sw_near(product[i][j], i == j ? 1.0 : 0.0, 1.0e-12, 1.0e-12));
        }
    }
}

static void test_metric(void) {
    sw_metric minkowski;
    double gamma[4][4][4];
    double x[4] = {0.0, 10.0, -2.0, 3.0};
    size_t mu;
    size_t a;
    size_t b;
    sw_geodesic_state state;
    memset(&minkowski, 0, sizeof(minkowski));
    minkowski.eval = sw_metric_minkowski_cartesian;
    SW_CHECK(sw_metric_christoffel(&minkowski, x, gamma));
    for (mu = 0U; mu < 4U; ++mu) {
        for (a = 0U; a < 4U; ++a) {
            for (b = 0U; b < 4U; ++b) {
                SW_CHECK(fabs(gamma[mu][a][b]) < 1.0e-10);
            }
        }
    }
    memset(&state, 0, sizeof(state));
    state.x[1] = 1.0;
    state.u[0] = 1.0;
    state.u[1] = 0.1;
    SW_CHECK(sw_geodesic_rk4_step(&minkowski, &state, 10.0));
    SW_CHECK(sw_near(state.x[0], 10.0, 1.0e-12, 1.0e-12));
    SW_CHECK(sw_near(state.x[1], 2.0, 1.0e-12, 1.0e-12));
    SW_CHECK(sw_near(sw_metric_interval(&minkowski, state.x, state.u), -0.99, 1.0e-12, 1.0e-12));
}


static void test_additional_metrics(void) {
    sw_schwarzschild_metric schwarzschild_params = {10.0 * SW_SOLAR_MASS};
    sw_kerr_metric kerr_params = {10.0 * SW_SOLAR_MASS, 0.0};
    sw_flrw_metric flrw_params = {1.0, 0.0};
    const double x_spherical[4] = {0.0, 1.0e9, 1.1, 0.3};
    const double x_cartesian[4] = {100.0, 1.0, 2.0, 3.0};
    double schwarzschild[4][4];
    double kerr[4][4];
    double flrw[4][4];
    double minkowski[4][4];
    size_t i;
    size_t j;

    sw_metric_schwarzschild(x_spherical, schwarzschild, &schwarzschild_params);
    sw_metric_kerr_boyer_lindquist(x_spherical, kerr, &kerr_params);
    sw_metric_flrw_de_sitter(x_cartesian, flrw, &flrw_params);
    sw_metric_minkowski_cartesian(x_cartesian, minkowski, NULL);
    for (i = 0U; i < 4U; ++i) {
        for (j = 0U; j < 4U; ++j) {
            SW_CHECK(sw_near(kerr[i][j], schwarzschild[i][j], 1.0e-13, 1.0e-8));
            SW_CHECK(sw_near(flrw[i][j], minkowski[i][j], 1.0e-14, 1.0e-14));
        }
    }
    {
        sw_kerr_metric spinning = {10.0 * SW_SOLAR_MASS, 0.8};
        sw_metric metric;
        sw_geodesic_state circular;
        sw_geodesic_state derivative;
        const double mass_length = SW_G * spinning.mass_kg / (SW_C * SW_C);
        memset(&metric, 0, sizeof(metric));
        metric.eval = sw_metric_kerr_boyer_lindquist;
        metric.ctx = &spinning;
        metric.derivative_step[0] = 1.0;
        metric.derivative_step[1] = 1.0e-5 * mass_length;
        metric.derivative_step[2] = 1.0e-6;
        metric.derivative_step[3] = 1.0e-6;
        SW_CHECK(sw_kerr_equatorial_circular_state(
            &spinning, 10.0 * mass_length, 1, &circular
        ));
        SW_CHECK(sw_near(
            sw_metric_interval(&metric, circular.x, circular.u), -1.0, 1.0e-12, 1.0e-12
        ));
        SW_CHECK(sw_geodesic_rhs(&metric, &circular, &derivative));
        SW_CHECK(fabs(derivative.u[1]) < 1.0e-10);
    }
}


static void test_gravity_grid(void) {
    const size_t nodes = 17U;
    const double length = 1.0e9;
    const double spacing = length / (double)(nodes - 1U);
    const double amplitude = 1.0e6;
    const double wave_number = SW_PI / length;
    sw_gravity_grid grid;
    size_t iterations = 0U;
    double residual = INFINITY;
    double error2 = 0.0;
    double exact2 = 0.0;
    size_t i;
    size_t j;
    size_t k;

    SW_CHECK(sw_gravity_grid_init(
        &grid, nodes, nodes, nodes, spacing, sw_v3(0.0, 0.0, 0.0)
    ));
    for (k = 0U; k < nodes; ++k) {
        for (j = 0U; j < nodes; ++j) {
            for (i = 0U; i < nodes; ++i) {
                const sw_vec3 position = sw_gravity_grid_position(&grid, i, j, k);
                const double exact = -amplitude
                    * sin(wave_number * position.x)
                    * sin(wave_number * position.y)
                    * sin(wave_number * position.z);
                const size_t index = sw_gravity_grid_index(&grid, i, j, k);
                grid.density_kg_m3[index] = -3.0 * wave_number * wave_number * exact
                    / (4.0 * SW_PI * SW_G);
            }
        }
    }
    SW_CHECK(sw_gravity_grid_solve_sor(
        &grid, 1.85, 1.0e-9, 10000U, &iterations, &residual
    ));
    SW_CHECK(iterations < 10000U);
    SW_CHECK(residual < 1.0e-9);
    for (k = 1U; k + 1U < nodes; ++k) {
        for (j = 1U; j + 1U < nodes; ++j) {
            for (i = 1U; i + 1U < nodes; ++i) {
                const sw_vec3 position = sw_gravity_grid_position(&grid, i, j, k);
                const double exact = -amplitude
                    * sin(wave_number * position.x)
                    * sin(wave_number * position.y)
                    * sin(wave_number * position.z);
                const double numeric = grid.potential_m2_s2[
                    sw_gravity_grid_index(&grid, i, j, k)
                ];
                const double error = numeric - exact;
                error2 += error * error;
                exact2 += exact * exact;
            }
        }
    }
    SW_CHECK(sqrt(error2 / exact2) < 0.004);
    {
        const sw_vec3 point = sw_v3(0.5 * length, 0.5 * length, 0.5 * length);
        double sampled = 0.0;
        sw_vec3 acceleration;
        double metric[4][4];
        SW_CHECK(sw_gravity_grid_sample_potential(&grid, point, &sampled));
        SW_CHECK(sw_near(sampled, -amplitude, 0.004, 1.0));
        SW_CHECK(sw_gravity_grid_sample_acceleration(&grid, point, &acceleration));
        SW_CHECK(sw_v3_norm(acceleration) < 1.0e-10);
        SW_CHECK(sw_gravity_grid_weak_field_metric(&grid, point, metric));
        SW_CHECK(metric[0][0] > -1.0);
        SW_CHECK(metric[1][1] > 1.0);
    }
    sw_gravity_grid_destroy(&grid);
}

static void test_parker_environment(void) {
    sw_solar_wind_model model;
    sw_environment_sample near_sample;
    sw_environment_sample one_au;
    sw_environment_sample far_sample;
    sw_solar_wind_default(&model);
    model.stream_fraction = 0.0;
    model.turbulence_mode_count = 0U;
    SW_CHECK(sw_environment_evaluate(&model, 0.0, sw_v3(0.1 * SW_AU, 0.0, 0.0), &near_sample));
    SW_CHECK(sw_environment_evaluate(&model, 0.0, sw_v3(SW_AU, 0.0, 0.0), &one_au));
    SW_CHECK(sw_environment_evaluate(&model, 0.0, sw_v3(5.0 * SW_AU, 0.0, 0.0), &far_sample));
    SW_CHECK(sw_v3_norm(near_sample.wind_velocity_m_s) < sw_v3_norm(one_au.wind_velocity_m_s));
    SW_CHECK(sw_v3_norm(one_au.wind_velocity_m_s) < sw_v3_norm(far_sample.wind_velocity_m_s));
    SW_CHECK(one_au.proton_number_density_m3 > far_sample.proton_number_density_m3);
    SW_CHECK(sw_near(one_au.proton_number_density_m3, 5.0e6, 1.0e-10, 1.0e-3));
    SW_CHECK(one_au.debye_length_m > 0.0);
    SW_CHECK(sw_v3_norm(one_au.magnetic_field_t) > 0.0);
    SW_CHECK(sw_v3_norm(one_au.electric_field_v_m) > 0.0);
}

typedef struct {
    double mu;
} orbit_context;

static sw_vec3 orbit_acceleration(double time, sw_vec3 position, void *context) {
    const orbit_context *orbit = (const orbit_context *)context;
    const double r = sw_v3_norm(position);
    (void)time;
    return sw_v3_scale(position, -orbit->mu / (r * r * r));
}

static void test_verlet_orbit(void) {
    orbit_context context = {SW_SOLAR_MU};
    sw_vec3 position = sw_v3(SW_AU, 0.0, 0.0);
    sw_vec3 velocity = sw_v3(0.0, sqrt(SW_SOLAR_MU / SW_AU), 0.0);
    const double initial_energy = sw_specific_orbital_energy(context.mu, position, velocity);
    const double dt = 900.0;
    size_t i;
    for (i = 0U; i < 35064U; ++i) {
        sw_velocity_verlet_step(orbit_acceleration, &context, (double)i * dt, dt, &position, &velocity);
    }
    {
        const double final_energy = sw_specific_orbital_energy(context.mu, position, velocity);
        SW_CHECK(fabs((final_energy - initial_energy) / initial_energy) < 1.0e-7);
        SW_CHECK(fabs(sw_v3_norm(position) / SW_AU - 1.0) < 2.0e-4);
    }
}


static void test_nbody(void) {
    sw_nbody_system system;
    const double radius = SW_AU;
    const double total_mass = SW_SOLAR_MASS + SW_EARTH_MASS;
    const double relative_speed = sqrt(SW_G * total_mass / radius);
    const double sun_radius_from_com = radius * SW_EARTH_MASS / total_mass;
    const double earth_radius_from_com = radius * SW_SOLAR_MASS / total_mass;
    const double step_s = 6.0 * 3600.0;
    const size_t steps = (size_t)llround(SW_YEAR / step_s);
    double initial_energy;
    sw_vec3 initial_momentum;
    size_t i;

    memset(&system, 0, sizeof(system));
    system.count = 2U;
    system.softening_m = 0.0;
    system.objects[0].id = 0U;
    system.objects[0].mass_kg = SW_SOLAR_MASS;
    system.objects[0].radius_m = SW_SOLAR_RADIUS;
    system.objects[0].position_m = sw_v3(-sun_radius_from_com, 0.0, 0.0);
    system.objects[0].velocity_m_s = sw_v3(
        0.0,
        -relative_speed * SW_EARTH_MASS / total_mass,
        0.0
    );
    system.objects[1].id = 1U;
    system.objects[1].mass_kg = SW_EARTH_MASS;
    system.objects[1].radius_m = SW_EARTH_RADIUS;
    system.objects[1].position_m = sw_v3(earth_radius_from_com, 0.0, 0.0);
    system.objects[1].velocity_m_s = sw_v3(
        0.0,
        relative_speed * SW_SOLAR_MASS / total_mass,
        0.0
    );

    initial_energy = sw_nbody_total_energy_j(&system);
    initial_momentum = sw_nbody_total_momentum_kg_m_s(&system);
    SW_CHECK(sw_v3_norm(initial_momentum) < 1.0e15);
    for (i = 0U; i < steps; ++i) {
        SW_CHECK(sw_nbody_leapfrog_step(&system, step_s));
    }
    SW_CHECK(fabs((sw_nbody_total_energy_j(&system) - initial_energy) / initial_energy) < 1.0e-8);
    SW_CHECK(sw_v3_norm(sw_nbody_total_momentum_kg_m_s(&system)) < 1.0e16);
    SW_CHECK(sw_v3_norm(sw_nbody_center_of_mass_m(&system)) < 1.0e-4);
    SW_CHECK(fabs(sw_v3_norm(sw_v3_sub(
        system.objects[1].position_m,
        system.objects[0].position_m
    )) / SW_AU - 1.0) < 5.0e-5);
}

static void test_boris(void) {
    sw_particle particle;
    const double magnetic = 5.0e-9;
    const double omega = SW_QE * magnetic / SW_MP;
    const double dt = 0.02 / omega;
    const double initial_speed = 400000.0;
    size_t i;
    memset(&particle, 0, sizeof(particle));
    particle.velocity_m_s = sw_v3(initial_speed, 0.0, 0.0);
    particle.charge_c = SW_QE;
    particle.mass_kg = SW_MP;
    for (i = 0U; i < 10000U; ++i) {
        sw_boris_push(&particle, sw_v3(0.0, 0.0, 0.0), sw_v3(0.0, 0.0, magnetic), dt);
    }
    SW_CHECK(sw_near(sw_v3_norm(particle.velocity_m_s), initial_speed, 1.0e-12, 1.0e-6));

    particle.velocity_m_s = sw_v3(0.99 * SW_C, 0.0, 0.0);
    for (i = 0U; i < 100U; ++i) {
        sw_relativistic_boris_push(
            &particle,
            sw_v3(0.0, 0.0, 0.0),
            sw_v3(0.0, 0.0, magnetic),
            dt
        );
    }
    SW_CHECK(sw_v3_norm(particle.velocity_m_s) < SW_C);
    SW_CHECK(sw_near(sw_v3_norm(particle.velocity_m_s), 0.99 * SW_C, 1.0e-12, 1.0e-3));
}

static void test_pic(void) {
    sw_pic1d pic;
    size_t i;
    SW_CHECK(sw_pic1d_init(&pic, 32U, 512U, 1000.0, 1.0e-7));
    sw_pic1d_quiet_start(&pic, -SW_QE, SW_ME, 1.0e6, 1.0e5, 0.01, 1U);
    sw_pic1d_deposit_charge(&pic);
    sw_pic1d_solve_poisson_spectral(&pic);
    SW_CHECK(sw_pic1d_plasma_frequency_rad_s(&pic) > 0.0);
    for (i = 0U; i < 100U; ++i) {
        sw_pic1d_step(&pic);
    }
    SW_CHECK(isfinite(sw_pic1d_kinetic_energy_j(&pic)));
    SW_CHECK(isfinite(sw_pic1d_field_energy_j(&pic)));
    sw_pic1d_destroy(&pic);
}

static void test_fdtd(void) {
    sw_fdtd1d grid;
    const size_t cells = 256U;
    const double length = 1000.0;
    const double dt = 0.5 * (length / (double)cells) / SW_C;
    double initial;
    size_t i;
    SW_CHECK(sw_fdtd1d_init(&grid, cells, length, dt));
    sw_fdtd1d_gaussian_right_wave(&grid, 250.0, 30.0, 1.0);
    initial = sw_fdtd1d_energy_j_m2(&grid);
    for (i = 0U; i < 500U; ++i) {
        sw_fdtd1d_step(&grid);
    }
    {
        const double final = sw_fdtd1d_energy_j_m2(&grid);
        SW_CHECK(fabs((final - initial) / initial) < 0.03);
    }
    sw_fdtd1d_destroy(&grid);
}

static void test_mhd(void) {
    sw_mhd1d grid;
    double initial_mass;
    double time = 0.0;
    size_t i;
    SW_CHECK(sw_mhd1d_init(&grid, 256U, 1.0, 2.0, 1.0, 0.75));
    sw_mhd1d_init_brio_wu(&grid);
    initial_mass = sw_mhd1d_total_mass(&grid);
    while (time < 0.05) {
        double dt = sw_mhd1d_cfl_step(&grid, 0.3);
        if (time + dt > 0.05) {
            dt = 0.05 - time;
        }
        SW_CHECK(sw_mhd1d_step(&grid, dt));
        time += dt;
    }
    for (i = 0U; i < grid.cell_count; ++i) {
        const sw_mhd_primitive p = sw_mhd1d_get_cell(&grid, i);
        SW_CHECK(p.density > 0.0);
        SW_CHECK(p.pressure > 0.0);
    }
    SW_CHECK(fabs((sw_mhd1d_total_mass(&grid) - initial_mass) / initial_mass) < 1.0e-10);
    sw_mhd1d_destroy(&grid);
}

static bool oscillator_rhs(double time, const double *state, double *derivative, size_t dimension, void *context) {
    const double omega = *(const double *)context;
    (void)time;
    if (dimension != 2U) {
        return false;
    }
    derivative[0] = state[1];
    derivative[1] = -omega * omega * state[0];
    return true;
}

static void test_integrators(void) {
    double y[2] = {1.0, 0.0};
    const double omega = 2.0;
    sw_ode_workspace workspace;
    size_t i;
    SW_CHECK(sw_ode_workspace_init(&workspace, 2U, 9U));
    for (i = 0U; i < 1000U; ++i) {
        SW_CHECK(sw_rk4_step(oscillator_rhs, (void *)&omega, (double)i * 0.001, 0.001, y, 2U, &workspace));
    }
    SW_CHECK(sw_near(y[0], cos(2.0), 1.0e-10, 1.0e-10));
    SW_CHECK(sw_near(y[1], -2.0 * sin(2.0), 1.0e-10, 1.0e-10));
    {
        double adaptive[2] = {1.0, 0.0};
        double accepted = 0.0;
        double suggested = 0.0;
        unsigned rejected = 0U;
        SW_CHECK(sw_dopri54_step(
            oscillator_rhs,
            (void *)&omega,
            0.0,
            5.0,
            1.0e-12,
            1.0e-12,
            adaptive,
            2U,
            &workspace,
            &accepted,
            &suggested,
            &rejected
        ));
        SW_CHECK(rejected > 0U);
        SW_CHECK(accepted < 5.0);
        SW_CHECK(suggested > 0.0);
    }
    sw_ode_workspace_destroy(&workspace);
}

static void test_force_and_simulation(void) {
    sw_simulation_config config;
    sw_simulation simulation;
    sw_simulation_diagnostics diagnostics;
    size_t i;
    sw_simulation_config_default(&config);
    config.duration_s = 600.0;
    config.initial_step_s = 1.0;
    config.max_step_s = 1.0;
    config.output_every_steps = 10U;
    SW_CHECK(sw_simulation_init(&simulation, &config));
    SW_CHECK(sw_simulation_diagnostics_evaluate(&simulation, &diagnostics));
    SW_CHECK(diagnostics.environment.proton_number_density_m3 > 0.0);
    SW_CHECK(sw_v3_norm(diagnostics.forces.radiation_force_n) > 0.0);
    SW_CHECK(sw_v3_norm(diagnostics.forces.electric_sail_force_n) > 0.0);
    for (i = 0U; i < 600U; ++i) {
        SW_CHECK(sw_simulation_step(&simulation));
    }
    SW_CHECK(sw_v3_isfinite(simulation.state.position_m));
    SW_CHECK(sw_q_norm(simulation.state.attitude_body_to_inertial) > 0.999999);
    SW_CHECK(sw_q_norm(simulation.state.attitude_body_to_inertial) < 1.000001);
    sw_simulation_destroy(&simulation);
}

int main(void) {
    test_vectors_and_matrix();
    test_metric();
    test_additional_metrics();
    test_gravity_grid();
    test_parker_environment();
    test_verlet_orbit();
    test_nbody();
    test_boris();
    test_pic();
    test_fdtd();
    test_mhd();
    test_integrators();
    test_force_and_simulation();
    if (sw_tests_failed == 0U) {
        printf("PASS: %u checks\n", sw_tests_run);
        return 0;
    }
    fprintf(stderr, "FAILED: %u of %u checks\n", sw_tests_failed, sw_tests_run);
    return 1;
}
