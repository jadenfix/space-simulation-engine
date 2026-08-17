#ifndef SPACEWIND_SIMULATION_H
#define SPACEWIND_SIMULATION_H

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#include "spacewind/fields.h"
#include "spacewind/integrator.h"
#include "spacewind/spacecraft.h"

typedef enum {
    SW_INTEGRATOR_RK4 = 0,
    SW_INTEGRATOR_DOPRI54 = 1
} sw_integrator_kind;

typedef struct {
    sw_solar_wind_model wind;
    sw_gravity_system gravity;
    sw_vehicle vehicle;
    sw_controller controller;
    sw_spacecraft_state initial_state;

    sw_integrator_kind integrator;
    double start_time_s;
    double duration_s;
    double initial_step_s;
    double min_step_s;
    double max_step_s;
    double absolute_tolerance;
    double relative_tolerance;
    size_t output_every_steps;

    bool enable_gravity;
    bool enable_j2;
    bool enable_1pn;
    bool enable_radiation_pressure;
    bool enable_electric_sail;
    bool enable_magnetic_sail;
    bool enable_lorentz_force;

    uint64_t seed;
    uint64_t source_config_fnv1a64;
    bool source_config_hash_available;
} sw_simulation_config;

typedef struct {
    sw_simulation_config config;
    sw_spacecraft_state state;
    sw_controller controller;
    double time_s;
    double next_step_s;
    uint64_t accepted_steps;
    uint64_t rejected_steps;
    sw_ode_workspace workspace;
} sw_simulation;

typedef struct {
    sw_environment_sample environment;
    sw_control_output control;
    sw_force_breakdown forces;
    double specific_orbital_energy_j_kg;
    sw_vec3 specific_angular_momentum_m2_s;
    double kinetic_energy_j;
    double potential_energy_j;
    double total_mechanical_energy_j;
} sw_simulation_diagnostics;

void sw_simulation_config_default(sw_simulation_config *config);
bool sw_simulation_init(sw_simulation *simulation, const sw_simulation_config *config);
void sw_simulation_destroy(sw_simulation *simulation);
bool sw_simulation_step(sw_simulation *simulation);
bool sw_simulation_diagnostics_evaluate(
    sw_simulation *simulation,
    sw_simulation_diagnostics *diagnostics
);
bool sw_simulation_run_csv(
    sw_simulation *simulation,
    FILE *csv,
    FILE *receipt_json
);

bool sw_simulation_config_load_file(
    const char *path,
    sw_simulation_config *config,
    char *error_buffer,
    size_t error_buffer_size
);

uint64_t sw_fnv1a64_file(const char *path, bool *ok);

#endif
