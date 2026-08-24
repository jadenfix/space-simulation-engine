#include "spacewind/simulation.h"
#include "spacewind/constants.h"
#include "spacewind/metric.h"
#include "spacewind/version.h"

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define SW_SPACECRAFT_ODE_DIMENSION 14U

typedef struct {
    sw_simulation *simulation;
    sw_control_output frozen_control;
} sw_dynamics_context;

static double sw_clamp_sim(double x, double lo, double hi) {
    return fmax(lo, fmin(hi, x));
}

static void sw_state_pack(const sw_spacecraft_state *state, double y[SW_SPACECRAFT_ODE_DIMENSION]) {
    y[0] = state->position_m.x;
    y[1] = state->position_m.y;
    y[2] = state->position_m.z;
    y[3] = state->velocity_m_s.x;
    y[4] = state->velocity_m_s.y;
    y[5] = state->velocity_m_s.z;
    y[6] = state->attitude_body_to_inertial.w;
    y[7] = state->attitude_body_to_inertial.x;
    y[8] = state->attitude_body_to_inertial.y;
    y[9] = state->attitude_body_to_inertial.z;
    y[10] = state->angular_velocity_body_rad_s.x;
    y[11] = state->angular_velocity_body_rad_s.y;
    y[12] = state->angular_velocity_body_rad_s.z;
    y[13] = state->proper_time_s;
}

static void sw_state_unpack(const double y[SW_SPACECRAFT_ODE_DIMENSION], sw_spacecraft_state *state) {
    state->position_m = sw_v3(y[0], y[1], y[2]);
    state->velocity_m_s = sw_v3(y[3], y[4], y[5]);
    state->attitude_body_to_inertial = sw_q_normalize(sw_q(y[6], y[7], y[8], y[9]));
    state->angular_velocity_body_rad_s = sw_v3(y[10], y[11], y[12]);
    state->proper_time_s = y[13];
}

void sw_simulation_config_default(sw_simulation_config *config) {
    sw_body sun;
    memset(config, 0, sizeof(*config));
    sw_solar_wind_default(&config->wind);
    sw_solar_wind_init_turbulence(
        &config->wind,
        42U,
        8U,
        0.5e-9,
        0.005 * SW_AU,
        0.1 * SW_AU
    );
    sw_vehicle_default(&config->vehicle);
    sw_controller_default(&config->controller);
    sw_spacecraft_state_default(&config->initial_state);

    memset(&sun, 0, sizeof(sun));
    sun.mass_kg = SW_SOLAR_MASS;
    sun.gravitational_mu = SW_SOLAR_MU;
    sun.radius_m = SW_SOLAR_RADIUS;
    sun.spin_axis = sw_v3(0.0, 0.0, 1.0);
    config->gravity.count = 1U;
    config->gravity.bodies[0] = sun;

    config->integrator = SW_INTEGRATOR_RK4;
    config->start_time_s = 0.0;
    config->duration_s = 30.0 * SW_DAY;
    config->initial_step_s = 60.0;
    config->min_step_s = 0.01;
    config->max_step_s = 3600.0;
    config->absolute_tolerance = 1.0e-7;
    config->relative_tolerance = 1.0e-10;
    config->output_every_steps = 60U;
    config->enable_gravity = true;
    config->enable_j2 = false;
    config->enable_1pn = true;
    config->enable_radiation_pressure = true;
    config->enable_electric_sail = true;
    config->enable_magnetic_sail = false;
    config->enable_lorentz_force = false;
    config->seed = 42U;
}

bool sw_simulation_init(sw_simulation *simulation, const sw_simulation_config *config) {
    if (simulation == NULL || config == NULL || config->vehicle.mass_kg <= 0.0
        || config->duration_s < 0.0 || config->initial_step_s <= 0.0) {
        return false;
    }
    memset(simulation, 0, sizeof(*simulation));
    simulation->config = *config;
    simulation->state = config->initial_state;
    simulation->state.attitude_body_to_inertial = sw_q_normalize(
        simulation->state.attitude_body_to_inertial
    );
    simulation->controller = config->controller;
    simulation->time_s = config->start_time_s;
    simulation->next_step_s = config->initial_step_s;
    if (!sw_ode_workspace_init(&simulation->workspace, SW_SPACECRAFT_ODE_DIMENSION, 9U)) {
        return false;
    }
    return true;
}

void sw_simulation_destroy(sw_simulation *simulation) {
    if (simulation != NULL) {
        sw_ode_workspace_destroy(&simulation->workspace);
    }
}

static bool sw_force_evaluate(
    sw_simulation *simulation,
    double time_s,
    const sw_spacecraft_state *state,
    const sw_control_output *control,
    sw_environment_sample *environment,
    sw_force_breakdown *forces
) {
    const sw_simulation_config *config = &simulation->config;
    sw_vec3 gravity_acceleration;
    sw_vec3 pn_acceleration;
    sw_vec3 actual_normal;
    sw_vec3 physical_force_body;
    sw_vec3 physical_torque_body;

    memset(forces, 0, sizeof(*forces));
    if (!sw_environment_evaluate(&config->wind, time_s, state->position_m, environment)) {
        return false;
    }

    actual_normal = sw_q_rotate(
        state->attitude_body_to_inertial,
        sw_v3(1.0, 0.0, 0.0)
    );

    if (config->enable_gravity) {
        gravity_acceleration = config->enable_j2
            ? sw_gravity_acceleration_with_j2(&config->gravity, state->position_m)
            : sw_gravity_acceleration(&config->gravity, state->position_m);
        forces->gravity_force_n = sw_v3_scale(gravity_acceleration, config->vehicle.mass_kg);
    }

    if (config->enable_1pn && config->gravity.count > 0U) {
        const sw_body *central = &config->gravity.bodies[0];
        const sw_vec3 relative_position = sw_v3_sub(state->position_m, central->position_m);
        const sw_vec3 relative_velocity = sw_v3_sub(state->velocity_m_s, central->velocity_m_s);
        pn_acceleration = sw_schwarzschild_1pn_acceleration(
            central->gravitational_mu,
            relative_position,
            relative_velocity
        );
        forces->post_newtonian_force_n = sw_v3_scale(pn_acceleration, config->vehicle.mass_kg);
    }

    if (config->enable_radiation_pressure) {
        forces->radiation_force_n = sw_solar_radiation_force(
            &config->vehicle,
            state,
            environment,
            actual_normal
        );
    }
    if (config->enable_electric_sail) {
        forces->electric_sail_force_n = sw_electric_sail_force(
            &config->vehicle,
            state,
            environment,
            control,
            &forces->electric_effective_radius_m
        );
    }
    if (config->enable_magnetic_sail) {
        forces->magnetic_sail_force_n = sw_magnetic_sail_force(
            &config->vehicle,
            state,
            environment,
            control,
            &forces->magnetic_effective_radius_m
        );
    }
    if (config->enable_lorentz_force) {
        forces->lorentz_force_n = sw_lorentz_force(&config->vehicle, state, environment);
    }

    forces->total_force_n = sw_v3_add(
        sw_v3_add(forces->gravity_force_n, forces->post_newtonian_force_n),
        sw_v3_add(
            sw_v3_add(forces->radiation_force_n, forces->electric_sail_force_n),
            sw_v3_add(forces->magnetic_sail_force_n, forces->lorentz_force_n)
        )
    );
    forces->relative_wind_speed_m_s = sw_v3_norm(sw_v3_sub(
        environment->wind_velocity_m_s,
        state->velocity_m_s
    ));

    physical_force_body = sw_q_rotate(
        sw_q_conj(state->attitude_body_to_inertial),
        forces->radiation_force_n
    );
    physical_torque_body = sw_v3_cross(
        config->vehicle.center_of_pressure_body_m,
        physical_force_body
    );
    forces->attitude_torque_body_n_m = sw_v3_add(
        physical_torque_body,
        sw_attitude_control_torque_body(
            &config->vehicle,
            state,
            control->desired_sail_normal
        )
    );
    return sw_v3_isfinite(forces->total_force_n);
}

static bool sw_spacecraft_rhs(
    double time,
    const double *y,
    double *dy,
    size_t dimension,
    void *context
) {
    sw_dynamics_context *dynamics = (sw_dynamics_context *)context;
    sw_spacecraft_state state;
    sw_environment_sample environment;
    sw_force_breakdown forces;
    sw_vec3 acceleration;
    sw_quat qdot;
    sw_vec3 angular_acceleration;
    double proper_time_rate = 1.0;

    if (dimension != SW_SPACECRAFT_ODE_DIMENSION) {
        return false;
    }
    sw_state_unpack(y, &state);
    if (!sw_force_evaluate(
        dynamics->simulation,
        time,
        &state,
        &dynamics->frozen_control,
        &environment,
        &forces
    )) {
        return false;
    }

    acceleration = sw_v3_scale(
        forces.total_force_n,
        1.0 / dynamics->simulation->config.vehicle.mass_kg
    );
    qdot = sw_q_derivative_body_rate(
        state.attitude_body_to_inertial,
        state.angular_velocity_body_rad_s
    );
    angular_acceleration = sw_rigid_body_angular_acceleration(
        &dynamics->simulation->config.vehicle,
        state.angular_velocity_body_rad_s,
        forces.attitude_torque_body_n_m
    );
    if (dynamics->simulation->config.gravity.count > 0U) {
        const sw_body *central = &dynamics->simulation->config.gravity.bodies[0];
        proper_time_rate = sw_weak_field_proper_time_rate(
            central->gravitational_mu,
            sw_v3_sub(state.position_m, central->position_m),
            sw_v3_sub(state.velocity_m_s, central->velocity_m_s)
        );
    }

    dy[0] = state.velocity_m_s.x;
    dy[1] = state.velocity_m_s.y;
    dy[2] = state.velocity_m_s.z;
    dy[3] = acceleration.x;
    dy[4] = acceleration.y;
    dy[5] = acceleration.z;
    dy[6] = qdot.w;
    dy[7] = qdot.x;
    dy[8] = qdot.y;
    dy[9] = qdot.z;
    dy[10] = angular_acceleration.x;
    dy[11] = angular_acceleration.y;
    dy[12] = angular_acceleration.z;
    dy[13] = proper_time_rate;
    return true;
}

bool sw_simulation_step(sw_simulation *simulation) {
    sw_environment_sample environment;
    sw_control_output control;
    sw_dynamics_context context;
    double y[SW_SPACECRAFT_ODE_DIMENSION];
    double remaining;
    double requested_step;
    bool ok;

    if (simulation == NULL) {
        return false;
    }
    remaining = simulation->config.start_time_s + simulation->config.duration_s - simulation->time_s;
    if (remaining <= 0.0) {
        return true;
    }
    if (!sw_environment_evaluate(
        &simulation->config.wind,
        simulation->time_s,
        simulation->state.position_m,
        &environment
    )) {
        return false;
    }
    sw_controller_evaluate(&simulation->controller, &simulation->state, &environment, &control);
    context.simulation = simulation;
    context.frozen_control = control;
    sw_state_pack(&simulation->state, y);

    requested_step = fmin(simulation->next_step_s, remaining);
    requested_step = sw_clamp_sim(
        requested_step,
        simulation->config.min_step_s,
        simulation->config.max_step_s
    );
    requested_step = fmin(requested_step, remaining);

    if (simulation->config.integrator == SW_INTEGRATOR_DOPRI54) {
        double accepted = 0.0;
        double suggested = requested_step;
        unsigned rejected_attempts = 0U;
        ok = sw_dopri54_step(
            sw_spacecraft_rhs,
            &context,
            simulation->time_s,
            requested_step,
            simulation->config.absolute_tolerance,
            simulation->config.relative_tolerance,
            y,
            SW_SPACECRAFT_ODE_DIMENSION,
            &simulation->workspace,
            &accepted,
            &suggested,
            &rejected_attempts
        );
        simulation->rejected_steps += (uint64_t)rejected_attempts;
        if (!ok) {
            return false;
        }
        simulation->time_s += accepted;
        simulation->next_step_s = sw_clamp_sim(
            suggested,
            simulation->config.min_step_s,
            simulation->config.max_step_s
        );
    } else {
        ok = sw_rk4_step(
            sw_spacecraft_rhs,
            &context,
            simulation->time_s,
            requested_step,
            y,
            SW_SPACECRAFT_ODE_DIMENSION,
            &simulation->workspace
        );
        if (!ok) {
            return false;
        }
        simulation->time_s += requested_step;
        simulation->next_step_s = simulation->config.initial_step_s;
    }

    sw_state_unpack(y, &simulation->state);
    simulation->accepted_steps += 1U;
    return sw_v3_isfinite(simulation->state.position_m)
        && sw_v3_isfinite(simulation->state.velocity_m_s);
}

bool sw_simulation_diagnostics_evaluate(
    sw_simulation *simulation,
    sw_simulation_diagnostics *diagnostics
) {
    sw_controller controller_copy;
    double central_mu = 0.0;
    if (simulation == NULL || diagnostics == NULL) {
        return false;
    }
    memset(diagnostics, 0, sizeof(*diagnostics));
    if (!sw_environment_evaluate(
        &simulation->config.wind,
        simulation->time_s,
        simulation->state.position_m,
        &diagnostics->environment
    )) {
        return false;
    }
    controller_copy = simulation->controller;
    sw_controller_evaluate(
        &controller_copy,
        &simulation->state,
        &diagnostics->environment,
        &diagnostics->control
    );
    if (!sw_force_evaluate(
        simulation,
        simulation->time_s,
        &simulation->state,
        &diagnostics->control,
        &diagnostics->environment,
        &diagnostics->forces
    )) {
        return false;
    }
    if (simulation->config.gravity.count > 0U) {
        central_mu = simulation->config.gravity.bodies[0].gravitational_mu;
    }
    diagnostics->specific_orbital_energy_j_kg = sw_specific_orbital_energy(
        central_mu,
        simulation->state.position_m,
        simulation->state.velocity_m_s
    );
    diagnostics->specific_angular_momentum_m2_s = sw_specific_angular_momentum(
        simulation->state.position_m,
        simulation->state.velocity_m_s
    );
    diagnostics->kinetic_energy_j = 0.5 * simulation->config.vehicle.mass_kg
        * sw_v3_norm2(simulation->state.velocity_m_s);
    diagnostics->potential_energy_j = central_mu > 0.0
        ? -simulation->config.vehicle.mass_kg * central_mu / sw_v3_norm(simulation->state.position_m)
        : 0.0;
    diagnostics->total_mechanical_energy_j = diagnostics->kinetic_energy_j
        + diagnostics->potential_energy_j;
    return true;
}

static void sw_csv_header(FILE *csv) {
    fprintf(csv,
        "time_s,proper_time_s,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,r_m,speed_m_s,"
        "wind_vx_m_s,wind_vy_m_s,wind_vz_m_s,number_density_m3,temperature_p_k,"
        "bx_t,by_t,bz_t,ex_v_m,ey_v_m,ez_v_m,photon_flux_w_m2,dynamic_pressure_pa,"
        "fx_n,fy_n,fz_n,frad_n,fesail_n,fmagsail_n,florentz_n,"
        "electric_radius_m,magnetic_radius_m,specific_energy_j_kg,mechanical_energy_j\n"
    );
}

static void sw_csv_row(FILE *csv, const sw_simulation *simulation, const sw_simulation_diagnostics *d) {
    fprintf(csv,
        "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
        "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
        "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
        simulation->time_s,
        simulation->state.proper_time_s,
        simulation->state.position_m.x,
        simulation->state.position_m.y,
        simulation->state.position_m.z,
        simulation->state.velocity_m_s.x,
        simulation->state.velocity_m_s.y,
        simulation->state.velocity_m_s.z,
        sw_v3_norm(simulation->state.position_m),
        sw_v3_norm(simulation->state.velocity_m_s),
        d->environment.wind_velocity_m_s.x,
        d->environment.wind_velocity_m_s.y,
        d->environment.wind_velocity_m_s.z,
        d->environment.proton_number_density_m3,
        d->environment.proton_temperature_k,
        d->environment.magnetic_field_t.x,
        d->environment.magnetic_field_t.y,
        d->environment.magnetic_field_t.z,
        d->environment.electric_field_v_m.x,
        d->environment.electric_field_v_m.y,
        d->environment.electric_field_v_m.z,
        d->environment.photon_flux_w_m2,
        d->environment.dynamic_pressure_pa,
        d->forces.total_force_n.x,
        d->forces.total_force_n.y,
        d->forces.total_force_n.z,
        sw_v3_norm(d->forces.radiation_force_n),
        sw_v3_norm(d->forces.electric_sail_force_n),
        sw_v3_norm(d->forces.magnetic_sail_force_n),
        sw_v3_norm(d->forces.lorentz_force_n),
        d->forces.electric_effective_radius_m,
        d->forces.magnetic_effective_radius_m,
        d->specific_orbital_energy_j_kg,
        d->total_mechanical_energy_j
    );
}

bool sw_simulation_run_csv(
    sw_simulation *simulation,
    FILE *csv,
    FILE *receipt_json
) {
    sw_simulation_diagnostics diagnostics;
    double end_time;
    if (csv == NULL || simulation == NULL) {
        return false;
    }
    end_time = simulation->config.start_time_s + simulation->config.duration_s;
    sw_csv_header(csv);
    if (!sw_simulation_diagnostics_evaluate(simulation, &diagnostics)) {
        return false;
    }
    sw_csv_row(csv, simulation, &diagnostics);
    while (simulation->time_s < end_time) {
        if (!sw_simulation_step(simulation)) {
            return false;
        }
        {
            const uint64_t output_stride = (uint64_t)(simulation->config.output_every_steps > 0U
                ? simulation->config.output_every_steps : 1U);
            if (simulation->accepted_steps % output_stride == 0U
                || simulation->time_s >= end_time) {
                if (!sw_simulation_diagnostics_evaluate(simulation, &diagnostics)) {
                    return false;
                }
                sw_csv_row(csv, simulation, &diagnostics);
            }
        }
    }
    if (receipt_json != NULL) {
        fprintf(receipt_json,
            "{\n"
            "  \"schema\": \"spacewind.receipt.v3\",\n"
            "  \"engine_version\": \"%s\",\n"
            "  \"build\": {\n"
            "    \"compiler_id\": \"%s\",\n"
            "    \"compiler_version\": \"%s\",\n"
            "    \"c_standard\": %ld,\n"
            "    \"sizeof_double\": %zu,\n"
            "    \"flt_radix\": %d,\n"
            "    \"dbl_mant_dig\": %d,\n"
            "    \"dbl_max_exp\": %d,\n"
            "    \"flt_eval_method\": %d,\n"
            "    \"fast_math\": %s\n"
            "  },\n"
            "  \"completed\": true,\n"
            "  \"config_fnv1a64\": ",
            SPACEWIND_VERSION_STRING,
            SPACEWIND_COMPILER_ID,
            SPACEWIND_COMPILER_VERSION,
            (long)__STDC_VERSION__,
            sizeof(double),
            FLT_RADIX,
            DBL_MANT_DIG,
            DBL_MAX_EXP,
            FLT_EVAL_METHOD,
            SPACEWIND_FAST_MATH_JSON
        );
        if (simulation->config.source_config_hash_available) {
            fprintf(receipt_json, "\"%016" PRIx64 "\",\n",
                simulation->config.source_config_fnv1a64);
        } else {
            fprintf(receipt_json, "null,\n");
        }
        fprintf(receipt_json,
            "  \"seed\": %" PRIu64 ",\n"
            "  \"integrator\": \"%s\",\n"
            "  \"controller\": \"%s\",\n"
            "  \"start_time_s\": %.17g,\n"
            "  \"end_time_s\": %.17g,\n"
            "  \"accepted_steps\": %" PRIu64 ",\n"
            "  \"rejected_trial_steps\": %" PRIu64 ",\n"
            "  \"final_position_m\": [%.17g, %.17g, %.17g],\n"
            "  \"final_velocity_m_s\": [%.17g, %.17g, %.17g],\n"
            "  \"final_radius_m\": %.17g,\n"
            "  \"final_speed_m_s\": %.17g,\n"
            "  \"proper_time_s\": %.17g,\n"
            "  \"enabled_models\": {\n"
            "    \"gravity\": %s,\n"
            "    \"j2\": %s,\n"
            "    \"schwarzschild_1pn\": %s,\n"
            "    \"solar_radiation_pressure\": %s,\n"
            "    \"parameterized_electric_sail\": %s,\n"
            "    \"parameterized_magnetic_sail\": %s,\n"
            "    \"lorentz_force\": %s,\n"
            "    \"synthetic_shear_layer\": %s,\n"
            "    \"synthetic_cme_shell\": %s\n"
            "  },\n"
            "  \"nonclaim\": \"Field-sail coefficients and effective radii are phenomenological inputs, not flight-validated performance.\"\n"
            "}\n",
            simulation->config.seed,
            simulation->config.integrator == SW_INTEGRATOR_DOPRI54 ? "dopri54" : "rk4",
            simulation->controller.mode == SW_CONTROL_DYNAMIC_SOAR ? "dynamic_soar"
                : (simulation->controller.mode == SW_CONTROL_SUN_CONE ? "sun_cone" : "fixed"),
            simulation->config.start_time_s,
            simulation->time_s,
            simulation->accepted_steps,
            simulation->rejected_steps,
            simulation->state.position_m.x,
            simulation->state.position_m.y,
            simulation->state.position_m.z,
            simulation->state.velocity_m_s.x,
            simulation->state.velocity_m_s.y,
            simulation->state.velocity_m_s.z,
            sw_v3_norm(simulation->state.position_m),
            sw_v3_norm(simulation->state.velocity_m_s),
            simulation->state.proper_time_s,
            simulation->config.enable_gravity ? "true" : "false",
            simulation->config.enable_j2 ? "true" : "false",
            simulation->config.enable_1pn ? "true" : "false",
            simulation->config.enable_radiation_pressure ? "true" : "false",
            simulation->config.enable_electric_sail ? "true" : "false",
            simulation->config.enable_magnetic_sail ? "true" : "false",
            simulation->config.enable_lorentz_force ? "true" : "false",
            simulation->config.wind.shear_enabled ? "true" : "false",
            simulation->config.wind.cme_enabled ? "true" : "false"
        );
    }
    return true;
}

static char *sw_trim(char *text) {
    char *end;
    while (*text != '\0' && isspace((unsigned char)*text)) {
        ++text;
    }
    if (*text == '\0') {
        return text;
    }
    end = text + strlen(text) - 1;
    while (end > text && isspace((unsigned char)*end)) {
        *end-- = '\0';
    }
    return text;
}

static bool sw_parse_double(const char *value, double *out) {
    char *end = NULL;
    errno = 0;
    *out = strtod(value, &end);
    return errno == 0 && end != value && *sw_trim(end) == '\0' && isfinite(*out);
}

static bool sw_parse_u64(const char *value, uint64_t *out) {
    char *end = NULL;
    unsigned long long parsed;
    if (value == NULL || out == NULL || value[0] == '-') {
        return false;
    }
    errno = 0;
    parsed = strtoull(value, &end, 10);
    if (errno != 0 || end == value || *sw_trim(end) != '\0') {
        return false;
    }
    *out = (uint64_t)parsed;
    return true;
}

static bool sw_parse_bool(const char *value, bool *out) {
    if (strcmp(value, "true") == 0 || strcmp(value, "1") == 0 || strcmp(value, "yes") == 0) {
        *out = true;
        return true;
    }
    if (strcmp(value, "false") == 0 || strcmp(value, "0") == 0 || strcmp(value, "no") == 0) {
        *out = false;
        return true;
    }
    return false;
}

static bool sw_parse_vec3(const char *value, sw_vec3 *out) {
    char buffer[256];
    char *first;
    char *second;
    double x;
    double y;
    double z;
    if (strlen(value) >= sizeof(buffer)) {
        return false;
    }
    strcpy(buffer, value);
    first = strchr(buffer, ',');
    if (first == NULL) {
        return false;
    }
    *first++ = '\0';
    second = strchr(first, ',');
    if (second == NULL) {
        return false;
    }
    *second++ = '\0';
    if (!sw_parse_double(sw_trim(buffer), &x)
        || !sw_parse_double(sw_trim(first), &y)
        || !sw_parse_double(sw_trim(second), &z)) {
        return false;
    }
    *out = sw_v3(x, y, z);
    return true;
}

static bool sw_apply_config_value(sw_simulation_config *c, const char *key, const char *value) {
#define SW_DBL(name, field) if (strcmp(key, name) == 0) return sw_parse_double(value, &(field))
#define SW_BOOL(name, field) if (strcmp(key, name) == 0) return sw_parse_bool(value, &(field))
#define SW_VEC(name, field) if (strcmp(key, name) == 0) return sw_parse_vec3(value, &(field))
    SW_DBL("start_time_s", c->start_time_s);
    SW_DBL("duration_s", c->duration_s);
    SW_DBL("step_s", c->initial_step_s);
    SW_DBL("min_step_s", c->min_step_s);
    SW_DBL("max_step_s", c->max_step_s);
    SW_DBL("absolute_tolerance", c->absolute_tolerance);
    SW_DBL("relative_tolerance", c->relative_tolerance);
    SW_DBL("mass_kg", c->vehicle.mass_kg);
    SW_DBL("sail_area_m2", c->vehicle.area_m2);
    SW_DBL("sail_reflectivity", c->vehicle.reflectivity);
    SW_DBL("sail_absorptivity", c->vehicle.absorptivity);
    SW_DBL("spacecraft_charge_c", c->vehicle.charge_c);
    SW_DBL("tether_length_m", c->vehicle.total_tether_length_m);
    SW_DBL("tether_voltage_v", c->vehicle.tether_voltage_v);
    SW_DBL("electric_sail_cd", c->vehicle.electric_sail_cd);
    SW_DBL("electric_sail_cl", c->vehicle.electric_sail_cl);
    SW_DBL("electric_sail_max_radius_m", c->vehicle.electric_sail_max_radius_m);
    SW_DBL("magnetic_dipole_a_m2", c->vehicle.magnetic_dipole_moment_a_m2);
    SW_DBL("magnetic_sail_cd", c->vehicle.magnetic_sail_cd);
    SW_DBL("magnetic_sail_cl", c->vehicle.magnetic_sail_cl);
    SW_DBL("magnetic_sail_max_radius_m", c->vehicle.magnetic_sail_max_radius_m);
    SW_VEC("position_m", c->initial_state.position_m);
    SW_VEC("velocity_m_s", c->initial_state.velocity_m_s);
    SW_VEC("fixed_sail_normal", c->controller.fixed_normal_inertial);
    if (strcmp(key, "sun_cone_angle_deg") == 0) {
        double degrees;
        if (!sw_parse_double(value, &degrees)) {
            return false;
        }
        c->controller.sun_cone_angle_rad = degrees * SW_PI / 180.0;
        return true;
    }
    SW_DBL("soar_turn_distance_m", c->controller.soar_turn_distance_m);
    SW_DBL("wind_proton_temperature_k", c->wind.proton_temperature_k);
    SW_DBL("wind_electron_temperature_k", c->wind.electron_temperature_k);
    SW_DBL("wind_density_1au_m3", c->wind.reference_number_density_m3);
    SW_DBL("wind_radial_b_1au_t", c->wind.reference_radial_field_t);
    SW_DBL("wind_stream_fraction", c->wind.stream_fraction);
    SW_BOOL("shear_enabled", c->wind.shear_enabled);
    SW_VEC("shear_normal", c->wind.shear_normal);
    SW_DBL("shear_offset_m", c->wind.shear_offset_m);
    SW_DBL("shear_width_m", c->wind.shear_width_m);
    SW_DBL("shear_delta_speed_m_s", c->wind.shear_delta_speed_m_s);
    SW_VEC("shear_flow_direction", c->wind.shear_flow_direction);
    SW_BOOL("cme_enabled", c->wind.cme_enabled);
    SW_DBL("cme_launch_time_s", c->wind.cme_launch_time_s);
    SW_DBL("cme_speed_m_s", c->wind.cme_speed_m_s);
    SW_DBL("cme_width_m", c->wind.cme_width_m);
    SW_DBL("cme_density_multiplier", c->wind.cme_density_multiplier);
    SW_DBL("cme_field_multiplier", c->wind.cme_field_multiplier);
    SW_DBL("cme_speed_increment_m_s", c->wind.cme_speed_increment_m_s);
    SW_BOOL("enable_gravity", c->enable_gravity);
    SW_BOOL("enable_j2", c->enable_j2);
    SW_BOOL("enable_1pn", c->enable_1pn);
    SW_BOOL("enable_radiation_pressure", c->enable_radiation_pressure);
    SW_BOOL("enable_electric_sail", c->enable_electric_sail);
    SW_BOOL("enable_magnetic_sail", c->enable_magnetic_sail);
    SW_BOOL("enable_lorentz_force", c->enable_lorentz_force);
    if (strcmp(key, "integrator") == 0) {
        if (strcmp(value, "rk4") == 0) {
            c->integrator = SW_INTEGRATOR_RK4;
            return true;
        }
        if (strcmp(value, "dopri54") == 0) {
            c->integrator = SW_INTEGRATOR_DOPRI54;
            return true;
        }
        return false;
    }
    if (strcmp(key, "controller") == 0) {
        if (strcmp(value, "fixed") == 0) {
            c->controller.mode = SW_CONTROL_FIXED;
            return true;
        }
        if (strcmp(value, "sun_cone") == 0) {
            c->controller.mode = SW_CONTROL_SUN_CONE;
            return true;
        }
        if (strcmp(value, "dynamic_soar") == 0) {
            c->controller.mode = SW_CONTROL_DYNAMIC_SOAR;
            return true;
        }
        return false;
    }
    if (strcmp(key, "output_every_steps") == 0) {
        uint64_t parsed;
        if (!sw_parse_u64(value, &parsed) || parsed > (uint64_t)SIZE_MAX) {
            return false;
        }
        c->output_every_steps = (size_t)parsed;
        return true;
    }
    if (strcmp(key, "seed") == 0) {
        return sw_parse_u64(value, &c->seed);
    }
    if (strcmp(key, "wind_stream_arms") == 0) {
        uint64_t parsed;
        if (!sw_parse_u64(value, &parsed) || parsed > (uint64_t)UINT_MAX) {
            return false;
        }
        c->wind.stream_arms = (unsigned)parsed;
        return true;
    }
#undef SW_DBL
#undef SW_BOOL
#undef SW_VEC
    return false;
}

bool sw_simulation_config_load_file(
    const char *path,
    sw_simulation_config *config,
    char *error_buffer,
    size_t error_buffer_size
) {
    FILE *file;
    char line[1024];
    unsigned line_number = 0U;
    if (path == NULL || config == NULL || error_buffer == NULL || error_buffer_size == 0U) {
        return false;
    }
    sw_simulation_config_default(config);
    file = fopen(path, "r");
    if (file == NULL) {
        snprintf(error_buffer, error_buffer_size, "cannot open %s", path);
        return false;
    }
    while (fgets(line, sizeof(line), file) != NULL) {
        char *text;
        char *equals;
        char *key;
        char *value;
        ++line_number;
        text = sw_trim(line);
        if (*text == '\0' || *text == '#' || *text == ';') {
            continue;
        }
        equals = strchr(text, '=');
        if (equals == NULL) {
            snprintf(error_buffer, error_buffer_size, "%s:%u: expected key=value", path, line_number);
            fclose(file);
            return false;
        }
        *equals = '\0';
        key = sw_trim(text);
        value = sw_trim(equals + 1);
        {
            char *comment = strchr(value, '#');
            if (comment != NULL) {
                *comment = '\0';
                value = sw_trim(value);
            }
        }
        if (!sw_apply_config_value(config, key, value)) {
            snprintf(error_buffer, error_buffer_size, "%s:%u: invalid or unknown key '%s'", path, line_number, key);
            fclose(file);
            return false;
        }
    }
    fclose(file);

    if (config->initial_step_s <= 0.0 || config->min_step_s <= 0.0
        || config->max_step_s < config->min_step_s || config->vehicle.mass_kg <= 0.0) {
        snprintf(error_buffer, error_buffer_size, "%s: inconsistent numerical or vehicle parameters", path);
        return false;
    }
    sw_solar_wind_init_turbulence(
        &config->wind,
        config->seed,
        8U,
        0.5e-9,
        0.005 * SW_AU,
        0.1 * SW_AU
    );
    config->source_config_fnv1a64 = sw_fnv1a64_file(
        path, &config->source_config_hash_available
    );
    return true;
}

uint64_t sw_fnv1a64_file(const char *path, bool *ok) {
    FILE *file;
    uint64_t hash = UINT64_C(14695981039346656037);
    int byte;
    if (ok != NULL) {
        *ok = false;
    }
    if (path == NULL) {
        return 0U;
    }
    file = fopen(path, "rb");
    if (file == NULL) {
        return 0U;
    }
    while ((byte = fgetc(file)) != EOF) {
        hash ^= (uint64_t)(unsigned char)byte;
        hash *= UINT64_C(1099511628211);
    }
    if (ok != NULL) {
        *ok = ferror(file) == 0;
    }
    fclose(file);
    return hash;
}
