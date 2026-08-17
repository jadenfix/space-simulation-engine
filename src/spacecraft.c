#include "spacewind/spacecraft.h"
#include "spacewind/constants.h"

#include <float.h>
#include <math.h>
#include <string.h>

static double sw_clamp_local(double x, double lo, double hi) {
    return fmax(lo, fmin(hi, x));
}

void sw_vehicle_default(sw_vehicle *vehicle) {
    memset(vehicle, 0, sizeof(*vehicle));
    vehicle->mass_kg = 12.0;
    vehicle->area_m2 = 100.0;
    vehicle->reflectivity = 0.88;
    vehicle->absorptivity = 0.12;
    vehicle->charge_c = 0.0;
    vehicle->total_tether_length_m = 20000.0;
    vehicle->tether_voltage_v = 20000.0;
    vehicle->electric_sail_cd = 1.0;
    vehicle->electric_sail_cl = 0.15;
    vehicle->electric_sail_max_radius_m = 200.0;
    vehicle->magnetic_dipole_moment_a_m2 = 0.0;
    vehicle->magnetic_sail_cd = 1.0;
    vehicle->magnetic_sail_cl = 0.1;
    vehicle->magnetic_sail_max_radius_m = 5000.0;
    vehicle->center_of_pressure_body_m = sw_v3(0.0, 0.0, 0.1);
    vehicle->inertia_kg_m2 = sw_v3(8.0, 8.0, 2.0);
    vehicle->attitude_kp = 1.0e-4;
    vehicle->attitude_kd = 1.0e-2;
}

void sw_spacecraft_state_default(sw_spacecraft_state *state) {
    memset(state, 0, sizeof(*state));
    state->position_m = sw_v3(SW_AU, 0.0, 0.0);
    state->velocity_m_s = sw_v3(0.0, sqrt(SW_SOLAR_MU / SW_AU), 0.0);
    state->attitude_body_to_inertial = sw_q_identity();
}

void sw_controller_default(sw_controller *controller) {
    memset(controller, 0, sizeof(*controller));
    controller->mode = SW_CONTROL_SUN_CONE;
    controller->fixed_normal_inertial = sw_v3(1.0, 0.0, 0.0);
    controller->sun_cone_angle_rad = 25.0 * SW_PI / 180.0;
    controller->crossing_direction = 1;
    controller->soar_turn_distance_m = 3.0e8;
}

static sw_vec3 sw_unit_perpendicular(sw_vec3 direction, sw_vec3 preferred) {
    sw_vec3 perpendicular = sw_v3_reject(preferred, direction);
    if (sw_v3_norm2(perpendicular) < 1.0e-20) {
        const sw_vec3 fallback = fabs(direction.z) < 0.8
            ? sw_v3(0.0, 0.0, 1.0)
            : sw_v3(0.0, 1.0, 0.0);
        perpendicular = sw_v3_reject(fallback, direction);
    }
    return sw_v3_normalize(perpendicular);
}

void sw_controller_evaluate(
    sw_controller *controller,
    const sw_spacecraft_state *state,
    const sw_environment_sample *environment,
    sw_control_output *output
) {
    const sw_vec3 radial = sw_v3_normalize(state->position_m);
    const sw_vec3 relative_wind = sw_v3_sub(
        environment->wind_velocity_m_s,
        state->velocity_m_s
    );
    const sw_vec3 wind_hat = sw_v3_normalize(relative_wind);
    memset(output, 0, sizeof(*output));
    output->electric_voltage_scale = 1.0;
    output->magnetic_moment_scale = 1.0;

    switch (controller->mode) {
        case SW_CONTROL_FIXED:
            output->desired_sail_normal = sw_v3_normalize(controller->fixed_normal_inertial);
            output->desired_lift_direction = sw_unit_perpendicular(
                wind_hat,
                environment->shear_normal
            );
            break;

        case SW_CONTROL_DYNAMIC_SOAR: {
            sw_vec3 lift_basis;
            if (!controller->initialized) {
                const double normal_velocity = sw_v3_dot(
                    state->velocity_m_s,
                    environment->shear_normal
                );
                controller->crossing_direction = normal_velocity >= 0.0 ? 1 : -1;
                controller->initialized = true;
            }
            if (fabs(environment->shear_coordinate_m) >= controller->soar_turn_distance_m) {
                controller->crossing_direction = environment->shear_coordinate_m > 0.0 ? -1 : 1;
            }
            lift_basis = sw_unit_perpendicular(wind_hat, environment->shear_normal);
            if (sw_v3_dot(lift_basis, environment->shear_normal) * (double)controller->crossing_direction < 0.0) {
                lift_basis = sw_v3_scale(lift_basis, -1.0);
            }
            output->desired_lift_direction = lift_basis;
            output->desired_sail_normal = sw_v3_normalize(sw_v3_add(radial, sw_v3_scale(lift_basis, 0.25)));
            break;
        }

        case SW_CONTROL_SUN_CONE:
        default: {
            const sw_vec3 prograde = sw_unit_perpendicular(radial, state->velocity_m_s);
            output->desired_sail_normal = sw_v3_normalize(sw_v3_add(
                sw_v3_scale(radial, cos(controller->sun_cone_angle_rad)),
                sw_v3_scale(prograde, sin(controller->sun_cone_angle_rad))
            ));
            output->desired_lift_direction = sw_unit_perpendicular(wind_hat, prograde);
            break;
        }
    }
}

sw_vec3 sw_solar_radiation_force(
    const sw_vehicle *vehicle,
    const sw_spacecraft_state *state,
    const sw_environment_sample *environment,
    sw_vec3 sail_normal_inertial
) {
    sw_vec3 normal = sw_v3_normalize(sail_normal_inertial);
    const sw_vec3 ray_direction = sw_v3_normalize(state->position_m);
    double cosine;
    double pressure;
    sw_vec3 absorption;
    sw_vec3 reflection;

    cosine = sw_v3_dot(normal, ray_direction);
    if (cosine < 0.0) {
        normal = sw_v3_scale(normal, -1.0);
        cosine = -cosine;
    }
    cosine = sw_clamp_local(cosine, 0.0, 1.0);
    pressure = environment->photon_flux_w_m2 / SW_C;
    absorption = sw_v3_scale(
        ray_direction,
        pressure * vehicle->area_m2 * vehicle->absorptivity * cosine
    );
    reflection = sw_v3_scale(
        normal,
        pressure * vehicle->area_m2 * 2.0 * vehicle->reflectivity * cosine * cosine
    );
    return sw_v3_add(absorption, reflection);
}

static sw_vec3 sw_pressure_force(
    sw_vec3 relative_wind,
    double dynamic_pressure_pa,
    double effective_area_m2,
    double cd,
    double cl,
    sw_vec3 requested_lift_direction
) {
    const sw_vec3 drag_direction = sw_v3_normalize(relative_wind);
    sw_vec3 lift_direction = sw_v3_reject(requested_lift_direction, drag_direction);
    lift_direction = sw_v3_normalize(lift_direction);
    return sw_v3_scale(
        sw_v3_add(
            sw_v3_scale(drag_direction, cd),
            sw_v3_scale(lift_direction, cl)
        ),
        dynamic_pressure_pa * effective_area_m2
    );
}

sw_vec3 sw_electric_sail_force(
    const sw_vehicle *vehicle,
    const sw_spacecraft_state *state,
    const sw_environment_sample *environment,
    const sw_control_output *control,
    double *effective_radius_m
) {
    const sw_vec3 relative_wind = sw_v3_sub(
        environment->wind_velocity_m_s,
        state->velocity_m_s
    );
    const double speed2 = sw_v3_norm2(relative_wind);
    const double dynamic_pressure = environment->mass_density_kg_m3 * speed2;
    const double electron_thermal_energy = SW_KB * fmax(environment->electron_temperature_k, 1.0);
    const double voltage_energy = SW_QE * fabs(vehicle->tether_voltage_v * control->electric_voltage_scale);
    const double sheath_multiplier = sqrt(log1p(voltage_energy / electron_thermal_energy));
    const double radius = sw_clamp_local(
        environment->debye_length_m * sheath_multiplier,
        0.0,
        fmax(vehicle->electric_sail_max_radius_m, 0.0)
    );
    const double effective_area = 2.0 * vehicle->total_tether_length_m * radius;
    if (effective_radius_m != NULL) {
        *effective_radius_m = radius;
    }
    if (effective_area <= 0.0 || dynamic_pressure <= 0.0) {
        return sw_v3(0.0, 0.0, 0.0);
    }
    return sw_pressure_force(
        relative_wind,
        dynamic_pressure,
        effective_area,
        vehicle->electric_sail_cd,
        vehicle->electric_sail_cl,
        control->desired_lift_direction
    );
}

sw_vec3 sw_magnetic_sail_force(
    const sw_vehicle *vehicle,
    const sw_spacecraft_state *state,
    const sw_environment_sample *environment,
    const sw_control_output *control,
    double *effective_radius_m
) {
    const sw_vec3 relative_wind = sw_v3_sub(
        environment->wind_velocity_m_s,
        state->velocity_m_s
    );
    const double dynamic_pressure = environment->mass_density_kg_m3 * sw_v3_norm2(relative_wind);
    const double moment = fabs(
        vehicle->magnetic_dipole_moment_a_m2 * control->magnetic_moment_scale
    );
    double radius = 0.0;
    double effective_area;

    if (moment > 0.0 && dynamic_pressure > DBL_MIN) {
        radius = pow(
            SW_MU0 * moment * moment / (8.0 * SW_PI * SW_PI * dynamic_pressure),
            1.0 / 6.0
        );
        radius = sw_clamp_local(radius, 0.0, fmax(vehicle->magnetic_sail_max_radius_m, 0.0));
    }
    if (effective_radius_m != NULL) {
        *effective_radius_m = radius;
    }
    effective_area = SW_PI * radius * radius;
    if (effective_area <= 0.0) {
        return sw_v3(0.0, 0.0, 0.0);
    }
    return sw_pressure_force(
        relative_wind,
        dynamic_pressure,
        effective_area,
        vehicle->magnetic_sail_cd,
        vehicle->magnetic_sail_cl,
        control->desired_lift_direction
    );
}

sw_vec3 sw_lorentz_force(
    const sw_vehicle *vehicle,
    const sw_spacecraft_state *state,
    const sw_environment_sample *environment
) {
    const sw_vec3 motional = sw_v3_cross(state->velocity_m_s, environment->magnetic_field_t);
    return sw_v3_scale(
        sw_v3_add(environment->electric_field_v_m, motional),
        vehicle->charge_c
    );
}

sw_vec3 sw_attitude_control_torque_body(
    const sw_vehicle *vehicle,
    const sw_spacecraft_state *state,
    sw_vec3 desired_normal_inertial
) {
    const sw_vec3 actual_normal_inertial = sw_q_rotate(
        state->attitude_body_to_inertial,
        sw_v3(1.0, 0.0, 0.0)
    );
    const sw_vec3 error_axis_inertial = sw_v3_cross(
        actual_normal_inertial,
        sw_v3_normalize(desired_normal_inertial)
    );
    const sw_vec3 error_axis_body = sw_q_rotate(
        sw_q_conj(state->attitude_body_to_inertial),
        error_axis_inertial
    );
    return sw_v3_sub(
        sw_v3_scale(error_axis_body, vehicle->attitude_kp),
        sw_v3_scale(state->angular_velocity_body_rad_s, vehicle->attitude_kd)
    );
}

sw_vec3 sw_rigid_body_angular_acceleration(
    const sw_vehicle *vehicle,
    sw_vec3 omega_body,
    sw_vec3 torque_body
) {
    const sw_vec3 inertia_omega = sw_v3(
        vehicle->inertia_kg_m2.x * omega_body.x,
        vehicle->inertia_kg_m2.y * omega_body.y,
        vehicle->inertia_kg_m2.z * omega_body.z
    );
    const sw_vec3 gyroscopic = sw_v3_cross(omega_body, inertia_omega);
    const sw_vec3 net = sw_v3_sub(torque_body, gyroscopic);
    return sw_v3(
        vehicle->inertia_kg_m2.x > DBL_MIN ? net.x / vehicle->inertia_kg_m2.x : 0.0,
        vehicle->inertia_kg_m2.y > DBL_MIN ? net.y / vehicle->inertia_kg_m2.y : 0.0,
        vehicle->inertia_kg_m2.z > DBL_MIN ? net.z / vehicle->inertia_kg_m2.z : 0.0
    );
}
