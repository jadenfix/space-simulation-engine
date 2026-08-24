#ifndef SPACEWIND_SPACECRAFT_H
#define SPACEWIND_SPACECRAFT_H

#include <stdbool.h>

#include "spacewind/fields.h"
#include "spacewind/math3.h"

typedef enum {
    SW_CONTROL_FIXED = 0,
    SW_CONTROL_SUN_CONE = 1,
    SW_CONTROL_DYNAMIC_SOAR = 2
} sw_control_mode;

typedef struct {
    double mass_kg;
    double area_m2;
    double reflectivity;
    double absorptivity;
    double charge_c;

    double total_tether_length_m;
    double tether_voltage_v;
    double electric_sail_cd;
    double electric_sail_cl;
    double electric_sail_max_radius_m;

    double magnetic_dipole_moment_a_m2;
    double magnetic_sail_cd;
    double magnetic_sail_cl;
    double magnetic_sail_max_radius_m;

    sw_vec3 center_of_pressure_body_m;
    sw_vec3 inertia_kg_m2;
    double attitude_kp;
    double attitude_kd;
} sw_vehicle;

typedef struct {
    sw_vec3 position_m;
    sw_vec3 velocity_m_s;
    sw_quat attitude_body_to_inertial;
    sw_vec3 angular_velocity_body_rad_s;
    double proper_time_s;
} sw_spacecraft_state;

typedef struct {
    sw_control_mode mode;
    sw_vec3 fixed_normal_inertial;
    double sun_cone_angle_rad;
    int crossing_direction;
    double soar_turn_distance_m;
    bool initialized;
} sw_controller;

typedef struct {
    sw_vec3 desired_sail_normal;
    sw_vec3 desired_lift_direction;
    double electric_voltage_scale;
    double magnetic_moment_scale;
} sw_control_output;

typedef struct {
    sw_vec3 gravity_force_n;
    sw_vec3 post_newtonian_force_n;
    sw_vec3 radiation_force_n;
    sw_vec3 electric_sail_force_n;
    sw_vec3 magnetic_sail_force_n;
    sw_vec3 lorentz_force_n;
    sw_vec3 total_force_n;
    sw_vec3 attitude_torque_body_n_m;
    double electric_effective_radius_m;
    double magnetic_effective_radius_m;
    double relative_wind_speed_m_s;
} sw_force_breakdown;

void sw_vehicle_default(sw_vehicle *vehicle);
void sw_spacecraft_state_default(sw_spacecraft_state *state);
void sw_controller_default(sw_controller *controller);

void sw_controller_evaluate(
    sw_controller *controller,
    const sw_spacecraft_state *state,
    const sw_environment_sample *environment,
    sw_control_output *output
);

sw_vec3 sw_solar_radiation_force(
    const sw_vehicle *vehicle,
    const sw_spacecraft_state *state,
    const sw_environment_sample *environment,
    sw_vec3 sail_normal_inertial
);

sw_vec3 sw_electric_sail_force(
    const sw_vehicle *vehicle,
    const sw_spacecraft_state *state,
    const sw_environment_sample *environment,
    const sw_control_output *control,
    double *effective_radius_m
);

sw_vec3 sw_magnetic_sail_force(
    const sw_vehicle *vehicle,
    const sw_spacecraft_state *state,
    const sw_environment_sample *environment,
    const sw_control_output *control,
    double *effective_radius_m
);

sw_vec3 sw_lorentz_force(
    const sw_vehicle *vehicle,
    const sw_spacecraft_state *state,
    const sw_environment_sample *environment
);

sw_vec3 sw_attitude_control_torque_body(
    const sw_vehicle *vehicle,
    const sw_spacecraft_state *state,
    sw_vec3 desired_normal_inertial
);

sw_vec3 sw_rigid_body_angular_acceleration(
    const sw_vehicle *vehicle,
    sw_vec3 omega_body,
    sw_vec3 torque_body
);

#endif
