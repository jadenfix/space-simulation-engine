#ifndef SPACEWIND_FIELDS_H
#define SPACEWIND_FIELDS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "spacewind/math3.h"

#define SW_MAX_TURBULENCE_MODES 16U
#define SW_MAX_BODIES 16U

typedef struct {
    sw_vec3 wave_vector;
    sw_vec3 polarization;
    double angular_frequency;
    double phase;
    double magnetic_amplitude_t;
    double velocity_sign;
} sw_turbulence_mode;

typedef struct {
    double central_mu;
    double proton_temperature_k;
    double electron_temperature_k;
    double reference_radius_m;
    double reference_number_density_m3;
    double reference_radial_field_t;
    double rotation_rate_rad_s;
    double source_surface_radius_m;

    double stream_fraction;
    unsigned stream_arms;
    double stream_phase_rad;
    double stream_latitude_width_rad;

    bool shear_enabled;
    sw_vec3 shear_normal;
    double shear_offset_m;
    double shear_width_m;
    double shear_delta_speed_m_s;
    sw_vec3 shear_flow_direction;

    bool cme_enabled;
    double cme_launch_time_s;
    double cme_launch_radius_m;
    double cme_speed_m_s;
    double cme_width_m;
    double cme_density_multiplier;
    double cme_field_multiplier;
    double cme_speed_increment_m_s;

    size_t turbulence_mode_count;
    sw_turbulence_mode turbulence_modes[SW_MAX_TURBULENCE_MODES];
} sw_solar_wind_model;

typedef struct {
    double mass_kg;
    double gravitational_mu;
    double radius_m;
    sw_vec3 position_m;
    sw_vec3 velocity_m_s;
    double j2;
    sw_vec3 spin_axis;
} sw_body;

typedef struct {
    size_t count;
    sw_body bodies[SW_MAX_BODIES];
} sw_gravity_system;

typedef struct {
    double radius_m;
    double photon_flux_w_m2;
    double proton_number_density_m3;
    double mass_density_kg_m3;
    double proton_temperature_k;
    double electron_temperature_k;
    double thermal_pressure_pa;
    double dynamic_pressure_pa;
    double magnetic_pressure_pa;
    double debye_length_m;
    double alfven_speed_m_s;
    double sound_speed_m_s;
    double fast_magnetosonic_speed_m_s;
    sw_vec3 wind_velocity_m_s;
    sw_vec3 electric_field_v_m;
    sw_vec3 magnetic_field_t;
    sw_vec3 shear_normal;
    double shear_coordinate_m;
} sw_environment_sample;

void sw_solar_wind_default(sw_solar_wind_model *model);
void sw_solar_wind_init_turbulence(
    sw_solar_wind_model *model,
    uint64_t seed,
    size_t mode_count,
    double rms_magnetic_field_t,
    double min_wavelength_m,
    double max_wavelength_m
);

double sw_parker_sound_speed(const sw_solar_wind_model *model);
double sw_parker_critical_radius(const sw_solar_wind_model *model);
double sw_parker_speed(const sw_solar_wind_model *model, double radius_m);

bool sw_environment_evaluate(
    const sw_solar_wind_model *model,
    double time_s,
    sw_vec3 position_m,
    sw_environment_sample *sample
);

sw_vec3 sw_gravity_acceleration(
    const sw_gravity_system *system,
    sw_vec3 position_m
);

sw_vec3 sw_gravity_acceleration_with_j2(
    const sw_gravity_system *system,
    sw_vec3 position_m
);

double sw_specific_orbital_energy(double central_mu, sw_vec3 position_m, sw_vec3 velocity_m_s);
sw_vec3 sw_specific_angular_momentum(sw_vec3 position_m, sw_vec3 velocity_m_s);

#endif
