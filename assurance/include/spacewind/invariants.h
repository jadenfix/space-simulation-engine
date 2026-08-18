#ifndef SPACEWIND_INVARIANTS_H
#define SPACEWIND_INVARIANTS_H

#include "spacewind/assurance.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double t;
    double x;
    double y;
    double z;
} swa_four_vector;

typedef struct {
    swa_vec3 electric_V_m;
    swa_vec3 magnetic_T;
    swa_vec3 velocity_mps;
    double charge_C;
} swa_em_state;

typedef struct {
    swa_vec3 lorentz_force_N;
    swa_vec3 electric_force_N;
    swa_vec3 magnetic_force_N;
    swa_vec3 poynting_W_m2;
    double field_energy_density_J_m3;
    double invariant_B2_minus_E2_over_c2_T2;
    double invariant_E_dot_B_over_c;
    double total_particle_power_W;
    double electric_particle_power_W;
    double magnetic_particle_power_W;
    double magnetic_power_relative;
    int finite;
    int magnetic_work_free;
} swa_em_audit;

typedef struct {
    double radius_m;
    double number_density_m3;
    double radial_speed_mps;
    double radial_magnetic_field_T;
    swa_vec3 bulk_velocity_mps;
    swa_vec3 magnetic_field_T;
    swa_vec3 electric_field_V_m;
} swa_parker_sample;

typedef struct {
    double mass_flux_relative_error;
    double magnetic_flux_relative_error;
    double motional_field_relative_error_a;
    double motional_field_relative_error_b;
    double electric_dot_magnetic_relative_a;
    double electric_dot_magnetic_relative_b;
    double electric_dot_velocity_relative_a;
    double electric_dot_velocity_relative_b;
    int finite;
    int passes;
} swa_parker_audit;

typedef struct {
    double requested_force_N;
    double upper_bound_force_N;
    double utilization;
    int finite;
    int passes;
} swa_force_bound;

typedef struct {
    double specific_energy_J_kg;
    swa_vec3 specific_angular_momentum_m2_s;
    swa_vec3 eccentricity_vector;
    double eccentricity;
    double semimajor_axis_m;
    double periapsis_m;
    double apoapsis_m;
    int bound;
    int finite;
} swa_orbit_invariants;

typedef struct {
    double total_mass_kg;
    swa_vec3 center_of_mass_m;
    swa_vec3 total_momentum_Ns;
    swa_vec3 angular_momentum_kg_m2_s;
    double kinetic_energy_J;
    double potential_energy_J;
    double total_energy_J;
    int finite;
} swa_nbody_invariants;

double swa_minkowski_norm(swa_four_vector v);
int swa_four_velocity_from_three_velocity(swa_vec3 velocity_mps, swa_four_vector *out);
int swa_audit_four_velocity(swa_four_vector four_velocity_mps, double relative_tolerance);
int swa_audit_em_state(const swa_em_state *state, double relative_tolerance,
                       swa_em_audit *out);
int swa_audit_parker_pair(const swa_parker_sample *a, const swa_parker_sample *b,
                          double relative_tolerance, swa_parker_audit *out);
int swa_audit_photon_force(double irradiance_W_m2, double area_m2,
                           double reflectivity, double requested_force_N,
                           double relative_tolerance, swa_force_bound *out);
int swa_audit_plasma_momentum_force(double mass_density_kg_m3,
                                    double relative_speed_mps,
                                    double interaction_area_m2,
                                    double momentum_multiplier,
                                    double requested_force_N,
                                    double relative_tolerance,
                                    swa_force_bound *out);
int swa_compute_orbit_invariants(double gravitational_parameter_m3_s2,
                                 swa_vec3 position_m, swa_vec3 velocity_mps,
                                 swa_orbit_invariants *out);
int swa_compute_nbody_invariants(size_t count, const double *mass_kg,
                                 const swa_vec3 *position_m,
                                 const swa_vec3 *velocity_mps,
                                 double gravitational_constant,
                                 double softening_m,
                                 swa_nbody_invariants *out);

#ifdef __cplusplus
}
#endif

#endif
