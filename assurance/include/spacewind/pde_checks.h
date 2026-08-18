#ifndef SPACEWIND_PDE_CHECKS_H
#define SPACEWIND_PDE_CHECKS_H

#include "spacewind/assurance.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t resolution;
    double spacing;
    double l2_error;
    double linf_error;
    int finite;
} swa_discretization_error;

typedef struct {
    size_t resolution;
    size_t steps;
    double cfl;
    double l2_electric_error_V_m;
    double linf_electric_error_V_m;
    double initial_electric_rms_V_m;
    double final_electric_rms_V_m;
    double rms_amplitude_relative_error;
    int finite;
    int stable;
} swa_maxwell_wave_result;

typedef struct {
    size_t steps_per_gyroperiod;
    double speed_initial_mps;
    double speed_final_mps;
    double speed_relative_error;
    double velocity_closure_error_mps;
    double position_closure_error_m;
    int finite;
} swa_boris_orbit_result;

int swa_poisson_manufactured_error(size_t intervals,
                                   swa_discretization_error *out);
int swa_continuity_manufactured_error(size_t intervals,
                                      double cfl,
                                      swa_discretization_error *out);
int swa_vacuum_maxwell_period(size_t cells,
                              double cfl,
                              swa_maxwell_wave_result *out);
int swa_uniform_boris_gyroperiod(size_t steps_per_period,
                                 double charge_to_mass_C_kg,
                                 double magnetic_field_T,
                                 double initial_speed_mps,
                                 swa_boris_orbit_result *out);
int swa_parker_implicit_residual(double radius_over_critical,
                                 double speed_over_sound,
                                 double *residual);

#ifdef __cplusplus
}
#endif

#endif
