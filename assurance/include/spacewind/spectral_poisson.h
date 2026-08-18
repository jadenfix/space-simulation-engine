#ifndef SPACEWIND_SPECTRAL_POISSON_H
#define SPACEWIND_SPECTRAL_POISSON_H

#include "spacewind/assurance.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t samples;
    double length_m;
    double permittivity_F_m;
    const double *charge_density_C_m3;
    const double *potential_V;
    const double *electric_field_V_m;
} swa_spectral_poisson_input;

typedef struct {
    double operator_relative_tolerance;
    double charge_absolute_tolerance_C_m3;
    double electric_absolute_tolerance_V_m;
    double energy_relative_tolerance;
    double energy_absolute_tolerance_J_m2;
    double neutrality_relative_tolerance;
    double neutrality_absolute_tolerance_C_m3;
    double harmonic_field_relative_tolerance;
    double harmonic_field_absolute_tolerance_V_m;
    double nyquist_charge_power_fraction_tolerance;
} swa_spectral_poisson_tolerances;

typedef struct {
    size_t samples;
    size_t resolved_nonzero_modes;
    int nyquist_present;
    double charge_mean_C_m3;
    double charge_rms_C_m3;
    double electric_mean_V_m;
    double electric_rms_V_m;
    double charge_zero_mode_relative;
    double harmonic_field_relative;
    double nyquist_charge_rms_C_m3;
    double nyquist_charge_power_fraction;
    double nyquist_gradient_residual_V_m;
    double nyquist_gauss_residual_C_m3;
    double poisson_residual_rms_C_m3;
    double poisson_residual_max_spectral_C_m3;
    double poisson_relative_rms;
    double poisson_relative_maximum;
    double gradient_residual_rms_V_m;
    double gradient_residual_max_spectral_V_m;
    double gradient_relative_rms;
    double gradient_relative_maximum;
    double gauss_residual_rms_C_m3;
    double gauss_residual_max_spectral_C_m3;
    double gauss_relative_rms;
    double gauss_relative_maximum;
    double field_energy_sample_J_m2;
    double field_energy_parseval_J_m2;
    double charge_potential_energy_J_m2;
    double parseval_energy_relative_error;
    double electrostatic_energy_relative_error;
    int finite;
    int neutrality_passes;
    int harmonic_field_passes;
    int nyquist_passes;
    int poisson_passes;
    int gradient_passes;
    int gauss_passes;
    int parseval_passes;
    int electrostatic_energy_passes;
    int passes;
} swa_spectral_poisson_result;

void swa_default_spectral_poisson_tolerances(
    swa_spectral_poisson_tolerances *out
);

int swa_audit_periodic_spectral_poisson(
    const swa_spectral_poisson_input *input,
    const swa_spectral_poisson_tolerances *tolerances,
    swa_spectral_poisson_result *out
);

int swa_write_spectral_poisson_receipt(
    FILE *fp,
    const swa_spectral_poisson_input *input,
    const swa_spectral_poisson_tolerances *tolerances,
    const swa_spectral_poisson_result *result
);

#ifdef __cplusplus
}
#endif

#endif
