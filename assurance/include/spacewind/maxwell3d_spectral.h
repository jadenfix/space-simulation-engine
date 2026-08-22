#ifndef SPACEWIND_MAXWELL3D_SPECTRAL_H
#define SPACEWIND_MAXWELL3D_SPECTRAL_H

#include "spacewind/assurance.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double imaginary_relative_tolerance;
    double electric_divergence_relative_tolerance;
    double magnetic_divergence_relative_tolerance;
    double energy_relative_tolerance;
} swa_maxwell3d_spectral_limits;

typedef struct {
    size_t mode_count;
    double duration_s;
    double maximum_modified_wavenumber_rad_m;
    double maximum_angular_frequency_rad_s;

    double initial_energy_J;
    double final_energy_J;
    double energy_residual_J;
    double energy_relative;

    double initial_electric_divergence_rms_V_m2;
    double final_electric_divergence_rms_V_m2;
    double initial_electric_divergence_relative;
    double final_electric_divergence_relative;

    double initial_magnetic_divergence_rms_T_m;
    double final_magnetic_divergence_rms_T_m;
    double initial_magnetic_divergence_relative;
    double final_magnetic_divergence_relative;

    double maximum_reconstruction_imaginary;
    double maximum_reconstruction_real;
    double reconstruction_imaginary_relative;

    int finite;
    int reconstruction_passes;
    int electric_divergence_passes;
    int magnetic_divergence_passes;
    int energy_passes;
    int passes;
} swa_maxwell3d_spectral_result;

/*
 * Advance a periodic, source-free, staggered Yee field exactly in time for
 * the semi-discrete centered-difference Maxwell operator. The component
 * locations are:
 *
 *   Ex: (i,     j+1/2, k+1/2)   Bx: (i+1/2, j,     k)
 *   Ey: (i+1/2, j,     k+1/2)   By: (i,     j+1/2, k)
 *   Ez: (i+1/2, j+1/2, k)       Bz: (i,     j,     k+1/2)
 *
 * A direct phase-aware DFT is used deliberately as a small-grid reference;
 * no FFT library or production solver implementation is shared. Negative
 * duration_s is supported for time-reversal checks. Inputs and corresponding
 * outputs may alias for in-place advancement.
 *
 * The function returns nonzero when the computation executes and populates
 * result. Scientific acceptance is reported separately in result->passes.
 */
int swa_maxwell3d_spectral_advance(
    size_t nx,
    size_t ny,
    size_t nz,
    double length_x_m,
    double length_y_m,
    double length_z_m,
    double duration_s,
    const double *ex_initial_V_m,
    const double *ey_initial_V_m,
    const double *ez_initial_V_m,
    const double *bx_initial_T,
    const double *by_initial_T,
    const double *bz_initial_T,
    double *ex_final_V_m,
    double *ey_final_V_m,
    double *ez_final_V_m,
    double *bx_final_T,
    double *by_final_T,
    double *bz_final_T,
    const swa_maxwell3d_spectral_limits *limits,
    swa_maxwell3d_spectral_result *result
);

int swa_maxwell3d_spectral_write_receipt(
    FILE *fp,
    size_t nx,
    size_t ny,
    size_t nz,
    double length_x_m,
    double length_y_m,
    double length_z_m,
    const swa_maxwell3d_spectral_limits *limits,
    const swa_maxwell3d_spectral_result *result
);

#ifdef __cplusplus
}
#endif

#endif
