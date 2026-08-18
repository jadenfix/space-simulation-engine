#ifndef SPACEWIND_MAXWELL1D_REFERENCE_H
#define SPACEWIND_MAXWELL1D_REFERENCE_H

#include "spacewind/assurance.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t samples;
    double length_m;
    double cross_section_area_m2;
    double *electric_y_V_m;
    double *electric_z_V_m;
    double *magnetic_y_T;
    double *magnetic_z_T;
} swa_maxwell1d_state;

typedef struct {
    double field_energy_initial_J;
    double field_energy_final_J;
    double right_moving_energy_initial_J;
    double right_moving_energy_final_J;
    double left_moving_energy_initial_J;
    double left_moving_energy_final_J;
    double field_momentum_initial_Ns;
    double field_momentum_final_Ns;
    double current_work_J;
    double energy_residual_J;
    double energy_relative_residual;
    double vacuum_momentum_residual_Ns;
    double vacuum_momentum_relative_residual;
    double spectral_imaginary_maximum;
    int vacuum_step;
    int finite;
    int energy_closes;
    int vacuum_momentum_closes;
    int passes;
} swa_maxwell1d_step_result;

int swa_maxwell1d_init(
    swa_maxwell1d_state *state,
    size_t samples,
    double length_m,
    double cross_section_area_m2
);

void swa_maxwell1d_destroy(swa_maxwell1d_state *state);

int swa_maxwell1d_copy(
    const swa_maxwell1d_state *source,
    swa_maxwell1d_state *destination
);

int swa_maxwell1d_state_is_finite(const swa_maxwell1d_state *state);

double swa_maxwell1d_field_energy_J(const swa_maxwell1d_state *state);
double swa_maxwell1d_field_momentum_Ns(const swa_maxwell1d_state *state);

int swa_maxwell1d_propagate_vacuum(
    swa_maxwell1d_state *state,
    double dt_s,
    double *spectral_imaginary_maximum
);

int swa_maxwell1d_advance(
    swa_maxwell1d_state *state,
    const double *current_y_A_m2,
    const double *current_z_A_m2,
    double dt_s,
    double relative_tolerance,
    double absolute_energy_tolerance_J,
    double absolute_momentum_tolerance_Ns,
    swa_maxwell1d_step_result *out
);

int swa_maxwell1d_audit_energy_balance(
    double field_energy_initial_J,
    double field_energy_final_J,
    double current_work_J,
    double relative_tolerance,
    double absolute_tolerance_J,
    double *residual_J,
    double *relative_residual,
    int *passes
);

int swa_write_maxwell1d_reference_receipt(
    FILE *fp,
    const swa_maxwell1d_state *state,
    const swa_maxwell1d_step_result *result
);

#ifdef __cplusplus
}
#endif

#endif
