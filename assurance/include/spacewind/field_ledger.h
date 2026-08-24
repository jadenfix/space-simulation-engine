#ifndef SPACEWIND_FIELD_LEDGER_H
#define SPACEWIND_FIELD_LEDGER_H

#include "spacewind/assurance.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t nx;
    size_t ny;
    double dx_m;
    double dy_m;
    double dt_s;
    const double *rho_initial_C_m3;
    const double *rho_final_C_m3;
    const double *jx_face_A_m2;
    const double *jy_face_A_m2;
} swa_charge_continuity_input;

typedef struct {
    double initial_charge_C_m;
    double final_charge_C_m;
    double outward_boundary_current_A_m;
    double residual_mean_absolute_A_m3;
    double residual_rms_A_m3;
    double residual_maximum_A_m3;
    double residual_relative_rms;
    double residual_relative_maximum;
    double residual_integral_A_m;
    double global_charge_residual_C_m;
    double global_charge_relative;
    double divergence_theorem_error_C_m;
    int finite;
    int local_passes;
    int global_passes;
    int passes;
} swa_charge_continuity_result;

typedef struct {
    size_t nx;
    size_t ny;
    double dx_m;
    double dy_m;
    const double *x_face_field;
    const double *y_face_field;
    const double *cell_source;
} swa_staggered_scalar_constraint_input;

typedef struct {
    double residual_mean_absolute;
    double residual_rms;
    double residual_maximum;
    double residual_relative_rms;
    double residual_relative_maximum;
    double residual_signed_integral;
    int finite;
    int passes;
} swa_field_constraint_result;

typedef struct {
    double field_energy_initial_J;
    double field_energy_final_J;
    double particle_work_J;
    double outward_poynting_energy_J;
    double impressed_source_energy_J;
} swa_poynting_ledger;

typedef struct {
    double residual_J;
    double relative_residual;
    int finite;
    int passes;
} swa_poynting_ledger_result;

typedef struct {
    swa_vec3 field_momentum_initial_Ns;
    swa_vec3 field_momentum_final_Ns;
    swa_vec3 mechanical_impulse_Ns;
    swa_vec3 outward_maxwell_impulse_Ns;
    swa_vec3 impressed_external_impulse_Ns;
} swa_field_momentum_ledger;

typedef struct {
    swa_vec3 residual_Ns;
    double residual_norm_Ns;
    double relative_residual;
    int finite;
    int passes;
} swa_field_momentum_result;

typedef struct {
    double energy_density_J_m3;
    swa_vec3 poynting_flux_W_m2;
    swa_vec3 momentum_density_Ns_m3;
    double invariant_B2_minus_E2_over_c2_T2;
    double invariant_E_dot_B_over_c;
    double dominant_energy_ratio;
    int finite;
} swa_local_em_result;

typedef struct {
    swa_vec3 unit_normal;
    swa_vec3 traction_N_m2;
    double normal_pressure_Pa;
    int finite;
} swa_maxwell_traction_result;

int swa_audit_charge_continuity(
    const swa_charge_continuity_input *input,
    double relative_tolerance,
    double local_absolute_tolerance_A_m3,
    double global_absolute_tolerance_C_m,
    swa_charge_continuity_result *out
);

int swa_audit_gauss_law(
    const swa_staggered_scalar_constraint_input *input,
    double relative_tolerance,
    double absolute_tolerance_C_m3,
    swa_field_constraint_result *out
);

int swa_audit_magnetic_divergence(
    const swa_staggered_scalar_constraint_input *input,
    double relative_tolerance,
    double absolute_tolerance_T_m,
    swa_field_constraint_result *out
);

int swa_audit_poynting_ledger(
    const swa_poynting_ledger *ledger,
    double relative_tolerance,
    double absolute_tolerance_J,
    swa_poynting_ledger_result *out
);

int swa_audit_field_momentum(
    const swa_field_momentum_ledger *ledger,
    double relative_tolerance,
    double absolute_tolerance_Ns,
    swa_field_momentum_result *out
);

int swa_local_electromagnetic_state(
    swa_vec3 electric_field_V_m,
    swa_vec3 magnetic_field_T,
    swa_local_em_result *out
);

int swa_maxwell_traction(
    swa_vec3 electric_field_V_m,
    swa_vec3 magnetic_field_T,
    swa_vec3 surface_normal,
    swa_maxwell_traction_result *out
);

int swa_write_field_ledger_receipt(
    FILE *fp,
    const swa_charge_continuity_result *continuity,
    const swa_field_constraint_result *gauss,
    const swa_field_constraint_result *magnetic_divergence,
    const swa_poynting_ledger_result *energy,
    const swa_field_momentum_result *momentum
);

#ifdef __cplusplus
}
#endif

#endif
