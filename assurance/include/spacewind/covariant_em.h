#ifndef SPACEWIND_COVARIANT_EM_H
#define SPACEWIND_COVARIANT_EM_H

#include "spacewind/assurance.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double component[4];
} swa_four;

typedef struct {
    double component[4][4];
} swa_tensor4;

typedef struct {
    swa_vec3 electric_field_V_m;
    swa_vec3 magnetic_field_T;
} swa_em_fields;

typedef struct {
    double invariant_B2_minus_E2_over_c2_T2;
    double invariant_E_dot_B_over_c;
    double energy_density_J_m3;
    swa_vec3 poynting_flux_W_m2;
    double dominant_energy_ratio;
    double tensor_trace_J_m3;
    int finite;
    int dominant_energy_condition;
} swa_covariant_em_diagnostics;

typedef struct {
    swa_four field_initial_Ns;
    swa_four field_final_Ns;
    swa_four matter_initial_Ns;
    swa_four matter_final_Ns;
    swa_four external_impulse_Ns;
    swa_four momentum_in_Ns;
    swa_four momentum_out_Ns;
} swa_four_momentum_ledger;

typedef struct {
    swa_four residual_Ns;
    double residual_norm_Ns;
    double relative_residual;
    int finite;
    int passes;
} swa_four_momentum_ledger_result;

swa_four swa_four_make(double time_component, double x, double y, double z);
swa_four swa_four_add(swa_four a, swa_four b);
swa_four swa_four_sub(swa_four a, swa_four b);
swa_four swa_four_scale(swa_four a, double scale);
double swa_four_minkowski_dot(swa_four a, swa_four b);
double swa_four_euclidean_norm(swa_four a);
int swa_four_is_finite(swa_four value);

int swa_lorentz_boost_matrix(
    swa_vec3 frame_velocity_mps,
    swa_tensor4 *out
);

int swa_lorentz_transform_four(
    swa_vec3 frame_velocity_mps,
    swa_four input,
    swa_four *out
);

int swa_lorentz_transform_tensor(
    swa_vec3 frame_velocity_mps,
    const swa_tensor4 *input,
    swa_tensor4 *out
);

int swa_em_field_tensor(
    const swa_em_fields *fields,
    swa_tensor4 *out
);

int swa_em_stress_energy_tensor(
    const swa_em_fields *fields,
    swa_tensor4 *out,
    swa_covariant_em_diagnostics *diagnostics
);

int swa_lorentz_transform_em_fields(
    swa_vec3 frame_velocity_mps,
    const swa_em_fields *input,
    swa_em_fields *out
);

int swa_four_current(
    double charge_density_C_m3,
    swa_vec3 current_density_A_m2,
    swa_four *out
);

int swa_lorentz_force_density(
    const swa_em_fields *fields,
    double charge_density_C_m3,
    swa_vec3 current_density_A_m2,
    swa_four *out_force_density_N_m3
);

int swa_particle_four_momentum(
    double rest_mass_kg,
    swa_vec3 velocity_mps,
    swa_four *out_momentum_Ns
);

int swa_audit_particle_mass_shell(
    double rest_mass_kg,
    swa_four momentum_Ns,
    double relative_tolerance,
    double absolute_tolerance_kg2_m2_s2,
    double *out_residual_kg2_m2_s2,
    double *out_relative_residual,
    int *out_passes
);

int swa_audit_four_momentum_ledger(
    const swa_four_momentum_ledger *ledger,
    double relative_tolerance,
    double absolute_tolerance_Ns,
    swa_four_momentum_ledger_result *out
);

int swa_write_covariant_em_receipt(
    FILE *fp,
    size_t checks,
    size_t failures,
    uint64_t deterministic_hash,
    double worst_field_invariant_relative,
    double worst_tensor_covariance_relative,
    double worst_four_force_covariance_relative,
    double worst_mass_shell_relative
);

#ifdef __cplusplus
}
#endif

#endif
