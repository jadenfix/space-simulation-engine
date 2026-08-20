
#ifndef SPACEWIND_EM_PIC3D_H
#define SPACEWIND_EM_PIC3D_H

#include "spacewind/assurance.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    swa_vec3 position_unwrapped_m;
    swa_vec3 velocity_mps;
    double charge_C;
    double mass_kg;
    double macro_weight;
} swa_em_pic3d_particle;

typedef struct {
    size_t nx;
    size_t ny;
    size_t nz;
    size_t cell_count;
    size_t particle_count;
    double length_x_m;
    double length_y_m;
    double length_z_m;
    double dt_s;
    swa_em_pic3d_particle *particles;

    /* Yee electric-field / current faces:
       Ex/Jx: (i, j+1/2, k+1/2)
       Ey/Jy: (i+1/2, j, k+1/2)
       Ez/Jz: (i+1/2, j+1/2, k) */
    double *ex_V_m;
    double *ey_V_m;
    double *ez_V_m;

    /* Yee magnetic-field edges:
       Bx: (i+1/2, j, k)
       By: (i, j+1/2, k)
       Bz: (i, j, k+1/2) */
    double *bx_T;
    double *by_T;
    double *bz_T;
} swa_em_pic3d;

typedef struct {
    double continuity_relative_tolerance;
    double gauss_relative_tolerance;
    double magnetic_divergence_relative_tolerance;
    double local_current_difference_relative_tolerance;
    double spectral_oracle_relative_tolerance;
    double spectral_oracle_curl_relative_tolerance;
    double energy_relative_tolerance;
    double particle_work_relative_tolerance;
    double field_work_relative_tolerance;
    double momentum_relative_tolerance;
    double nonlinear_relative_tolerance;
    double maximum_courant;
    double maximum_particle_cells_per_step;
    size_t maximum_nonlinear_iterations;
} swa_em_pic3d_limits;

typedef struct {
    double raw_current_l2_A_m2;
    double correction_l2_A_m2;
    double correction_relative;
    double correction_curl_max_A_m3;
    double correction_curl_relative;
    double continuity_before_max_A_m3;
    double continuity_before_relative;
    double continuity_after_max_A_m3;
    double continuity_after_rms_A_m3;
    double continuity_after_relative;
    double zero_mode_A_m3;
    int finite;
    int zero_mode_passes;
    int continuity_passes;
    int correction_curl_passes;
    int passes;
} swa_em_pic3d_projection_result;

typedef struct {
    double courant_number;
    double maximum_particle_displacement_cells;
    size_t maximum_path_segments;
    size_t nonlinear_iterations;
    double nonlinear_relative_residual;

    double initial_charge_C;
    double final_charge_C;
    double charge_change_C;

    double continuity_residual_max_A_m3;
    double continuity_residual_rms_A_m3;
    double continuity_relative_max;
    double initial_gauss_residual_max_C_m3;
    double final_gauss_residual_max_C_m3;
    double initial_gauss_relative_max;
    double final_gauss_relative_max;
    double initial_magnetic_divergence_max_T_m;
    double final_magnetic_divergence_max_T_m;
    double initial_magnetic_divergence_relative;
    double final_magnetic_divergence_relative;

    double raw_current_l2_A_m2;
    double current_correction_l2_A_m2;
    double current_correction_relative;
    double local_current_difference_relative;
    double spectral_oracle_difference_relative;
    double current_projection_curl_relative;
    double current_projection_zero_mode_A_m3;

    swa_vec3 mean_current_A_m2;
    swa_vec3 transport_mean_current_A_m2;
    swa_vec3 mean_current_residual_A_m2;

    double initial_particle_energy_J;
    double final_particle_energy_J;
    double initial_field_energy_J;
    double final_field_energy_J;
    double particle_energy_change_J;
    double field_energy_change_J;
    double raw_current_work_J;
    double corrected_current_work_J;
    double particle_midpoint_current_work_J;
    double field_midpoint_current_work_J;
    double projection_work_correction_J;
    double particle_work_residual_J;
    double field_work_residual_J;
    double total_energy_residual_J;
    double total_energy_relative;
    double particle_work_relative;
    double field_work_relative;

    swa_vec3 initial_particle_momentum_Ns;
    swa_vec3 final_particle_momentum_Ns;
    swa_vec3 initial_field_momentum_Ns;
    swa_vec3 final_field_momentum_Ns;
    swa_vec3 total_momentum_residual_Ns;
    double total_momentum_relative;

    double maximum_speed_fraction_c;
    int finite;
    int courant_passes;
    int particle_courant_passes;
    int charge_closure_passes;
    int continuity_passes;
    int initial_gauss_passes;
    int final_gauss_passes;
    int initial_magnetic_divergence_passes;
    int final_magnetic_divergence_passes;
    int current_projection_passes;
    int local_current_passes;
    int spectral_oracle_passes;
    int mean_current_passes;
    int nonlinear_passes;
    int subluminal_passes;
    int energy_passes;
    int particle_work_passes;
    int field_work_passes;
    int momentum_passes;
    int passes;
} swa_em_pic3d_result;

int swa_em_pic3d_init(
    swa_em_pic3d *state,
    size_t nx,
    size_t ny,
    size_t nz,
    size_t particle_count,
    double length_x_m,
    double length_y_m,
    double length_z_m,
    double dt_s
);

void swa_em_pic3d_destroy(swa_em_pic3d *state);

int swa_em_pic3d_state_is_finite(const swa_em_pic3d *state);

int swa_em_pic3d_set_particle(
    swa_em_pic3d *state,
    size_t index,
    swa_vec3 position_unwrapped_m,
    swa_vec3 velocity_mps,
    double charge_C,
    double mass_kg,
    double macro_weight
);

int swa_em_pic3d_higuera_cary_push(
    swa_vec3 velocity_initial_mps,
    double charge_C,
    double mass_kg,
    swa_vec3 electric_field_V_m,
    swa_vec3 magnetic_field_T,
    double dt_s,
    swa_vec3 *velocity_final_mps
);

int swa_em_pic3d_deposit_charge_density(
    const swa_em_pic3d *state,
    const swa_em_pic3d_particle *particles,
    size_t particle_count,
    double *rho_C_m3
);

int swa_em_pic3d_deposit_charge_conserving_current(
    const swa_em_pic3d *state,
    const swa_em_pic3d_particle *initial_particles,
    const swa_em_pic3d_particle *final_particles,
    size_t particle_count,
    double *jx_A_m2,
    double *jy_A_m2,
    double *jz_A_m2,
    swa_vec3 *transport_mean_current_A_m2
);

int swa_em_pic3d_deposit_transport_current(
    const swa_em_pic3d *state,
    const swa_em_pic3d_particle *initial_particles,
    const swa_em_pic3d_particle *final_particles,
    size_t particle_count,
    double *raw_jx_A_m2,
    double *raw_jy_A_m2,
    double *raw_jz_A_m2,
    size_t *maximum_path_segments,
    swa_vec3 *transport_mean_current_A_m2
);

int swa_em_pic3d_project_periodic_current(
    const swa_em_pic3d *state,
    const double *rho_initial_C_m3,
    const double *rho_final_C_m3,
    const double *raw_jx_A_m2,
    const double *raw_jy_A_m2,
    const double *raw_jz_A_m2,
    double relative_tolerance,
    double curl_relative_tolerance,
    double *corrected_jx_A_m2,
    double *corrected_jy_A_m2,
    double *corrected_jz_A_m2,
    swa_em_pic3d_projection_result *out
);

int swa_em_pic3d_initialize_gauss_field(
    swa_em_pic3d *state,
    swa_vec3 requested_mean_electric_field_V_m,
    double absolute_charge_tolerance_C
);

double swa_em_pic3d_field_energy_J(const swa_em_pic3d *state);
swa_vec3 swa_em_pic3d_field_momentum_Ns(const swa_em_pic3d *state);

double swa_em_pic3d_electric_divergence_max_C_m3(
    const swa_em_pic3d *state,
    const double *rho_C_m3,
    double *relative
);

double swa_em_pic3d_magnetic_divergence_max_T_m(
    const swa_em_pic3d *state,
    double *relative
);

int swa_em_pic3d_step(
    swa_em_pic3d *state,
    const swa_em_pic3d_limits *limits,
    swa_em_pic3d_result *out
);

int swa_em_pic3d_write_receipt(
    FILE *fp,
    const swa_em_pic3d *state,
    const swa_em_pic3d_limits *limits,
    const swa_em_pic3d_result *result
);

#ifdef __cplusplus
}
#endif

#endif
