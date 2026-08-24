#ifndef SPACEWIND_EM_PIC1D_H
#define SPACEWIND_EM_PIC1D_H

#include "spacewind/assurance.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double position_unwrapped_m;
    swa_vec3 velocity_mps;
    double charge_C;
    double mass_kg;
    double macro_weight;
} swa_em_pic1d_particle;

typedef struct {
    size_t cells;
    size_t particle_count;
    double domain_length_m;
    double cross_section_area_m2;
    double dt_s;
    double bx_uniform_T;
    swa_em_pic1d_particle *particles;
    double *ex_face_V_m;
    double *ey_face_V_m;
    double *ez_face_V_m;
    double *by_cell_T;
    double *bz_cell_T;
} swa_em_pic1d;

typedef struct {
    double continuity_relative_tolerance;
    double gauss_relative_tolerance;
    double magnetic_divergence_relative_tolerance;
    double energy_relative_tolerance;
    double momentum_relative_tolerance;
    double maximum_courant;
    double maximum_particle_cells_per_step;
} swa_em_pic1d_limits;

typedef struct {
    double courant_number;
    double maximum_particle_displacement_cells;
    size_t maximum_path_segments;
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
    double magnetic_divergence_residual_T_m;
    double mean_jx_A_m2;
    double mean_jy_A_m2;
    double mean_jz_A_m2;
    double transverse_current_partition_relative_max;
    double initial_mean_ex_V_m;
    double final_mean_ex_V_m;
    double harmonic_ampere_residual_V_m;
    double harmonic_ampere_relative;
    double initial_particle_energy_J;
    double final_particle_energy_J;
    double initial_field_energy_J;
    double final_field_energy_J;
    double particle_energy_change_J;
    double field_energy_change_J;
    double midpoint_current_work_J;
    double particle_work_residual_J;
    double field_work_residual_J;
    double total_energy_residual_J;
    double total_energy_relative;
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
    int magnetic_divergence_passes;
    int transverse_current_partition_passes;
    int harmonic_ampere_passes;
    int subluminal_passes;
    int energy_passes;
    int momentum_passes;
    int passes;
} swa_em_pic1d_result;

int swa_em_pic1d_init(
    swa_em_pic1d *state,
    size_t cells,
    size_t particle_count,
    double domain_length_m,
    double cross_section_area_m2,
    double dt_s
);

void swa_em_pic1d_destroy(swa_em_pic1d *state);

int swa_em_pic1d_set_particle(
    swa_em_pic1d *state,
    size_t index,
    double position_unwrapped_m,
    swa_vec3 velocity_mps,
    double charge_C,
    double mass_kg,
    double macro_weight
);

int swa_em_pic1d_initialize_gauss_field(
    swa_em_pic1d *state,
    double requested_mean_ex_V_m,
    double absolute_charge_tolerance_C
);

int swa_em_pic1d_step(
    swa_em_pic1d *state,
    const swa_em_pic1d_limits *limits,
    swa_em_pic1d_result *out
);

int swa_em_pic1d_write_receipt(
    FILE *fp,
    const swa_em_pic1d *state,
    const swa_em_pic1d_limits *limits,
    const swa_em_pic1d_result *result
);

#ifdef __cplusplus
}
#endif

#endif
