#ifndef SPACEWIND_CHARGE_DEPOSITION_H
#define SPACEWIND_CHARGE_DEPOSITION_H

#include "spacewind/field_ledger.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t cells;
    size_t particles;
    double domain_length_m;
    double cross_section_area_m2;
    double dt_s;
    const double *initial_position_unwrapped_m;
    const double *final_position_unwrapped_m;
    const double *particle_charge_C;
    double continuity_relative_tolerance;
    double continuity_absolute_tolerance_A_m3;
    double continuity_global_absolute_tolerance_C_m;
    double charge_absolute_tolerance_C;
    double mean_current_relative_tolerance;
    double mean_current_absolute_tolerance_A_m2;
} swa_periodic_cic_transport_input;

typedef struct {
    double *rho_initial_C_m3;
    double *rho_final_C_m3;
    size_t rho_capacity;
    double *current_face_A_m2;
    size_t current_face_capacity;
} swa_periodic_cic_transport_buffers;

typedef struct {
    double particle_charge_sum_C;
    double deposited_initial_charge_C;
    double deposited_final_charge_C;
    double initial_charge_error_C;
    double final_charge_error_C;
    double transported_charge_turns_C;
    double target_mean_current_A_m2;
    double deposited_mean_current_A_m2;
    double mean_current_error_A_m2;
    double maximum_particle_displacement_cells;
    swa_charge_continuity_result continuity;
    int charge_closure_passes;
    int mean_current_passes;
    int finite;
    int passes;
} swa_periodic_cic_transport_result;

int swa_deposit_periodic_cic_transport(
    const swa_periodic_cic_transport_input *input,
    const swa_periodic_cic_transport_buffers *buffers,
    swa_periodic_cic_transport_result *out
);

int swa_write_periodic_cic_transport_receipt(
    FILE *fp,
    const swa_periodic_cic_transport_input *input,
    const swa_periodic_cic_transport_result *result
);

#ifdef __cplusplus
}
#endif

#endif
