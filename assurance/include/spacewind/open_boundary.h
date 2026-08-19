#ifndef SPACEWIND_OPEN_BOUNDARY_H
#define SPACEWIND_OPEN_BOUNDARY_H

#include "spacewind/covariant_em.h"
#include "spacewind/field_ledger.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    swa_vec3 unit_normal;
    swa_four outward_four_momentum_flux_N_m2;
    double outward_power_W_m2;
    swa_vec3 outward_momentum_flux_N_m2;
    int finite;
} swa_em_boundary_flux;

typedef struct {
    swa_four outward_Ns;
    swa_four inward_Ns;
    double normal_velocity_mps;
    int is_outward;
    int is_inward;
    int finite;
} swa_particle_boundary_transport;

int swa_em_boundary_flux_density(
    const swa_em_fields *fields,
    swa_vec3 outward_normal,
    swa_em_boundary_flux *out
);

int swa_integrate_em_boundary_flux(
    const swa_em_boundary_flux *flux,
    double area_m2,
    double duration_s,
    swa_four *out_four_impulse_Ns
);

int swa_particle_boundary_crossing(
    double rest_mass_kg,
    swa_vec3 velocity_mps,
    double macro_weight,
    swa_vec3 outward_normal,
    swa_particle_boundary_transport *out
);

int swa_audit_open_boundary_four_momentum(
    swa_four field_initial_Ns,
    swa_four field_final_Ns,
    swa_four matter_initial_Ns,
    swa_four matter_final_Ns,
    swa_four external_impulse_Ns,
    swa_four field_outward_Ns,
    swa_four particle_outward_Ns,
    swa_four particle_inward_Ns,
    double relative_tolerance,
    double absolute_tolerance_Ns,
    swa_four_momentum_ledger_result *out
);

int swa_write_open_boundary_receipt(
    FILE *fp,
    const swa_em_boundary_flux *right_wave,
    const swa_em_boundary_flux *left_wave,
    double worst_traction_relative,
    double worst_control_volume_relative,
    size_t checks,
    size_t failures,
    uint64_t deterministic_hash
);

#ifdef __cplusplus
}
#endif

#endif
