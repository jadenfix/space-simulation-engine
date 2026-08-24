#ifndef SPACEWIND_FOUR_MOMENTUM_DELTA_H
#define SPACEWIND_FOUR_MOMENTUM_DELTA_H

#include "spacewind/covariant_em.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    swa_four field_change_Ns;
    swa_four matter_change_Ns;
    swa_four external_impulse_Ns;
    swa_four momentum_in_Ns;
    swa_four momentum_out_Ns;
} swa_four_momentum_delta_ledger;

typedef struct {
    swa_four residual_Ns;
    double residual_norm_Ns;
    double physical_scale_Ns;
    double relative_residual;
    int finite;
    int passes;
} swa_four_momentum_delta_result;

int swa_audit_four_momentum_delta_ledger(
    const swa_four_momentum_delta_ledger *ledger,
    double physical_scale_floor_Ns,
    double relative_tolerance,
    double absolute_tolerance_Ns,
    swa_four_momentum_delta_result *out
);

int swa_lorentz_transform_four_momentum_delta_ledger(
    swa_vec3 frame_velocity_mps,
    const swa_four_momentum_delta_ledger *input,
    swa_four_momentum_delta_ledger *output
);

#ifdef __cplusplus
}
#endif

#endif
