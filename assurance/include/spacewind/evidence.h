#ifndef SPACEWIND_EVIDENCE_H
#define SPACEWIND_EVIDENCE_H

#include "spacewind/assurance.h"
#include "spacewind/experiment.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    SWA_CLAIM_NONE = 0,
    SWA_CLAIM_SOFTWARE_FINITE_CHECK = 1,
    SWA_CLAIM_CONVERGED_MODEL_RESULT = 2,
    SWA_CLAIM_CHAMBER_LOCAL_FORCE = 3,
    SWA_CLAIM_CLOSED_CYCLE_CANDIDATE = 4,
    SWA_CLAIM_FLIGHT_PROPULSION = 5
} swa_claim_tier;

enum {
    SWA_BLOCK_DIMENSIONS = UINT64_C(1) << 0,
    SWA_BLOCK_CONSERVATION = UINT64_C(1) << 1,
    SWA_BLOCK_DETERMINISM = UINT64_C(1) << 2,
    SWA_BLOCK_MEMORY_SAFETY = UINT64_C(1) << 3,
    SWA_BLOCK_CONVERGENCE = UINT64_C(1) << 4,
    SWA_BLOCK_INDEPENDENT_IMPLEMENTATION = UINT64_C(1) << 5,
    SWA_BLOCK_SIMILARITY = UINT64_C(1) << 6,
    SWA_BLOCK_POWER = UINT64_C(1) << 7,
    SWA_BLOCK_THERMAL = UINT64_C(1) << 8,
    SWA_BLOCK_STRUCTURE = UINT64_C(1) << 9,
    SWA_BLOCK_FAULT_CONTAINMENT = UINT64_C(1) << 10,
    SWA_BLOCK_REVERSAL_EXPERIMENT = UINT64_C(1) << 11,
    SWA_BLOCK_INDEPENDENT_REPLICATION = UINT64_C(1) << 12,
    SWA_BLOCK_CLOSED_CYCLE = UINT64_C(1) << 13,
    SWA_BLOCK_ROBUST_UNCERTAINTY = UINT64_C(1) << 14,
    SWA_BLOCK_DEPLOYMENT = UINT64_C(1) << 15,
    SWA_BLOCK_FLIGHT_EVIDENCE = UINT64_C(1) << 16
};

typedef struct {
    int dimensional_contract_passes;
    int conservation_ledgers_pass;
    int deterministic_replay_passes;
    int memory_and_undefined_behavior_checks_pass;
    int convergence_study_passes;
    int independent_implementation_agrees;
    int similarity_contract_passes;
    double available_to_required_power_ratio;
    double thermal_margin_fraction;
    int structural_envelope_passes;
    int fault_containment_passes;
    const swa_reversal_result *reversal_experiment;
    const swa_replication_result *independent_replication;
    const swa_cycle_result *closed_cycle;
    swa_interval mission_net_gain_interval_J;
    double deployment_survival_lower_bound;
    double required_deployment_survival_lower_bound;
    int flight_demonstration_passes;
} swa_claim_gate_input;

typedef struct {
    swa_claim_tier highest_tier;
    uint64_t blockers;
    int software_finite_check_ready;
    int converged_model_ready;
    int chamber_local_force_ready;
    int closed_cycle_candidate_ready;
    int flight_propulsion_ready;
    int propulsion_claim_ready;
} swa_claim_gate_result;

const char *swa_claim_tier_name(swa_claim_tier tier);
int swa_evaluate_claim_gate(const swa_claim_gate_input *input,
                            swa_claim_gate_result *out);
int swa_write_claim_gate_receipt(FILE *fp,
                                 const swa_claim_gate_input *input,
                                 const swa_claim_gate_result *result);

#ifdef __cplusplus
}
#endif

#endif
