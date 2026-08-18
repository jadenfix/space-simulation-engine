#include "spacewind/evidence.h"

#include <math.h>
#include <string.h>

const char *swa_claim_tier_name(swa_claim_tier tier) {
    switch (tier) {
        case SWA_CLAIM_SOFTWARE_FINITE_CHECK: return "software_finite_check";
        case SWA_CLAIM_CONVERGED_MODEL_RESULT: return "converged_model_result";
        case SWA_CLAIM_CHAMBER_LOCAL_FORCE: return "chamber_local_force";
        case SWA_CLAIM_CLOSED_CYCLE_CANDIDATE: return "closed_cycle_candidate";
        case SWA_CLAIM_FLIGHT_PROPULSION: return "flight_propulsion";
        case SWA_CLAIM_NONE:
        default: return "none";
    }
}

int swa_evaluate_claim_gate(const swa_claim_gate_input *i,
                            swa_claim_gate_result *o) {
    int robust_mission_gain;
    if (i == NULL || o == NULL ||
        !isfinite(i->available_to_required_power_ratio) ||
        !isfinite(i->thermal_margin_fraction) ||
        !isfinite(i->mission_net_gain_interval_J.lo) ||
        !isfinite(i->mission_net_gain_interval_J.hi) ||
        !isfinite(i->deployment_survival_lower_bound) ||
        !isfinite(i->required_deployment_survival_lower_bound)) {
        return 0;
    }
    memset(o, 0, sizeof(*o));
    if (!i->dimensional_contract_passes) o->blockers |= SWA_BLOCK_DIMENSIONS;
    if (!i->conservation_ledgers_pass) o->blockers |= SWA_BLOCK_CONSERVATION;
    if (!i->deterministic_replay_passes) o->blockers |= SWA_BLOCK_DETERMINISM;
    if (!i->memory_and_undefined_behavior_checks_pass) o->blockers |= SWA_BLOCK_MEMORY_SAFETY;
    if (!i->convergence_study_passes) o->blockers |= SWA_BLOCK_CONVERGENCE;
    if (!i->independent_implementation_agrees) o->blockers |= SWA_BLOCK_INDEPENDENT_IMPLEMENTATION;
    if (!i->similarity_contract_passes) o->blockers |= SWA_BLOCK_SIMILARITY;
    if (!(i->available_to_required_power_ratio >= 1.0)) o->blockers |= SWA_BLOCK_POWER;
    if (!(i->thermal_margin_fraction > 0.0)) o->blockers |= SWA_BLOCK_THERMAL;
    if (!i->structural_envelope_passes) o->blockers |= SWA_BLOCK_STRUCTURE;
    if (!i->fault_containment_passes) o->blockers |= SWA_BLOCK_FAULT_CONTAINMENT;
    if (i->reversal_experiment == NULL || !i->reversal_experiment->passes)
        o->blockers |= SWA_BLOCK_REVERSAL_EXPERIMENT;
    if (i->independent_replication == NULL || !i->independent_replication->passes)
        o->blockers |= SWA_BLOCK_INDEPENDENT_REPLICATION;
    if (i->closed_cycle == NULL || !i->closed_cycle->numerical_cycle_passes)
        o->blockers |= SWA_BLOCK_CLOSED_CYCLE;
    robust_mission_gain = i->mission_net_gain_interval_J.lo > 0.0;
    if (!robust_mission_gain || i->closed_cycle == NULL ||
        !i->closed_cycle->robust_positive_gain)
        o->blockers |= SWA_BLOCK_ROBUST_UNCERTAINTY;
    if (!(i->required_deployment_survival_lower_bound >= 0.0 &&
          i->required_deployment_survival_lower_bound <= 1.0) ||
        !(i->deployment_survival_lower_bound >=
          i->required_deployment_survival_lower_bound))
        o->blockers |= SWA_BLOCK_DEPLOYMENT;
    if (!i->flight_demonstration_passes) o->blockers |= SWA_BLOCK_FLIGHT_EVIDENCE;

    o->software_finite_check_ready =
        (o->blockers & (SWA_BLOCK_DIMENSIONS | SWA_BLOCK_CONSERVATION |
                        SWA_BLOCK_DETERMINISM | SWA_BLOCK_MEMORY_SAFETY)) == 0U;
    o->converged_model_ready = o->software_finite_check_ready &&
        (o->blockers & (SWA_BLOCK_CONVERGENCE |
                        SWA_BLOCK_INDEPENDENT_IMPLEMENTATION)) == 0U;
    o->chamber_local_force_ready = o->converged_model_ready &&
        (o->blockers & (SWA_BLOCK_SIMILARITY | SWA_BLOCK_POWER |
                        SWA_BLOCK_THERMAL | SWA_BLOCK_STRUCTURE |
                        SWA_BLOCK_FAULT_CONTAINMENT |
                        SWA_BLOCK_REVERSAL_EXPERIMENT |
                        SWA_BLOCK_INDEPENDENT_REPLICATION)) == 0U;
    o->closed_cycle_candidate_ready = o->chamber_local_force_ready &&
        (o->blockers & (SWA_BLOCK_CLOSED_CYCLE |
                        SWA_BLOCK_ROBUST_UNCERTAINTY)) == 0U;
    o->flight_propulsion_ready = o->closed_cycle_candidate_ready &&
        (o->blockers & (SWA_BLOCK_DEPLOYMENT |
                        SWA_BLOCK_FLIGHT_EVIDENCE)) == 0U;
    o->propulsion_claim_ready = o->flight_propulsion_ready;

    if (o->flight_propulsion_ready) o->highest_tier = SWA_CLAIM_FLIGHT_PROPULSION;
    else if (o->closed_cycle_candidate_ready) o->highest_tier = SWA_CLAIM_CLOSED_CYCLE_CANDIDATE;
    else if (o->chamber_local_force_ready) o->highest_tier = SWA_CLAIM_CHAMBER_LOCAL_FORCE;
    else if (o->converged_model_ready) o->highest_tier = SWA_CLAIM_CONVERGED_MODEL_RESULT;
    else if (o->software_finite_check_ready) o->highest_tier = SWA_CLAIM_SOFTWARE_FINITE_CHECK;
    else o->highest_tier = SWA_CLAIM_NONE;
    return 1;
}

int swa_write_claim_gate_receipt(FILE *fp,
                                 const swa_claim_gate_input *i,
                                 const swa_claim_gate_result *r) {
    if (fp == NULL || i == NULL || r == NULL) {
        return 0;
    }
    return fprintf(fp,
        "{\n  \"schema\": \"spacewind.claim-gate/v1\",\n"
        "  \"highest_tier\": \"%s\",\n"
        "  \"blockers_hex\": \"%016llx\",\n"
        "  \"software_finite_check_ready\": %s,\n"
        "  \"converged_model_ready\": %s,\n"
        "  \"chamber_local_force_ready\": %s,\n"
        "  \"closed_cycle_candidate_ready\": %s,\n"
        "  \"flight_propulsion_ready\": %s,\n"
        "  \"propulsion_claim_ready\": %s,\n"
        "  \"power_margin_ratio\": %.17g,\n"
        "  \"thermal_margin_fraction\": %.17g,\n"
        "  \"mission_net_gain_interval_J\": [%.17g, %.17g],\n"
        "  \"deployment_survival_lower_bound\": %.17g,\n"
        "  \"required_deployment_survival_lower_bound\": %.17g,\n"
        "  \"nonclaim\": \"a lower-tier result must not be worded as a higher-tier physical claim\"\n}\n",
        swa_claim_tier_name(r->highest_tier),
        (unsigned long long)r->blockers,
        r->software_finite_check_ready ? "true" : "false",
        r->converged_model_ready ? "true" : "false",
        r->chamber_local_force_ready ? "true" : "false",
        r->closed_cycle_candidate_ready ? "true" : "false",
        r->flight_propulsion_ready ? "true" : "false",
        r->propulsion_claim_ready ? "true" : "false",
        i->available_to_required_power_ratio,
        i->thermal_margin_fraction,
        i->mission_net_gain_interval_J.lo,
        i->mission_net_gain_interval_J.hi,
        i->deployment_survival_lower_bound,
        i->required_deployment_survival_lower_bound) > 0;
}
