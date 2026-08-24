#include "spacewind/evidence.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t receipt_hash = UINT64_C(14695981039346656037);

static void check_true(int condition, const char *name) {
    ++checks;
    receipt_hash ^= swa_fnv1a64(name, strlen(name));
    receipt_hash *= UINT64_C(1099511628211);
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}

static swa_claim_gate_input empty_gate(void) {
    swa_claim_gate_input input;
    memset(&input, 0, sizeof(input));
    input.mission_net_gain_interval_J = swa_interval_make(-1.0, 1.0);
    input.required_deployment_survival_lower_bound = 0.99;
    return input;
}

static void fill_software(swa_claim_gate_input *i) {
    i->dimensional_contract_passes = 1;
    i->conservation_ledgers_pass = 1;
    i->field_particle_ledgers_pass = 1;
    i->deterministic_replay_passes = 1;
    i->memory_and_undefined_behavior_checks_pass = 1;
}

static void fill_model(swa_claim_gate_input *i) {
    fill_software(i);
    i->convergence_study_passes = 1;
    i->independent_implementation_agrees = 1;
}

static void fill_chamber(swa_claim_gate_input *i,
                         swa_reversal_result *reversal,
                         swa_replication_result *replication) {
    fill_model(i);
    i->similarity_contract_passes = 1;
    i->available_to_required_power_ratio = 1.5;
    i->thermal_margin_fraction = 0.25;
    i->structural_envelope_passes = 1;
    i->fault_containment_passes = 1;
    reversal->finite = 1;
    reversal->passes = 1;
    reversal->reversed_effect_mean_N = 1.0;
    reversal->effect_interval_N = swa_interval_make(0.8, 1.2);
    replication->passing_replicates = 2U;
    replication->same_sign = 1;
    replication->interval_overlap = 1;
    replication->pooled_effect_N = 1.0;
    replication->relative_heterogeneity = 0.1;
    replication->passes = 1;
    i->reversal_experiment = reversal;
    i->independent_replication = replication;
}

static void fill_cycle(swa_claim_gate_input *i, swa_cycle_result *cycle) {
    cycle->net_gain_J = 15.0;
    cycle->net_gain_interval_J = swa_interval_make(10.0, 20.0);
    cycle->numerical_cycle_passes = 1;
    cycle->robust_positive_gain = 1;
    cycle->physical_promotion_passes = 1;
    i->closed_cycle = cycle;
    i->mission_net_gain_interval_J = swa_interval_make(10.0, 20.0);
}

static void test_tier_progression(void) {
    swa_claim_gate_input i = empty_gate();
    swa_claim_gate_result r;
    swa_reversal_result reversal;
    swa_replication_result replication;
    swa_cycle_result cycle;
    memset(&reversal, 0, sizeof(reversal));
    memset(&replication, 0, sizeof(replication));
    memset(&cycle, 0, sizeof(cycle));

    check_true(swa_evaluate_claim_gate(&i, &r), "empty claim gate evaluates");
    check_true(r.highest_tier == SWA_CLAIM_NONE, "empty evidence yields no claim");
    check_true(!r.propulsion_claim_ready, "empty evidence cannot claim propulsion");

    fill_software(&i);
    check_true(swa_evaluate_claim_gate(&i, &r), "software gate evaluates");
    check_true(r.highest_tier == SWA_CLAIM_SOFTWARE_FINITE_CHECK,
               "software evidence stops at finite-check tier");

    fill_model(&i);
    check_true(swa_evaluate_claim_gate(&i, &r), "model gate evaluates");
    check_true(r.highest_tier == SWA_CLAIM_CONVERGED_MODEL_RESULT,
               "convergence and independent implementation permit model tier");

    fill_chamber(&i, &reversal, &replication);
    check_true(swa_evaluate_claim_gate(&i, &r), "chamber gate evaluates");
    check_true(r.highest_tier == SWA_CLAIM_CHAMBER_LOCAL_FORCE,
               "reversal and replication permit local-force tier");
    check_true(!r.propulsion_claim_ready,
               "local chamber force is not a propulsion claim");

    fill_cycle(&i, &cycle);
    check_true(swa_evaluate_claim_gate(&i, &r), "cycle gate evaluates");
    check_true(r.highest_tier == SWA_CLAIM_CLOSED_CYCLE_CANDIDATE,
               "robust closed cycle permits candidate tier");
    check_true(!r.flight_propulsion_ready,
               "closed-cycle candidate still lacks flight qualification");

    i.deployment_survival_lower_bound = 0.995;
    i.flight_demonstration_passes = 1;
    check_true(swa_evaluate_claim_gate(&i, &r), "flight gate evaluates");
    check_true(r.highest_tier == SWA_CLAIM_FLIGHT_PROPULSION,
               "flight plus deployment evidence permits top tier");
    check_true(r.propulsion_claim_ready,
               "propulsion claim only ready at top tier");
}

static void test_fail_closed_blockers(void) {
    swa_claim_gate_input i = empty_gate();
    swa_claim_gate_result r;
    swa_reversal_result reversal;
    swa_replication_result replication;
    swa_cycle_result cycle;
    memset(&reversal, 0, sizeof(reversal));
    memset(&replication, 0, sizeof(replication));
    memset(&cycle, 0, sizeof(cycle));
    fill_chamber(&i, &reversal, &replication);
    fill_cycle(&i, &cycle);
    i.deployment_survival_lower_bound = 0.999;
    i.flight_demonstration_passes = 1;

    i.available_to_required_power_ratio = 0.99;
    check_true(swa_evaluate_claim_gate(&i, &r) &&
               (r.blockers & SWA_BLOCK_POWER) != 0U,
               "power deficit blocks chamber promotion");
    check_true(r.highest_tier == SWA_CLAIM_CONVERGED_MODEL_RESULT,
               "power deficit falls back to model tier");
    i.available_to_required_power_ratio = 1.5;

    i.mission_net_gain_interval_J = swa_interval_make(-0.1, 20.0);
    check_true(swa_evaluate_claim_gate(&i, &r) &&
               (r.blockers & SWA_BLOCK_ROBUST_UNCERTAINTY) != 0U,
               "uncertainty crossing zero blocks closed cycle");
    check_true(r.highest_tier == SWA_CLAIM_CHAMBER_LOCAL_FORCE,
               "uncertain cycle falls back to chamber tier");
    i.mission_net_gain_interval_J = swa_interval_make(10.0, 20.0);

    replication.passes = 0;
    check_true(swa_evaluate_claim_gate(&i, &r) &&
               (r.blockers & SWA_BLOCK_INDEPENDENT_REPLICATION) != 0U,
               "failed replication blocks chamber claim");
    replication.passes = 1;

    replication.passing_replicates = 1U;
    check_true(swa_evaluate_claim_gate(&i, &r) &&
               (r.blockers & SWA_BLOCK_INDEPENDENT_REPLICATION) != 0U,
               "one replicate is not independent replication");
    replication.passing_replicates = 2U;

    cycle.physical_promotion_passes = 0;
    check_true(swa_evaluate_claim_gate(&i, &r) &&
               (r.blockers & SWA_BLOCK_CLOSED_CYCLE) != 0U,
               "numerical cycle without physical promotion is blocked");
    cycle.physical_promotion_passes = 1;

    cycle.net_gain_J = 30.0;
    check_true(swa_evaluate_claim_gate(&i, &r) &&
               (r.blockers & SWA_BLOCK_CLOSED_CYCLE) != 0U,
               "cycle nominal gain outside its interval is rejected");
    cycle.net_gain_J = 15.0;

    i.field_particle_ledgers_pass = 0;
    check_true(swa_evaluate_claim_gate(&i, &r) &&
               (r.blockers & SWA_BLOCK_FIELD_PARTICLE_LEDGER) != 0U,
               "missing field-particle ledger blocks software promotion");
    check_true(r.highest_tier == SWA_CLAIM_NONE,
               "missing field-particle ledger yields no promotable claim");
    i.field_particle_ledgers_pass = 1;

    reversal.reversed_effect_mean_N = NAN;
    check_true(swa_evaluate_claim_gate(&i, &r) &&
               (r.blockers & SWA_BLOCK_REVERSAL_EXPERIMENT) != 0U,
               "nonfinite reversal evidence fails closed");
    reversal.reversed_effect_mean_N = 1.0;

    i.dimensional_contract_passes = 0;
    check_true(swa_evaluate_claim_gate(&i, &r) &&
               (r.blockers & SWA_BLOCK_DIMENSIONS) != 0U,
               "dimension failure blocks even software tier");
    check_true(r.highest_tier == SWA_CLAIM_NONE,
               "dimension failure yields no promotable claim");
}

static void test_invalid_inputs(void) {
    swa_claim_gate_input i = empty_gate();
    swa_claim_gate_result r;
    i.mission_net_gain_interval_J.lo = 2.0;
    i.mission_net_gain_interval_J.hi = 1.0;
    check_true(!swa_evaluate_claim_gate(&i, &r),
               "inverted mission interval rejected as invalid input");
    i = empty_gate();
    i.available_to_required_power_ratio = NAN;
    check_true(!swa_evaluate_claim_gate(&i, &r),
               "nonfinite power ratio rejected as invalid input");
}

static void test_current_project_receipt(void) {
    swa_claim_gate_input current = empty_gate();
    swa_claim_gate_result result;
    FILE *fp;
    fill_model(&current);
    current.similarity_contract_passes = 0;
    current.available_to_required_power_ratio = 0.0;
    current.thermal_margin_fraction = 0.0;
    current.structural_envelope_passes = 0;
    current.fault_containment_passes = 1;
    current.deployment_survival_lower_bound = 0.0;
    current.flight_demonstration_passes = 0;
    check_true(swa_evaluate_claim_gate(&current, &result),
               "current project claim gate evaluates");
    check_true(result.highest_tier == SWA_CLAIM_CONVERGED_MODEL_RESULT,
               "current project receipt remains a model result");
    check_true(!result.propulsion_claim_ready,
               "current project receipt explicitly denies propulsion readiness");
    fp = fopen("output/current_claim_gate.json", "w");
    check_true(fp != NULL, "open current claim gate receipt");
    if (fp != NULL) {
        check_true(swa_write_claim_gate_receipt(fp, &current, &result),
                   "write current claim gate receipt");
        (void)fclose(fp);
    }
}

int main(void) {
    FILE *fp;
    test_tier_progression();
    test_fail_closed_blockers();
    test_invalid_inputs();
    test_current_project_receipt();
    fp = fopen("output/evidence_assurance_receipt.json", "w");
    if (fp != NULL) {
        (void)swa_write_assurance_receipt(fp, checks, failures, receipt_hash);
        (void)fclose(fp);
    } else {
        ++failures;
    }
    printf("spacewind evidence gates: %zu checks, %zu failures, hash=%016llx\n",
           checks, failures, (unsigned long long)receipt_hash);
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
