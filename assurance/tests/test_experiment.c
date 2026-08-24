#include "spacewind/experiment.h"

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

static void fill_clean_experiment(double *plus, double *minus, double *zero,
                                  size_t n, double effect, double bias,
                                  double phase) {
    size_t i;
    for (i = 0U; i < n; ++i) {
        const double q1 = swa_halton(i + 1U, 2U) - 0.5;
        const double q2 = swa_halton(i + 1U, 3U) - 0.5;
        const double q3 = swa_halton(i + 1U, 5U) - 0.5;
        plus[i] = bias + effect + 0.04 * q1 + 0.01 * sin((double)i + phase);
        minus[i] = bias - effect + 0.04 * q2 - 0.01 * cos((double)i + phase);
        zero[i] = bias + 0.02 * q3;
    }
}

static swa_reversal_experiment config(const double *plus, const double *minus,
                                      const double *zero, size_t n) {
    swa_reversal_experiment e;
    e.positive_command_force_N = plus;
    e.negative_command_force_N = minus;
    e.zero_command_force_N = zero;
    e.count = n;
    e.alpha = 0.01;
    e.minimum_effect_N = 1.0;
    e.maximum_common_mode_bias_N = 0.1;
    e.maximum_reversal_asymmetry = 0.05;
    e.maximum_lag1_autocorrelation = 0.70;
    e.maximum_zero_drift_N_per_sample = 0.01;
    return e;
}

static void test_clean_reversal(void) {
    enum { N = 12 };
    double plus[N], minus[N], zero[N];
    swa_reversal_experiment e;
    swa_reversal_result r;
    FILE *fp;
    fill_clean_experiment(plus, minus, zero, N, 10.0, 2.0, 0.0);
    e = config(plus, minus, zero, N);
    check_true(swa_audit_reversal_experiment(&e, &r), "clean reversal audit executes");
    check_true(r.finite, "clean reversal outputs finite");
    check_true(r.effect_interval_N.lo > 1.0, "clean reversal lower bound exceeds minimum");
    check_true(r.exact_or_deterministic_sign_flip_p <= e.alpha,
               "clean reversal exact sign-flip test significant");
    check_true(r.common_mode_controlled, "clean reversal common mode controlled");
    check_true(r.reversal_symmetric, "clean reversal symmetry controlled");
    check_true(r.autocorrelation_controlled, "clean reversal autocorrelation controlled");
    check_true(r.drift_controlled, "clean reversal drift controlled");
    check_true(r.passes, "clean reversal passes all gates");
    fp = fopen("output/reversal_receipt.json", "w");
    check_true(fp != NULL, "open reversal receipt");
    if (fp != NULL) {
        check_true(swa_write_reversal_receipt(fp, &e, &r), "write reversal receipt");
        (void)fclose(fp);
    }
}

static void test_null_effect_rejected(void) {
    enum { N = 12 };
    double plus[N], minus[N], zero[N];
    swa_reversal_experiment e;
    swa_reversal_result r;
    fill_clean_experiment(plus, minus, zero, N, 0.0, 2.0, 0.3);
    e = config(plus, minus, zero, N);
    check_true(swa_audit_reversal_experiment(&e, &r), "null reversal audit executes");
    check_true(!r.effect_above_minimum, "null effect fails minimum effect gate");
    check_true(!r.passes, "null effect cannot pass reversal protocol");
}

static void test_common_mode_rejected(void) {
    enum { N = 12 };
    double plus[N], minus[N], zero[N];
    swa_reversal_experiment e;
    swa_reversal_result r;
    size_t i;
    fill_clean_experiment(plus, minus, zero, N, 10.0, 2.0, 0.6);
    for (i = 0U; i < N; ++i) zero[i] -= 1.0;
    e = config(plus, minus, zero, N);
    check_true(swa_audit_reversal_experiment(&e, &r), "biased reversal audit executes");
    check_true(!r.common_mode_controlled, "common-mode force bias detected");
    check_true(!r.passes, "common-mode biased experiment rejected");
}

static void test_asymmetry_rejected(void) {
    enum { N = 12 };
    double plus[N], minus[N], zero[N];
    swa_reversal_experiment e;
    swa_reversal_result r;
    size_t i;
    fill_clean_experiment(plus, minus, zero, N, 10.0, 2.0, 0.9);
    for (i = 0U; i < N; ++i) minus[i] += 5.0;
    e = config(plus, minus, zero, N);
    check_true(swa_audit_reversal_experiment(&e, &r), "asymmetric reversal audit executes");
    check_true(!r.reversal_symmetric, "command reversal asymmetry detected");
    check_true(!r.passes, "asymmetric experiment rejected");
}

static void test_drift_rejected(void) {
    enum { N = 12 };
    double plus[N], minus[N], zero[N];
    swa_reversal_experiment e;
    swa_reversal_result r;
    size_t i;
    fill_clean_experiment(plus, minus, zero, N, 10.0, 2.0, 1.2);
    for (i = 0U; i < N; ++i) zero[i] += 0.2 * (double)i;
    e = config(plus, minus, zero, N);
    e.maximum_common_mode_bias_N = 10.0;
    check_true(swa_audit_reversal_experiment(&e, &r), "drifting reversal audit executes");
    check_true(!r.drift_controlled, "zero-control drift detected");
    check_true(!r.passes, "drifting experiment rejected");
}

static void test_replication(void) {
    enum { N = 12 };
    double p1[N], m1[N], z1[N], p2[N], m2[N], z2[N], p3[N], m3[N], z3[N];
    swa_reversal_experiment e1, e2, e3;
    swa_reversal_result results[3];
    swa_replication_audit audit;
    swa_replication_result replication;
    fill_clean_experiment(p1, m1, z1, N, 10.0, 2.0, 0.0);
    fill_clean_experiment(p2, m2, z2, N, 10.15, -1.0, 0.4);
    fill_clean_experiment(p3, m3, z3, N, -10.0, 0.5, 0.8);
    e1 = config(p1, m1, z1, N);
    e2 = config(p2, m2, z2, N);
    e3 = config(p3, m3, z3, N);
    check_true(swa_audit_reversal_experiment(&e1, &results[0]) && results[0].passes,
               "first independent reversal passes");
    check_true(swa_audit_reversal_experiment(&e2, &results[1]) && results[1].passes,
               "second independent reversal passes");
    audit.replicates = results;
    audit.replicate_count = 2U;
    audit.maximum_relative_heterogeneity = 0.10;
    audit.minimum_passing_replicates = 2U;
    check_true(swa_audit_independent_replication(&audit, &replication),
               "independent replication audit executes");
    check_true(replication.passes, "two consistent independent reversals replicate");
    check_true(swa_audit_reversal_experiment(&e3, &results[2]),
               "opposite-sign reversal audit executes");
    audit.replicate_count = 3U;
    check_true(swa_audit_independent_replication(&audit, &replication),
               "inconsistent replication audit executes");
    check_true(!replication.same_sign && !replication.passes,
               "opposite-sign replication is rejected");
}

static void test_large_deterministic_sign_flip(void) {
    enum { N = 32 };
    double plus[N], minus[N], zero[N];
    swa_reversal_experiment e;
    swa_reversal_result a, b;
    fill_clean_experiment(plus, minus, zero, N, 5.0, 0.2, 0.7);
    e = config(plus, minus, zero, N);
    e.maximum_lag1_autocorrelation = 0.95;
    check_true(swa_audit_reversal_experiment(&e, &a),
               "large deterministic sign-flip audit executes");
    check_true(swa_audit_reversal_experiment(&e, &b),
               "large deterministic sign-flip replay executes");
    check_true(a.exact_or_deterministic_sign_flip_p ==
               b.exact_or_deterministic_sign_flip_p,
               "deterministic sign-flip p-value replays exactly");
}

int main(void) {
    FILE *fp;
    test_clean_reversal();
    test_null_effect_rejected();
    test_common_mode_rejected();
    test_asymmetry_rejected();
    test_drift_rejected();
    test_replication();
    test_large_deterministic_sign_flip();
    fp = fopen("output/experiment_assurance_receipt.json", "w");
    if (fp != NULL) {
        (void)swa_write_assurance_receipt(fp, checks, failures, receipt_hash);
        (void)fclose(fp);
    } else {
        ++failures;
    }
    printf("spacewind experiment assurance: %zu checks, %zu failures, hash=%016llx\n",
           checks, failures, (unsigned long long)receipt_hash);
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
