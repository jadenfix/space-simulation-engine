#include "spacewind/experiment.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static double max2(double a, double b) { return a > b ? a : b; }
static double min2(double a, double b) { return a < b ? a : b; }
static double clamp(double x, double lo, double hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}
static uint64_t splitmix64(uint64_t x) {
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30U)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27U)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31U);
}

static int finite_series(const double *x, size_t n) {
    size_t i;
    if (x == NULL || n == 0U) {
        return 0;
    }
    for (i = 0U; i < n; ++i) {
        if (!isfinite(x[i])) {
            return 0;
        }
    }
    return 1;
}

static double mean_series(const double *x, size_t n) {
    swa_kahan sum = {0.0, 0.0};
    size_t i;
    for (i = 0U; i < n; ++i) {
        swa_kahan_add(&sum, x[i]);
    }
    return sum.sum / (double)n;
}

static double sample_variance(const double *x, size_t n, double mean) {
    swa_kahan sum = {0.0, 0.0};
    size_t i;
    if (n < 2U) {
        return 0.0;
    }
    for (i = 0U; i < n; ++i) {
        const double d = x[i] - mean;
        swa_kahan_add(&sum, d * d);
    }
    return sum.sum / (double)(n - 1U);
}

static double lag1_autocorrelation(const double *x, size_t n, double mean) {
    swa_kahan numerator = {0.0, 0.0};
    swa_kahan denominator = {0.0, 0.0};
    size_t i;
    if (n < 3U) {
        return 0.0;
    }
    for (i = 0U; i < n; ++i) {
        const double d = x[i] - mean;
        swa_kahan_add(&denominator, d * d);
        if (i + 1U < n) {
            swa_kahan_add(&numerator, d * (x[i + 1U] - mean));
        }
    }
    if (!(denominator.sum > 0.0)) {
        return 0.0;
    }
    return clamp(numerator.sum / denominator.sum, -0.999, 0.999);
}

static double linear_slope(const double *x, size_t n) {
    const double center = 0.5 * (double)(n - 1U);
    swa_kahan numerator = {0.0, 0.0};
    swa_kahan denominator = {0.0, 0.0};
    size_t i;
    if (n < 2U) {
        return 0.0;
    }
    for (i = 0U; i < n; ++i) {
        const double t = (double)i - center;
        swa_kahan_add(&numerator, t * x[i]);
        swa_kahan_add(&denominator, t * t);
    }
    return denominator.sum > 0.0 ? numerator.sum / denominator.sum : 0.0;
}

static double sign_flip_p_value(const double *paired, size_t n, double observed_abs_mean) {
    uint64_t extreme = 0U;
    uint64_t total = 0U;
    size_t permutation;
    const double comparison_slack = 16.0 * DBL_EPSILON * max2(observed_abs_mean, 1.0);
    if (n <= SWA_EXACT_SIGN_FLIP_MAX) {
        const uint64_t combinations = UINT64_C(1) << n;
        uint64_t mask;
        for (mask = 0U; mask < combinations; ++mask) {
            swa_kahan sum = {0.0, 0.0};
            size_t j;
            for (j = 0U; j < n; ++j) {
                const double sign = ((mask >> j) & UINT64_C(1)) != 0U ? 1.0 : -1.0;
                swa_kahan_add(&sum, sign * paired[j]);
            }
            if (fabs(sum.sum / (double)n) + comparison_slack >= observed_abs_mean) {
                ++extreme;
            }
        }
        total = combinations;
    } else {
        for (permutation = 0U; permutation < SWA_DETERMINISTIC_SIGN_FLIP_SAMPLES; ++permutation) {
            swa_kahan sum = {0.0, 0.0};
            size_t j;
            for (j = 0U; j < n; ++j) {
                const uint64_t bits = splitmix64((uint64_t)permutation * UINT64_C(0x9e3779b97f4a7c15) +
                                                 (uint64_t)j * UINT64_C(0xd1b54a32d192ed03));
                const double sign = (bits & UINT64_C(1)) != 0U ? 1.0 : -1.0;
                swa_kahan_add(&sum, sign * paired[j]);
            }
            if (fabs(sum.sum / (double)n) + comparison_slack >= observed_abs_mean) {
                ++extreme;
            }
        }
        total = SWA_DETERMINISTIC_SIGN_FLIP_SAMPLES;
    }
    return ((double)extreme + 1.0) / ((double)total + 1.0);
}

int swa_audit_reversal_experiment(const swa_reversal_experiment *e,
                                  swa_reversal_result *o) {
    double *paired;
    double paired_min = INFINITY;
    double paired_max = -INFINITY;
    double variance;
    double log_term;
    double positive_effect;
    double negative_effect;
    size_t i;
    if (e == NULL || o == NULL || e->count < 3U ||
        !finite_series(e->positive_command_force_N, e->count) ||
        !finite_series(e->negative_command_force_N, e->count) ||
        !finite_series(e->zero_command_force_N, e->count) ||
        !(e->alpha > 0.0 && e->alpha < 1.0) ||
        !(e->minimum_effect_N >= 0.0) ||
        !(e->maximum_common_mode_bias_N >= 0.0) ||
        !(e->maximum_reversal_asymmetry >= 0.0) ||
        !(e->maximum_lag1_autocorrelation >= 0.0 &&
          e->maximum_lag1_autocorrelation < 1.0) ||
        !(e->maximum_zero_drift_N_per_sample >= 0.0)) {
        return 0;
    }
    memset(o, 0, sizeof(*o));
    paired = (double *)malloc(e->count * sizeof(*paired));
    if (paired == NULL) {
        return 0;
    }
    o->count = e->count;
    o->positive_mean_N = mean_series(e->positive_command_force_N, e->count);
    o->negative_mean_N = mean_series(e->negative_command_force_N, e->count);
    o->zero_mean_N = mean_series(e->zero_command_force_N, e->count);
    for (i = 0U; i < e->count; ++i) {
        paired[i] = 0.5 * (e->positive_command_force_N[i] -
                           e->negative_command_force_N[i]);
        if (paired[i] < paired_min) {
            paired_min = paired[i];
        }
        if (paired[i] > paired_max) {
            paired_max = paired[i];
        }
    }
    o->reversed_effect_mean_N = mean_series(paired, e->count);
    o->common_mode_bias_N = 0.5 * (o->positive_mean_N + o->negative_mean_N) -
                            o->zero_mean_N;
    positive_effect = o->positive_mean_N - o->zero_mean_N;
    negative_effect = o->zero_mean_N - o->negative_mean_N;
    o->reversal_asymmetry = fabs(positive_effect - negative_effect) /
                            max2(0.5 * (fabs(positive_effect) + fabs(negative_effect)), DBL_MIN);
    variance = sample_variance(paired, e->count, o->reversed_effect_mean_N);
    o->paired_standard_deviation_N = sqrt(max2(variance, 0.0));
    o->lag1_autocorrelation = lag1_autocorrelation(paired, e->count,
                                                   o->reversed_effect_mean_N);
    o->effective_sample_size = clamp((double)e->count *
        (1.0 - o->lag1_autocorrelation) / (1.0 + o->lag1_autocorrelation),
        1.0, (double)e->count);
    o->standard_error_N = o->paired_standard_deviation_N /
                          sqrt(o->effective_sample_size);
    log_term = log(3.0 / e->alpha);
    o->empirical_bernstein_radius_N =
        sqrt(2.0 * variance * log_term / o->effective_sample_size) +
        3.0 * (paired_max - paired_min) * log_term / o->effective_sample_size;
    o->effect_interval_N = swa_interval_make(
        o->reversed_effect_mean_N - o->empirical_bernstein_radius_N,
        o->reversed_effect_mean_N + o->empirical_bernstein_radius_N);
    o->exact_or_deterministic_sign_flip_p = sign_flip_p_value(
        paired, e->count, fabs(o->reversed_effect_mean_N));
    o->zero_drift_N_per_sample = linear_slope(e->zero_command_force_N, e->count);
    o->finite = isfinite(o->positive_mean_N) && isfinite(o->negative_mean_N) &&
                isfinite(o->zero_mean_N) && isfinite(o->reversed_effect_mean_N) &&
                isfinite(o->common_mode_bias_N) && isfinite(o->reversal_asymmetry) &&
                isfinite(o->paired_standard_deviation_N) &&
                isfinite(o->empirical_bernstein_radius_N) &&
                isfinite(o->exact_or_deterministic_sign_flip_p) &&
                isfinite(o->lag1_autocorrelation) &&
                isfinite(o->effective_sample_size) &&
                isfinite(o->zero_drift_N_per_sample);
    o->effect_above_minimum = o->effect_interval_N.lo > e->minimum_effect_N;
    o->common_mode_controlled = fabs(o->common_mode_bias_N) <=
                                e->maximum_common_mode_bias_N;
    o->reversal_symmetric = o->reversal_asymmetry <=
                            e->maximum_reversal_asymmetry;
    o->autocorrelation_controlled = fabs(o->lag1_autocorrelation) <=
                                    e->maximum_lag1_autocorrelation;
    o->drift_controlled = fabs(o->zero_drift_N_per_sample) <=
                          e->maximum_zero_drift_N_per_sample;
    o->randomization_significant = o->exact_or_deterministic_sign_flip_p <=
                                   e->alpha;
    o->passes = o->finite && o->effect_above_minimum &&
                o->common_mode_controlled && o->reversal_symmetric &&
                o->autocorrelation_controlled && o->drift_controlled &&
                o->randomization_significant;
    free(paired);
    return 1;
}

int swa_audit_independent_replication(const swa_replication_audit *a,
                                      swa_replication_result *o) {
    double minimum_effect = INFINITY;
    double maximum_effect = -INFINITY;
    double maximum_lower = -INFINITY;
    double minimum_upper = INFINITY;
    swa_kahan pooled = {0.0, 0.0};
    int sign = 0;
    size_t i;
    if (a == NULL || o == NULL || a->replicates == NULL ||
        a->replicate_count == 0U ||
        !(a->maximum_relative_heterogeneity >= 0.0) ||
        a->minimum_passing_replicates == 0U) {
        return 0;
    }
    memset(o, 0, sizeof(*o));
    o->same_sign = 1;
    for (i = 0U; i < a->replicate_count; ++i) {
        const swa_reversal_result *r = &a->replicates[i];
        const int current_sign = r->reversed_effect_mean_N > 0.0 ? 1 :
                                 (r->reversed_effect_mean_N < 0.0 ? -1 : 0);
        if (!r->finite) {
            o->same_sign = 0;
            continue;
        }
        if (r->passes) {
            ++o->passing_replicates;
        }
        if (sign == 0) {
            sign = current_sign;
        } else if (current_sign == 0 || current_sign != sign) {
            o->same_sign = 0;
        }
        if (r->reversed_effect_mean_N < minimum_effect) {
            minimum_effect = r->reversed_effect_mean_N;
        }
        if (r->reversed_effect_mean_N > maximum_effect) {
            maximum_effect = r->reversed_effect_mean_N;
        }
        if (r->effect_interval_N.lo > maximum_lower) {
            maximum_lower = r->effect_interval_N.lo;
        }
        if (r->effect_interval_N.hi < minimum_upper) {
            minimum_upper = r->effect_interval_N.hi;
        }
        swa_kahan_add(&pooled, r->reversed_effect_mean_N);
    }
    o->pooled_effect_N = pooled.sum / (double)a->replicate_count;
    o->minimum_lower_bound_N = maximum_lower;
    o->maximum_upper_bound_N = minimum_upper;
    o->relative_heterogeneity = (maximum_effect - minimum_effect) /
                                max2(fabs(o->pooled_effect_N), DBL_MIN);
    o->interval_overlap = maximum_lower <= minimum_upper;
    o->passes = o->passing_replicates >= a->minimum_passing_replicates &&
                o->same_sign && o->interval_overlap &&
                o->relative_heterogeneity <= a->maximum_relative_heterogeneity;
    return 1;
}

int swa_write_reversal_receipt(FILE *fp,
                               const swa_reversal_experiment *e,
                               const swa_reversal_result *r) {
    if (fp == NULL || e == NULL || r == NULL) {
        return 0;
    }
    return fprintf(fp,
        "{\n  \"schema\": \"spacewind.force-reversal/v1\",\n"
        "  \"sample_count\": %zu,\n  \"positive_mean_N\": %.17g,\n"
        "  \"negative_mean_N\": %.17g,\n  \"zero_mean_N\": %.17g,\n"
        "  \"reversed_effect_mean_N\": %.17g,\n"
        "  \"effect_interval_N\": [%.17g, %.17g],\n"
        "  \"common_mode_bias_N\": %.17g,\n"
        "  \"reversal_asymmetry\": %.17g,\n"
        "  \"sign_flip_p\": %.17g,\n"
        "  \"lag1_autocorrelation\": %.17g,\n"
        "  \"effective_sample_size\": %.17g,\n"
        "  \"zero_drift_N_per_sample\": %.17g,\n"
        "  \"passes\": %s,\n"
        "  \"assumptions\": [\"bounded paired observations\", \"exchangeable signs under the null\", \"declared drift and autocorrelation limits\"],\n"
        "  \"nonclaim\": \"a passing synthetic or chamber reversal audit is not a flight-qualified propulsion result\"\n}\n",
        r->count, r->positive_mean_N, r->negative_mean_N, r->zero_mean_N,
        r->reversed_effect_mean_N, r->effect_interval_N.lo,
        r->effect_interval_N.hi, r->common_mode_bias_N,
        r->reversal_asymmetry, r->exact_or_deterministic_sign_flip_p,
        r->lag1_autocorrelation, r->effective_sample_size,
        r->zero_drift_N_per_sample, r->passes ? "true" : "false") > 0;
}
