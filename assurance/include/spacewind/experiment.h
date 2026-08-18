#ifndef SPACEWIND_EXPERIMENT_H
#define SPACEWIND_EXPERIMENT_H

#include "spacewind/assurance.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SWA_EXACT_SIGN_FLIP_MAX 20U
#define SWA_DETERMINISTIC_SIGN_FLIP_SAMPLES 65536U

typedef struct {
    const double *positive_command_force_N;
    const double *negative_command_force_N;
    const double *zero_command_force_N;
    size_t count;
    double alpha;
    double minimum_effect_N;
    double maximum_common_mode_bias_N;
    double maximum_reversal_asymmetry;
    double maximum_lag1_autocorrelation;
    double maximum_zero_drift_N_per_sample;
} swa_reversal_experiment;

typedef struct {
    size_t count;
    double positive_mean_N;
    double negative_mean_N;
    double zero_mean_N;
    double reversed_effect_mean_N;
    double common_mode_bias_N;
    double reversal_asymmetry;
    double paired_standard_deviation_N;
    double standard_error_N;
    double empirical_bernstein_radius_N;
    swa_interval effect_interval_N;
    double exact_or_deterministic_sign_flip_p;
    double lag1_autocorrelation;
    double effective_sample_size;
    double zero_drift_N_per_sample;
    int finite;
    int effect_above_minimum;
    int common_mode_controlled;
    int reversal_symmetric;
    int autocorrelation_controlled;
    int drift_controlled;
    int randomization_significant;
    int passes;
} swa_reversal_result;

typedef struct {
    const swa_reversal_result *replicates;
    size_t replicate_count;
    double maximum_relative_heterogeneity;
    size_t minimum_passing_replicates;
} swa_replication_audit;

typedef struct {
    size_t passing_replicates;
    double pooled_effect_N;
    double minimum_lower_bound_N;
    double maximum_upper_bound_N;
    double relative_heterogeneity;
    int same_sign;
    int interval_overlap;
    int passes;
} swa_replication_result;

int swa_audit_reversal_experiment(const swa_reversal_experiment *experiment,
                                  swa_reversal_result *out);
int swa_audit_independent_replication(const swa_replication_audit *audit,
                                      swa_replication_result *out);
int swa_write_reversal_receipt(FILE *fp,
                               const swa_reversal_experiment *experiment,
                               const swa_reversal_result *result);

#ifdef __cplusplus
}
#endif

#endif
