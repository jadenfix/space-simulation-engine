#ifndef SPACEWIND_ASSURANCE_H
#define SPACEWIND_ASSURANCE_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SWA_PI 3.141592653589793238462643383279502884
#define SWA_C 299792458.0
#define SWA_EPS0 8.8541878128e-12
#define SWA_MU0 1.25663706212e-6
#define SWA_KB 1.380649e-23
#define SWA_QE 1.602176634e-19
#define SWA_ME 9.1093837139e-31
#define SWA_MP 1.67262192369e-27
#define SWA_G 6.67430e-11
#define SWA_AU 149597870700.0
#define SWA_SOLAR_IRRADIANCE_1AU 1361.0
#define SWA_DIM_COUNT 7
#define SWA_SIMILARITY_COUNT 12

typedef struct { double x, y, z; } swa_vec3;

typedef struct {
    int8_t exponent[SWA_DIM_COUNT];
} swa_dimension;

typedef struct {
    double value;
    swa_dimension dimension;
} swa_quantity;

typedef struct {
    double lo;
    double hi;
} swa_interval;

typedef struct {
    double sum;
    double correction;
} swa_kahan;

typedef struct {
    size_t count;
    double mean;
    double m2;
    double minimum;
    double maximum;
} swa_online_stats;

typedef struct {
    double mass_initial_kg;
    double mass_final_kg;
    double mass_in_kg;
    double mass_out_kg;
    double charge_initial_C;
    double charge_final_C;
    double charge_in_C;
    double charge_out_C;
    swa_vec3 momentum_initial_Ns;
    swa_vec3 momentum_final_Ns;
    swa_vec3 impulse_external_Ns;
    swa_vec3 momentum_in_Ns;
    swa_vec3 momentum_out_Ns;
    double energy_initial_J;
    double energy_final_J;
    double work_source_J;
    double heat_in_J;
    double energy_in_J;
    double energy_out_J;
} swa_ledger;

typedef struct {
    double mass_residual_kg;
    double charge_residual_C;
    swa_vec3 momentum_residual_Ns;
    double momentum_residual_norm_Ns;
    double energy_residual_J;
    double mass_relative;
    double charge_relative;
    double momentum_relative;
    double energy_relative;
    int finite;
    int passes;
} swa_ledger_result;

typedef struct {
    double ion_density_m3;
    double electron_density_m3;
    double ion_temperature_K;
    double electron_temperature_K;
    double ion_mass_kg;
    double ion_charge_C;
    double neutral_density_m3;
    double collision_cross_section_m2;
    double electrode_radius_m;
    double interaction_length_m;
    double electrode_voltage_V;
    double circuit_capacitance_F;
    double available_power_W;
    double collection_current_A;
    swa_vec3 bulk_velocity_mps;
    swa_vec3 magnetic_field_T;
} swa_plasma_state;

typedef struct {
    double debye_length_m;
    double ion_plasma_frequency_rad_s;
    double electron_plasma_frequency_rad_s;
    double ion_inertial_length_m;
    double ion_gyrofrequency_rad_s;
    double ion_thermal_speed_mps;
    double ion_gyroradius_m;
    double sound_speed_mps;
    double alfven_speed_mps;
    double fast_magnetosonic_speed_mps;
    double plasma_beta;
    double sonic_mach;
    double alfven_mach;
    double fast_mach;
    double knudsen;
    double normalized_voltage;
    double radius_to_debye;
    double gyroradius_to_length;
    double inertial_to_length;
    double circuit_to_flow_time;
    double power_margin;
    int finite;
} swa_plasma_regime;

typedef struct {
    double value[SWA_SIMILARITY_COUNT];
    double tolerance_log[SWA_SIMILARITY_COUNT];
    uint32_t evidence_level;
} swa_similarity_signature;

typedef struct {
    double weighted_log_distance;
    double maximum_log_mismatch;
    size_t worst_index;
    uint32_t minimum_evidence_level;
    int inside;
} swa_similarity_result;

typedef struct {
    double observed_order;
    double extrapolated_value;
    double fine_error_estimate;
    double gci_fine;
    double asymptotic_ratio;
    int monotone;
    int finite;
    int passes;
} swa_convergence_result;

typedef struct {
    swa_vec3 position_m;
    swa_vec3 velocity_mps;
    double stored_energy_J;
    double controller_phase_rad;
} swa_cycle_state;

typedef struct {
    swa_cycle_state start;
    swa_cycle_state end;
    double craft_work_J;
    double source_flow_work_J;
    double relative_dissipation_J;
    double actuator_energy_J;
    double thermal_loss_J;
    double numerical_residual_J;
    swa_interval uncertain_net_gain_J;
    double position_scale_m;
    double velocity_scale_mps;
    double energy_scale_J;
    double phase_scale_rad;
    uint32_t evidence_level;
    uint32_t required_evidence_level;
} swa_cycle_audit;

typedef struct {
    double normalized_state_closure;
    double moving_medium_identity_relative;
    double energy_identity_relative;
    double net_gain_J;
    swa_interval net_gain_interval_J;
    int state_closed;
    int identity_closed;
    int evidence_sufficient;
    int robust_positive_gain;
    int numerical_cycle_passes;
    int physical_promotion_passes;
} swa_cycle_result;

typedef struct {
    double sensor_bias_fraction;
    double dropout_fraction;
    double tether_cut_fraction;
    double arc_energy_J;
    double temperature_K;
    double temperature_limit_K;
    double bus_voltage_V;
    double minimum_bus_voltage_V;
    double angular_rate_rad_s;
    double maximum_angular_rate_rad_s;
} swa_fault_state;

typedef struct {
    int sensor_degraded;
    int communications_degraded;
    int tether_degraded;
    int arc_detected;
    int thermal_trip;
    int undervoltage;
    int attitude_trip;
    int safe_mode_required;
    double severity;
} swa_fault_result;

swa_vec3 swa_v3(double x, double y, double z);
swa_vec3 swa_vadd(swa_vec3 a, swa_vec3 b);
swa_vec3 swa_vsub(swa_vec3 a, swa_vec3 b);
swa_vec3 swa_vscale(swa_vec3 a, double s);
double swa_vdot(swa_vec3 a, swa_vec3 b);
swa_vec3 swa_vcross(swa_vec3 a, swa_vec3 b);
double swa_vnorm(swa_vec3 a);
swa_vec3 swa_rotate_z(swa_vec3 a, double angle_rad);

swa_dimension swa_dim(int8_t mass, int8_t length, int8_t time, int8_t current,
                      int8_t temperature, int8_t amount, int8_t luminous);
int swa_dim_equal(swa_dimension a, swa_dimension b);
swa_dimension swa_dim_mul(swa_dimension a, swa_dimension b);
swa_dimension swa_dim_div(swa_dimension a, swa_dimension b);
swa_quantity swa_quantity_make(double value, swa_dimension dimension);
int swa_quantity_add(swa_quantity a, swa_quantity b, swa_quantity *out);
swa_quantity swa_quantity_mul(swa_quantity a, swa_quantity b);
swa_quantity swa_quantity_div(swa_quantity a, swa_quantity b);

swa_interval swa_interval_make(double a, double b);
swa_interval swa_interval_add(swa_interval a, swa_interval b);
swa_interval swa_interval_sub(swa_interval a, swa_interval b);
swa_interval swa_interval_mul(swa_interval a, swa_interval b);
int swa_interval_div(swa_interval a, swa_interval b, swa_interval *out);
int swa_interval_sqrt(swa_interval a, swa_interval *out);
int swa_interval_contains(swa_interval a, double x);
double swa_interval_width(swa_interval a);

void swa_kahan_add(swa_kahan *accumulator, double value);
void swa_stats_reset(swa_online_stats *stats);
void swa_stats_push(swa_online_stats *stats, double value);
double swa_stats_variance(const swa_online_stats *stats);
double swa_halton(size_t index, unsigned base);
uint64_t swa_fnv1a64(const void *data, size_t size);

int swa_audit_ledger(const swa_ledger *ledger, double relative_tolerance,
                     double absolute_tolerance, swa_ledger_result *out);
int swa_compute_plasma_regime(const swa_plasma_state *state, swa_plasma_regime *out);
void swa_make_similarity_signature(const swa_plasma_regime *regime,
                                   uint32_t evidence_level,
                                   swa_similarity_signature *out);
int swa_compare_similarity(const swa_similarity_signature *target,
                           const swa_similarity_signature *candidate,
                           uint32_t required_evidence_level,
                           swa_similarity_result *out);
int swa_richardson_three_level(double coarse, double medium, double fine,
                               double refinement_ratio, double expected_order,
                               double order_tolerance, swa_convergence_result *out);
int swa_audit_cycle(const swa_cycle_audit *audit, double state_tolerance,
                    double identity_tolerance, swa_cycle_result *out);
int swa_assess_faults(const swa_fault_state *state, swa_fault_result *out);
int swa_write_cycle_receipt(FILE *fp, const swa_cycle_audit *audit,
                            const swa_cycle_result *result);
int swa_write_assurance_receipt(FILE *fp, size_t checks, size_t failures,
                                uint64_t deterministic_hash);

#ifdef __cplusplus
}
#endif

#endif
