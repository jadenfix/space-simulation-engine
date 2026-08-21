#include "spacewind/em_pic3d.h"

#include <float.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MODE_NX 8U
#define MODE_NY 1U
#define MODE_NZ 1U
#define MODE_PARTICLES (2U * MODE_NX * MODE_NY * MODE_NZ)
#define MODE_MAX_CROSSINGS 16U
#define MODE_FIELD_AMPLITUDE_V_M 1000.0
#define MODE_PERIODS 4.5
#define MODE_MAX_NONLINEAR_ITERATIONS 32U
/* Instantaneous current normalization is weakest at mode zero crossings. */
#define MODE_SPECTRAL_ORACLE_ENVELOPE 1.0e-7

typedef struct {
    double density_m3;
    double negative_mass_kg;
    double positive_mass_kg;
    double dt_s;
    double continuum_omega_rad_s;
    double midpoint_omega_rad_s;
    double measured_omega_rad_s;
    double initial_mean_e_V_m;
    double first_positive_peak_V_m;
    double last_positive_peak_V_m;
    double peak_relative_change;
    double maximum_total_energy_relative_drift;
    double maximum_step_energy_relative;
    double maximum_particle_work_relative;
    double maximum_field_work_relative;
    double maximum_continuity_relative;
    double maximum_gauss_relative;
    double maximum_spectral_oracle_relative;
    double maximum_momentum_relative;
    size_t steps;
    size_t crossings;
    size_t positive_peaks;
    int finite;
    int passes;
} mode_result;

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);

static double max2(double a, double b) {
    return a > b ? a : b;
}

static double max3(double a, double b, double c) {
    return max2(max2(a, b), c);
}

static double relative_scalar(double a, double b) {
    return fabs(a - b) / max3(fabs(a), fabs(b), DBL_MIN);
}

static void check_true(int condition, const char *name) {
    ++checks;
    transcript_hash ^= swa_fnv1a64(name, strlen(name));
    transcript_hash *= UINT64_C(1099511628211);
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}

static double plasma_omega_rad_s(
    double density_m3,
    double negative_mass_kg,
    double positive_mass_kg
) {
    const double inverse_mass =
        1.0 / negative_mass_kg + 1.0 / positive_mass_kg;
    return sqrt(
        density_m3 * SWA_QE * SWA_QE * inverse_mass / SWA_EPS0
    );
}

static double midpoint_omega_rad_s(
    double continuum_omega_rad_s,
    double dt_s
) {
    return 2.0 * atan(0.5 * continuum_omega_rad_s * dt_s) / dt_s;
}

static swa_em_pic3d_limits mode_limits(void) {
    swa_em_pic3d_limits value;
    memset(&value, 0, sizeof(value));
    value.continuity_relative_tolerance = 2.0e-10;
    value.gauss_relative_tolerance = 2.0e-10;
    value.magnetic_divergence_relative_tolerance = 2.0e-11;
    value.local_current_difference_relative_tolerance = 0.90;
    value.spectral_oracle_relative_tolerance =
        MODE_SPECTRAL_ORACLE_ENVELOPE;
    value.spectral_oracle_curl_relative_tolerance = 5.0e-10;
    value.energy_relative_tolerance = 5.0e-8;
    value.particle_work_relative_tolerance = 2.0e-4;
    value.field_work_relative_tolerance = 5.0e-8;
    value.momentum_relative_tolerance = 5.0e-6;
    value.nonlinear_relative_tolerance = 2.0e-10;
    value.maximum_courant = 0.25;
    value.maximum_particle_cells_per_step = 0.20;
    value.maximum_nonlinear_iterations = MODE_MAX_NONLINEAR_ITERATIONS;
    return value;
}

static double stable_gamma_minus_one(double speed2_m2_s2) {
