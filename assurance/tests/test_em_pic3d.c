#include "spacewind/em_pic3d.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);
static double worst_pusher_reversal = 0.0;
static double worst_pusher_rotation = 0.0;
static double worst_projection_continuity = 0.0;
static double worst_projection_curl = 0.0;
static double worst_plane_wave_error = 0.0;
static double worst_gauss_relative = 0.0;
static double worst_div_b_relative = 0.0;
static double worst_coupled_energy = 0.0;
static double worst_coupled_momentum = 0.0;
static double worst_coupled_projection = 0.0;
static double worst_spectral_oracle_difference = 0.0;
static swa_em_pic3d_result receipt_result;

static double max2(double a, double b) { return a > b ? a : b; }
static double max3(double a, double b, double c) {
    return max2(max2(a, b), c);
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

static size_t index3(
    const swa_em_pic3d *state,
    size_t i,
    size_t j,
    size_t k
) {
    return (k * state->ny + j) * state->nx + i;
}

static double relative_scalar(double a, double b) {
    return fabs(a - b) / max3(fabs(a), fabs(b), DBL_MIN);
}

static double relative_vec(swa_vec3 a, swa_vec3 b) {
    return swa_vnorm(swa_vsub(a, b)) /
        max3(swa_vnorm(a), swa_vnorm(b), DBL_MIN);
}

static swa_vec3 cycle_xyz(swa_vec3 value) {
    return swa_v3(value.y, value.z, value.x);
}

static swa_vec3 negate_vec(swa_vec3 value) {
    return swa_vscale(value, -1.0);
}

static swa_vec3 deterministic_direction(
    size_t i,
    unsigned a,
    unsigned b,
    unsigned c
) {
    swa_vec3 direction = swa_v3(
        2.0 * swa_halton(i, a) - 1.0,
        2.0 * swa_halton(i, b) - 1.0,
        2.0 * swa_halton(i, c) - 1.0
    );
    const double norm = swa_vnorm(direction);
    if (!(norm > 1.0e-15)) {
        return swa_v3(1.0, 0.0, 0.0);
    }
    return swa_vscale(direction, 1.0 / norm);
}

static swa_em_pic3d_limits permissive_limits(void) {
    swa_em_pic3d_limits limits;
    memset(&limits, 0, sizeof(limits));
    limits.continuity_relative_tolerance = 1.0e-8;
    limits.gauss_relative_tolerance = 1.0e-8;
    limits.magnetic_divergence_relative_tolerance = 1.0e-8;
    limits.local_current_difference_relative_tolerance = 1.0;
    limits.spectral_oracle_relative_tolerance = 1.0e-8;
    limits.spectral_oracle_curl_relative_tolerance = 1.0e-8;
    limits.energy_relative_tolerance = 1.0;
    limits.particle_work_relative_tolerance = 1.0;
    limits.field_work_relative_tolerance = 1.0;
    limits.momentum_relative_tolerance = 1.0;
    limits.nonlinear_relative_tolerance = 1.0e-9;
    limits.maximum_courant = 0.95;
    limits.maximum_particle_cells_per_step = 0.95;
    limits.maximum_nonlinear_iterations = 12U;
    return limits;
}

static swa_em_pic3d_limits strict_coupled_limits(void) {
    swa_em_pic3d_limits limits = permissive_limits();
    limits.continuity_relative_tolerance = 2.0e-10;
    limits.gauss_relative_tolerance = 2.0e-11;
    limits.magnetic_divergence_relative_tolerance = 2.0e-11;
    limits.local_current_difference_relative_tolerance = 0.85;
    limits.spectral_oracle_relative_tolerance = 2.0e-10;
    limits.spectral_oracle_curl_relative_tolerance = 2.0e-10;
    limits.energy_relative_tolerance = 2.0e-10;
    limits.particle_work_relative_tolerance = 5.0e-6;
    limits.field_work_relative_tolerance = 2.0e-10;
    limits.momentum_relative_tolerance = 2.0e-4;
    limits.nonlinear_relative_tolerance = 2.0e-10;
    limits.maximum_courant = 0.8;
    limits.maximum_particle_cells_per_step = 0.8;
    limits.maximum_nonlinear_iterations = 16U;
    return limits;
}

static void test_invalid_contracts(void) {
    swa_em_pic3d state;
    swa_em_pic3d_result result;
    swa_em_pic3d_limits limits = permissive_limits();
    swa_vec3 output;
    memset(&state, 0, sizeof(state));
    check_true(!swa_em_pic3d_init(
                   NULL, 4U, 4U, 4U, 0U,
                   1.0, 1.0, 1.0, 1.0e-10),
               "null 3D3V state rejected");
    check_true(!swa_em_pic3d_init(
                   &state, 0U, 4U, 4U, 0U,
                   1.0, 1.0, 1.0, 1.0e-10),
               "zero 3D3V dimension rejected");
    check_true(!swa_em_pic3d_init(
                   &state, 4U, 4U, 4U, 0U,
                   -1.0, 1.0, 1.0, 1.0e-10),
               "negative 3D3V length rejected");
    check_true(!swa_em_pic3d_init(
                   &state, 4U, 4U, 4U, 0U,
                   1.0, 1.0, 1.0, 0.0),
               "zero 3D3V timestep rejected");
    check_true(!swa_em_pic3d_higuera_cary_push(
                   swa_v3(SWA_C, 0.0, 0.0),
                   1.0, 1.0,
                   swa_v3(0.0, 0.0, 0.0),
                   swa_v3(0.0, 0.0, 0.0),
                   1.0, &output),
               "luminal Higuera-Cary input rejected");
    check_true(!swa_em_pic3d_higuera_cary_push(
                   swa_v3(0.0, 0.0, 0.0),
                   1.0, 0.0,
                   swa_v3(0.0, 0.0, 0.0),
                   swa_v3(0.0, 0.0, 0.0),
                   1.0, &output),
               "zero pusher mass rejected");
    check_true(swa_em_pic3d_init(
                   &state, 3U, 3U, 3U, 1U,
                   1.0, 1.0, 1.0, 1.0e-10),
               "valid invalid-input test state initializes");
    check_true(!swa_em_pic3d_set_particle(
                   &state, 0U,
                   swa_v3(0.1, 0.2, 0.3),
                   swa_v3(SWA_C, 0.0, 0.0),
                   1.0, 1.0, 1.0),
               "luminal 3D3V particle rejected");
    check_true(!swa_em_pic3d_set_particle(
                   &state, 0U,
                   swa_v3(0.1, 0.2, 0.3),
                   swa_v3(0.0, 0.0, 0.0),
                   1.0, -1.0, 1.0),
               "negative 3D3V particle mass rejected");
    check_true(!swa_em_pic3d_step(&state, &limits, &result),
               "unconfigured particle state rejected before step");
    swa_em_pic3d_destroy(&state);
}

static swa_vec3 proper_velocity_from_velocity(swa_vec3 velocity) {
    const double speed2 = swa_vdot(velocity, velocity);
    const double gamma = 1.0 / sqrt(
        1.0 - speed2 / (SWA_C * SWA_C)
    );
    return swa_vscale(velocity, gamma);
}

static swa_vec3 velocity_from_proper_velocity(swa_vec3 proper) {
    const double gamma = sqrt(
        1.0 + swa_vdot(proper, proper) / (SWA_C * SWA_C)
    );
    return swa_vscale(proper, 1.0 / gamma);
}

static swa_vec3 uniform_lorentz_rhs(
    swa_vec3 proper,
    double charge_to_mass,
    swa_vec3 electric,
    swa_vec3 magnetic
) {
    const swa_vec3 velocity = velocity_from_proper_velocity(proper);
    return swa_vscale(
        swa_vadd(electric, swa_vcross(velocity, magnetic)),
        charge_to_mass
    );
}

static swa_vec3 rk4_uniform_lorentz(
    swa_vec3 velocity_initial,
    double charge_to_mass,
    swa_vec3 electric,
    swa_vec3 magnetic,
    double dt,
    size_t substeps
) {
    swa_vec3 proper = proper_velocity_from_velocity(velocity_initial);
    const double h = dt / (double)substeps;
    size_t step;
    for (step = 0U; step < substeps; ++step) {
        const swa_vec3 k1 = uniform_lorentz_rhs(
            proper, charge_to_mass, electric, magnetic
        );
        const swa_vec3 k2 = uniform_lorentz_rhs(
            swa_vadd(proper, swa_vscale(k1, 0.5 * h)),
            charge_to_mass, electric, magnetic
        );
        const swa_vec3 k3 = uniform_lorentz_rhs(
            swa_vadd(proper, swa_vscale(k2, 0.5 * h)),
            charge_to_mass, electric, magnetic
        );
        const swa_vec3 k4 = uniform_lorentz_rhs(
            swa_vadd(proper, swa_vscale(k3, h)),
            charge_to_mass, electric, magnetic
        );
        proper = swa_vadd(
            proper,
            swa_vscale(
                swa_vadd(
                    swa_vadd(k1, swa_vscale(k2, 2.0)),
                    swa_vadd(swa_vscale(k3, 2.0), k4)
                ),
                h / 6.0
            )
        );
    }
    return velocity_from_proper_velocity(proper);
}

static void test_higuera_cary_pusher(void) {
    size_t i;
    for (i = 1U; i <= 1200U; ++i) {
        const swa_vec3 direction = deterministic_direction(
            i, 2U, 3U, 5U
        );
        const swa_vec3 electric = swa_vscale(
            deterministic_direction(i, 7U, 11U, 13U),
            1.0e3 + 5.0e6 * swa_halton(i, 17U)
        );
        const swa_vec3 magnetic = swa_vscale(
            deterministic_direction(i, 19U, 23U, 29U),
            1.0e-7 + 5.0e-3 * swa_halton(i, 31U)
        );
        const swa_vec3 velocity = swa_vscale(
            direction, 0.75 * SWA_C * swa_halton(i, 37U)
        );
        const double charge =
            (i % 2U == 0U ? 1.0 : -1.0) * 1.602176634e-19;
        const double mass = 9.1093837139e-31 +
            1.0e-26 * swa_halton(i, 41U);
        const double dt = 1.0e-14 +
            2.0e-11 * swa_halton(i, 43U);
        swa_vec3 advanced;
        swa_vec3 reversed;
        swa_vec3 rotated;
        swa_vec3 expected_rotated;
        double reversal_error;
        double rotation_error;
        check_true(swa_em_pic3d_higuera_cary_push(
                       velocity, charge, mass,
                       electric, magnetic, dt, &advanced),
                   "Higuera-Cary deterministic push executes");
        check_true(swa_em_pic3d_higuera_cary_push(
                       negate_vec(advanced), charge, mass,
                       electric, negate_vec(magnetic), dt, &reversed),
                   "Higuera-Cary time-reversed push executes");
        reversal_error = relative_vec(reversed, negate_vec(velocity));
        worst_pusher_reversal = max2(
            worst_pusher_reversal, reversal_error
        );
        check_true(reversal_error < 2.0e-12,
                   "Higuera-Cary push is time reversible");

        check_true(swa_em_pic3d_higuera_cary_push(
                       cycle_xyz(velocity), charge, mass,
                       cycle_xyz(electric), cycle_xyz(magnetic),
                       dt, &rotated),
                   "cyclically rotated Higuera-Cary push executes");
        expected_rotated = cycle_xyz(advanced);
        rotation_error = relative_vec(rotated, expected_rotated);
        worst_pusher_rotation = max2(
            worst_pusher_rotation, rotation_error
        );
        check_true(rotation_error < 2.0e-13,
                   "Higuera-Cary push is cyclic-rotation covariant");
    }

    for (i = 1U; i <= 500U; ++i) {
        const swa_vec3 velocity = swa_vscale(
            deterministic_direction(i, 2U, 5U, 7U),
            0.8 * SWA_C * swa_halton(i, 11U)
        );
        const swa_vec3 electric = swa_vscale(
            deterministic_direction(i, 13U, 17U, 19U),
            2.0e5 * swa_halton(i, 23U)
        );
        const double charge = 2.0e-17;
        const double mass = 3.0e-24;
        const double dt = 1.0e-9 * swa_halton(i, 29U);
        const double speed2 = swa_vdot(velocity, velocity);
        const double gamma = 1.0 / sqrt(
            1.0 - speed2 / (SWA_C * SWA_C)
        );
        const swa_vec3 proper_initial = swa_vscale(velocity, gamma);
        const swa_vec3 proper_expected = swa_vadd(
            proper_initial,
            swa_vscale(electric, charge * dt / mass)
        );
        const double gamma_expected = sqrt(
            1.0 + swa_vdot(proper_expected, proper_expected) /
            (SWA_C * SWA_C)
        );
        const swa_vec3 velocity_expected = swa_vscale(
            proper_expected, 1.0 / gamma_expected
        );
        swa_vec3 actual;
        check_true(swa_em_pic3d_higuera_cary_push(
                       velocity, charge, mass, electric,
                       swa_v3(0.0, 0.0, 0.0), dt, &actual),
                   "pure-electric Higuera-Cary push executes");
        check_true(relative_vec(actual, velocity_expected) < 3.0e-14,
                   "pure-electric proper momentum is analytic");
    }

    {
        size_t case_index;
        size_t refined_cases = 0U;
        for (case_index = 1U; case_index <= 240U; ++case_index) {
            const swa_vec3 velocity = swa_vscale(
                deterministic_direction(case_index, 2U, 3U, 5U),
                0.55 * SWA_C * swa_halton(case_index, 7U)
            );
            const swa_vec3 electric = swa_vscale(
                deterministic_direction(case_index, 11U, 13U, 17U),
                0.015 * SWA_C *
                (0.2 + swa_halton(case_index, 19U))
            );
            const swa_vec3 magnetic = swa_vscale(
                deterministic_direction(case_index, 23U, 29U, 31U),
                0.2 + 0.8 * swa_halton(case_index, 37U)
            );
            const double dt = 0.03 +
                0.07 * swa_halton(case_index, 41U);
            const swa_vec3 reference = rk4_uniform_lorentz(
                velocity, 1.0, electric, magnetic, dt, 128U
            );
            swa_vec3 coarse;
            swa_vec3 half;
            swa_vec3 fine;
            double coarse_error;
            double fine_error;
            check_true(swa_em_pic3d_higuera_cary_push(
                           velocity, 1.0, 1.0,
                           electric, magnetic, dt, &coarse),
                       "coarse Higuera-Cary versus RK4 push executes");
            check_true(swa_em_pic3d_higuera_cary_push(
                           velocity, 1.0, 1.0,
                           electric, magnetic, 0.5 * dt, &half),
                       "first refined Higuera-Cary half-step executes");
            check_true(swa_em_pic3d_higuera_cary_push(
                           half, 1.0, 1.0,
                           electric, magnetic, 0.5 * dt, &fine),
                       "second refined Higuera-Cary half-step executes");
            coarse_error = relative_vec(coarse, reference);
            fine_error = relative_vec(fine, reference);
            check_true(coarse_error < 4.0e-4,
                       "Higuera-Cary agrees with independent RK4 reference");
            if (coarse_error > 2.0e-12) {
                ++refined_cases;
                check_true(fine_error < 0.40 * coarse_error,
                           "Higuera-Cary error decreases at second order");
            }
        }
        check_true(refined_cases > 200U,
                   "RK4 pusher comparison exercises resolved truncation error");
    }

    {
        swa_vec3 velocity = swa_v3(
            0.12 * SWA_C, -0.04 * SWA_C, 0.03 * SWA_C
        );
        const double initial_speed = swa_vnorm(velocity);
        const swa_vec3 magnetic = swa_v3(0.0, 0.0, 0.02);
        size_t step;
        for (step = 0U; step < 50000U; ++step) {
            swa_vec3 next;
            check_true(swa_em_pic3d_higuera_cary_push(
                           velocity, 1.602176634e-19,
                           9.1093837139e-31,
                           swa_v3(0.0, 0.0, 0.0), magnetic,
                           5.0e-15, &next),
                       "long pure-magnetic orbit step executes");
            velocity = next;
        }
        check_true(fabs(swa_vnorm(velocity) - initial_speed) /
                   initial_speed < 2.0e-12,
                   "pure-magnetic Higuera-Cary orbit preserves speed");
    }

    {
        swa_vec3 velocity = swa_v3(0.0, 0.0, 0.0);
        const swa_vec3 electric = swa_v3(0.05 * SWA_C, 0.0, 0.0);
        const swa_vec3 magnetic = swa_v3(0.0, 0.0, 1.0);
        const double dt = 0.01;
        const size_t steps = 20000U;
        swa_vec3 displacement = swa_v3(0.0, 0.0, 0.0);
        size_t step;
        for (step = 0U; step < steps; ++step) {
            swa_vec3 next;
            check_true(swa_em_pic3d_higuera_cary_push(
                           velocity, 1.0, 1.0,
                           electric, magnetic, dt, &next),
                       "crossed-field drift step executes");
            displacement = swa_vadd(
                displacement,
                swa_vscale(swa_vadd(velocity, next), 0.5 * dt)
            );
            velocity = next;
        }
        {
            const swa_vec3 average = swa_vscale(
                displacement, 1.0 / (dt * (double)steps)
            );
            const swa_vec3 expected = swa_v3(
                0.0, -0.05 * SWA_C, 0.0
            );
            check_true(fabs(average.y - expected.y) /
                       fabs(expected.y) < 7.0e-3,
                       "Higuera-Cary resolves relativistic E cross B drift");
            check_true(fabs(average.x) < 8.0e-3 * fabs(expected.y),
                       "crossed-field gyro average has bounded transverse bias");
        }
    }
}

static void configure_neutral_particles(
    swa_em_pic3d *state,
    double mass_kg
) {
    size_t p;
    for (p = 0U; p < state->particle_count; ++p) {
        const double d = (double)p;
        const double sign = p % 2U == 0U ? 1.0 : -1.0;
        check_true(swa_em_pic3d_set_particle(
                       state, p,
                       swa_v3(
                           0.10 + 0.11 * d,
                           0.08 + 0.13 * d,
                           0.06 + 0.17 * d
                       ),
                       swa_v3(
                           sign * 0.005 * SWA_C,
                           (0.001 + 0.0001 * d) * SWA_C,
                           -sign * 0.0007 * SWA_C
                       ),
                       sign * 1.0e-16,
                       mass_kg,
                       1.0),
                   "configure neutral 3D3V particle");
    }
}

static void test_gauss_initialization(void) {
    swa_em_pic3d state;
    double *rho;
    double gauss_relative;
    double mean_x = 0.0;
    double mean_y = 0.0;
    double mean_z = 0.0;
    size_t index;
    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, 5U, 4U, 3U, 6U,
                   1.0, 1.2, 1.4, 1.0e-11),
               "Gauss initialization grid initializes");
    configure_neutral_particles(&state, 2.0e-27);
    check_true(swa_em_pic3d_initialize_gauss_field(
                   &state, swa_v3(3.0e4, -2.0e4, 1.0e4),
                   1.0e-25),
               "periodic neutral Gauss field initializes");
    rho = calloc(state.cell_count, sizeof(double));
    check_true(rho != NULL, "allocate Gauss audit charge");
    if (rho != NULL) {
        check_true(swa_em_pic3d_deposit_charge_density(
                       &state, state.particles,
                       state.particle_count, rho),
                   "deposit charge for Gauss audit");
        check_true(isfinite(
                       swa_em_pic3d_electric_divergence_max_C_m3(
                           &state, rho, &gauss_relative)),
                   "Gauss residual diagnostic executes");
        worst_gauss_relative = max2(
            worst_gauss_relative, gauss_relative
        );
        check_true(gauss_relative < 2.0e-12,
                   "periodic spectral Gauss initialization closes");
        for (index = 0U; index < state.cell_count; ++index) {
            mean_x += state.ex_V_m[index];
            mean_y += state.ey_V_m[index];
            mean_z += state.ez_V_m[index];
        }
        mean_x /= (double)state.cell_count;
        mean_y /= (double)state.cell_count;
        mean_z /= (double)state.cell_count;
        check_true(relative_scalar(mean_x, 3.0e4) < 2.0e-14,
                   "Gauss solve preserves requested mean Ex");
        check_true(relative_scalar(mean_y, -2.0e4) < 2.0e-14,
                   "Gauss solve preserves requested mean Ey");
        check_true(relative_scalar(mean_z, 1.0e4) < 2.0e-14,
                   "Gauss solve preserves requested mean Ez");
    }
    free(rho);
    swa_em_pic3d_destroy(&state);

    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, 4U, 4U, 4U, 1U,
                   1.0, 1.0, 1.0, 1.0e-11),
               "nonneutral Gauss grid initializes");
    check_true(swa_em_pic3d_set_particle(
                   &state, 0U,
                   swa_v3(0.2, 0.3, 0.4),
                   swa_v3(0.0, 0.0, 0.0),
                   1.0e-15, 1.0e-20, 1.0),
               "nonneutral Gauss particle configures");
    check_true(!swa_em_pic3d_initialize_gauss_field(
                   &state, swa_v3(0.0, 0.0, 0.0), 1.0e-20),
               "periodic nonneutral Gauss solve rejected");
    swa_em_pic3d_destroy(&state);
}

static swa_vec3 array_mean_current(
    const swa_em_pic3d *state,
    const double *jx,
    const double *jy,
    const double *jz
) {
    swa_kahan x = {0.0, 0.0};
    swa_kahan y = {0.0, 0.0};
    swa_kahan z = {0.0, 0.0};
    size_t index;
    for (index = 0U; index < state->cell_count; ++index) {
        swa_kahan_add(&x, jx[index]);
        swa_kahan_add(&y, jy[index]);
        swa_kahan_add(&z, jz[index]);
    }
    return swa_v3(
        x.sum / (double)state->cell_count,
        y.sum / (double)state->cell_count,
        z.sum / (double)state->cell_count
    );
}

static void test_transport_and_projection(void) {
    swa_em_pic3d state;
    swa_em_pic3d_particle initial[4];
    swa_em_pic3d_particle final[4];
    swa_em_pic3d_projection_result projection;
    const size_t arrays = 8U;
    double *memory;
    double *rho_initial;
    double *rho_final;
    double *raw_jx;
    double *raw_jy;
    double *raw_jz;
    double *jx;
    double *jy;
    double *jz;
    size_t maximum_segments;
    swa_vec3 transport_mean;
    swa_vec3 raw_mean;
    swa_vec3 corrected_mean;
    size_t p;
    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, 4U, 5U, 3U, 4U,
                   2.0, 3.0, 4.0, 0.2),
               "transport projection grid initializes");
    for (p = 0U; p < 4U; ++p) {
        const double sign = p % 2U == 0U ? 1.0 : -1.0;
        initial[p].position_unwrapped_m = swa_v3(
            0.1 + 0.4 * (double)p,
            0.2 + 0.5 * (double)p,
            0.3 + 0.6 * (double)p
        );
        initial[p].velocity_mps = swa_v3(0.0, 0.0, 0.0);
        initial[p].charge_C = sign * 2.0e-12;
        initial[p].mass_kg = 1.0e-9;
        initial[p].macro_weight = 1.0 + 0.25 * (double)p;
        final[p] = initial[p];
    }
    final[0].position_unwrapped_m = swa_vadd(
        initial[0].position_unwrapped_m,
        swa_v3(2.0 * 2.0 + 0.37, -1.2, 1.1)
    );
    final[1].position_unwrapped_m = swa_vadd(
        initial[1].position_unwrapped_m,
        swa_v3(-1.5, 3.0 * 3.0 + 0.21, -0.8)
    );
    final[2].position_unwrapped_m = swa_vadd(
        initial[2].position_unwrapped_m,
        swa_v3(0.9, -1.1, -2.0 * 4.0 - 0.33)
    );
    final[3].position_unwrapped_m = swa_vadd(
        initial[3].position_unwrapped_m,
        swa_v3(-0.7, 1.3, 0.6)
    );
    memory = calloc(arrays * state.cell_count, sizeof(double));
    check_true(memory != NULL, "allocate projection work arrays");
    if (memory != NULL) {
        rho_initial = memory;
        rho_final = rho_initial + state.cell_count;
        raw_jx = rho_final + state.cell_count;
        raw_jy = raw_jx + state.cell_count;
        raw_jz = raw_jy + state.cell_count;
        jx = raw_jz + state.cell_count;
        jy = jx + state.cell_count;
        jz = jy + state.cell_count;
        check_true(swa_em_pic3d_deposit_charge_density(
                       &state, initial, 4U, rho_initial),
                   "deposit initial transport charge");
        check_true(swa_em_pic3d_deposit_charge_density(
                       &state, final, 4U, rho_final),
                   "deposit final transport charge");
        check_true(swa_em_pic3d_deposit_transport_current(
                       &state, initial, final, 4U,
                       raw_jx, raw_jy, raw_jz,
                       &maximum_segments, &transport_mean),
                   "path-integrated 3D transport current deposits");
        raw_mean = array_mean_current(
            &state, raw_jx, raw_jy, raw_jz
        );
        check_true(relative_vec(raw_mean, transport_mean) < 5.0e-14,
                   "path current preserves complete unwrapped mean transport");
        check_true(maximum_segments > 8U,
                   "multi-cell and multi-wrap path segmentation executes");
        check_true(swa_em_pic3d_project_periodic_current(
                       &state, rho_initial, rho_final,
                       raw_jx, raw_jy, raw_jz,
                       2.0e-11, 2.0e-11,
                       jx, jy, jz, &projection),
                   "periodic minimum-norm current projection executes");
        corrected_mean = array_mean_current(&state, jx, jy, jz);
        worst_projection_continuity = max2(
            worst_projection_continuity,
            projection.continuity_after_relative
        );
        worst_projection_curl = max2(
            worst_projection_curl,
            projection.correction_curl_relative
        );
        check_true(projection.passes,
                   "periodic current projection passes constraints");
        check_true(projection.continuity_after_relative < 2.0e-12,
                   "projected current satisfies 3D discrete continuity");
        check_true(projection.correction_curl_relative < 2.0e-12,
                   "minimum-norm projection correction is curl free");
        check_true(relative_vec(corrected_mean, transport_mean) < 5.0e-13,
                   "longitudinal projection preserves harmonic mean current");
        {
            swa_kahan correction_dot_solenoidal = {0.0, 0.0};
            swa_kahan correction2 = {0.0, 0.0};
            swa_kahan solenoidal2 = {0.0, 0.0};
            size_t i;
            size_t j;
            size_t k;
            for (k = 0U; k < state.nz; ++k) {
                const double z = 2.0 * SWA_PI * (double)k /
                                 (double)state.nz;
                for (j = 0U; j < state.ny; ++j) {
                    const double y = 2.0 * SWA_PI * (double)j /
                                     (double)state.ny;
                    for (i = 0U; i < state.nx; ++i) {
                        const double x = 2.0 * SWA_PI * (double)i /
                                         (double)state.nx;
                        const size_t grid = index3(&state, i, j, k);
                        const double kx = sin(y) + 0.2 * cos(2.0 * z);
                        const double ky = sin(z) - 0.3 * cos(2.0 * x);
                        const double kz = sin(x) + 0.4 * cos(2.0 * y);
                        const double cx = jx[grid] - raw_jx[grid];
                        const double cy = jy[grid] - raw_jy[grid];
                        const double cz = jz[grid] - raw_jz[grid];
                        swa_kahan_add(
                            &correction_dot_solenoidal,
                            cx * kx + cy * ky + cz * kz
                        );
                        swa_kahan_add(
                            &correction2,
                            cx * cx + cy * cy + cz * cz
                        );
                        swa_kahan_add(
                            &solenoidal2,
                            kx * kx + ky * ky + kz * kz
                        );
                    }
                }
            }
            check_true(
                fabs(correction_dot_solenoidal.sum) /
                max2(sqrt(correction2.sum * solenoidal2.sum), DBL_MIN) <
                    2.0e-13,
                "current correction is orthogonal to divergence-free currents"
            );
            check_true(
                correction2.sum + 0.09 * solenoidal2.sum >=
                    correction2.sum,
                "minimum-norm projection cannot improve along a solenoidal perturbation"
            );
        }

        check_true(swa_em_pic3d_deposit_charge_conserving_current(
                       &state, initial, final, 4U,
                       jx, jy, jz, &transport_mean),
                   "symmetric zigzag charge-conserving current deposits");
        corrected_mean = array_mean_current(&state, jx, jy, jz);
        check_true(relative_vec(corrected_mean, transport_mean) < 5.0e-13,
                   "local charge-conserving current preserves unwrapped mean");
        check_true(swa_em_pic3d_project_periodic_current(
                       &state, rho_initial, rho_final,
                       jx, jy, jz,
                       2.0e-11, 2.0e-11,
                       raw_jx, raw_jy, raw_jz, &projection),
                   "spectral oracle audits local charge-conserving current");
        check_true(projection.continuity_before_relative < 2.0e-12,
                   "local zigzag current closes continuity before projection");
        check_true(projection.correction_relative < 2.0e-12,
                   "spectral oracle leaves local charge-conserving current unchanged");

        rho_final[0] += 1.0e-3;
        check_true(swa_em_pic3d_project_periodic_current(
                       &state, rho_initial, rho_final,
                       raw_jx, raw_jy, raw_jz,
                       2.0e-11, 2.0e-11,
                       jx, jy, jz, &projection),
                   "nonzero-mode-corrupted projection returns audit");
        check_true(!projection.zero_mode_passes && !projection.passes,
                   "global charge inconsistency fails current projection");
    }
    free(memory);
    swa_em_pic3d_destroy(&state);
}


static void test_local_current_symmetries(void) {
    enum { PARTICLES = 6 };
    swa_em_pic3d state;
    swa_em_pic3d_particle initial[PARTICLES];
    swa_em_pic3d_particle final[PARTICLES];
    swa_em_pic3d_particle reversed_initial[PARTICLES];
    swa_em_pic3d_particle reversed_final[PARTICLES];
    swa_em_pic3d_particle translated_initial[PARTICLES];
    swa_em_pic3d_particle translated_final[PARTICLES];
    double *memory;
    double *forward_x;
    double *forward_y;
    double *forward_z;
    double *reverse_x;
    double *reverse_y;
    double *reverse_z;
    double *translated_x;
    double *translated_y;
    double *translated_z;
    swa_vec3 mean_forward;
    swa_vec3 mean_reverse;
    swa_vec3 mean_translated;
    const size_t shift_x = 1U;
    const size_t shift_y = 2U;
    const size_t shift_z = 1U;
    size_t p;
    size_t i;
    size_t j;
    size_t k;
    double reversal_error = 0.0;
    double translation_error = 0.0;
    double scale = DBL_MIN;

    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, 5U, 4U, 3U, PARTICLES,
                   2.5, 2.0, 1.5, 0.125),
               "local-current symmetry grid initializes");
    for (p = 0U; p < PARTICLES; ++p) {
        const double sign = p % 2U == 0U ? 1.0 : -1.0;
        initial[p].position_unwrapped_m = swa_v3(
            0.07 + 0.31 * (double)p,
            0.11 + 0.23 * (double)p,
            0.13 + 0.19 * (double)p
        );
        initial[p].velocity_mps = swa_v3(0.0, 0.0, 0.0);
        initial[p].charge_C = sign * (1.0 + 0.1 * (double)p) * 1.0e-12;
        initial[p].mass_kg = 1.0e-9;
        initial[p].macro_weight = 1.0 + 0.2 * (double)p;
        final[p] = initial[p];
        final[p].position_unwrapped_m = swa_vadd(
            initial[p].position_unwrapped_m,
            swa_v3(
                (0.23 + 0.07 * (double)p) * (p % 3U == 0U ? 1.0 : -1.0),
                (0.17 + 0.05 * (double)p) * (p % 3U == 1U ? 1.0 : -1.0),
                (0.19 + 0.03 * (double)p) * (p % 3U == 2U ? 1.0 : -1.0)
            )
        );
        translated_initial[p] = initial[p];
        translated_final[p] = final[p];
        translated_initial[p].position_unwrapped_m = swa_vadd(
            translated_initial[p].position_unwrapped_m,
            swa_v3(
                (double)shift_x * state.length_x_m / (double)state.nx,
                (double)shift_y * state.length_y_m / (double)state.ny,
                (double)shift_z * state.length_z_m / (double)state.nz
            )
        );
        translated_final[p].position_unwrapped_m = swa_vadd(
            translated_final[p].position_unwrapped_m,
            swa_v3(
                (double)shift_x * state.length_x_m / (double)state.nx,
                (double)shift_y * state.length_y_m / (double)state.ny,
                (double)shift_z * state.length_z_m / (double)state.nz
            )
        );
    }
    for (p = 0U; p < PARTICLES; ++p) {
        reversed_initial[p] = final[PARTICLES - 1U - p];
        reversed_final[p] = initial[PARTICLES - 1U - p];
    }
    memory = calloc(9U * state.cell_count, sizeof(double));
    check_true(memory != NULL, "allocate local-current symmetry arrays");
    if (memory != NULL) {
        forward_x = memory;
        forward_y = forward_x + state.cell_count;
        forward_z = forward_y + state.cell_count;
        reverse_x = forward_z + state.cell_count;
        reverse_y = reverse_x + state.cell_count;
        reverse_z = reverse_y + state.cell_count;
        translated_x = reverse_z + state.cell_count;
        translated_y = translated_x + state.cell_count;
        translated_z = translated_y + state.cell_count;
        check_true(swa_em_pic3d_deposit_charge_conserving_current(
                       &state, initial, final, PARTICLES,
                       forward_x, forward_y, forward_z, &mean_forward),
                   "forward local charge-conserving current deposits");
        check_true(swa_em_pic3d_deposit_charge_conserving_current(
                       &state, reversed_initial, reversed_final, PARTICLES,
                       reverse_x, reverse_y, reverse_z, &mean_reverse),
                   "reversed and reordered local current deposits");
        check_true(swa_em_pic3d_deposit_charge_conserving_current(
                       &state, translated_initial, translated_final, PARTICLES,
                       translated_x, translated_y, translated_z,
                       &mean_translated),
                   "integer-cell translated local current deposits");
        for (k = 0U; k < state.nz; ++k) {
            const size_t ks = (k + shift_z) % state.nz;
            for (j = 0U; j < state.ny; ++j) {
                const size_t js = (j + shift_y) % state.ny;
                for (i = 0U; i < state.nx; ++i) {
                    const size_t is = (i + shift_x) % state.nx;
                    const size_t original = index3(&state, i, j, k);
                    const size_t shifted = index3(&state, is, js, ks);
                    reversal_error = max3(
                        reversal_error,
                        fabs(forward_x[original] + reverse_x[original]),
                        max2(
                            fabs(forward_y[original] + reverse_y[original]),
                            fabs(forward_z[original] + reverse_z[original])
                        )
                    );
                    translation_error = max3(
                        translation_error,
                        fabs(forward_x[original] - translated_x[shifted]),
                        max2(
                            fabs(forward_y[original] - translated_y[shifted]),
                            fabs(forward_z[original] - translated_z[shifted])
                        )
                    );
                    scale = max3(
                        scale,
                        fabs(forward_x[original]),
                        max2(fabs(forward_y[original]),
                             fabs(forward_z[original]))
                    );
                }
            }
        }
        check_true(reversal_error / scale < 2.0e-13,
                   "local current is antisymmetric under trajectory reversal and particle reordering");
        check_true(translation_error / scale < 2.0e-13,
                   "local current is covariant under integer-cell periodic translation");
        check_true(relative_vec(mean_reverse, negate_vec(mean_forward)) <
                   2.0e-13,
                   "reversed local current negates complete harmonic transport");
        check_true(relative_vec(mean_translated, mean_forward) < 2.0e-13,
                   "periodic translation preserves local-current harmonic mode");
    }
    free(memory);
    swa_em_pic3d_destroy(&state);
}

static void fill_divergence_free_fields(swa_em_pic3d *state) {
    size_t i;
    size_t j;
    size_t k;
    for (k = 0U; k < state->nz; ++k) {
        const double z = 2.0 * SWA_PI * (double)k /
                         (double)state->nz;
        for (j = 0U; j < state->ny; ++j) {
            const double y = 2.0 * SWA_PI * (double)j /
                             (double)state->ny;
            for (i = 0U; i < state->nx; ++i) {
                const double x = 2.0 * SWA_PI * (double)i /
                                 (double)state->nx;
                const size_t index = index3(state, i, j, k);
                state->ex_V_m[index] =
                    2.0e3 * (sin(2.0 * y) + 0.3 * cos(3.0 * z));
                state->ey_V_m[index] =
                    1.5e3 * (sin(3.0 * z) - 0.2 * cos(2.0 * x));
                state->ez_V_m[index] =
                    1.2e3 * (sin(2.0 * x) + 0.4 * cos(3.0 * y));
                state->bx_T[index] =
                    2.0e-6 * (sin(y) + 0.3 * cos(2.0 * z));
                state->by_T[index] =
                    1.5e-6 * (sin(z) - 0.2 * cos(2.0 * x));
                state->bz_T[index] =
                    1.2e-6 * (sin(x) + 0.4 * cos(2.0 * y));
            }
        }
    }
}

static void test_discrete_maxwell_identities(void) {
    swa_em_pic3d state;
    swa_em_pic3d_result result;
    swa_em_pic3d_limits limits = permissive_limits();
    double div_b_relative;
    double *rho;
    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, 7U, 6U, 5U, 0U,
                   1.0, 1.1, 1.2, 1.0e-11),
               "3D Maxwell identity grid initializes");
    fill_divergence_free_fields(&state);
    rho = calloc(state.cell_count, sizeof(double));
    check_true(rho != NULL, "allocate Maxwell identity charge");
    if (rho != NULL) {
        check_true(swa_em_pic3d_magnetic_divergence_max_T_m(
                       &state, &div_b_relative) >= 0.0,
                   "initial magnetic divergence audit executes");
        check_true(div_b_relative < 5.0e-15,
                   "manufactured 3D magnetic field is divergence free");
        check_true(swa_em_pic3d_step(&state, &limits, &result),
                   "3D vacuum Maxwell identity step executes");
        worst_gauss_relative = max2(
            worst_gauss_relative, result.final_gauss_relative_max
        );
        worst_div_b_relative = max2(
            worst_div_b_relative,
            result.final_magnetic_divergence_relative
        );
        check_true(result.continuity_passes,
                   "vacuum continuity is exact");
        check_true(result.final_gauss_relative_max < 2.0e-13,
                   "discrete divergence of Ampere curl vanishes");
        check_true(result.final_magnetic_divergence_relative < 2.0e-13,
                   "discrete divergence of Faraday curl vanishes");
        check_true(result.passes,
                   "divergence-free 3D vacuum state passes all gates");

        state.bx_T[index3(&state, 0U, 0U, 0U)] += 1.0e-4;
        check_true(swa_em_pic3d_magnetic_divergence_max_T_m(
                       &state, &div_b_relative) > 0.0,
                   "corrupted 3D magnetic field divergence executes");
        check_true(div_b_relative > 1.0e-3,
                   "localized magnetic monopole corruption is detected");
    }
    free(rho);
    swa_em_pic3d_destroy(&state);
}

static void initialize_plane_wave(
    swa_em_pic3d *state,
    int axis,
    unsigned mode,
    double amplitude,
    int magnetic_sign
) {
    const double dx = state->length_x_m / (double)state->nx;
    const double dy = state->length_y_m / (double)state->ny;
    const double dz = state->length_z_m / (double)state->nz;
    const double length = axis == 0 ? state->length_x_m :
        (axis == 1 ? state->length_y_m : state->length_z_m);
    size_t i;
    size_t j;
    size_t k;
    for (k = 0U; k < state->nz; ++k) {
        for (j = 0U; j < state->ny; ++j) {
            for (i = 0U; i < state->nx; ++i) {
                const size_t index = index3(state, i, j, k);
                if (axis == 0) {
                    const double electric_x = ((double)i + 0.5) * dx;
                    const double magnetic_x = (double)i * dx;
                    state->ey_V_m[index] = amplitude * sin(
                        2.0 * SWA_PI * (double)mode *
                        electric_x / length
                    );
                    state->bz_T[index] =
                        (double)magnetic_sign * amplitude / SWA_C * sin(
                            2.0 * SWA_PI * (double)mode *
                            magnetic_x / length
                        );
                } else if (axis == 1) {
                    const double electric_y = ((double)j + 0.5) * dy;
                    const double magnetic_y = (double)j * dy;
                    state->ez_V_m[index] = amplitude * sin(
                        2.0 * SWA_PI * (double)mode *
                        electric_y / length
                    );
                    state->bx_T[index] =
                        (double)magnetic_sign * amplitude / SWA_C * sin(
                            2.0 * SWA_PI * (double)mode *
                            magnetic_y / length
                        );
                } else {
                    const double electric_z = ((double)k + 0.5) * dz;
                    const double magnetic_z = (double)k * dz;
                    state->ex_V_m[index] = amplitude * sin(
                        2.0 * SWA_PI * (double)mode *
                        electric_z / length
                    );
                    state->by_T[index] =
                        (double)magnetic_sign * amplitude / SWA_C * sin(
                            2.0 * SWA_PI * (double)mode *
                            magnetic_z / length
                        );
                }
            }
        }
    }
}

static double plane_wave_error(
    const swa_em_pic3d *state,
    int axis,
    unsigned mode,
    double amplitude,
    double time_s
) {
    const double dx = state->length_x_m / (double)state->nx;
    const double dy = state->length_y_m / (double)state->ny;
    const double dz = state->length_z_m / (double)state->nz;
    const double length = axis == 0 ? state->length_x_m :
        (axis == 1 ? state->length_y_m : state->length_z_m);
    swa_kahan error2 = {0.0, 0.0};
    swa_kahan exact2 = {0.0, 0.0};
    size_t i;
    size_t j;
    size_t k;
    for (k = 0U; k < state->nz; ++k) {
        for (j = 0U; j < state->ny; ++j) {
            for (i = 0U; i < state->nx; ++i) {
                const size_t index = index3(state, i, j, k);
                double coordinate;
                double actual;
                double exact;
                if (axis == 0) {
                    coordinate = ((double)i + 0.5) * dx;
                    actual = state->ey_V_m[index];
                } else if (axis == 1) {
                    coordinate = ((double)j + 0.5) * dy;
                    actual = state->ez_V_m[index];
                } else {
                    coordinate = ((double)k + 0.5) * dz;
                    actual = state->ex_V_m[index];
                }
                exact = amplitude * sin(
                    2.0 * SWA_PI * (double)mode *
                    (coordinate - SWA_C * time_s) / length
                );
                swa_kahan_add(&error2, (actual - exact) * (actual - exact));
                swa_kahan_add(&exact2, exact * exact);
            }
        }
    }
    return sqrt(error2.sum / max2(exact2.sum, DBL_MIN));
}

static double run_reduced_plane_wave(size_t cells) {
    swa_em_pic3d state;
    swa_em_pic3d_result result;
    swa_em_pic3d_limits limits = permissive_limits();
    const double length = 1.0;
    const double dx = length / (double)cells;
    const double dt = 0.4 * dx / SWA_C;
    const size_t steps = cells / 2U;
    size_t step;
    double error;
    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, cells, 1U, 1U, 0U,
                   length, 1000.0, 1000.0, dt),
               "reduced 3D plane-wave grid initializes");
    initialize_plane_wave(&state, 0, 1U, 1.0, 1);
    for (step = 0U; step < steps; ++step) {
        check_true(swa_em_pic3d_step(&state, &limits, &result),
                   "reduced 3D plane-wave step executes");
    }
    error = plane_wave_error(
        &state, 0, 1U, 1.0, (double)steps * dt
    );
    worst_plane_wave_error = max2(worst_plane_wave_error, error);
    check_true(result.final_gauss_relative_max < 1.0e-14,
               "reduced plane wave preserves Gauss constraint");
    check_true(result.final_magnetic_divergence_relative < 1.0e-14,
               "reduced plane wave preserves magnetic divergence");
    swa_em_pic3d_destroy(&state);
    return error;
}

static void test_plane_wave_convergence_and_symmetry(void) {
    const double error_16 = run_reduced_plane_wave(16U);
    const double error_32 = run_reduced_plane_wave(32U);
    const double error_64 = run_reduced_plane_wave(64U);
    const double order_16_32 = log(error_16 / error_32) / log(2.0);
    const double order_32_64 = log(error_32 / error_64) / log(2.0);
    double directional_error[3];
    int axis;
    check_true(error_64 < error_32 && error_32 < error_16,
               "3D Yee plane-wave error decreases monotonically");
    check_true(order_16_32 > 1.90 && order_16_32 < 2.10,
               "coarse-to-medium 3D Yee convergence is second order");
    check_true(order_32_64 > 1.90 && order_32_64 < 2.10,
               "medium-to-fine 3D Yee convergence is second order");

    for (axis = 0; axis < 3; ++axis) {
        swa_em_pic3d state;
        swa_em_pic3d_result result;
        swa_em_pic3d_limits limits = permissive_limits();
        const size_t cells = 12U;
        const double dx = 1.0 / (double)cells;
        const double dt = 0.30 * dx / (SWA_C * sqrt(3.0));
        const size_t steps = 9U;
        size_t step;
        memset(&state, 0, sizeof(state));
        check_true(swa_em_pic3d_init(
                       &state, cells, cells, cells, 0U,
                       1.0, 1.0, 1.0, dt),
                   "cubic directional plane-wave grid initializes");
        initialize_plane_wave(&state, axis, 2U, 2.5e4, 1);
        for (step = 0U; step < steps; ++step) {
            check_true(swa_em_pic3d_step(&state, &limits, &result),
                       "cubic directional plane-wave step executes");
        }
        directional_error[axis] = plane_wave_error(
            &state, axis, 2U, 2.5e4, (double)steps * dt
        );
        check_true(result.passes,
                   "cubic directional plane wave passes gates");
        swa_em_pic3d_destroy(&state);
    }
    check_true(max3(
                   directional_error[0],
                   directional_error[1],
                   directional_error[2]) /
               max2(fmin(directional_error[0],
                         fmin(directional_error[1],
                              directional_error[2])), DBL_MIN) < 1.05,
               "3D Yee propagation is cyclic-axis symmetric");

    {
        swa_em_pic3d correct;
        swa_em_pic3d wrong;
        swa_em_pic3d_result result;
        swa_em_pic3d_limits limits = permissive_limits();
        const size_t cells = 24U;
        const double dx = 1.0 / (double)cells;
        const double dt = 0.35 * dx / SWA_C;
        const size_t steps = 7U;
        size_t step;
        double correct_error;
        double wrong_error;
        memset(&correct, 0, sizeof(correct));
        memset(&wrong, 0, sizeof(wrong));
        check_true(swa_em_pic3d_init(
                       &correct, cells, 1U, 1U, 0U,
                       1.0, 1000.0, 1000.0, dt),
                   "correct-sign plane-wave grid initializes");
        check_true(swa_em_pic3d_init(
                       &wrong, cells, 1U, 1U, 0U,
                       1.0, 1000.0, 1000.0, dt),
                   "wrong-sign plane-wave grid initializes");
        initialize_plane_wave(&correct, 0, 3U, 1.0, 1);
        initialize_plane_wave(&wrong, 0, 3U, 1.0, -1);
        for (step = 0U; step < steps; ++step) {
            check_true(swa_em_pic3d_step(&correct, &limits, &result),
                       "correct-sign plane-wave step executes");
            check_true(swa_em_pic3d_step(&wrong, &limits, &result),
                       "wrong-sign plane-wave step executes");
        }
        correct_error = plane_wave_error(
            &correct, 0, 3U, 1.0, (double)steps * dt
        );
        wrong_error = plane_wave_error(
            &wrong, 0, 3U, 1.0, (double)steps * dt
        );
        check_true(wrong_error > 8.0 * correct_error,
                   "wrong 3D magnetic polarization is rejected analytically");
        swa_em_pic3d_destroy(&correct);
        swa_em_pic3d_destroy(&wrong);
    }
}

static void test_vacuum_time_reversal(void) {
    swa_em_pic3d state;
    swa_em_pic3d_result result;
    swa_em_pic3d_limits limits = permissive_limits();
    double *initial;
    double *initial_ex;
    double *initial_ey;
    double *initial_ez;
    double *initial_bx;
    double *initial_by;
    double *initial_bz;
    double difference2 = 0.0;
    double scale2 = 0.0;
    size_t index;
    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, 6U, 5U, 4U, 0U,
                   1.0, 1.2, 1.4, 1.0e-11),
               "vacuum reversal grid initializes");
    fill_divergence_free_fields(&state);
    initial = calloc(6U * state.cell_count, sizeof(double));
    check_true(initial != NULL, "allocate vacuum reversal snapshot");
    if (initial != NULL) {
        initial_ex = initial;
        initial_ey = initial_ex + state.cell_count;
        initial_ez = initial_ey + state.cell_count;
        initial_bx = initial_ez + state.cell_count;
        initial_by = initial_bx + state.cell_count;
        initial_bz = initial_by + state.cell_count;
        memcpy(initial_ex, state.ex_V_m,
               state.cell_count * sizeof(double));
        memcpy(initial_ey, state.ey_V_m,
               state.cell_count * sizeof(double));
        memcpy(initial_ez, state.ez_V_m,
               state.cell_count * sizeof(double));
        memcpy(initial_bx, state.bx_T,
               state.cell_count * sizeof(double));
        memcpy(initial_by, state.by_T,
               state.cell_count * sizeof(double));
        memcpy(initial_bz, state.bz_T,
               state.cell_count * sizeof(double));
        check_true(swa_em_pic3d_step(&state, &limits, &result),
                   "forward vacuum reversal step executes");
        for (index = 0U; index < state.cell_count; ++index) {
            state.bx_T[index] = -state.bx_T[index];
            state.by_T[index] = -state.by_T[index];
            state.bz_T[index] = -state.bz_T[index];
        }
        check_true(swa_em_pic3d_step(&state, &limits, &result),
                   "time-reversed vacuum step executes");
        for (index = 0U; index < state.cell_count; ++index) {
            const double values[6] = {
                state.ex_V_m[index] - initial_ex[index],
                state.ey_V_m[index] - initial_ey[index],
                state.ez_V_m[index] - initial_ez[index],
                state.bx_T[index] + initial_bx[index],
                state.by_T[index] + initial_by[index],
                state.bz_T[index] + initial_bz[index]
            };
            const double scales[6] = {
                initial_ex[index], initial_ey[index], initial_ez[index],
                initial_bx[index], initial_by[index], initial_bz[index]
            };
            size_t component;
            for (component = 0U; component < 6U; ++component) {
                difference2 += values[component] * values[component];
                scale2 += scales[component] * scales[component];
            }
        }
        check_true(sqrt(difference2 / max2(scale2, DBL_MIN)) < 2.0e-12,
                   "symmetric 3D Maxwell update is time reversible");
    }
    free(initial);
    swa_em_pic3d_destroy(&state);
}

static void cycle_state_into(
    const swa_em_pic3d *source,
    swa_em_pic3d *destination
) {
    size_t i;
    size_t j;
    size_t k;
    size_t p;
    for (k = 0U; k < source->nz; ++k) {
        for (j = 0U; j < source->ny; ++j) {
            for (i = 0U; i < source->nx; ++i) {
                const size_t old_index = index3(source, i, j, k);
                const size_t new_index = index3(destination, j, k, i);
                destination->ex_V_m[new_index] = source->ey_V_m[old_index];
                destination->ey_V_m[new_index] = source->ez_V_m[old_index];
                destination->ez_V_m[new_index] = source->ex_V_m[old_index];
                destination->bx_T[new_index] = source->by_T[old_index];
                destination->by_T[new_index] = source->bz_T[old_index];
                destination->bz_T[new_index] = source->bx_T[old_index];
            }
        }
    }
    for (p = 0U; p < source->particle_count; ++p) {
        destination->particles[p] = source->particles[p];
        destination->particles[p].position_unwrapped_m = cycle_xyz(
            source->particles[p].position_unwrapped_m
        );
        destination->particles[p].velocity_mps = cycle_xyz(
            source->particles[p].velocity_mps
        );
    }
}

static double cyclic_state_relative(
    const swa_em_pic3d *source,
    const swa_em_pic3d *cycled
) {
    swa_kahan difference2 = {0.0, 0.0};
    swa_kahan scale2 = {0.0, 0.0};
    size_t i;
    size_t j;
    size_t k;
    for (k = 0U; k < source->nz; ++k) {
        for (j = 0U; j < source->ny; ++j) {
            for (i = 0U; i < source->nx; ++i) {
                const size_t old_index = index3(source, i, j, k);
                const size_t new_index = index3(cycled, j, k, i);
                const double expected[6] = {
                    source->ey_V_m[old_index],
                    source->ez_V_m[old_index],
                    source->ex_V_m[old_index],
                    source->by_T[old_index],
                    source->bz_T[old_index],
                    source->bx_T[old_index]
                };
                const double actual[6] = {
                    cycled->ex_V_m[new_index],
                    cycled->ey_V_m[new_index],
                    cycled->ez_V_m[new_index],
                    cycled->bx_T[new_index],
                    cycled->by_T[new_index],
                    cycled->bz_T[new_index]
                };
                size_t component;
                for (component = 0U; component < 6U; ++component) {
                    const double delta = actual[component] - expected[component];
                    swa_kahan_add(&difference2, delta * delta);
                    swa_kahan_add(&scale2, expected[component] * expected[component]);
                }
            }
        }
    }
    return sqrt(difference2.sum / max2(scale2.sum, DBL_MIN));
}

static void test_coupled_cyclic_covariance(void) {
    swa_em_pic3d original;
    swa_em_pic3d cycled;
    swa_em_pic3d_result original_result;
    swa_em_pic3d_result cycled_result;
    swa_em_pic3d_limits limits = permissive_limits();
    const size_t cells = 4U;
    const double dx = 1.0 / (double)cells;
    const double dt = 0.06 * dx / (SWA_C * sqrt(3.0));
    size_t p;
    size_t index;
    memset(&original, 0, sizeof(original));
    memset(&cycled, 0, sizeof(cycled));
    check_true(swa_em_pic3d_init(
                   &original, cells, cells, cells, 6U,
                   1.0, 1.0, 1.0, dt),
               "cyclic-covariance original grid initializes");
    check_true(swa_em_pic3d_init(
                   &cycled, cells, cells, cells, 6U,
                   1.0, 1.0, 1.0, dt),
               "cyclic-covariance rotated grid initializes");
    for (p = 0U; p < original.particle_count; ++p) {
        const double d = (double)p;
        const double sign = p % 2U == 0U ? 1.0 : -1.0;
        check_true(swa_em_pic3d_set_particle(
                       &original, p,
                       swa_v3(
                           0.09 + 0.13 * d,
                           0.07 + 0.11 * d,
                           0.05 + 0.09 * d
                       ),
                       swa_v3(
                           sign * (0.004 + 0.0002 * d) * SWA_C,
                           (0.0012 - 0.0001 * d) * SWA_C,
                           -sign * (0.0008 + 0.00005 * d) * SWA_C
                       ),
                       sign * 8.0e-17,
                       3.0e-27,
                       1.0),
                   "cyclic-covariance particle configures");
    }
    check_true(swa_em_pic3d_initialize_gauss_field(
                   &original, swa_v3(1.7e5, -0.9e5, 0.6e5),
                   1.0e-25),
               "cyclic-covariance Gauss field initializes");
    for (index = 0U; index < original.cell_count; ++index) {
        original.bx_T[index] = 1.1e-3;
        original.by_T[index] = -0.7e-3;
        original.bz_T[index] = 0.4e-3;
    }
    cycle_state_into(&original, &cycled);
    check_true(swa_em_pic3d_step(
                   &original, &limits, &original_result),
               "cyclic-covariance original step executes");
    check_true(swa_em_pic3d_step(
                   &cycled, &limits, &cycled_result),
               "cyclic-covariance rotated step executes");
    check_true(cyclic_state_relative(&original, &cycled) < 3.0e-12,
               "full 3D3V step is cyclic-axis covariant");
    for (p = 0U; p < original.particle_count; ++p) {
        check_true(relative_vec(
                       cycled.particles[p].position_unwrapped_m,
                       cycle_xyz(original.particles[p].position_unwrapped_m)
                   ) < 3.0e-13,
                   "cycled 3D3V particle position matches");
        check_true(relative_vec(
                       cycled.particles[p].velocity_mps,
                       cycle_xyz(original.particles[p].velocity_mps)
                   ) < 3.0e-12,
                   "cycled 3D3V particle velocity matches");
    }
    check_true(original_result.energy_passes &&
               cycled_result.energy_passes,
               "cycled 3D3V energy audits both pass physical gates");
    check_true(relative_scalar(
                   original_result.current_correction_relative,
                   cycled_result.current_correction_relative) < 2.0e-3,
               "cycled raw-versus-local current diagnostic is axis stable");
    swa_em_pic3d_destroy(&original);
    swa_em_pic3d_destroy(&cycled);
}

static double coupled_dt_case(double courant_factor) {
    swa_em_pic3d state;
    swa_em_pic3d_result result;
    swa_em_pic3d_limits limits = permissive_limits();
    const double dx = 1.0 / 4.0;
    const double dy = 1.2 / 5.0;
    const double dz = 1.4 / 6.0;
    const double dt = courant_factor / (
        SWA_C * sqrt(
            1.0 / (dx * dx) +
            1.0 / (dy * dy) +
            1.0 / (dz * dz)
        )
    );
    size_t index;
    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, 4U, 5U, 6U, 8U,
                   1.0, 1.2, 1.4, dt),
               "coupled refinement grid initializes");
    configure_neutral_particles(&state, 2.0e-27);
    check_true(swa_em_pic3d_initialize_gauss_field(
                   &state, swa_v3(2.0e5, -1.0e5, 0.5e5),
                   1.0e-25),
               "coupled refinement Gauss field initializes");
    for (index = 0U; index < state.cell_count; ++index) {
        state.bx_T[index] = 2.0e-3;
        state.by_T[index] = -1.0e-3;
        state.bz_T[index] = 0.5e-3;
    }
    check_true(swa_em_pic3d_step(&state, &limits, &result),
               "coupled refinement step executes");
    swa_em_pic3d_destroy(&state);
    return result.particle_work_relative;
}

static void test_coupled_3d3v_reference(void) {
    swa_em_pic3d state;
    swa_em_pic3d_result result;
    swa_em_pic3d_limits limits = strict_coupled_limits();
    const double dx = 1.0 / 4.0;
    const double dy = 1.2 / 5.0;
    const double dz = 1.4 / 6.0;
    const double dt = 0.08 / (
        SWA_C * sqrt(
            1.0 / (dx * dx) +
            1.0 / (dy * dy) +
            1.0 / (dz * dz)
        )
    );
    size_t index;
    size_t step;
    memset(&state, 0, sizeof(state));
    memset(&receipt_result, 0, sizeof(receipt_result));
    check_true(swa_em_pic3d_init(
                   &state, 4U, 5U, 6U, 8U,
                   1.0, 1.2, 1.4, dt),
               "coupled periodic 3D3V grid initializes");
    configure_neutral_particles(&state, 2.0e-27);
    check_true(swa_em_pic3d_initialize_gauss_field(
                   &state, swa_v3(2.0e5, -1.0e5, 0.5e5),
                   1.0e-25),
               "coupled periodic 3D3V Gauss field initializes");
    for (index = 0U; index < state.cell_count; ++index) {
        state.bx_T[index] = 2.0e-3;
        state.by_T[index] = -1.0e-3;
        state.bz_T[index] = 0.5e-3;
    }
    for (step = 0U; step < 12U; ++step) {
        check_true(swa_em_pic3d_step(&state, &limits, &result),
                   "coupled periodic 3D3V step executes");
        check_true(result.passes,
                   "coupled periodic 3D3V step passes declared gates");
        check_true(result.continuity_passes &&
                   result.initial_gauss_passes &&
                   result.final_gauss_passes,
                   "coupled 3D3V charge and Gauss constraints close");
        check_true(result.initial_magnetic_divergence_passes &&
                   result.final_magnetic_divergence_passes,
                   "coupled 3D3V magnetic divergence closes");
        check_true(result.local_current_passes &&
                   result.spectral_oracle_passes &&
                   result.current_projection_passes,
                   "coupled 3D3V local current and independent spectral oracle pass");
        check_true(result.mean_current_passes,
                   "coupled 3D3V mean current preserves unwrapped transport");
        check_true(result.nonlinear_passes && result.subluminal_passes,
                   "coupled 3D3V midpoint solve remains converged and subluminal");
        worst_coupled_energy = max2(
            worst_coupled_energy, result.total_energy_relative
        );
        worst_coupled_momentum = max2(
            worst_coupled_momentum, result.total_momentum_relative
        );
        worst_coupled_projection = max2(
            worst_coupled_projection,
            result.local_current_difference_relative
        );
        worst_spectral_oracle_difference = max2(
            worst_spectral_oracle_difference,
            result.spectral_oracle_difference_relative
        );
        receipt_result = result;
    }
    {
        FILE *step_fp = fopen("output/em_pic3d_step_receipt.json", "w");
        check_true(step_fp != NULL, "open 3D3V solver-step receipt");
        if (step_fp != NULL) {
            check_true(swa_em_pic3d_write_receipt(
                           step_fp, &state, &limits, &receipt_result),
                       "write 3D3V solver-step receipt");
            (void)fclose(step_fp);
        }
    }
    swa_em_pic3d_destroy(&state);

    {
        const double coarse = coupled_dt_case(0.16);
        const double medium = coupled_dt_case(0.08);
        const double fine = coupled_dt_case(0.04);
        check_true(isfinite(coarse) && isfinite(medium) && isfinite(fine),
                   "coupled timestep-refinement defects are finite");
        check_true(medium < 0.35 * coarse,
                   "halving timestep strongly reduces particle-work defect");
        check_true(fine < 0.35 * medium,
                   "second timestep halving strongly reduces particle-work defect");
    }
}

static void test_fail_closed_gates(void) {
    swa_em_pic3d state;
    swa_em_pic3d_result result;
    swa_em_pic3d_limits limits = permissive_limits();
    size_t index;
    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, 4U, 4U, 4U, 0U,
                   1.0, 1.0, 1.0,
                   1.1 * (1.0 / 4.0) /
                   (SWA_C * sqrt(3.0))),
               "over-CFL 3D grid initializes");
    limits.maximum_courant = 0.8;
    check_true(swa_em_pic3d_step(&state, &limits, &result),
               "over-CFL 3D step returns finite audit");
    check_true(!result.courant_passes && !result.passes,
               "over-CFL 3D Maxwell state fails closed");
    swa_em_pic3d_destroy(&state);

    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, 4U, 5U, 6U, 8U,
                   1.0, 1.2, 1.4, 3.71026e-11),
               "projection-promotion failure grid initializes");
    configure_neutral_particles(&state, 2.0e-27);
    check_true(swa_em_pic3d_initialize_gauss_field(
                   &state, swa_v3(2.0e5, -1.0e5, 0.5e5),
                   1.0e-25),
               "projection-promotion failure Gauss field initializes");
    for (index = 0U; index < state.cell_count; ++index) {
        state.bx_T[index] = 2.0e-3;
        state.by_T[index] = -1.0e-3;
        state.bz_T[index] = 0.5e-3;
    }
    limits = permissive_limits();
    limits.local_current_difference_relative_tolerance = 0.10;
    check_true(swa_em_pic3d_step(&state, &limits, &result),
               "strict local-current promotion audit executes");
    check_true(result.continuity_passes &&
               !result.current_projection_passes && !result.passes,
               "large path-current versus local-current discrepancy blocks promotion");
    swa_em_pic3d_destroy(&state);

    memset(&state, 0, sizeof(state));
    check_true(swa_em_pic3d_init(
                   &state, 5U, 4U, 3U, 6U,
                   1.0, 1.2, 1.4, 1.0e-11),
               "Gauss-corruption failure grid initializes");
    configure_neutral_particles(&state, 2.0e-27);
    check_true(swa_em_pic3d_initialize_gauss_field(
                   &state, swa_v3(0.0, 0.0, 0.0), 1.0e-25),
               "Gauss-corruption reference initializes");
    state.ex_V_m[index3(&state, 0U, 0U, 0U)] += 1.0e6;
    limits = permissive_limits();
    check_true(swa_em_pic3d_step(&state, &limits, &result),
               "Gauss-corrupted step returns audit");
    check_true(!result.initial_gauss_passes &&
               !result.final_gauss_passes && !result.passes,
               "initial Gauss corruption remains visible and fails closed");
    swa_em_pic3d_destroy(&state);
}


static void capture_field_state(
    const swa_em_pic3d *state,
    double *snapshot
) {
    const size_t bytes = state->cell_count * sizeof(double);
    memcpy(snapshot, state->ex_V_m, bytes);
    memcpy(snapshot + state->cell_count, state->ey_V_m, bytes);
    memcpy(snapshot + 2U * state->cell_count, state->ez_V_m, bytes);
    memcpy(snapshot + 3U * state->cell_count, state->bx_T, bytes);
    memcpy(snapshot + 4U * state->cell_count, state->by_T, bytes);
    memcpy(snapshot + 5U * state->cell_count, state->bz_T, bytes);
}

static int field_state_matches_snapshot(
    const swa_em_pic3d *state,
    const double *snapshot
) {
    const size_t bytes = state->cell_count * sizeof(double);
    return memcmp(state->ex_V_m, snapshot, bytes) == 0 &&
           memcmp(state->ey_V_m, snapshot + state->cell_count, bytes) == 0 &&
           memcmp(state->ez_V_m, snapshot + 2U * state->cell_count, bytes) == 0 &&
           memcmp(state->bx_T, snapshot + 3U * state->cell_count, bytes) == 0 &&
           memcmp(state->by_T, snapshot + 4U * state->cell_count, bytes) == 0 &&
           memcmp(state->bz_T, snapshot + 5U * state->cell_count, bytes) == 0;
}

static void test_transactional_step_acceptance(void) {
    swa_em_pic3d state;
    swa_em_pic3d_result result;
    swa_em_pic3d_limits limits = permissive_limits();
    double *snapshot = NULL;
    const size_t rejected_cells = 4U;
    const double rejected_dx = 1.0 / (double)rejected_cells;
    const double rejected_dt = 1.05 * rejected_dx /
        (SWA_C * sqrt(3.0));
    size_t index;

    memset(&state, 0, sizeof(state));
    check_true(
        swa_em_pic3d_init(
            &state,
            rejected_cells, rejected_cells, rejected_cells, 0U,
            1.0, 1.0, 1.0, rejected_dt
        ),
        "transactional rejected-step grid initializes"
    );
    fill_divergence_free_fields(&state);
    snapshot = calloc(6U * state.cell_count, sizeof(double));
    check_true(snapshot != NULL,
               "allocate rejected-step field snapshot");
    if (snapshot != NULL) {
        capture_field_state(&state, snapshot);
        limits.maximum_courant = 0.8;
        check_true(
            swa_em_pic3d_step(&state, &limits, &result),
            "rejected transactional step returns an audit receipt"
        );
        check_true(!result.passes && !result.state_committed,
                   "failed scientific gates do not commit candidate state");
        check_true(field_state_matches_snapshot(&state, snapshot),
                   "rejected 3D3V step preserves every field bit");
    }
    free(snapshot);
    swa_em_pic3d_destroy(&state);

    memset(&state, 0, sizeof(state));
    limits = permissive_limits();
    check_true(
        swa_em_pic3d_init(
            &state, 17U, 1U, 1U, 0U,
            1.0, 1000.0, 1000.0,
            0.20 * (1.0 / 17.0) / SWA_C
        ),
        "transactional accepted-step grid initializes"
    );
    initialize_plane_wave(&state, 0, 2U, 1.0, 1);
    snapshot = calloc(6U * state.cell_count, sizeof(double));
    check_true(snapshot != NULL,
               "allocate accepted-step field snapshot");
    if (snapshot != NULL) {
        capture_field_state(&state, snapshot);
        check_true(
            swa_em_pic3d_step(&state, &limits, &result),
            "accepted transactional step executes"
        );
        check_true(result.passes && result.state_committed,
                   "passing 3D3V step commits exactly one candidate state");
        check_true(!field_state_matches_snapshot(&state, snapshot),
                   "accepted 3D3V step advances the field state");
        for (index = 0U; index < state.cell_count; ++index) {
            check_true(isfinite(state.ex_V_m[index]) &&
                       isfinite(state.ey_V_m[index]) &&
                       isfinite(state.ez_V_m[index]) &&
                       isfinite(state.bx_T[index]) &&
                       isfinite(state.by_T[index]) &&
                       isfinite(state.bz_T[index]),
                       "committed transactional state remains finite");
        }
    }
    free(snapshot);
    swa_em_pic3d_destroy(&state);
}

int main(void) {
    FILE *fp;
    test_invalid_contracts();
    test_higuera_cary_pusher();
    test_gauss_initialization();
    test_transport_and_projection();
    test_local_current_symmetries();
    test_discrete_maxwell_identities();
    test_plane_wave_convergence_and_symmetry();
    test_vacuum_time_reversal();
    test_coupled_cyclic_covariance();
    test_coupled_3d3v_reference();
    test_fail_closed_gates();
    test_transactional_step_acceptance();

    fp = fopen("output/em_pic3d_reference_receipt.json", "w");
    check_true(fp != NULL, "open 3D3V reference receipt");
    if (fp != NULL) {
        check_true(fprintf(
            fp,
            "{\n"
            "  \"schema\": \"spacewind.em-pic3d-assurance/v1\",\n"
            "  \"checks\": %zu,\n"
            "  \"failures\": %zu,\n"
            "  \"deterministic_hash\": \"%016llx\",\n"
            "  \"worst_pusher_reversal\": %.17g,\n"
            "  \"worst_pusher_rotation\": %.17g,\n"
            "  \"worst_projection_continuity\": %.17g,\n"
            "  \"worst_projection_curl\": %.17g,\n"
            "  \"worst_plane_wave_error\": %.17g,\n"
            "  \"worst_gauss_relative\": %.17g,\n"
            "  \"worst_div_b_relative\": %.17g,\n"
            "  \"worst_coupled_energy\": %.17g,\n"
            "  \"worst_coupled_momentum\": %.17g,\n"
            "  \"worst_local_current_difference\": %.17g,\n"
            "  \"worst_spectral_oracle_difference\": %.17g,\n"
            "  \"passes\": %s,\n"
            "  \"nonclaim\": \"the bounded periodic 3D3V reference uses a symmetric coordinate-split local charge-conserving current and an independent global spectral oracle; passing checks do not establish a final Esirkepov production deposition, open-boundary plasma-wing force, chamber validity, or propulsion\"\n"
            "}\n",
            checks,
            failures,
            (unsigned long long)transcript_hash,
            worst_pusher_reversal,
            worst_pusher_rotation,
            worst_projection_continuity,
            worst_projection_curl,
            worst_plane_wave_error,
            worst_gauss_relative,
            worst_div_b_relative,
            worst_coupled_energy,
            worst_coupled_momentum,
            worst_coupled_projection,
            worst_spectral_oracle_difference,
            failures == 0U ? "true" : "false"
        ) > 0, "write 3D3V reference receipt");
        (void)fclose(fp);
    }

    printf(
        "spacewind 3D3V EM PIC assurance: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
