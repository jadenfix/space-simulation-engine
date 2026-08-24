#include "spacewind/em_pic1d.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);

static void check_true(int condition, const char *name) {
    ++checks;
    transcript_hash ^= swa_fnv1a64(name, strlen(name));
    transcript_hash *= UINT64_C(1099511628211);
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}

static swa_em_pic1d_limits limits(double energy, double momentum) {
    swa_em_pic1d_limits result;
    result.continuity_relative_tolerance = 2.0e-11;
    result.gauss_relative_tolerance = 2.0e-11;
    result.magnetic_divergence_relative_tolerance = 1.0e-15;
    result.energy_relative_tolerance = energy;
    result.momentum_relative_tolerance = momentum;
    result.maximum_courant = 0.8;
    result.maximum_particle_cells_per_step = 0.8;
    return result;
}

static void initialize_right_wave(swa_em_pic1d *state,
                                  unsigned mode,
                                  double amplitude_V_m) {
    const double dx = state->domain_length_m / (double)state->cells;
    const double wave_number = 2.0 * SWA_PI * (double)mode /
                               state->domain_length_m;
    size_t i;
    for (i = 0U; i < state->cells; ++i) {
        const double x_face = (double)i * dx;
        const double x_cell = ((double)i + 0.5) * dx;
        state->ey_face_V_m[i] = amplitude_V_m * sin(wave_number * x_face);
        state->bz_cell_T[i] = amplitude_V_m / SWA_C *
                              sin(wave_number * x_cell);
    }
}

static double total_energy(const swa_em_pic1d *state) {
    const double dx = state->domain_length_m / (double)state->cells;
    const double volume = dx * state->cross_section_area_m2;
    double energy = 0.0;
    size_t i;
    for (i = 0U; i < state->cells; ++i) {
        const double e2 = state->ex_face_V_m[i] * state->ex_face_V_m[i] +
                          state->ey_face_V_m[i] * state->ey_face_V_m[i] +
                          state->ez_face_V_m[i] * state->ez_face_V_m[i];
        const double b2 = state->bx_uniform_T * state->bx_uniform_T +
                          state->by_cell_T[i] * state->by_cell_T[i] +
                          state->bz_cell_T[i] * state->bz_cell_T[i];
        energy += 0.5 * volume * (SWA_EPS0 * e2 + b2 / SWA_MU0);
    }
    for (i = 0U; i < state->particle_count; ++i) {
        const double speed2 = swa_vdot(state->particles[i].velocity_mps,
                                       state->particles[i].velocity_mps);
        const double gamma = 1.0 / sqrt(1.0 -
            speed2 / (SWA_C * SWA_C));
        energy += state->particles[i].macro_weight *
                  state->particles[i].mass_kg * SWA_C * SWA_C *
                  (gamma - 1.0);
    }
    return energy;
}

static double vacuum_one_step_error(double courant) {
    swa_em_pic1d state;
    swa_em_pic1d_result result;
    swa_em_pic1d_limits gate = limits(1.0, 1.0);
    const size_t cells = 64U;
    const double length = 1.0;
    const double dx = length / (double)cells;
    const double dt = courant * dx / SWA_C;
    if (!swa_em_pic1d_init(&state, cells, 0U, length, 1.0, dt)) {
        return INFINITY;
    }
    initialize_right_wave(&state, 3U, 2.0e4);
    if (!swa_em_pic1d_step(&state, &gate, &result)) {
        swa_em_pic1d_destroy(&state);
        return INFINITY;
    }
    swa_em_pic1d_destroy(&state);
    return result.total_energy_relative;
}

static void test_vacuum_wave(void) {
    swa_em_pic1d state;
    swa_em_pic1d_result result;
    swa_em_pic1d_limits gate = limits(2.0e-3, 2.0e-3);
    const size_t cells = 96U;
    const double length = 3.0;
    const double dx = length / (double)cells;
    const double dt = 0.35 * dx / SWA_C;
    const double initial_energy_expected_floor = 1.0e-10;
    double initial_energy;
    double final_energy;
    double relative_change;
    size_t step;
    FILE *fp;

    check_true(swa_em_pic1d_init(&state, cells, 0U,
                                 length, 2.0, dt),
               "initialize vacuum electromagnetic grid");
    initialize_right_wave(&state, 4U, 1.5e4);
    initial_energy = total_energy(&state);
    check_true(initial_energy > initial_energy_expected_floor,
               "vacuum wave carries field energy");
    for (step = 0U; step < 500U; ++step) {
        check_true(swa_em_pic1d_step(&state, &gate, &result),
                   "vacuum electromagnetic step executes");
        check_true(result.finite,
                   "vacuum electromagnetic receipt finite");
        check_true(result.continuity_passes &&
                   result.initial_gauss_passes &&
                   result.final_gauss_passes,
                   "vacuum charge and Gauss constraints close");
        check_true(result.magnetic_divergence_passes,
                   "one-dimensional magnetic divergence closes");
        check_true(result.energy_passes,
                   "vacuum one-step energy residual bounded");
    }
    final_energy = total_energy(&state);
    relative_change = fabs(final_energy - initial_energy) / initial_energy;
    check_true(relative_change < 3.0e-3,
               "vacuum wave long-run energy remains bounded");
    fp = fopen("output/em_pic1d_vacuum_receipt.json", "w");
    check_true(fp != NULL, "open electromagnetic PIC receipt");
    if (fp != NULL) {
        check_true(swa_em_pic1d_write_receipt(fp, &state, &gate, &result),
                   "write electromagnetic PIC receipt");
        (void)fclose(fp);
    }
    swa_em_pic1d_destroy(&state);
}

static void test_vacuum_timestep_refinement(void) {
    const double coarse = vacuum_one_step_error(0.6);
    const double medium = vacuum_one_step_error(0.3);
    const double fine = vacuum_one_step_error(0.15);
    check_true(isfinite(coarse) && isfinite(medium) && isfinite(fine),
               "vacuum refinement errors finite");
    check_true(medium < coarse,
               "halving timestep reduces vacuum energy defect");
    check_true(fine < medium,
               "second timestep halving reduces vacuum energy defect");
}

static void test_neutral_relativistic_plasma(void) {
    enum { PARTICLES = 8 };
    swa_em_pic1d state;
    swa_em_pic1d_result result;
    swa_em_pic1d_limits gate = limits(6.0e-2, 2.5e-1);
    const size_t cells = 64U;
    const double length = 2.0;
    const double dx = length / (double)cells;
    const double dt = 0.12 * dx / SWA_C;
    const double charge = 2.0e-18;
    const double mass = 5.0e-22;
    size_t i;
    size_t step;
    double worst_energy = 0.0;
    double worst_momentum = 0.0;

    check_true(swa_em_pic1d_init(&state, cells, PARTICLES,
                                 length, 1.5, dt),
               "initialize neutral 1D3V plasma");
    state.bx_uniform_T = 2.0e-4;
    for (i = 0U; i < PARTICLES; ++i) {
        const double sign = (i % 2U) == 0U ? 1.0 : -1.0;
        const double x = length * ((double)i + 0.35) /
                         (double)PARTICLES;
        const swa_vec3 velocity = swa_v3(
            sign * 0.04 * SWA_C,
            (0.01 + 0.001 * (double)i) * SWA_C,
            sign * 0.006 * SWA_C
        );
        check_true(swa_em_pic1d_set_particle(
                       &state, i, x, velocity,
                       sign * charge, mass, 1.0),
                   "configure neutral plasma particle");
    }
    check_true(swa_em_pic1d_initialize_gauss_field(
                   &state, 0.0, 1.0e-28),
               "initialize neutral plasma Gauss field");
    initialize_right_wave(&state, 2U, 50.0);
    for (step = 0U; step < 120U; ++step) {
        check_true(swa_em_pic1d_step(&state, &gate, &result),
                   "neutral electromagnetic PIC step executes");
        check_true(result.continuity_passes,
                   "neutral plasma continuity closes");
        check_true(result.initial_gauss_passes &&
                   result.final_gauss_passes,
                   "neutral plasma Gauss law propagates");
        check_true(result.subluminal_passes,
                   "relativistic Boris update remains subluminal");
        check_true(result.particle_courant_passes,
                   "particle transport remains spatially resolved");
        if (result.total_energy_relative > worst_energy) {
            worst_energy = result.total_energy_relative;
        }
        if (result.total_momentum_relative > worst_momentum) {
            worst_momentum = result.total_momentum_relative;
        }
    }
    check_true(worst_energy < gate.energy_relative_tolerance,
               "coupled plasma energy defect remains inside declared gate");
    check_true(worst_momentum < gate.momentum_relative_tolerance,
               "coupled plasma momentum defect remains inside declared gate");
    swa_em_pic1d_destroy(&state);
}


static void test_harmonic_current_and_node_crossing(void) {
    swa_em_pic1d state;
    swa_em_pic1d_result result;
    swa_em_pic1d_limits gate = limits(1.0, 1.0);
    const size_t cells = 32U;
    const double length = 1.0;
    const double dx = length / (double)cells;
    const double dt = 0.5 * dx / SWA_C;
    const double start = 0.90 * dx;
    const double speed = 0.60 * SWA_C;
    const double charge = 1.0e-18;
    const double mass = 1.0e-21;

    check_true(swa_em_pic1d_init(&state, cells, 2U,
                                 length, 1.0, dt),
               "initialize harmonic-current crossing case");
    check_true(swa_em_pic1d_set_particle(
                   &state, 0U, start,
                   swa_v3(speed, 0.02 * SWA_C, 0.01 * SWA_C),
                   charge, mass, 1.0),
               "set positive crossing particle");
    check_true(swa_em_pic1d_set_particle(
                   &state, 1U, start,
                   swa_v3(-speed, -0.02 * SWA_C, -0.01 * SWA_C),
                   -charge, mass, 1.0),
               "set negative crossing particle");
    check_true(swa_em_pic1d_initialize_gauss_field(
                   &state, 0.0, 1.0e-28),
               "initialize harmonic-current Gauss field");
    check_true(swa_em_pic1d_step(&state, &gate, &result),
               "harmonic-current electromagnetic step executes");
    check_true(result.maximum_path_segments >= 2U,
               "node-crossing path is split across grid faces");
    check_true(fabs(result.mean_jx_A_m2) > 0.0,
               "opposite charge and velocity pair carries harmonic current");
    check_true(result.harmonic_ampere_passes,
               "mean electric field follows harmonic Ampere update");
    check_true(result.transverse_current_partition_passes,
               "transverse path weights preserve deposited current sum");
    check_true(result.continuity_passes && result.final_gauss_passes,
               "node-crossing current preserves continuity and Gauss law");
    swa_em_pic1d_destroy(&state);
}

static void test_fail_closed_inputs(void) {
    swa_em_pic1d state;
    swa_em_pic1d_result result;
    swa_em_pic1d_limits gate = limits(1.0, 1.0);
    const size_t cells = 16U;
    const double length = 1.0;
    const double dx = length / (double)cells;

    check_true(!swa_em_pic1d_init(NULL, cells, 0U,
                                  length, 1.0, 1.0e-12),
               "null electromagnetic PIC state rejected");
    check_true(swa_em_pic1d_init(&state, cells, 1U,
                                 length, 1.0,
                                 1.1 * dx / SWA_C),
               "initialize deliberately over-CFL state");
    check_true(swa_em_pic1d_set_particle(
                   &state, 0U, 0.25,
                   swa_v3(0.0, 0.0, 0.0),
                   1.0e-19, 1.0e-20, 1.0),
               "configure over-CFL particle");
    check_true(swa_em_pic1d_step(&state, &gate, &result),
               "over-CFL audit returns receipt");
    check_true(!result.courant_passes && !result.passes,
               "over-CFL state fails closed");
    swa_em_pic1d_destroy(&state);

    check_true(swa_em_pic1d_init(&state, cells, 1U,
                                 length, 1.0,
                                 0.1 * dx / SWA_C),
               "initialize nonneutral state");
    check_true(swa_em_pic1d_set_particle(
                   &state, 0U, 0.2,
                   swa_v3(0.0, 0.0, 0.0),
                   1.0e-15, 1.0e-20, 1.0),
               "configure nonneutral particle");
    check_true(!swa_em_pic1d_initialize_gauss_field(
                   &state, 0.0, 1.0e-20),
               "periodic nonneutral Gauss initialization rejected");
    check_true(!swa_em_pic1d_set_particle(
                   &state, 0U, 0.2,
                   swa_v3(SWA_C, 0.0, 0.0),
                   1.0e-15, 1.0e-20, 1.0),
               "luminal particle input rejected");
    swa_em_pic1d_destroy(&state);
}

static void test_gauss_corruption_detected(void) {
    swa_em_pic1d state;
    swa_em_pic1d_result result;
    swa_em_pic1d_limits gate = limits(1.0, 1.0);
    const size_t cells = 32U;
    const double length = 1.0;
    const double dx = length / (double)cells;
    const double dt = 0.08 * dx / SWA_C;

    check_true(swa_em_pic1d_init(&state, cells, 2U,
                                 length, 1.0, dt),
               "initialize Gauss corruption case");
    check_true(swa_em_pic1d_set_particle(
                   &state, 0U, 0.3,
                   swa_v3(1.0e4, 0.0, 0.0),
                   1.0e-16, 1.0e-19, 1.0),
               "set positive Gauss particle");
    check_true(swa_em_pic1d_set_particle(
                   &state, 1U, 0.7,
                   swa_v3(-1.0e4, 0.0, 0.0),
                   -1.0e-16, 1.0e-19, 1.0),
               "set negative Gauss particle");
    check_true(swa_em_pic1d_initialize_gauss_field(
                   &state, 0.0, 1.0e-25),
               "initialize valid Gauss field before corruption");
    state.ex_face_V_m[7U] += 100.0;
    check_true(swa_em_pic1d_step(&state, &gate, &result),
               "corrupted Gauss step returns receipt");
    check_true(!result.initial_gauss_passes,
               "initial Gauss corruption detected");
    check_true(!result.final_gauss_passes,
               "Gauss corruption propagates as failure");
    check_true(!result.passes,
               "corrupted Gauss state cannot pass");
    swa_em_pic1d_destroy(&state);
}

int main(void) {
    test_vacuum_wave();
    test_vacuum_timestep_refinement();
    test_neutral_relativistic_plasma();
    test_harmonic_current_and_node_crossing();
    test_fail_closed_inputs();
    test_gauss_corruption_detected();
    printf("spacewind 1D3V EM PIC assurance: %zu checks, %zu failures, hash=%016llx\n",
           checks, failures, (unsigned long long)transcript_hash);
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
