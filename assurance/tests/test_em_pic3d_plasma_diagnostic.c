#include "spacewind/em_pic3d.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NX 8U
#define NY 1U
#define NZ 1U
#define PARTICLES (2U * NX * NY * NZ)

static swa_em_pic3d_limits limits(void) {
    swa_em_pic3d_limits value;
    memset(&value, 0, sizeof(value));
    value.continuity_relative_tolerance = 2.0e-10;
    value.gauss_relative_tolerance = 2.0e-10;
    value.magnetic_divergence_relative_tolerance = 2.0e-11;
    value.local_current_difference_relative_tolerance = 0.90;
    value.spectral_oracle_relative_tolerance = 5.0e-10;
    value.spectral_oracle_curl_relative_tolerance = 5.0e-10;
    value.energy_relative_tolerance = 5.0e-8;
    value.particle_work_relative_tolerance = 2.0e-4;
    value.field_work_relative_tolerance = 5.0e-8;
    value.momentum_relative_tolerance = 5.0e-6;
    value.nonlinear_relative_tolerance = 2.0e-10;
    value.maximum_courant = 0.25;
    value.maximum_particle_cells_per_step = 0.20;
    value.maximum_nonlinear_iterations = 24U;
    return value;
}

static int configure(
    swa_em_pic3d *state,
    double number_density_m3,
    double electron_mass_kg
) {
    const double dx = state->length_x_m / (double)state->nx;
    const double dy = state->length_y_m / (double)state->ny;
    const double dz = state->length_z_m / (double)state->nz;
    const double macro_weight = number_density_m3 * dx * dy * dz;
    size_t cell;
    if (!isfinite(macro_weight) || !(macro_weight > 0.0)) {
        return 0;
    }
    for (cell = 0U; cell < state->cell_count; ++cell) {
        const size_t i = cell % state->nx;
        const size_t j = (cell / state->nx) % state->ny;
        const size_t k = cell / (state->nx * state->ny);
        const swa_vec3 position = swa_v3(
            ((double)i + 0.5) * dx,
            ((double)j + 0.5) * dy,
            ((double)k + 0.5) * dz
        );
        if (!swa_em_pic3d_set_particle(
                state, 2U * cell, position,
                swa_v3(0.0, 0.0, 0.0),
                -SWA_QE, electron_mass_kg, macro_weight) ||
            !swa_em_pic3d_set_particle(
                state, 2U * cell + 1U, position,
                swa_v3(0.0, 0.0, 0.0),
                SWA_QE, SWA_MP, macro_weight)) {
            return 0;
        }
    }
    return 1;
}

static void print_result(const swa_em_pic3d_result *r, int executed) {
    fprintf(stderr,
        "DIAG exec=%d finite=%d passes=%d committed=%d "
        "flags[cfl=%d pcfl=%d q=%d cont=%d g0=%d g1=%d d0=%d d1=%d "
        "proj=%d local=%d spectral=%d mean=%d nonlinear=%d sublum=%d "
        "energy=%d pwork=%d fwork=%d mom=%d]\n",
        executed, r->finite, r->passes, r->state_committed,
        r->courant_passes, r->particle_courant_passes,
        r->charge_closure_passes, r->continuity_passes,
        r->initial_gauss_passes, r->final_gauss_passes,
        r->initial_magnetic_divergence_passes,
        r->final_magnetic_divergence_passes,
        r->current_projection_passes, r->local_current_passes,
        r->spectral_oracle_passes, r->mean_current_passes,
        r->nonlinear_passes, r->subluminal_passes,
        r->energy_passes, r->particle_work_passes,
        r->field_work_passes, r->momentum_passes);
    fprintf(stderr,
        "DIAG metrics[cfl=%.17g disp=%.17g iterations=%zu nonlinear=%.17g "
        "cont=%.17g g0=%.17g g1=%.17g d0=%.17g d1=%.17g "
        "local=%.17g spectral=%.17g mean=(%.17g,%.17g,%.17g) "
        "energy=%.17g pwork=%.17g fwork=%.17g mom=%.17g speed=%.17g]\n",
        r->courant_number, r->maximum_particle_displacement_cells,
        r->nonlinear_iterations, r->nonlinear_relative_residual,
        r->continuity_relative_max,
        r->initial_gauss_relative_max, r->final_gauss_relative_max,
        r->initial_magnetic_divergence_relative,
        r->final_magnetic_divergence_relative,
        r->local_current_difference_relative,
        r->spectral_oracle_difference_relative,
        r->mean_current_residual_A_m2.x,
        r->mean_current_residual_A_m2.y,
        r->mean_current_residual_A_m2.z,
        r->total_energy_relative, r->particle_work_relative,
        r->field_work_relative, r->total_momentum_relative,
        r->maximum_speed_fraction_c);
    fprintf(stderr,
        "DIAG charge[initial=%.17g final=%.17g change=%.17g] "
        "continuity[max=%.17g rms=%.17g zero=%.17g] "
        "current[raw_l2=%.17g correction_l2=%.17g correction_rel=%.17g "
        "curl_rel=%.17g mean=(%.17g,%.17g,%.17g) "
        "transport=(%.17g,%.17g,%.17g)]\n",
        r->initial_charge_C, r->final_charge_C, r->charge_change_C,
        r->continuity_residual_max_A_m3,
        r->continuity_residual_rms_A_m3,
        r->current_projection_zero_mode_A_m3,
        r->raw_current_l2_A_m2,
        r->current_correction_l2_A_m2,
        r->current_correction_relative,
        r->current_projection_curl_relative,
        r->mean_current_A_m2.x,
        r->mean_current_A_m2.y,
        r->mean_current_A_m2.z,
        r->transport_mean_current_A_m2.x,
        r->transport_mean_current_A_m2.y,
        r->transport_mean_current_A_m2.z);
    fprintf(stderr,
        "DIAG work[particle_dE=%.17g field_dE=%.17g particle_JE=%.17g "
        "field_JE=%.17g projection_work=%.17g particle_res=%.17g "
        "field_res=%.17g total_res=%.17g]\n",
        r->particle_energy_change_J, r->field_energy_change_J,
        r->particle_midpoint_current_work_J,
        r->field_midpoint_current_work_J,
        r->projection_work_correction_J,
        r->particle_work_residual_J, r->field_work_residual_J,
        r->total_energy_residual_J);
    fprintf(stderr,
        "DIAG momentum[particle0=(%.17g,%.17g,%.17g) "
        "particle1=(%.17g,%.17g,%.17g) field0=(%.17g,%.17g,%.17g) "
        "field1=(%.17g,%.17g,%.17g) residual=(%.17g,%.17g,%.17g)]\n",
        r->initial_particle_momentum_Ns.x,
        r->initial_particle_momentum_Ns.y,
        r->initial_particle_momentum_Ns.z,
        r->final_particle_momentum_Ns.x,
        r->final_particle_momentum_Ns.y,
        r->final_particle_momentum_Ns.z,
        r->initial_field_momentum_Ns.x,
        r->initial_field_momentum_Ns.y,
        r->initial_field_momentum_Ns.z,
        r->final_field_momentum_Ns.x,
        r->final_field_momentum_Ns.y,
        r->final_field_momentum_Ns.z,
        r->total_momentum_residual_Ns.x,
        r->total_momentum_residual_Ns.y,
        r->total_momentum_residual_Ns.z);
}

int main(void) {
    swa_em_pic3d state;
    swa_em_pic3d_result result;
    swa_em_pic3d_limits gate = limits();
    const double lx = 1.0;
    const double ly = 1.0;
    const double lz = 1.0;
    const double dx = lx / (double)NX;
    const double dy = ly / (double)NY;
    const double dz = lz / (double)NZ;
    const double inverse_spacing2 =
        1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz);
    const double dt = 0.12 / (SWA_C * sqrt(inverse_spacing2));
    const double target_omega = 2.0 * SWA_PI / (1024.0 * dt);
    const double density = target_omega * target_omega * SWA_EPS0 /
        (SWA_QE * SWA_QE * (1.0 / SWA_ME + 1.0 / SWA_MP));
    const double absolute_charge =
        2.0 * (double)(NX * NY * NZ) * density * dx * dy * dz * SWA_QE;
    const double charge_tolerance =
        fmax(1.0e-13 * absolute_charge, 1.0e-30);
    int initialized;
    int configured;
    int gauss;
    int executed;

    memset(&state, 0, sizeof(state));
    memset(&result, 0, sizeof(result));
    initialized = swa_em_pic3d_init(
        &state, NX, NY, NZ, PARTICLES,
        lx, ly, lz, dt
    );
    fprintf(stderr,
            "DIAG setup dt=%.17g omega=%.17g density=%.17g init=%d\n",
            dt, target_omega, density, initialized);
    if (!initialized) {
        return EXIT_FAILURE;
    }
    configured = configure(&state, density, SWA_ME);
    fprintf(stderr, "DIAG configured=%d state_finite=%d\n",
            configured, swa_em_pic3d_state_is_finite(&state));
    if (!configured) {
        swa_em_pic3d_destroy(&state);
        return EXIT_FAILURE;
    }
    gauss = swa_em_pic3d_initialize_gauss_field(
        &state, swa_v3(1000.0, 0.0, 0.0), charge_tolerance
    );
    fprintf(stderr,
            "DIAG gauss=%d state_finite=%d abs_charge=%.17g tol=%.17g\n",
            gauss, swa_em_pic3d_state_is_finite(&state),
            absolute_charge, charge_tolerance);
    if (!gauss) {
        swa_em_pic3d_destroy(&state);
        return EXIT_FAILURE;
    }
    executed = swa_em_pic3d_step(&state, &gate, &result);
    print_result(&result, executed);
    swa_em_pic3d_destroy(&state);
    return executed && result.passes && result.state_committed ?
        EXIT_SUCCESS : EXIT_FAILURE;
}
