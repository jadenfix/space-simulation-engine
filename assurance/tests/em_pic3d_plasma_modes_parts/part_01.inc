    const double beta2 = speed2_m2_s2 / (SWA_C * SWA_C);
    const double inverse_gamma = sqrt(1.0 - beta2);
    return beta2 / (inverse_gamma * (1.0 + inverse_gamma));
}

static double particle_energy_J(const swa_em_pic3d *state) {
    swa_kahan energy = {0.0, 0.0};
    size_t particle;
    for (particle = 0U; particle < state->particle_count; ++particle) {
        const swa_em_pic3d_particle *item = &state->particles[particle];
        const double speed2 = swa_vdot(
            item->velocity_mps, item->velocity_mps
        );
        const double kinetic_factor = stable_gamma_minus_one(speed2);
        swa_kahan_add(
            &energy,
            item->macro_weight * item->mass_kg *
            SWA_C * SWA_C * kinetic_factor
        );
    }
    return energy.sum;
}

static double total_energy_J(const swa_em_pic3d *state) {
    return particle_energy_J(state) + swa_em_pic3d_field_energy_J(state);
}

static double mean_ex_V_m(const swa_em_pic3d *state) {
    swa_kahan sum = {0.0, 0.0};
    size_t index;
    for (index = 0U; index < state->cell_count; ++index) {
        swa_kahan_add(&sum, state->ex_V_m[index]);
    }
    return sum.sum / (double)state->cell_count;
}

static int configure_neutral_pair(
    swa_em_pic3d *state,
    double density_m3,
    double negative_mass_kg,
    double positive_mass_kg
) {
    const double dx = state->length_x_m / (double)state->nx;
    const double dy = state->length_y_m / (double)state->ny;
    const double dz = state->length_z_m / (double)state->nz;
    const double macro_weight = density_m3 * dx * dy * dz;
    size_t cell;
    if (!isfinite(macro_weight) || !(macro_weight > 0.0) ||
        !isfinite(negative_mass_kg) || !(negative_mass_kg > 0.0) ||
        !isfinite(positive_mass_kg) || !(positive_mass_kg > 0.0)) {
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
                -SWA_QE, negative_mass_kg, macro_weight
            ) ||
            !swa_em_pic3d_set_particle(
                state, 2U * cell + 1U, position,
                swa_v3(0.0, 0.0, 0.0),
                SWA_QE, positive_mass_kg, macro_weight
            )) {
            return 0;
        }
    }
    return 1;
}

static int same_sign_charge_control_rejected(
    double dt_s,
    double density_m3
) {
    swa_em_pic3d state;
    const double length_x_m = 1.0;
    const double length_y_m = 1.0;
    const double length_z_m = 1.0;
    const double dx = length_x_m / (double)MODE_NX;
    const double dy = length_y_m / (double)MODE_NY;
    const double dz = length_z_m / (double)MODE_NZ;
    const double macro_weight = density_m3 * dx * dy * dz;
    const double absolute_charge =
        2.0 * (double)(MODE_NX * MODE_NY * MODE_NZ) *
        macro_weight * SWA_QE;
    const double charge_tolerance =
        max2(1.0e-13 * absolute_charge, 1.0e-30);
    size_t cell;
    int configured = 1;
    int gauss_initialized;

    memset(&state, 0, sizeof(state));
    if (!swa_em_pic3d_init(
            &state,
            MODE_NX, MODE_NY, MODE_NZ, MODE_PARTICLES,
            length_x_m, length_y_m, length_z_m, dt_s
        )) {
        return 0;
    }
    for (cell = 0U; cell < state.cell_count; ++cell) {
        const size_t i = cell % state.nx;
        const size_t j = (cell / state.nx) % state.ny;
        const size_t k = cell / (state.nx * state.ny);
