static double particle_energy_J(
    const swa_em_pic1d_particle *particles,
    size_t count
) {
    swa_kahan energy = {0.0, 0.0};
    size_t i;
    for (i = 0U; i < count; ++i) {
        const double speed2 = swa_vdot(particles[i].velocity_mps,
                                       particles[i].velocity_mps);
        const double gamma = 1.0 / sqrt(1.0 -
            speed2 / (SWA_C * SWA_C));
        swa_kahan_add(&energy,
            particles[i].macro_weight * particles[i].mass_kg *
            SWA_C * SWA_C * (gamma - 1.0));
    }
    return energy.sum;
}

static swa_vec3 particle_momentum_Ns(
    const swa_em_pic1d_particle *particles,
    size_t count
) {
    swa_kahan px = {0.0, 0.0};
    swa_kahan py = {0.0, 0.0};
    swa_kahan pz = {0.0, 0.0};
    size_t i;
    for (i = 0U; i < count; ++i) {
        const double speed2 = swa_vdot(particles[i].velocity_mps,
                                       particles[i].velocity_mps);
        const double gamma = 1.0 / sqrt(1.0 -
            speed2 / (SWA_C * SWA_C));
        const double scale = particles[i].macro_weight *
                             particles[i].mass_kg * gamma;
        swa_kahan_add(&px, scale * particles[i].velocity_mps.x);
        swa_kahan_add(&py, scale * particles[i].velocity_mps.y);
        swa_kahan_add(&pz, scale * particles[i].velocity_mps.z);
    }
    return swa_v3(px.sum, py.sum, pz.sum);
}

static double field_energy_J(
    const swa_em_pic1d *state,
    const double *ex_face,
    const double *ey_face,
    const double *ez_face,
    const double *by_cell,
    const double *bz_cell
) {
    const double volume = state->cross_section_area_m2 *
        state->domain_length_m / (double)state->cells;
    swa_kahan energy = {0.0, 0.0};
    size_t i;
    for (i = 0U; i < state->cells; ++i) {
        const double e2 = ex_face[i] * ex_face[i] +
                          ey_face[i] * ey_face[i] +
                          ez_face[i] * ez_face[i];
        const double b2 = state->bx_uniform_T * state->bx_uniform_T +
                          by_cell[i] * by_cell[i] +
                          bz_cell[i] * bz_cell[i];
        swa_kahan_add(&energy,
            0.5 * volume * (SWA_EPS0 * e2 + b2 / SWA_MU0));
    }
    return energy.sum;
}

static swa_vec3 field_momentum_Ns(
    const swa_em_pic1d *state,
    const double *ex_face,
    const double *ey_face,
    const double *ez_face,
    const double *by_cell,
    const double *bz_cell
) {
    const double volume = state->cross_section_area_m2 *
        state->domain_length_m / (double)state->cells;
    swa_kahan px = {0.0, 0.0};
    swa_kahan py = {0.0, 0.0};
    swa_kahan pz = {0.0, 0.0};
    size_t i;
    for (i = 0U; i < state->cells; ++i) {
        const size_t previous = (i + state->cells - 1U) % state->cells;
        const swa_vec3 electric = swa_v3(ex_face[i], ey_face[i],
                                          ez_face[i]);
        const swa_vec3 magnetic = swa_v3(
            state->bx_uniform_T,
            0.5 * (by_cell[i] + by_cell[previous]),
            0.5 * (bz_cell[i] + bz_cell[previous])
        );
        const swa_vec3 density = swa_vscale(
            swa_vcross(electric, magnetic), SWA_EPS0
        );
        swa_kahan_add(&px, density.x * volume);
        swa_kahan_add(&py, density.y * volume);
        swa_kahan_add(&pz, density.z * volume);
    }
    return swa_v3(px.sum, py.sum, pz.sum);
}

static void gauss_metrics(
    const swa_em_pic1d *state,
    const double *ex_face,
    const double *rho_C_m3,
    double *maximum_abs,
    double *maximum_relative
) {
    const double dx_m = state->domain_length_m / (double)state->cells;
    double global_scale = DBL_MIN;
    size_t i;
    *maximum_abs = 0.0;
    *maximum_relative = 0.0;
    for (i = 0U; i < state->cells; ++i) {
        const size_t next = (i + 1U) % state->cells;
        const double divergence_charge = SWA_EPS0 *
            (ex_face[next] - ex_face[i]) / dx_m;
        global_scale = max3(global_scale, fabs(divergence_charge),
                            fabs(rho_C_m3[i]));
    }
    for (i = 0U; i < state->cells; ++i) {
        const size_t next = (i + 1U) % state->cells;
        const double divergence_charge = SWA_EPS0 *
            (ex_face[next] - ex_face[i]) / dx_m;
        const double residual = divergence_charge - rho_C_m3[i];
        *maximum_abs = max2(*maximum_abs, fabs(residual));
    }
    *maximum_relative = *maximum_abs / global_scale;
}

static void continuity_metrics(
    const swa_em_pic1d *state,
    const double *rho_initial_C_m3,
    const double *rho_final_C_m3,
    const double *jx_face_A_m2,
    double *maximum_abs,
    double *rms,
    double *maximum_relative
) {
    const double dx_m = state->domain_length_m / (double)state->cells;
    double global_scale = DBL_MIN;
    swa_kahan squared = {0.0, 0.0};
    size_t i;
    *maximum_abs = 0.0;
    for (i = 0U; i < state->cells; ++i) {
        const size_t next = (i + 1U) % state->cells;
        const double time_derivative =
            (rho_final_C_m3[i] - rho_initial_C_m3[i]) / state->dt_s;
        const double divergence =
            (jx_face_A_m2[next] - jx_face_A_m2[i]) / dx_m;
        global_scale = max3(global_scale, fabs(time_derivative),
                            fabs(divergence));
        *maximum_abs = max2(*maximum_abs,
                            fabs(time_derivative + divergence));
        swa_kahan_add(&squared,
                      (time_derivative + divergence) *
                      (time_derivative + divergence));
    }
    *maximum_relative = *maximum_abs / global_scale;
    *rms = sqrt(squared.sum / (double)state->cells);
}

