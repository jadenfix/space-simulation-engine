int swa_em_pic1d_step(
    swa_em_pic1d *state,
    const swa_em_pic1d_limits *limits,
    swa_em_pic1d_result *out
) {
    const size_t cells = state != NULL ? state->cells : 0U;
    const size_t particle_count = state != NULL ? state->particle_count : 0U;
    double dx_m;
    double *position_initial = NULL;
    double *position_final = NULL;
    double *rho_initial = NULL;
    double *rho_final = NULL;
    double *jx = NULL;
    double *jy = NULL;
    double *jz = NULL;
    double *ex_initial = NULL;
    double *ey_initial = NULL;
    double *ez_initial = NULL;
    double *by_initial = NULL;
    double *bz_initial = NULL;
    double *by_half = NULL;
    double *bz_half = NULL;
    double *ex_final = NULL;
    double *ey_final = NULL;
    double *ez_final = NULL;
    double *by_final = NULL;
    double *bz_final = NULL;
    swa_em_pic1d_particle *particle_initial = NULL;
    swa_em_pic1d_particle *particle_final = NULL;
    double total_charge_initial;
    double total_charge_final;
    double absolute_charge_initial;
    double absolute_charge_final;
    double target_mean_current;
    double charge_scale;
    double energy_scale;
    double momentum_scale;
    double current_work = 0.0;
    swa_kahan target_jy_sum = {0.0, 0.0};
    swa_kahan target_jz_sum = {0.0, 0.0};
    swa_kahan deposited_jy_sum = {0.0, 0.0};
    swa_kahan deposited_jz_sum = {0.0, 0.0};
    swa_kahan initial_ex_sum = {0.0, 0.0};
    swa_kahan final_ex_sum = {0.0, 0.0};
    size_t i;
    int ok = 0;

    if (!state_values_valid(state) || !limits_valid(limits) || out == NULL ||
        cells > SIZE_MAX / sizeof(double) ||
        particle_count > SIZE_MAX / sizeof(double) ||
        particle_count > SIZE_MAX / sizeof(swa_em_pic1d_particle)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    dx_m = state->domain_length_m / (double)cells;
    out->courant_number = SWA_C * state->dt_s / dx_m;
    out->courant_passes = out->courant_number <= limits->maximum_courant &&
                          out->courant_number <= 1.0;
    if (!out->courant_passes) {
        out->finite = isfinite(out->courant_number);
        return 1;
    }

#define ALLOC_DOUBLE(name, count) do { \
    name = calloc(count_or_one((count)), sizeof(double)); \
    if ((name) == NULL) goto cleanup; \
} while (0)
    ALLOC_DOUBLE(position_initial, particle_count);
    ALLOC_DOUBLE(position_final, particle_count);
    ALLOC_DOUBLE(rho_initial, cells);
    ALLOC_DOUBLE(rho_final, cells);
    ALLOC_DOUBLE(jx, cells);
    ALLOC_DOUBLE(jy, cells);
    ALLOC_DOUBLE(jz, cells);
    ALLOC_DOUBLE(ex_initial, cells);
    ALLOC_DOUBLE(ey_initial, cells);
    ALLOC_DOUBLE(ez_initial, cells);
    ALLOC_DOUBLE(by_initial, cells);
    ALLOC_DOUBLE(bz_initial, cells);
    ALLOC_DOUBLE(by_half, cells);
    ALLOC_DOUBLE(bz_half, cells);
    ALLOC_DOUBLE(ex_final, cells);
    ALLOC_DOUBLE(ey_final, cells);
    ALLOC_DOUBLE(ez_final, cells);
    ALLOC_DOUBLE(by_final, cells);
    ALLOC_DOUBLE(bz_final, cells);
#undef ALLOC_DOUBLE
    particle_initial = calloc(count_or_one(particle_count), sizeof(*particle_initial));
    particle_final = calloc(count_or_one(particle_count), sizeof(*particle_final));
    if (particle_initial == NULL || particle_final == NULL) {
        goto cleanup;
    }

    memcpy(ex_initial, state->ex_face_V_m, cells * sizeof(double));
    memcpy(ey_initial, state->ey_face_V_m, cells * sizeof(double));
    memcpy(ez_initial, state->ez_face_V_m, cells * sizeof(double));
    memcpy(by_initial, state->by_cell_T, cells * sizeof(double));
    memcpy(bz_initial, state->bz_cell_T, cells * sizeof(double));
    memcpy(particle_initial, state->particles,
           particle_count * sizeof(*particle_initial));
    memcpy(particle_final, state->particles,
           particle_count * sizeof(*particle_final));
    for (i = 0U; i < particle_count; ++i) {
        position_initial[i] = particle_initial[i].position_unwrapped_m;
    }
    if (!deposit_charge(state, position_initial, rho_initial,
                        &total_charge_initial, &absolute_charge_initial)) {
        goto cleanup;
    }

    for (i = 0U; i < cells; ++i) {
        const size_t next = (i + 1U) % cells;
        by_half[i] = by_initial[i] + 0.5 * state->dt_s *
            (ez_initial[next] - ez_initial[i]) / dx_m;
        bz_half[i] = bz_initial[i] - 0.5 * state->dt_s *
            (ey_initial[next] - ey_initial[i]) / dx_m;
    }

    out->maximum_speed_fraction_c = 0.0;
    out->maximum_particle_displacement_cells = 0.0;
    out->maximum_path_segments = 0U;
    for (i = 0U; i < particle_count; ++i) {
        swa_vec3 magnetic;
        const swa_vec3 electric = interpolate_fields(
            state, ex_initial, ey_initial, ez_initial,
            by_half, bz_half, position_initial[i], &magnetic
        );
        swa_vec3 velocity_final;
        const swa_vec3 velocity_initial = particle_initial[i].velocity_mps;
        const double charge_to_mass = particle_initial[i].charge_C /
                                      particle_initial[i].mass_kg;
        double displacement;
        double speed_fraction;
        size_t segments;
        if (!relativistic_boris(velocity_initial, electric, magnetic,
                                charge_to_mass, state->dt_s,
                                &velocity_final)) {
            goto cleanup;
        }
        particle_final[i].velocity_mps = velocity_final;
        displacement = 0.5 * state->dt_s *
            (velocity_initial.x + velocity_final.x);
        particle_final[i].position_unwrapped_m =
            position_initial[i] + displacement;
        position_final[i] = particle_final[i].position_unwrapped_m;
        speed_fraction = swa_vnorm(velocity_final) / SWA_C;
        out->maximum_speed_fraction_c = max2(
            out->maximum_speed_fraction_c, speed_fraction
        );
        out->maximum_particle_displacement_cells = max2(
            out->maximum_particle_displacement_cells,
            fabs(displacement) / dx_m
        );
        {
            const double jy_scale = particle_initial[i].charge_C *
                particle_initial[i].macro_weight *
                0.5 * (velocity_initial.y + velocity_final.y) /
                state->cross_section_area_m2;
            const double jz_scale = particle_initial[i].charge_C *
                particle_initial[i].macro_weight *
                0.5 * (velocity_initial.z + velocity_final.z) /
                state->cross_section_area_m2;
            swa_kahan_add(&target_jy_sum, jy_scale);
            swa_kahan_add(&target_jz_sum, jz_scale);
            if (!add_transverse_path_current(
                    state,
                    position_initial[i], position_final[i],
                    jy_scale, jz_scale,
                    jy, jz, &segments)) {
                goto cleanup;
            }
        }
        if (segments > out->maximum_path_segments) {
            out->maximum_path_segments = segments;
        }
    }
    out->particle_courant_passes =
        out->maximum_particle_displacement_cells <=
        limits->maximum_particle_cells_per_step;
    out->subluminal_passes = out->maximum_speed_fraction_c < 1.0;

    if (!deposit_charge(state, position_final, rho_final,
                        &total_charge_final, &absolute_charge_final) ||
        !deposit_longitudinal_current(
            state, position_initial, position_final,
            rho_initial, rho_final, jx, &target_mean_current)) {
        goto cleanup;
    }
    out->mean_jx_A_m2 = target_mean_current;
    for (i = 0U; i < cells; ++i) {
        swa_kahan_add(&deposited_jy_sum, jy[i]);
        swa_kahan_add(&deposited_jz_sum, jz[i]);
    }
    out->mean_jy_A_m2 = deposited_jy_sum.sum / (double)cells;
    out->mean_jz_A_m2 = deposited_jz_sum.sum / (double)cells;
    out->transverse_current_partition_relative_max = max2(
        fabs(deposited_jy_sum.sum - target_jy_sum.sum) /
            max3(fabs(deposited_jy_sum.sum), fabs(target_jy_sum.sum), DBL_MIN),
        fabs(deposited_jz_sum.sum - target_jz_sum.sum) /
            max3(fabs(deposited_jz_sum.sum), fabs(target_jz_sum.sum), DBL_MIN)
    );
    out->transverse_current_partition_passes =
        out->transverse_current_partition_relative_max <=
        limits->continuity_relative_tolerance;

    for (i = 0U; i < cells; ++i) {
        const size_t previous = (i + cells - 1U) % cells;
        ex_final[i] = ex_initial[i] -
                      state->dt_s * jx[i] / SWA_EPS0;
        ey_final[i] = ey_initial[i] -
            SWA_C * SWA_C * state->dt_s *
                (bz_half[i] - bz_half[previous]) / dx_m -
            state->dt_s * jy[i] / SWA_EPS0;
        ez_final[i] = ez_initial[i] +
            SWA_C * SWA_C * state->dt_s *
                (by_half[i] - by_half[previous]) / dx_m -
            state->dt_s * jz[i] / SWA_EPS0;
    }
    for (i = 0U; i < cells; ++i) {
        const size_t next = (i + 1U) % cells;
        by_final[i] = by_half[i] + 0.5 * state->dt_s *
            (ez_final[next] - ez_final[i]) / dx_m;
        bz_final[i] = bz_half[i] - 0.5 * state->dt_s *
            (ey_final[next] - ey_final[i]) / dx_m;
    }

    for (i = 0U; i < cells; ++i) {
        swa_kahan_add(&initial_ex_sum, ex_initial[i]);
        swa_kahan_add(&final_ex_sum, ex_final[i]);
    }
    out->initial_mean_ex_V_m = initial_ex_sum.sum / (double)cells;
    out->final_mean_ex_V_m = final_ex_sum.sum / (double)cells;
    out->harmonic_ampere_residual_V_m =
        out->final_mean_ex_V_m - out->initial_mean_ex_V_m +
        state->dt_s * out->mean_jx_A_m2 / SWA_EPS0;
    out->harmonic_ampere_relative =
        fabs(out->harmonic_ampere_residual_V_m) /
        max3(fabs(out->initial_mean_ex_V_m),
             fabs(out->final_mean_ex_V_m),
             max2(fabs(state->dt_s * out->mean_jx_A_m2 / SWA_EPS0),
                  DBL_MIN));
    out->harmonic_ampere_passes = out->harmonic_ampere_relative <=
                                  limits->continuity_relative_tolerance;

    out->initial_charge_C = total_charge_initial;
    out->final_charge_C = total_charge_final;
    out->charge_change_C = total_charge_final - total_charge_initial;
    charge_scale = max3(absolute_charge_initial, absolute_charge_final,
                        DBL_MIN);
    out->charge_closure_passes =
        fabs(out->charge_change_C) / charge_scale <=
        limits->continuity_relative_tolerance;
    continuity_metrics(state, rho_initial, rho_final, jx,
                       &out->continuity_residual_max_A_m3,
                       &out->continuity_residual_rms_A_m3,
                       &out->continuity_relative_max);
    out->continuity_passes = out->continuity_relative_max <=
                             limits->continuity_relative_tolerance;
    gauss_metrics(state, ex_initial, rho_initial,
                  &out->initial_gauss_residual_max_C_m3,
                  &out->initial_gauss_relative_max);
    gauss_metrics(state, ex_final, rho_final,
                  &out->final_gauss_residual_max_C_m3,
                  &out->final_gauss_relative_max);
    out->initial_gauss_passes = out->initial_gauss_relative_max <=
                                limits->gauss_relative_tolerance;
    out->final_gauss_passes = out->final_gauss_relative_max <=
                              limits->gauss_relative_tolerance;
    out->magnetic_divergence_residual_T_m = 0.0;
    out->magnetic_divergence_passes =
        out->magnetic_divergence_residual_T_m <=
        limits->magnetic_divergence_relative_tolerance;

    out->initial_particle_energy_J = particle_energy_J(
        particle_initial, particle_count
    );
    out->final_particle_energy_J = particle_energy_J(
        particle_final, particle_count
    );
    out->initial_field_energy_J = field_energy_J(
        state, ex_initial, ey_initial, ez_initial, by_initial, bz_initial
    );
    out->final_field_energy_J = field_energy_J(
        state, ex_final, ey_final, ez_final, by_final, bz_final
    );
    out->particle_energy_change_J = out->final_particle_energy_J -
                                    out->initial_particle_energy_J;
    out->field_energy_change_J = out->final_field_energy_J -
                                 out->initial_field_energy_J;
    for (i = 0U; i < cells; ++i) {
        current_work += state->dt_s * dx_m *
            state->cross_section_area_m2 *
            (jx[i] * 0.5 * (ex_initial[i] + ex_final[i]) +
             jy[i] * 0.5 * (ey_initial[i] + ey_final[i]) +
             jz[i] * 0.5 * (ez_initial[i] + ez_final[i]));
    }
    out->midpoint_current_work_J = current_work;
    out->particle_work_residual_J = out->particle_energy_change_J -
                                    current_work;
    out->field_work_residual_J = out->field_energy_change_J +
                                 current_work;
    out->total_energy_residual_J = out->particle_energy_change_J +
                                   out->field_energy_change_J;
    energy_scale = max3(
        fabs(out->initial_particle_energy_J + out->initial_field_energy_J),
        fabs(out->final_particle_energy_J + out->final_field_energy_J),
        max2(fabs(current_work), DBL_MIN)
    );
    out->total_energy_relative = fabs(out->total_energy_residual_J) /
                                 energy_scale;
    out->energy_passes = out->total_energy_relative <=
                         limits->energy_relative_tolerance;

    out->initial_particle_momentum_Ns = particle_momentum_Ns(
        particle_initial, particle_count
    );
    out->final_particle_momentum_Ns = particle_momentum_Ns(
        particle_final, particle_count
    );
    out->initial_field_momentum_Ns = field_momentum_Ns(
        state, ex_initial, ey_initial, ez_initial, by_initial, bz_initial
    );
    out->final_field_momentum_Ns = field_momentum_Ns(
        state, ex_final, ey_final, ez_final, by_final, bz_final
    );
    out->total_momentum_residual_Ns = swa_vsub(
        swa_vadd(out->final_particle_momentum_Ns,
                 out->final_field_momentum_Ns),
        swa_vadd(out->initial_particle_momentum_Ns,
                 out->initial_field_momentum_Ns)
    );
    momentum_scale = max3(
        swa_vnorm(swa_vadd(out->initial_particle_momentum_Ns,
                           out->initial_field_momentum_Ns)),
        swa_vnorm(swa_vadd(out->final_particle_momentum_Ns,
                           out->final_field_momentum_Ns)),
        max2(swa_vnorm(out->initial_particle_momentum_Ns), DBL_MIN)
    );
    out->total_momentum_relative =
        swa_vnorm(out->total_momentum_residual_Ns) / momentum_scale;
    out->momentum_passes = out->total_momentum_relative <=
                           limits->momentum_relative_tolerance;

    out->finite = isfinite(out->courant_number) &&
        isfinite(out->maximum_particle_displacement_cells) &&
        isfinite(out->initial_charge_C) &&
        isfinite(out->final_charge_C) &&
        isfinite(out->continuity_relative_max) &&
        isfinite(out->initial_gauss_relative_max) &&
        isfinite(out->final_gauss_relative_max) &&
        isfinite(out->transverse_current_partition_relative_max) &&
        isfinite(out->harmonic_ampere_relative) &&
        isfinite(out->initial_particle_energy_J) &&
        isfinite(out->final_particle_energy_J) &&
        isfinite(out->initial_field_energy_J) &&
        isfinite(out->final_field_energy_J) &&
        isfinite(out->midpoint_current_work_J) &&
        isfinite(out->total_energy_relative) &&
        finite_vec(out->total_momentum_residual_Ns) &&
        isfinite(out->total_momentum_relative) &&
        isfinite(out->maximum_speed_fraction_c);
    out->passes = out->finite && out->courant_passes &&
        out->particle_courant_passes && out->charge_closure_passes &&
        out->continuity_passes && out->initial_gauss_passes &&
        out->final_gauss_passes && out->magnetic_divergence_passes &&
        out->transverse_current_partition_passes &&
        out->harmonic_ampere_passes && out->subluminal_passes &&
        out->energy_passes &&
        out->momentum_passes;

    memcpy(state->particles, particle_final,
           particle_count * sizeof(*particle_final));
    memcpy(state->ex_face_V_m, ex_final, cells * sizeof(double));
    memcpy(state->ey_face_V_m, ey_final, cells * sizeof(double));
    memcpy(state->ez_face_V_m, ez_final, cells * sizeof(double));
    memcpy(state->by_cell_T, by_final, cells * sizeof(double));
    memcpy(state->bz_cell_T, bz_final, cells * sizeof(double));
    ok = 1;

cleanup:
    free(position_initial);
    free(position_final);
    free(rho_initial);
    free(rho_final);
    free(jx);
    free(jy);
    free(jz);
    free(ex_initial);
    free(ey_initial);
    free(ez_initial);
    free(by_initial);
    free(bz_initial);
    free(by_half);
    free(bz_half);
    free(ex_final);
    free(ey_final);
    free(ez_final);
    free(by_final);
    free(bz_final);
    free(particle_initial);
    free(particle_final);
    return ok;
}

