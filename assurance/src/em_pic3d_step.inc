
static void pic3d_faraday_update(
    const swa_em_pic3d *state,
    const double *bx_initial,
    const double *by_initial,
    const double *bz_initial,
    const double *ex,
    const double *ey,
    const double *ez,
    double duration_s,
    double *bx_final,
    double *by_final,
    double *bz_final
) {
    double *curl;
    double *curl_x;
    double *curl_y;
    double *curl_z;
    size_t index;
    if (state->cell_count > SIZE_MAX / (3U * sizeof(double))) {
        for (index = 0U; index < state->cell_count; ++index) {
            bx_final[index] = NAN;
            by_final[index] = NAN;
            bz_final[index] = NAN;
        }
        return;
    }
    curl = calloc(3U * state->cell_count, sizeof(double));
    if (curl == NULL) {
        for (index = 0U; index < state->cell_count; ++index) {
            bx_final[index] = NAN;
            by_final[index] = NAN;
            bz_final[index] = NAN;
        }
        return;
    }
    curl_x = curl;
    curl_y = curl_x + state->cell_count;
    curl_z = curl_y + state->cell_count;
    pic3d_curl_e_locations_to_b(
        state, ex, ey, ez, curl_x, curl_y, curl_z
    );
    for (index = 0U; index < state->cell_count; ++index) {
        bx_final[index] = bx_initial[index] - duration_s * curl_x[index];
        by_final[index] = by_initial[index] - duration_s * curl_y[index];
        bz_final[index] = bz_initial[index] - duration_s * curl_z[index];
    }
    free(curl);
}

static void pic3d_ampere_update(
    const swa_em_pic3d *state,
    const double *ex_initial,
    const double *ey_initial,
    const double *ez_initial,
    const double *bx_half,
    const double *by_half,
    const double *bz_half,
    const double *jx,
    const double *jy,
    const double *jz,
    double *ex_final,
    double *ey_final,
    double *ez_final
) {
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    const double c2dt = SWA_C * SWA_C * state->dt_s;
    const double current_factor = state->dt_s / SWA_EPS0;
    size_t i;
    size_t j;
    size_t k;
    for (k = 0U; k < state->nz; ++k) {
        const size_t kp = (k + 1U) % state->nz;
        for (j = 0U; j < state->ny; ++j) {
            const size_t jp = (j + 1U) % state->ny;
            for (i = 0U; i < state->nx; ++i) {
                const size_t ip = (i + 1U) % state->nx;
                const size_t center = pic3d_index(state, i, j, k);
                const double curl_x =
                    (bz_half[pic3d_index(state, i, jp, k)] -
                     bz_half[center]) / dy -
                    (by_half[pic3d_index(state, i, j, kp)] -
                     by_half[center]) / dz;
                const double curl_y =
                    (bx_half[pic3d_index(state, i, j, kp)] -
                     bx_half[center]) / dz -
                    (bz_half[pic3d_index(state, ip, j, k)] -
                     bz_half[center]) / dx;
                const double curl_z =
                    (by_half[pic3d_index(state, ip, j, k)] -
                     by_half[center]) / dx -
                    (bx_half[pic3d_index(state, i, jp, k)] -
                     bx_half[center]) / dy;
                ex_final[center] =
                    ex_initial[center] + c2dt * curl_x -
                    current_factor * jx[center];
                ey_final[center] =
                    ey_initial[center] + c2dt * curl_y -
                    current_factor * jy[center];
                ez_final[center] =
                    ez_initial[center] + c2dt * curl_z -
                    current_factor * jz[center];
            }
        }
    }
}

static int pic3d_push_particles(
    const swa_em_pic3d *state,
    const swa_em_pic3d_particle *initial,
    const double *ex_mid,
    const double *ey_mid,
    const double *ez_mid,
    const double *bx_half,
    const double *by_half,
    const double *bz_half,
    swa_em_pic3d_particle *final,
    double *maximum_displacement_cells,
    double *maximum_speed_fraction_c
) {
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    size_t p;
    *maximum_displacement_cells = 0.0;
    *maximum_speed_fraction_c = 0.0;
    for (p = 0U; p < state->particle_count; ++p) {
        const swa_vec3 e_initial = pic3d_interpolate_e(
            state, ex_mid, ey_mid, ez_mid,
            initial[p].position_unwrapped_m
        );
        const swa_vec3 b_initial = pic3d_interpolate_b(
            state, bx_half, by_half, bz_half,
            initial[p].position_unwrapped_m
        );
        swa_vec3 velocity_guess;
        swa_vec3 midpoint_position;
        swa_vec3 electric_midpoint;
        swa_vec3 magnetic_midpoint;
        swa_vec3 velocity_final;
        swa_vec3 displacement;
        double displacement_cells;
        double speed_fraction;
        if (!swa_em_pic3d_higuera_cary_push(
                initial[p].velocity_mps,
                initial[p].charge_C,
                initial[p].mass_kg,
                e_initial, b_initial, state->dt_s,
                &velocity_guess
            )) {
            return 0;
        }
        midpoint_position = swa_vadd(
            initial[p].position_unwrapped_m,
            swa_vscale(
                swa_vadd(initial[p].velocity_mps, velocity_guess),
                0.25 * state->dt_s
            )
        );
        electric_midpoint = pic3d_interpolate_e(
            state, ex_mid, ey_mid, ez_mid, midpoint_position
        );
        magnetic_midpoint = pic3d_interpolate_b(
            state, bx_half, by_half, bz_half, midpoint_position
        );
        if (!swa_em_pic3d_higuera_cary_push(
                initial[p].velocity_mps,
                initial[p].charge_C,
                initial[p].mass_kg,
                electric_midpoint, magnetic_midpoint,
                state->dt_s, &velocity_final
            )) {
            return 0;
        }
        displacement = swa_vscale(
            swa_vadd(initial[p].velocity_mps, velocity_final),
            0.5 * state->dt_s
        );
        final[p] = initial[p];
        final[p].velocity_mps = velocity_final;
        final[p].position_unwrapped_m = swa_vadd(
            initial[p].position_unwrapped_m, displacement
        );
        displacement_cells = pic3d_max3(
            fabs(displacement.x) / dx,
            fabs(displacement.y) / dy,
            fabs(displacement.z) / dz
        );
        speed_fraction = swa_vnorm(velocity_final) / SWA_C;
        *maximum_displacement_cells = pic3d_max2(
            *maximum_displacement_cells, displacement_cells
        );
        *maximum_speed_fraction_c = pic3d_max2(
            *maximum_speed_fraction_c, speed_fraction
        );
        if (!pic3d_particle_is_finite(&final[p])) {
            return 0;
        }
    }
    return 1;
}

static double pic3d_vector_arrays_relative_difference(
    const double *ax,
    const double *ay,
    const double *az,
    const double *bx,
    const double *by,
    const double *bz,
    size_t count
) {
    swa_kahan difference2 = {0.0, 0.0};
    swa_kahan a2 = {0.0, 0.0};
    swa_kahan b2 = {0.0, 0.0};
    size_t index;
    for (index = 0U; index < count; ++index) {
        const double dx = ax[index] - bx[index];
        const double dy = ay[index] - by[index];
        const double dz = az[index] - bz[index];
        swa_kahan_add(&difference2, dx * dx + dy * dy + dz * dz);
        swa_kahan_add(
            &a2,
            ax[index] * ax[index] +
            ay[index] * ay[index] +
            az[index] * az[index]
        );
        swa_kahan_add(
            &b2,
            bx[index] * bx[index] +
            by[index] * by[index] +
            bz[index] * bz[index]
        );
    }
    return sqrt(difference2.sum) /
           pic3d_max3(sqrt(a2.sum), sqrt(b2.sum), DBL_MIN);
}


static double pic3d_vector_arrays_l2(
    const double *x,
    const double *y,
    const double *z,
    size_t count
) {
    swa_kahan norm2 = {0.0, 0.0};
    size_t index;
    for (index = 0U; index < count; ++index) {
        swa_kahan_add(
            &norm2,
            x[index] * x[index] +
            y[index] * y[index] +
            z[index] * z[index]
        );
    }
    return sqrt(norm2.sum / (double)count);
}

static double pic3d_vector_arrays_difference_l2(
    const double *ax,
    const double *ay,
    const double *az,
    const double *bx,
    const double *by,
    const double *bz,
    size_t count
) {
    swa_kahan difference2 = {0.0, 0.0};
    size_t index;
    for (index = 0U; index < count; ++index) {
        const double dx = ax[index] - bx[index];
        const double dy = ay[index] - by[index];
        const double dz = az[index] - bz[index];
        swa_kahan_add(&difference2, dx * dx + dy * dy + dz * dz);
    }
    return sqrt(difference2.sum / (double)count);
}

static int pic3d_limits_valid(const swa_em_pic3d_limits *limits) {
    return limits != NULL &&
           isfinite(limits->continuity_relative_tolerance) &&
           limits->continuity_relative_tolerance >= 0.0 &&
           isfinite(limits->gauss_relative_tolerance) &&
           limits->gauss_relative_tolerance >= 0.0 &&
           isfinite(limits->magnetic_divergence_relative_tolerance) &&
           limits->magnetic_divergence_relative_tolerance >= 0.0 &&
           isfinite(limits->local_current_difference_relative_tolerance) &&
           limits->local_current_difference_relative_tolerance >= 0.0 &&
           isfinite(limits->spectral_oracle_relative_tolerance) &&
           limits->spectral_oracle_relative_tolerance >= 0.0 &&
           isfinite(limits->spectral_oracle_curl_relative_tolerance) &&
           limits->spectral_oracle_curl_relative_tolerance >= 0.0 &&
           isfinite(limits->energy_relative_tolerance) &&
           limits->energy_relative_tolerance >= 0.0 &&
           isfinite(limits->particle_work_relative_tolerance) &&
           limits->particle_work_relative_tolerance >= 0.0 &&
           isfinite(limits->field_work_relative_tolerance) &&
           limits->field_work_relative_tolerance >= 0.0 &&
           isfinite(limits->momentum_relative_tolerance) &&
           limits->momentum_relative_tolerance >= 0.0 &&
           isfinite(limits->nonlinear_relative_tolerance) &&
           limits->nonlinear_relative_tolerance >= 0.0 &&
           isfinite(limits->maximum_courant) &&
           limits->maximum_courant > 0.0 &&
           isfinite(limits->maximum_particle_cells_per_step) &&
           limits->maximum_particle_cells_per_step > 0.0 &&
           limits->maximum_nonlinear_iterations > 0U &&
           limits->maximum_nonlinear_iterations <= 64U;
}

int swa_em_pic3d_step(
    swa_em_pic3d *state,
    const swa_em_pic3d_limits *limits,
    swa_em_pic3d_result *out
) {
    const size_t count = state != NULL ? state->cell_count : 0U;
    const size_t doubles_needed = 26U;
    double *scratch;
    double *rho_initial;
    double *rho_final;
    double *raw_jx;
    double *raw_jy;
    double *raw_jz;
    double *jx;
    double *jy;
    double *jz;
    double *oracle_jx;
    double *oracle_jy;
    double *oracle_jz;
    double *bx_half;
    double *by_half;
    double *bz_half;
    double *ex_final;
    double *ey_final;
    double *ez_final;
    double *bx_final;
    double *by_final;
    double *bz_final;
    double *ex_mid;
    double *ey_mid;
    double *ez_mid;
    double *ex_mid_new;
    double *ey_mid_new;
    double *ez_mid_new;
    swa_em_pic3d_particle *particle_initial;
    swa_em_pic3d_particle *particle_final;
    swa_em_pic3d_projection_result projection;
    swa_vec3 transport_mean = {0.0, 0.0, 0.0};
    swa_vec3 local_transport_mean = {0.0, 0.0, 0.0};
    double initial_gauss_rms;
    double final_gauss_rms;
    double initial_divb_rms;
    double final_divb_rms;
    double charge_scale;
    double charge_l1_initial = 0.0;
    double charge_l1_final = 0.0;
    double energy_scale;
    double particle_work_scale;
    double field_work_scale;
    double momentum_scale;
    double current_mean_scale;
    size_t iteration;
    size_t index;
    int converged = 0;
    int ok = 0;
    if (!swa_em_pic3d_state_is_finite(state) ||
        !pic3d_limits_valid(limits) || out == NULL ||
        count > SIZE_MAX / doubles_needed ||
        doubles_needed * count > SIZE_MAX / sizeof(double) ||
        state->particle_count >
            SIZE_MAX / sizeof(swa_em_pic3d_particle)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    memset(&projection, 0, sizeof(projection));
    scratch = calloc(doubles_needed * count, sizeof(double));
    particle_initial = state->particle_count > 0U ?
        calloc(state->particle_count, sizeof(*particle_initial)) : NULL;
    particle_final = state->particle_count > 0U ?
        calloc(state->particle_count, sizeof(*particle_final)) : NULL;
    if (scratch == NULL ||
        (state->particle_count > 0U &&
         (particle_initial == NULL || particle_final == NULL))) {
        free(scratch);
        free(particle_initial);
        free(particle_final);
        return 0;
    }
    rho_initial = scratch;
    rho_final = rho_initial + count;
    raw_jx = rho_final + count;
    raw_jy = raw_jx + count;
    raw_jz = raw_jy + count;
    jx = raw_jz + count;
    jy = jx + count;
    jz = jy + count;
    oracle_jx = jz + count;
    oracle_jy = oracle_jx + count;
    oracle_jz = oracle_jy + count;
    bx_half = oracle_jz + count;
    by_half = bx_half + count;
    bz_half = by_half + count;
    ex_final = bz_half + count;
    ey_final = ex_final + count;
    ez_final = ey_final + count;
    bx_final = ez_final + count;
    by_final = bx_final + count;
    bz_final = by_final + count;
    ex_mid = bz_final + count;
    ey_mid = ex_mid + count;
    ez_mid = ey_mid + count;
    ex_mid_new = ez_mid + count;
    ey_mid_new = ex_mid_new + count;
    ez_mid_new = ey_mid_new + count;

    if (state->particle_count > 0U) {
        memcpy(
            particle_initial, state->particles,
            state->particle_count * sizeof(*particle_initial)
        );
    }
    if (!swa_em_pic3d_deposit_charge_density(
            state, state->particles, state->particle_count, rho_initial
        )) {
        goto cleanup;
    }
    {
        swa_kahan initial_charge = {0.0, 0.0};
        for (index = 0U; index < count; ++index) {
            const double cell_charge =
                rho_initial[index] * pic3d_cell_volume(state);
            swa_kahan_add(&initial_charge, cell_charge);
            charge_l1_initial += fabs(cell_charge);
        }
        out->initial_charge_C = initial_charge.sum;
    }
    pic3d_electric_constraint_metrics(
        state,
        state->ex_V_m, state->ey_V_m, state->ez_V_m,
        rho_initial,
        &out->initial_gauss_residual_max_C_m3,
        &initial_gauss_rms,
        &out->initial_gauss_relative_max
    );
    pic3d_magnetic_constraint_metrics(
        state, state->bx_T, state->by_T, state->bz_T,
        &out->initial_magnetic_divergence_max_T_m,
        &initial_divb_rms,
        &out->initial_magnetic_divergence_relative
    );

    pic3d_faraday_update(
        state,
        state->bx_T, state->by_T, state->bz_T,
        state->ex_V_m, state->ey_V_m, state->ez_V_m,
        0.5 * state->dt_s,
        bx_half, by_half, bz_half
    );
    memcpy(ex_mid, state->ex_V_m, count * sizeof(double));
    memcpy(ey_mid, state->ey_V_m, count * sizeof(double));
    memcpy(ez_mid, state->ez_V_m, count * sizeof(double));

    for (iteration = 0U;
         iteration < limits->maximum_nonlinear_iterations;
         ++iteration) {
        memset(raw_jx, 0, count * sizeof(double));
        memset(raw_jy, 0, count * sizeof(double));
        memset(raw_jz, 0, count * sizeof(double));
        memset(jx, 0, count * sizeof(double));
        memset(jy, 0, count * sizeof(double));
        memset(jz, 0, count * sizeof(double));
        memset(oracle_jx, 0, count * sizeof(double));
        memset(oracle_jy, 0, count * sizeof(double));
        memset(oracle_jz, 0, count * sizeof(double));
        memset(rho_final, 0, count * sizeof(double));

        if (state->particle_count > 0U) {
            if (!pic3d_push_particles(
                    state, particle_initial,
                    ex_mid, ey_mid, ez_mid,
                    bx_half, by_half, bz_half,
                    particle_final,
                    &out->maximum_particle_displacement_cells,
                    &out->maximum_speed_fraction_c
                ) ||
                !swa_em_pic3d_deposit_charge_density(
                    state, particle_final, state->particle_count, rho_final
                ) ||
                !swa_em_pic3d_deposit_transport_current(
                    state, particle_initial, particle_final,
                    state->particle_count,
                    raw_jx, raw_jy, raw_jz,
                    &out->maximum_path_segments,
                    &transport_mean
                ) ||
                !swa_em_pic3d_deposit_charge_conserving_current(
                    state, particle_initial, particle_final,
                    state->particle_count,
                    jx, jy, jz, &local_transport_mean
                ) ||
                !swa_em_pic3d_project_periodic_current(
                    state, rho_initial, rho_final,
                    jx, jy, jz,
                    limits->continuity_relative_tolerance,
                    limits->spectral_oracle_curl_relative_tolerance,
                    oracle_jx, oracle_jy, oracle_jz, &projection
                )) {
                goto cleanup;
            }
        } else {
            memcpy(rho_final, rho_initial, count * sizeof(double));
            memset(&projection, 0, sizeof(projection));
            projection.finite = 1;
            projection.zero_mode_passes = 1;
            projection.continuity_passes = 1;
            projection.correction_curl_passes = 1;
            projection.passes = 1;
            transport_mean = swa_v3(0.0, 0.0, 0.0);
            local_transport_mean = swa_v3(0.0, 0.0, 0.0);
        }

        pic3d_ampere_update(
            state,
            state->ex_V_m, state->ey_V_m, state->ez_V_m,
            bx_half, by_half, bz_half,
            jx, jy, jz,
            ex_final, ey_final, ez_final
        );
        for (index = 0U; index < count; ++index) {
            ex_mid_new[index] =
                0.5 * (state->ex_V_m[index] + ex_final[index]);
            ey_mid_new[index] =
                0.5 * (state->ey_V_m[index] + ey_final[index]);
            ez_mid_new[index] =
                0.5 * (state->ez_V_m[index] + ez_final[index]);
        }
        out->nonlinear_relative_residual =
            pic3d_vector_arrays_relative_difference(
                ex_mid, ey_mid, ez_mid,
                ex_mid_new, ey_mid_new, ez_mid_new, count
            );
        out->nonlinear_iterations = iteration + 1U;
        if (state->particle_count == 0U) {
            out->nonlinear_relative_residual = 0.0;
            converged = 1;
            break;
        }
        if (out->nonlinear_relative_residual <=
                limits->nonlinear_relative_tolerance) {
            converged = 1;
            break;
        }
        memcpy(ex_mid, ex_mid_new, count * sizeof(double));
        memcpy(ey_mid, ey_mid_new, count * sizeof(double));
        memcpy(ez_mid, ez_mid_new, count * sizeof(double));
    }

    pic3d_faraday_update(
        state,
        bx_half, by_half, bz_half,
        ex_final, ey_final, ez_final,
        0.5 * state->dt_s,
        bx_final, by_final, bz_final
    );

    {
        swa_kahan final_charge = {0.0, 0.0};
        for (index = 0U; index < count; ++index) {
            const double cell_charge =
                rho_final[index] * pic3d_cell_volume(state);
            swa_kahan_add(&final_charge, cell_charge);
            charge_l1_final += fabs(cell_charge);
        }
        out->final_charge_C = final_charge.sum;
    }
    out->charge_change_C =
        out->final_charge_C - out->initial_charge_C;
    pic3d_continuity_metrics(
        state, rho_initial, rho_final, jx, jy, jz,
        &out->continuity_residual_max_A_m3,
        &out->continuity_residual_rms_A_m3,
        &out->continuity_relative_max,
        &out->current_projection_zero_mode_A_m3
    );
    pic3d_electric_constraint_metrics(
        state,
        ex_final, ey_final, ez_final, rho_final,
        &out->final_gauss_residual_max_C_m3,
        &final_gauss_rms,
        &out->final_gauss_relative_max
    );
    pic3d_magnetic_constraint_metrics(
        state, bx_final, by_final, bz_final,
        &out->final_magnetic_divergence_max_T_m,
        &final_divb_rms,
        &out->final_magnetic_divergence_relative
    );

    out->raw_current_l2_A_m2 = pic3d_vector_arrays_l2(
        raw_jx, raw_jy, raw_jz, count
    );
    out->current_correction_l2_A_m2 =
        pic3d_vector_arrays_difference_l2(
            raw_jx, raw_jy, raw_jz,
            jx, jy, jz, count
        );
    out->local_current_difference_relative =
        pic3d_vector_arrays_relative_difference(
            raw_jx, raw_jy, raw_jz,
            jx, jy, jz, count
        );
    out->spectral_oracle_difference_relative =
        projection.correction_relative;
    out->current_correction_relative =
        out->local_current_difference_relative;
    out->current_projection_curl_relative =
        projection.correction_curl_relative;
    out->mean_current_A_m2 =
        pic3d_mean_current(state, jx, jy, jz);
    out->transport_mean_current_A_m2 = local_transport_mean;
    out->mean_current_residual_A_m2 = swa_vsub(
        out->mean_current_A_m2,
        out->transport_mean_current_A_m2
    );

    out->initial_particle_energy_J = pic3d_particle_energy_J(
        particle_initial, state->particle_count
    );
    out->final_particle_energy_J = pic3d_particle_energy_J(
        particle_final, state->particle_count
    );
    out->initial_field_energy_J = pic3d_field_array_energy_J(
        state,
        state->ex_V_m, state->ey_V_m, state->ez_V_m,
        state->bx_T, state->by_T, state->bz_T
    );
    out->final_field_energy_J = pic3d_field_array_energy_J(
        state,
        ex_final, ey_final, ez_final,
        bx_final, by_final, bz_final
    );
    out->particle_energy_change_J =
        out->final_particle_energy_J -
        out->initial_particle_energy_J;
    out->field_energy_change_J =
        out->final_field_energy_J -
        out->initial_field_energy_J;
    out->raw_current_work_J = pic3d_current_work_J(
        state, raw_jx, raw_jy, raw_jz,
        ex_mid, ey_mid, ez_mid
    );
    out->corrected_current_work_J = pic3d_current_work_J(
        state, jx, jy, jz,
        ex_mid_new, ey_mid_new, ez_mid_new
    );
    out->particle_midpoint_current_work_J =
        out->raw_current_work_J;
    out->field_midpoint_current_work_J =
        out->corrected_current_work_J;
    out->projection_work_correction_J =
        out->corrected_current_work_J -
        out->raw_current_work_J;
    out->particle_work_residual_J =
        out->particle_energy_change_J -
        out->raw_current_work_J;
    out->field_work_residual_J =
        out->field_energy_change_J +
        out->corrected_current_work_J;
    out->total_energy_residual_J =
        out->particle_energy_change_J +
        out->field_energy_change_J;

    energy_scale = pic3d_max4(
        fabs(out->initial_particle_energy_J +
             out->initial_field_energy_J),
        fabs(out->final_particle_energy_J +
             out->final_field_energy_J),
        fabs(out->particle_energy_change_J),
        pic3d_max2(fabs(out->field_energy_change_J), DBL_MIN)
    );
    out->total_energy_relative =
        fabs(out->total_energy_residual_J) / energy_scale;

    particle_work_scale = pic3d_max4(
        fabs(out->initial_particle_energy_J),
        fabs(out->final_particle_energy_J),
        fabs(out->particle_energy_change_J),
        pic3d_max2(fabs(out->raw_current_work_J), DBL_MIN)
    );
    out->particle_work_relative =
        fabs(out->particle_work_residual_J) / particle_work_scale;

    field_work_scale = pic3d_max4(
        fabs(out->initial_field_energy_J),
        fabs(out->final_field_energy_J),
        fabs(out->field_energy_change_J),
        pic3d_max2(fabs(out->corrected_current_work_J), DBL_MIN)
    );
    out->field_work_relative =
        fabs(out->field_work_residual_J) / field_work_scale;

    out->initial_particle_momentum_Ns =
        pic3d_particle_momentum_Ns(
            particle_initial, state->particle_count
        );
    out->final_particle_momentum_Ns =
        pic3d_particle_momentum_Ns(
            particle_final, state->particle_count
        );
    out->initial_field_momentum_Ns =
        pic3d_field_array_momentum_Ns(
            state,
            state->ex_V_m, state->ey_V_m, state->ez_V_m,
            state->bx_T, state->by_T, state->bz_T
        );
    out->final_field_momentum_Ns =
        pic3d_field_array_momentum_Ns(
            state,
            ex_final, ey_final, ez_final,
            bx_final, by_final, bz_final
        );
    out->total_momentum_residual_Ns = swa_vsub(
        swa_vadd(
            out->final_particle_momentum_Ns,
            out->final_field_momentum_Ns
        ),
        swa_vadd(
            out->initial_particle_momentum_Ns,
            out->initial_field_momentum_Ns
        )
    );
    momentum_scale = pic3d_max4(
        swa_vnorm(swa_vadd(
            out->initial_particle_momentum_Ns,
            out->initial_field_momentum_Ns
        )),
        swa_vnorm(swa_vadd(
            out->final_particle_momentum_Ns,
            out->final_field_momentum_Ns
        )),
        swa_vnorm(swa_vsub(
            out->final_particle_momentum_Ns,
            out->initial_particle_momentum_Ns
        )),
        pic3d_max2(
            swa_vnorm(swa_vsub(
                out->final_field_momentum_Ns,
                out->initial_field_momentum_Ns
            )),
            DBL_MIN
        )
    );
    out->total_momentum_relative =
        swa_vnorm(out->total_momentum_residual_Ns) / momentum_scale;

    out->courant_number = SWA_C * state->dt_s * sqrt(
        1.0 / (pic3d_dx(state) * pic3d_dx(state)) +
        1.0 / (pic3d_dy(state) * pic3d_dy(state)) +
        1.0 / (pic3d_dz(state) * pic3d_dz(state))
    );
    out->courant_passes =
        out->courant_number <= 1.0 + 128.0 * DBL_EPSILON &&
        out->courant_number <= limits->maximum_courant;
    out->particle_courant_passes =
        out->maximum_particle_displacement_cells <=
        limits->maximum_particle_cells_per_step;
    charge_scale = pic3d_max4(
        fabs(out->initial_charge_C),
        fabs(out->final_charge_C),
        charge_l1_initial,
        pic3d_max2(charge_l1_final, DBL_MIN)
    );
    out->charge_closure_passes =
        fabs(out->charge_change_C) <=
            4096.0 * DBL_EPSILON * charge_scale ||
        fabs(out->charge_change_C) / charge_scale <=
            limits->continuity_relative_tolerance;
    out->continuity_passes =
        projection.continuity_passes &&
        out->continuity_relative_max <=
            limits->continuity_relative_tolerance;
    out->initial_gauss_passes =
        out->initial_gauss_relative_max <=
            limits->gauss_relative_tolerance;
    out->final_gauss_passes =
        out->final_gauss_relative_max <=
            limits->gauss_relative_tolerance;
    out->initial_magnetic_divergence_passes =
        out->initial_magnetic_divergence_relative <=
            limits->magnetic_divergence_relative_tolerance;
    out->final_magnetic_divergence_passes =
        out->final_magnetic_divergence_relative <=
            limits->magnetic_divergence_relative_tolerance;
    out->local_current_passes =
        projection.continuity_before_relative <=
            limits->continuity_relative_tolerance &&
        out->local_current_difference_relative <=
            limits->local_current_difference_relative_tolerance;
    out->spectral_oracle_passes =
        projection.finite &&
        projection.continuity_passes &&
        projection.correction_curl_passes &&
        out->spectral_oracle_difference_relative <=
            pic3d_max2(
                limits->spectral_oracle_relative_tolerance,
                16384.0 * DBL_EPSILON
            );
    out->current_projection_passes =
        out->local_current_passes && out->spectral_oracle_passes;
    current_mean_scale = pic3d_max3(
        swa_vnorm(out->mean_current_A_m2),
        swa_vnorm(out->transport_mean_current_A_m2),
        DBL_MIN
    );
    out->mean_current_passes =
        swa_vnorm(out->mean_current_residual_A_m2) /
            current_mean_scale <=
                limits->continuity_relative_tolerance &&
        swa_vnorm(swa_vsub(transport_mean, local_transport_mean)) /
            pic3d_max3(
                swa_vnorm(transport_mean),
                swa_vnorm(local_transport_mean),
                DBL_MIN
            ) <= limits->continuity_relative_tolerance;
    out->nonlinear_passes = converged;
    out->subluminal_passes =
        out->maximum_speed_fraction_c < 1.0;
    out->energy_passes =
        out->total_energy_relative <=
            limits->energy_relative_tolerance;
    out->particle_work_passes =
        out->particle_work_relative <=
            limits->particle_work_relative_tolerance;
    out->field_work_passes =
        out->field_work_relative <=
            limits->field_work_relative_tolerance;
    out->momentum_passes =
        out->total_momentum_relative <=
            limits->momentum_relative_tolerance;

    out->finite =
        isfinite(out->courant_number) &&
        isfinite(out->maximum_particle_displacement_cells) &&
        isfinite(out->nonlinear_relative_residual) &&
        isfinite(out->initial_charge_C) &&
        isfinite(out->final_charge_C) &&
        isfinite(out->continuity_relative_max) &&
        isfinite(out->initial_gauss_relative_max) &&
        isfinite(out->final_gauss_relative_max) &&
        isfinite(out->initial_magnetic_divergence_relative) &&
        isfinite(out->final_magnetic_divergence_relative) &&
        isfinite(out->current_correction_relative) &&
        isfinite(out->local_current_difference_relative) &&
        isfinite(out->spectral_oracle_difference_relative) &&
        isfinite(out->current_projection_curl_relative) &&
        pic3d_finite_vec3(out->mean_current_A_m2) &&
        pic3d_finite_vec3(out->mean_current_residual_A_m2) &&
        isfinite(out->initial_particle_energy_J) &&
        isfinite(out->final_particle_energy_J) &&
        isfinite(out->initial_field_energy_J) &&
        isfinite(out->final_field_energy_J) &&
        isfinite(out->raw_current_work_J) &&
        isfinite(out->corrected_current_work_J) &&
        isfinite(out->particle_midpoint_current_work_J) &&
        isfinite(out->field_midpoint_current_work_J) &&
        isfinite(out->total_energy_relative) &&
        isfinite(out->particle_work_relative) &&
        isfinite(out->field_work_relative) &&
        pic3d_finite_vec3(out->total_momentum_residual_Ns) &&
        isfinite(out->total_momentum_relative) &&
        isfinite(out->maximum_speed_fraction_c) &&
        pic3d_finite_array(ex_final, count) &&
        pic3d_finite_array(ey_final, count) &&
        pic3d_finite_array(ez_final, count) &&
        pic3d_finite_array(bx_final, count) &&
        pic3d_finite_array(by_final, count) &&
        pic3d_finite_array(bz_final, count);

    out->passes = out->finite &&
        out->courant_passes &&
        out->particle_courant_passes &&
        out->charge_closure_passes &&
        out->continuity_passes &&
        out->initial_gauss_passes &&
        out->final_gauss_passes &&
        out->initial_magnetic_divergence_passes &&
        out->final_magnetic_divergence_passes &&
        out->current_projection_passes &&
        out->local_current_passes &&
        out->spectral_oracle_passes &&
        out->mean_current_passes &&
        out->nonlinear_passes &&
        out->subluminal_passes &&
        out->energy_passes &&
        out->particle_work_passes &&
        out->field_work_passes &&
        out->momentum_passes;

    if (!out->finite) {
        goto cleanup;
    }
    if (state->particle_count > 0U) {
        memcpy(
            state->particles, particle_final,
            state->particle_count * sizeof(*particle_final)
        );
    }
    memcpy(state->ex_V_m, ex_final, count * sizeof(double));
    memcpy(state->ey_V_m, ey_final, count * sizeof(double));
    memcpy(state->ez_V_m, ez_final, count * sizeof(double));
    memcpy(state->bx_T, bx_final, count * sizeof(double));
    memcpy(state->by_T, by_final, count * sizeof(double));
    memcpy(state->bz_T, bz_final, count * sizeof(double));
    ok = 1;

cleanup:
    free(scratch);
    free(particle_initial);
    free(particle_final);
    return ok;
}

int swa_em_pic3d_write_receipt(
    FILE *fp,
    const swa_em_pic3d *state,
    const swa_em_pic3d_limits *limits,
    const swa_em_pic3d_result *result
) {
    if (fp == NULL || state == NULL || limits == NULL || result == NULL) {
        return 0;
    }
    return fprintf(
        fp,
        "{\n"
        "  \"schema\": \"spacewind.em-pic3d-reference/v2\",\n"
        "  \"grid\": [%zu, %zu, %zu],\n"
        "  \"particles\": %zu,\n"
        "  \"length_m\": [%.17g, %.17g, %.17g],\n"
        "  \"dt_s\": %.17g,\n"
        "  \"courant\": %.17g,\n"
        "  \"maximum_particle_displacement_cells\": %.17g,\n"
        "  \"maximum_path_segments\": %zu,\n"
        "  \"nonlinear_iterations\": %zu,\n"
        "  \"nonlinear_relative_residual\": %.17g,\n"
        "  \"continuity_relative_max\": %.17g,\n"
        "  \"initial_gauss_relative_max\": %.17g,\n"
        "  \"final_gauss_relative_max\": %.17g,\n"
        "  \"initial_div_b_relative\": %.17g,\n"
        "  \"final_div_b_relative\": %.17g,\n"
        "  \"path_to_local_current_relative\": %.17g,\n"
        "  \"local_to_spectral_oracle_relative\": %.17g,\n"
        "  \"spectral_correction_curl_relative\": %.17g,\n"
        "  \"particle_midpoint_current_work_J\": %.17g,\n"
        "  \"field_midpoint_current_work_J\": %.17g,\n"
        "  \"particle_work_relative\": %.17g,\n"
        "  \"field_work_relative\": %.17g,\n"
        "  \"total_energy_relative\": %.17g,\n"
        "  \"total_momentum_relative\": %.17g,\n"
        "  \"maximum_speed_fraction_c\": %.17g,\n"
        "  \"local_current_passes\": %s,\n"
        "  \"spectral_oracle_passes\": %s,\n"
        "  \"passes\": %s,\n"
        "  \"current_method\": \"six-permutation symmetric coordinate-split charge-conserving face current with complete unwrapped harmonic transport\",\n"
        "  \"oracle_method\": \"independent direct-DFT minimum-norm periodic continuity projection\",\n"
        "  \"nonclaim\": \"this is a bounded periodic first-order-CIC 3D3V reference; the symmetric coordinate-split current is not yet claimed equivalent to a production Esirkepov deposition, and the model does not establish open-boundary plasma-wing force, material or sheath behavior, chamber validity, or propulsion\"\n"
        "}\n",
        state->nx, state->ny, state->nz,
        state->particle_count,
        state->length_x_m, state->length_y_m, state->length_z_m,
        state->dt_s,
        result->courant_number,
        result->maximum_particle_displacement_cells,
        result->maximum_path_segments,
        result->nonlinear_iterations,
        result->nonlinear_relative_residual,
        result->continuity_relative_max,
        result->initial_gauss_relative_max,
        result->final_gauss_relative_max,
        result->initial_magnetic_divergence_relative,
        result->final_magnetic_divergence_relative,
        result->local_current_difference_relative,
        result->spectral_oracle_difference_relative,
        result->current_projection_curl_relative,
        result->particle_midpoint_current_work_J,
        result->field_midpoint_current_work_J,
        result->particle_work_relative,
        result->field_work_relative,
        result->total_energy_relative,
        result->total_momentum_relative,
        result->maximum_speed_fraction_c,
        result->local_current_passes ? "true" : "false",
        result->spectral_oracle_passes ? "true" : "false",
        result->passes ? "true" : "false"
    ) > 0;
}
