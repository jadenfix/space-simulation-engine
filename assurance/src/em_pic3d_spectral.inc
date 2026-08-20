
static void pic3d_continuity_metrics(
    const swa_em_pic3d *state,
    const double *rho_initial,
    const double *rho_final,
    const double *jx,
    const double *jy,
    const double *jz,
    double *maximum,
    double *rms,
    double *relative,
    double *zero_mode
) {
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    swa_kahan squared = {0.0, 0.0};
    swa_kahan mean = {0.0, 0.0};
    double scale = DBL_MIN;
    size_t i;
    size_t j;
    size_t k;
    *maximum = 0.0;
    for (k = 0U; k < state->nz; ++k) {
        const size_t kp = (k + 1U) % state->nz;
        for (j = 0U; j < state->ny; ++j) {
            const size_t jp = (j + 1U) % state->ny;
            for (i = 0U; i < state->nx; ++i) {
                const size_t ip = (i + 1U) % state->nx;
                const size_t center = pic3d_index(state, i, j, k);
                const double time_derivative =
                    (rho_final[center] - rho_initial[center]) /
                    state->dt_s;
                const double divergence =
                    (jx[pic3d_index(state, ip, j, k)] - jx[center]) / dx +
                    (jy[pic3d_index(state, i, jp, k)] - jy[center]) / dy +
                    (jz[pic3d_index(state, i, j, kp)] - jz[center]) / dz;
                const double residual = time_derivative + divergence;
                scale = pic3d_max3(
                    scale, fabs(time_derivative), fabs(divergence)
                );
                *maximum = pic3d_max2(*maximum, fabs(residual));
                swa_kahan_add(&squared, residual * residual);
                swa_kahan_add(&mean, residual);
            }
        }
    }
    *rms = sqrt(squared.sum / (double)state->cell_count);
    *relative = *maximum / scale;
    *zero_mode = mean.sum / (double)state->cell_count;
}

static int pic3d_solve_periodic_laplacian(
    const swa_em_pic3d *state,
    const double *rhs,
    double *solution,
    double *zero_mode,
    double *maximum_imaginary
) {
    const size_t count = state->cell_count;
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    long double complex *coefficient;
    size_t kx;
    size_t ky;
    size_t kz;
    size_t i;
    size_t j;
    size_t k;
    double max_imaginary = 0.0;
    if (rhs == NULL || solution == NULL ||
        zero_mode == NULL || maximum_imaginary == NULL ||
        !pic3d_finite_array(rhs, count) ||
        count > SIZE_MAX / sizeof(*coefficient)) {
        return 0;
    }
    coefficient = calloc(count, sizeof(*coefficient));
    if (coefficient == NULL) {
        return 0;
    }
    *zero_mode = 0.0;
    for (kz = 0U; kz < state->nz; ++kz) {
        const long mode_z = kz <= state->nz / 2U ?
            (long)kz : (long)kz - (long)state->nz;
        const long double sz = sinl(
            (long double)SWA_PI * (long double)mode_z /
            (long double)state->nz
        );
        for (ky = 0U; ky < state->ny; ++ky) {
            const long mode_y = ky <= state->ny / 2U ?
                (long)ky : (long)ky - (long)state->ny;
            const long double sy = sinl(
                (long double)SWA_PI * (long double)mode_y /
                (long double)state->ny
            );
            for (kx = 0U; kx < state->nx; ++kx) {
                const long mode_x = kx <= state->nx / 2U ?
                    (long)kx : (long)kx - (long)state->nx;
                const long double sx = sinl(
                    (long double)SWA_PI * (long double)mode_x /
                    (long double)state->nx
                );
                long double complex sum = 0.0L + 0.0L * I;
                const long double eigenvalue =
                    -4.0L * sx * sx /
                        ((long double)dx * (long double)dx) -
                    4.0L * sy * sy /
                        ((long double)dy * (long double)dy) -
                    4.0L * sz * sz /
                        ((long double)dz * (long double)dz);
                for (k = 0U; k < state->nz; ++k) {
                    for (j = 0U; j < state->ny; ++j) {
                        for (i = 0U; i < state->nx; ++i) {
                            const long double angle =
                                -2.0L * (long double)SWA_PI * (
                                    (long double)mode_x * (long double)i /
                                        (long double)state->nx +
                                    (long double)mode_y * (long double)j /
                                        (long double)state->ny +
                                    (long double)mode_z * (long double)k /
                                        (long double)state->nz
                                );
                            sum += (long double)rhs[
                                pic3d_index(state, i, j, k)
                            ] * (cosl(angle) + I * sinl(angle));
                        }
                    }
                }
                if (mode_x == 0L && mode_y == 0L && mode_z == 0L) {
                    *zero_mode = (double)(creall(sum) /
                        (long double)count);
                    coefficient[pic3d_index(state, kx, ky, kz)] =
                        0.0L + 0.0L * I;
                } else {
                    if (!isfinite((double)eigenvalue) ||
                        !(eigenvalue < 0.0L)) {
                        free(coefficient);
                        return 0;
                    }
                    coefficient[pic3d_index(state, kx, ky, kz)] =
                        sum / eigenvalue;
                }
            }
        }
    }
    for (k = 0U; k < state->nz; ++k) {
        for (j = 0U; j < state->ny; ++j) {
            for (i = 0U; i < state->nx; ++i) {
                long double complex sum = 0.0L + 0.0L * I;
                for (kz = 0U; kz < state->nz; ++kz) {
                    const long mode_z = kz <= state->nz / 2U ?
                        (long)kz : (long)kz - (long)state->nz;
                    for (ky = 0U; ky < state->ny; ++ky) {
                        const long mode_y = ky <= state->ny / 2U ?
                            (long)ky : (long)ky - (long)state->ny;
                        for (kx = 0U; kx < state->nx; ++kx) {
                            const long mode_x = kx <= state->nx / 2U ?
                                (long)kx : (long)kx - (long)state->nx;
                            const long double angle =
                                2.0L * (long double)SWA_PI * (
                                    (long double)mode_x * (long double)i /
                                        (long double)state->nx +
                                    (long double)mode_y * (long double)j /
                                        (long double)state->ny +
                                    (long double)mode_z * (long double)k /
                                        (long double)state->nz
                                );
                            sum += coefficient[
                                pic3d_index(state, kx, ky, kz)
                            ] * (cosl(angle) + I * sinl(angle));
                        }
                    }
                }
                sum /= (long double)count;
                solution[pic3d_index(state, i, j, k)] =
                    (double)creall(sum);
                max_imaginary = pic3d_max2(
                    max_imaginary, fabs((double)cimagl(sum))
                );
            }
        }
    }
    free(coefficient);
    *maximum_imaginary = max_imaginary;
    return pic3d_finite_array(solution, count) &&
           isfinite(*zero_mode) && isfinite(*maximum_imaginary);
}

static void pic3d_gradient_cell_to_faces(
    const swa_em_pic3d *state,
    const double *scalar,
    double *gradient_x,
    double *gradient_y,
    double *gradient_z
) {
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    size_t i;
    size_t j;
    size_t k;
    for (k = 0U; k < state->nz; ++k) {
        const size_t km = (k + state->nz - 1U) % state->nz;
        for (j = 0U; j < state->ny; ++j) {
            const size_t jm = (j + state->ny - 1U) % state->ny;
            for (i = 0U; i < state->nx; ++i) {
                const size_t im = (i + state->nx - 1U) % state->nx;
                const size_t center = pic3d_index(state, i, j, k);
                gradient_x[center] =
                    (scalar[center] -
                     scalar[pic3d_index(state, im, j, k)]) / dx;
                gradient_y[center] =
                    (scalar[center] -
                     scalar[pic3d_index(state, i, jm, k)]) / dy;
                gradient_z[center] =
                    (scalar[center] -
                     scalar[pic3d_index(state, i, j, km)]) / dz;
            }
        }
    }
}

static void pic3d_curl_e_locations_to_b(
    const swa_em_pic3d *state,
    const double *ex,
    const double *ey,
    const double *ez,
    double *curl_x,
    double *curl_y,
    double *curl_z
) {
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    size_t i;
    size_t j;
    size_t k;
    for (k = 0U; k < state->nz; ++k) {
        const size_t km = (k + state->nz - 1U) % state->nz;
        for (j = 0U; j < state->ny; ++j) {
            const size_t jm = (j + state->ny - 1U) % state->ny;
            for (i = 0U; i < state->nx; ++i) {
                const size_t im = (i + state->nx - 1U) % state->nx;
                const size_t center = pic3d_index(state, i, j, k);
                curl_x[center] =
                    (ez[center] -
                     ez[pic3d_index(state, i, jm, k)]) / dy -
                    (ey[center] -
                     ey[pic3d_index(state, i, j, km)]) / dz;
                curl_y[center] =
                    (ex[center] -
                     ex[pic3d_index(state, i, j, km)]) / dz -
                    (ez[center] -
                     ez[pic3d_index(state, im, j, k)]) / dx;
                curl_z[center] =
                    (ey[center] -
                     ey[pic3d_index(state, im, j, k)]) / dx -
                    (ex[center] -
                     ex[pic3d_index(state, i, jm, k)]) / dy;
            }
        }
    }
}

int swa_em_pic3d_project_periodic_current(
    const swa_em_pic3d *state,
    const double *rho_initial_C_m3,
    const double *rho_final_C_m3,
    const double *raw_jx_A_m2,
    const double *raw_jy_A_m2,
    const double *raw_jz_A_m2,
    double relative_tolerance,
    double curl_relative_tolerance,
    double *corrected_jx_A_m2,
    double *corrected_jy_A_m2,
    double *corrected_jz_A_m2,
    swa_em_pic3d_projection_result *out
) {
    const size_t count = state != NULL ? state->cell_count : 0U;
    double *scratch;
    double *defect;
    double *potential;
    double *gradient_x;
    double *gradient_y;
    double *gradient_z;
    double *curl_x;
    double *curl_y;
    double *curl_z;
    double dummy_rms;
    double dummy_zero;
    double imaginary;
    swa_kahan raw2 = {0.0, 0.0};
    swa_kahan correction2 = {0.0, 0.0};
    double raw_scale;
    double curl_scale = DBL_MIN;
    const double minimum_spacing = state != NULL ?
        fmin(pic3d_dx(state), fmin(pic3d_dy(state), pic3d_dz(state))) :
        NAN;
    size_t index;
    int ok = 0;
    if (!pic3d_valid_shape(state) ||
        rho_initial_C_m3 == NULL || rho_final_C_m3 == NULL ||
        raw_jx_A_m2 == NULL || raw_jy_A_m2 == NULL ||
        raw_jz_A_m2 == NULL ||
        corrected_jx_A_m2 == NULL || corrected_jy_A_m2 == NULL ||
        corrected_jz_A_m2 == NULL || out == NULL ||
        !isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
        !isfinite(curl_relative_tolerance) ||
        curl_relative_tolerance < 0.0 ||
        !pic3d_finite_array(rho_initial_C_m3, count) ||
        !pic3d_finite_array(rho_final_C_m3, count) ||
        !pic3d_finite_array(raw_jx_A_m2, count) ||
        !pic3d_finite_array(raw_jy_A_m2, count) ||
        !pic3d_finite_array(raw_jz_A_m2, count) ||
        !isfinite(minimum_spacing) || !(minimum_spacing > 0.0) ||
        count > SIZE_MAX / (8U * sizeof(double))) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    scratch = calloc(8U * count, sizeof(double));
    if (scratch == NULL) {
        return 0;
    }
    defect = scratch;
    potential = defect + count;
    gradient_x = potential + count;
    gradient_y = gradient_x + count;
    gradient_z = gradient_y + count;
    curl_x = gradient_z + count;
    curl_y = curl_x + count;
    curl_z = curl_y + count;

    pic3d_continuity_metrics(
        state, rho_initial_C_m3, rho_final_C_m3,
        raw_jx_A_m2, raw_jy_A_m2, raw_jz_A_m2,
        &out->continuity_before_max_A_m3,
        &dummy_rms,
        &out->continuity_before_relative,
        &out->zero_mode_A_m3
    );
    {
        const double dx = pic3d_dx(state);
        const double dy = pic3d_dy(state);
        const double dz = pic3d_dz(state);
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
                    defect[center] =
                        (rho_final_C_m3[center] -
                         rho_initial_C_m3[center]) / state->dt_s +
                        (raw_jx_A_m2[
                            pic3d_index(state, ip, j, k)
                         ] - raw_jx_A_m2[center]) / dx +
                        (raw_jy_A_m2[
                            pic3d_index(state, i, jp, k)
                         ] - raw_jy_A_m2[center]) / dy +
                        (raw_jz_A_m2[
                            pic3d_index(state, i, j, kp)
                         ] - raw_jz_A_m2[center]) / dz;
                }
            }
        }
    }
    if (!pic3d_solve_periodic_laplacian(
            state, defect, potential, &dummy_zero, &imaginary
        )) {
        goto cleanup;
    }
    out->zero_mode_A_m3 = dummy_zero;
    pic3d_gradient_cell_to_faces(
        state, potential, gradient_x, gradient_y, gradient_z
    );
    for (index = 0U; index < count; ++index) {
        corrected_jx_A_m2[index] =
            raw_jx_A_m2[index] - gradient_x[index];
        corrected_jy_A_m2[index] =
            raw_jy_A_m2[index] - gradient_y[index];
        corrected_jz_A_m2[index] =
            raw_jz_A_m2[index] - gradient_z[index];
        swa_kahan_add(
            &raw2,
            raw_jx_A_m2[index] * raw_jx_A_m2[index] +
            raw_jy_A_m2[index] * raw_jy_A_m2[index] +
            raw_jz_A_m2[index] * raw_jz_A_m2[index]
        );
        swa_kahan_add(
            &correction2,
            gradient_x[index] * gradient_x[index] +
            gradient_y[index] * gradient_y[index] +
            gradient_z[index] * gradient_z[index]
        );
    }
    out->raw_current_l2_A_m2 =
        sqrt(raw2.sum / (double)count);
    out->correction_l2_A_m2 =
        sqrt(correction2.sum / (double)count);
    raw_scale = pic3d_max3(
        out->raw_current_l2_A_m2,
        out->correction_l2_A_m2,
        DBL_MIN
    );
    out->correction_relative =
        out->correction_l2_A_m2 / raw_scale;

    pic3d_continuity_metrics(
        state, rho_initial_C_m3, rho_final_C_m3,
        corrected_jx_A_m2, corrected_jy_A_m2, corrected_jz_A_m2,
        &out->continuity_after_max_A_m3,
        &out->continuity_after_rms_A_m3,
        &out->continuity_after_relative,
        &dummy_zero
    );
    pic3d_curl_e_locations_to_b(
        state, gradient_x, gradient_y, gradient_z,
        curl_x, curl_y, curl_z
    );
    for (index = 0U; index < count; ++index) {
        out->correction_curl_max_A_m3 = pic3d_max4(
            out->correction_curl_max_A_m3,
            fabs(curl_x[index]), fabs(curl_y[index]), fabs(curl_z[index])
        );
        curl_scale = pic3d_max4(
            curl_scale,
            fabs(gradient_x[index]) / minimum_spacing,
            fabs(gradient_y[index]) / minimum_spacing,
            fabs(gradient_z[index]) / minimum_spacing
        );
    }
    out->correction_curl_relative =
        out->correction_curl_max_A_m3 / curl_scale;
    {
        const double defect_scale = pic3d_max3(
            out->continuity_before_max_A_m3,
            out->continuity_after_max_A_m3,
            DBL_MIN
        );
        const double roundoff_floor =
            8192.0 * DBL_EPSILON * defect_scale;
        out->zero_mode_passes =
            fabs(out->zero_mode_A_m3) <= roundoff_floor ||
            fabs(out->zero_mode_A_m3) <=
                relative_tolerance * defect_scale;
        out->continuity_passes =
            out->continuity_after_max_A_m3 <= roundoff_floor ||
            out->continuity_after_relative <= relative_tolerance;
    }
    out->correction_curl_passes =
        out->correction_curl_max_A_m3 <=
            8192.0 * DBL_EPSILON * curl_scale ||
        out->correction_curl_relative <= curl_relative_tolerance;
    out->finite =
        pic3d_finite_array(corrected_jx_A_m2, count) &&
        pic3d_finite_array(corrected_jy_A_m2, count) &&
        pic3d_finite_array(corrected_jz_A_m2, count) &&
        isfinite(out->raw_current_l2_A_m2) &&
        isfinite(out->correction_l2_A_m2) &&
        isfinite(out->correction_relative) &&
        isfinite(out->continuity_before_relative) &&
        isfinite(out->continuity_after_relative) &&
        isfinite(out->zero_mode_A_m3) &&
        isfinite(out->correction_curl_relative) &&
        isfinite(imaginary);
    out->passes = out->finite && out->zero_mode_passes &&
                  out->continuity_passes &&
                  out->correction_curl_passes;
    ok = 1;

cleanup:
    free(scratch);
    return ok;
}

int swa_em_pic3d_initialize_gauss_field(
    swa_em_pic3d *state,
    swa_vec3 requested_mean_electric_field_V_m,
    double absolute_charge_tolerance_C
) {
    const size_t count = state != NULL ? state->cell_count : 0U;
    const double volume = state != NULL ? pic3d_cell_volume(state) : NAN;
    double *scratch;
    double *rho;
    double *rhs;
    double *potential;
    double *gradient_x;
    double *gradient_y;
    double *gradient_z;
    swa_kahan total_charge = {0.0, 0.0};
    double zero_mode;
    double imaginary;
    size_t index;
    int ok = 0;
    if (!pic3d_valid_shape(state) ||
        !pic3d_finite_vec3(requested_mean_electric_field_V_m) ||
        !isfinite(absolute_charge_tolerance_C) ||
        absolute_charge_tolerance_C < 0.0 ||
        !isfinite(volume) || !(volume > 0.0) ||
        count > SIZE_MAX / (6U * sizeof(double))) {
        return 0;
    }
    scratch = calloc(6U * count, sizeof(double));
    if (scratch == NULL) {
        return 0;
    }
    rho = scratch;
    rhs = rho + count;
    potential = rhs + count;
    gradient_x = potential + count;
    gradient_y = gradient_x + count;
    gradient_z = gradient_y + count;
    if (!swa_em_pic3d_deposit_charge_density(
            state, state->particles, state->particle_count, rho
        )) {
        goto cleanup;
    }
    for (index = 0U; index < count; ++index) {
        swa_kahan_add(&total_charge, rho[index] * volume);
        rhs[index] = -rho[index] / SWA_EPS0;
    }
    if (fabs(total_charge.sum) > absolute_charge_tolerance_C) {
        goto cleanup;
    }
    if (!pic3d_solve_periodic_laplacian(
            state, rhs, potential, &zero_mode, &imaginary
        )) {
        goto cleanup;
    }
    pic3d_gradient_cell_to_faces(
        state, potential, gradient_x, gradient_y, gradient_z
    );
    for (index = 0U; index < count; ++index) {
        state->ex_V_m[index] =
            requested_mean_electric_field_V_m.x - gradient_x[index];
        state->ey_V_m[index] =
            requested_mean_electric_field_V_m.y - gradient_y[index];
        state->ez_V_m[index] =
            requested_mean_electric_field_V_m.z - gradient_z[index];
    }
    ok = swa_em_pic3d_state_is_finite(state) &&
         isfinite(zero_mode) && isfinite(imaginary);

cleanup:
    free(scratch);
    return ok;
}
