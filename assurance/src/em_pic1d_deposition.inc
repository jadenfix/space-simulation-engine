static int deposit_charge(
    const swa_em_pic1d *state,
    const double *position_unwrapped_m,
    double *rho_C_m3,
    double *total_charge_C,
    double *absolute_charge_C
) {
    const double dx_m = state->domain_length_m / (double)state->cells;
    const double inverse_volume =
        1.0 / (dx_m * state->cross_section_area_m2);
    double *correction;
    swa_kahan total = {0.0, 0.0};
    swa_kahan absolute = {0.0, 0.0};
    size_t particle_index;
    if (position_unwrapped_m == NULL || rho_C_m3 == NULL ||
        total_charge_C == NULL || absolute_charge_C == NULL ||
        !isfinite(inverse_volume)) {
        return 0;
    }
    correction = calloc(state->cells, sizeof(double));
    if (correction == NULL) {
        return 0;
    }
    memset(rho_C_m3, 0, state->cells * sizeof(double));
    for (particle_index = 0U;
         particle_index < state->particle_count;
         ++particle_index) {
        const swa_em_pic1d_particle *particle =
            &state->particles[particle_index];
        const double charge = particle->charge_C * particle->macro_weight;
        const double wrapped = wrap_periodic(
            position_unwrapped_m[particle_index], state->domain_length_m
        );
        const double coordinate = wrapped / dx_m - 0.5;
        const double lower_real = floor(coordinate);
        double upper_weight;
        double lower_weight;
        long long lower_signed;
        size_t lower;
        size_t upper;
        double density_scale;
        double y;
        double t;
        if (!isfinite(charge) || !isfinite(coordinate) ||
            lower_real < (double)LLONG_MIN ||
            lower_real > (double)LLONG_MAX) {
            free(correction);
            return 0;
        }
        upper_weight = coordinate - lower_real;
        lower_weight = 1.0 - upper_weight;
        lower_signed = (long long)lower_real;
        lower = periodic_index(lower_signed, state->cells);
        upper = (lower + 1U) % state->cells;
        density_scale = charge * inverse_volume;
        y = density_scale * lower_weight - correction[lower];
        t = rho_C_m3[lower] + y;
        correction[lower] = (t - rho_C_m3[lower]) - y;
        rho_C_m3[lower] = t;
        y = density_scale * upper_weight - correction[upper];
        t = rho_C_m3[upper] + y;
        correction[upper] = (t - rho_C_m3[upper]) - y;
        rho_C_m3[upper] = t;
        swa_kahan_add(&total, charge);
        swa_kahan_add(&absolute, fabs(charge));
    }
    free(correction);
    *total_charge_C = total.sum;
    *absolute_charge_C = absolute.sum;
    return finite_array(rho_C_m3, state->cells) &&
           isfinite(*total_charge_C) && isfinite(*absolute_charge_C);
}

int swa_em_pic1d_initialize_gauss_field(
    swa_em_pic1d *state,
    double requested_mean_ex_V_m,
    double absolute_charge_tolerance_C
) {
    double *position;
    double *rho;
    double total_charge;
    double absolute_charge;
    double dx_m;
    swa_kahan mean = {0.0, 0.0};
    size_t i;
    if (!state_values_valid(state) ||
        !isfinite(requested_mean_ex_V_m) ||
        !isfinite(absolute_charge_tolerance_C) ||
        absolute_charge_tolerance_C < 0.0) {
        return 0;
    }
    position = calloc(count_or_one(state->particle_count), sizeof(double));
    rho = calloc(state->cells, sizeof(double));
    if (position == NULL || rho == NULL) {
        free(position);
        free(rho);
        return 0;
    }
    for (i = 0U; i < state->particle_count; ++i) {
        position[i] = state->particles[i].position_unwrapped_m;
    }
    if (!deposit_charge(state, position, rho,
                        &total_charge, &absolute_charge) ||
        fabs(total_charge) > max2(absolute_charge_tolerance_C,
                                  32.0 * DBL_EPSILON * absolute_charge)) {
        free(position);
        free(rho);
        return 0;
    }
    dx_m = state->domain_length_m / (double)state->cells;
    state->ex_face_V_m[0] = 0.0;
    for (i = 0U; i + 1U < state->cells; ++i) {
        state->ex_face_V_m[i + 1U] = state->ex_face_V_m[i] +
            dx_m * rho[i] / SWA_EPS0;
    }
    for (i = 0U; i < state->cells; ++i) {
        swa_kahan_add(&mean, state->ex_face_V_m[i]);
    }
    {
        const double shift = requested_mean_ex_V_m -
                             mean.sum / (double)state->cells;
        for (i = 0U; i < state->cells; ++i) {
            state->ex_face_V_m[i] += shift;
        }
    }
    free(position);
    free(rho);
    return finite_array(state->ex_face_V_m, state->cells);
}

static swa_vec3 interpolate_fields(
    const swa_em_pic1d *state,
    const double *ex_face,
    const double *ey_face,
    const double *ez_face,
    const double *by_cell,
    const double *bz_cell,
    double position_m,
    swa_vec3 *magnetic_T
) {
    const double dx_m = state->domain_length_m / (double)state->cells;
    const double wrapped = wrap_periodic(position_m, state->domain_length_m);
    const double face_coordinate = wrapped / dx_m;
    const long long face_lower_signed = (long long)floor(face_coordinate);
    const size_t face_lower = periodic_index(face_lower_signed, state->cells);
    const size_t face_upper = (face_lower + 1U) % state->cells;
    const double face_fraction = face_coordinate - floor(face_coordinate);
    const double cell_coordinate = wrapped / dx_m - 0.5;
    const long long cell_lower_signed = (long long)floor(cell_coordinate);
    const size_t cell_lower = periodic_index(cell_lower_signed, state->cells);
    const size_t cell_upper = (cell_lower + 1U) % state->cells;
    const double cell_fraction = cell_coordinate - floor(cell_coordinate);
    const swa_vec3 electric = swa_v3(
        ex_face[face_lower] * (1.0 - face_fraction) +
            ex_face[face_upper] * face_fraction,
        ey_face[face_lower] * (1.0 - face_fraction) +
            ey_face[face_upper] * face_fraction,
        ez_face[face_lower] * (1.0 - face_fraction) +
            ez_face[face_upper] * face_fraction
    );
    *magnetic_T = swa_v3(
        state->bx_uniform_T,
        by_cell[cell_lower] * (1.0 - cell_fraction) +
            by_cell[cell_upper] * cell_fraction,
        bz_cell[cell_lower] * (1.0 - cell_fraction) +
            bz_cell[cell_upper] * cell_fraction
    );
    return electric;
}

static int relativistic_boris(
    swa_vec3 velocity_initial_mps,
    swa_vec3 electric_V_m,
    swa_vec3 magnetic_T,
    double charge_to_mass_C_kg,
    double dt_s,
    swa_vec3 *velocity_final_mps
) {
    const double c2 = SWA_C * SWA_C;
    const double speed2 = swa_vdot(velocity_initial_mps,
                                   velocity_initial_mps);
    double gamma_initial;
    double gamma_minus;
    double gamma_final;
    double factor;
    swa_vec3 u;
    swa_vec3 u_minus;
    swa_vec3 t;
    swa_vec3 s;
    swa_vec3 u_prime;
    swa_vec3 u_plus;
    swa_vec3 u_final;
    if (velocity_final_mps == NULL || !finite_vec(velocity_initial_mps) ||
        !finite_vec(electric_V_m) || !finite_vec(magnetic_T) ||
        !isfinite(charge_to_mass_C_kg) || !isfinite(dt_s) ||
        !(dt_s > 0.0) || !(speed2 >= 0.0) || speed2 >= c2) {
        return 0;
    }
    gamma_initial = 1.0 / sqrt(1.0 - speed2 / c2);
    u = swa_vscale(velocity_initial_mps, gamma_initial);
    factor = 0.5 * charge_to_mass_C_kg * dt_s;
    u_minus = swa_vadd(u, swa_vscale(electric_V_m, factor));
    gamma_minus = sqrt(1.0 + swa_vdot(u_minus, u_minus) / c2);
    t = swa_vscale(magnetic_T, factor / gamma_minus);
    s = swa_vscale(t, 2.0 / (1.0 + swa_vdot(t, t)));
    u_prime = swa_vadd(u_minus, swa_vcross(u_minus, t));
    u_plus = swa_vadd(u_minus, swa_vcross(u_prime, s));
    u_final = swa_vadd(u_plus, swa_vscale(electric_V_m, factor));
    gamma_final = sqrt(1.0 + swa_vdot(u_final, u_final) / c2);
    *velocity_final_mps = swa_vscale(u_final, 1.0 / gamma_final);
    return finite_vec(*velocity_final_mps) &&
           swa_vdot(*velocity_final_mps, *velocity_final_mps) < c2;
}

static int deposit_longitudinal_current(
    const swa_em_pic1d *state,
    const double *position_initial_m,
    const double *position_final_m,
    const double *rho_initial_C_m3,
    const double *rho_final_C_m3,
    double *jx_face_A_m2,
    double *target_mean_A_m2
) {
    const double dx_m = state->domain_length_m / (double)state->cells;
    swa_kahan raw_mean = {0.0, 0.0};
    swa_kahan transported = {0.0, 0.0};
    size_t i;
    if (position_initial_m == NULL || position_final_m == NULL ||
        rho_initial_C_m3 == NULL || rho_final_C_m3 == NULL ||
        jx_face_A_m2 == NULL || target_mean_A_m2 == NULL) {
        return 0;
    }
    jx_face_A_m2[0] = 0.0;
    for (i = 0U; i + 1U < state->cells; ++i) {
        jx_face_A_m2[i + 1U] = jx_face_A_m2[i] -
            dx_m * (rho_final_C_m3[i] - rho_initial_C_m3[i]) /
            state->dt_s;
    }
    for (i = 0U; i < state->particle_count; ++i) {
        const double charge = state->particles[i].charge_C *
                              state->particles[i].macro_weight;
        const double displacement = position_final_m[i] -
                                    position_initial_m[i];
        swa_kahan_add(&transported, charge * displacement);
    }
    *target_mean_A_m2 = transported.sum /
        (state->domain_length_m * state->cross_section_area_m2 *
         state->dt_s);
    for (i = 0U; i < state->cells; ++i) {
        swa_kahan_add(&raw_mean, jx_face_A_m2[i]);
    }
    {
        const double shift = *target_mean_A_m2 -
                             raw_mean.sum / (double)state->cells;
        for (i = 0U; i < state->cells; ++i) {
            jx_face_A_m2[i] += shift;
        }
    }
    return finite_array(jx_face_A_m2, state->cells) &&
           isfinite(*target_mean_A_m2);
}

static int add_transverse_path_current(
    const swa_em_pic1d *state,
    double x_initial_m,
    double x_final_m,
    double current_y_scale_A_m2,
    double current_z_scale_A_m2,
    double *jy_face_A_m2,
    double *jz_face_A_m2,
    size_t *segments_out
) {
    const double length_m = state->domain_length_m;
    const double dx_m = length_m / (double)state->cells;
    const double displacement_m = x_final_m - x_initial_m;
    double remainder_m;
    double x_m;
    double turns_real;
    long long turns;
    size_t segment_count = 0U;
    size_t i;
    if (state == NULL || state->cells < 2U ||
        state->cells > (size_t)LLONG_MAX ||
        !isfinite(length_m) || !(length_m > 0.0) ||
        !isfinite(dx_m) || !(dx_m > 0.0) ||
        !isfinite(x_initial_m) || !isfinite(x_final_m) ||
        !isfinite(current_y_scale_A_m2) ||
        !isfinite(current_z_scale_A_m2) ||
        jy_face_A_m2 == NULL || jz_face_A_m2 == NULL ||
        segments_out == NULL) {
        return 0;
    }
    if (fabs(displacement_m) <=
        32.0 * DBL_EPSILON * max3(fabs(x_initial_m), length_m, 1.0)) {
        const double wrapped = wrap_periodic(x_initial_m, length_m);
        const double coordinate = wrapped / dx_m;
        const long long lower_signed = (long long)floor(coordinate);
        const size_t lower = periodic_index(lower_signed, state->cells);
        const size_t upper = (lower + 1U) % state->cells;
        const double fraction = coordinate - floor(coordinate);
        jy_face_A_m2[lower] += current_y_scale_A_m2 * (1.0 - fraction);
        jy_face_A_m2[upper] += current_y_scale_A_m2 * fraction;
        jz_face_A_m2[lower] += current_z_scale_A_m2 * (1.0 - fraction);
        jz_face_A_m2[upper] += current_z_scale_A_m2 * fraction;
        *segments_out = 1U;
        return 1;
    }

    turns_real = trunc(displacement_m / length_m);
    if (turns_real < (double)LLONG_MIN || turns_real > (double)LLONG_MAX) {
        return 0;
    }
    turns = (long long)turns_real;
    remainder_m = displacement_m - (double)turns * length_m;
    if (turns != 0) {
        const double uniform_weight = (double)turns * dx_m /
                                      displacement_m;
        for (i = 0U; i < state->cells; ++i) {
            jy_face_A_m2[i] += current_y_scale_A_m2 * uniform_weight;
            jz_face_A_m2[i] += current_z_scale_A_m2 * uniform_weight;
        }
    }

    x_m = wrap_periodic(x_initial_m, length_m);
    while (fabs(remainder_m) >
           64.0 * DBL_EPSILON * max3(fabs(x_m), length_m, 1.0)) {
        long long cell_signed;
        double boundary_m;
        double segment_m;
        double cell_origin_m;
        double t0;
        double t1;
        double left_integral_m;
        double right_integral_m;
        double left_weight;
        double right_weight;
        size_t left;
        size_t right;
        if (segment_count > state->cells + 2U) {
            return 0;
        }
        if (remainder_m > 0.0) {
            cell_signed = (long long)floor(x_m / dx_m);
            boundary_m = ((double)cell_signed + 1.0) * dx_m;
            segment_m = remainder_m < boundary_m - x_m ?
                        remainder_m : boundary_m - x_m;
        } else {
            cell_signed = (long long)ceil(x_m / dx_m) - 1LL;
            boundary_m = (double)cell_signed * dx_m;
            segment_m = remainder_m > boundary_m - x_m ?
                        remainder_m : boundary_m - x_m;
        }
        if (segment_m == 0.0 || !isfinite(segment_m)) {
            return 0;
        }
        cell_origin_m = (double)cell_signed * dx_m;
        t0 = (x_m - cell_origin_m) / dx_m;
        t1 = (x_m + segment_m - cell_origin_m) / dx_m;
        left_integral_m = dx_m *
            ((t1 - 0.5 * t1 * t1) -
             (t0 - 0.5 * t0 * t0));
        right_integral_m = 0.5 * dx_m *
            (t1 * t1 - t0 * t0);
        left_weight = left_integral_m / displacement_m;
        right_weight = right_integral_m / displacement_m;
        left = periodic_index(cell_signed, state->cells);
        right = (left + 1U) % state->cells;
        jy_face_A_m2[left] += current_y_scale_A_m2 * left_weight;
        jy_face_A_m2[right] += current_y_scale_A_m2 * right_weight;
        jz_face_A_m2[left] += current_z_scale_A_m2 * left_weight;
        jz_face_A_m2[right] += current_z_scale_A_m2 * right_weight;
        x_m += segment_m;
        remainder_m -= segment_m;
        ++segment_count;
    }
    *segments_out = segment_count + (turns != 0 ? 1U : 0U);
    return finite_array(jy_face_A_m2, state->cells) &&
           finite_array(jz_face_A_m2, state->cells);
}

