
int swa_em_pic3d_deposit_charge_density(
    const swa_em_pic3d *state,
    const swa_em_pic3d_particle *particles,
    size_t particle_count,
    double *rho_C_m3
) {
    const double volume = state != NULL ? pic3d_cell_volume(state) : NAN;
    size_t p;
    if (!pic3d_valid_shape(state) || rho_C_m3 == NULL ||
        (particle_count > 0U && particles == NULL) ||
        particle_count > state->particle_count ||
        !isfinite(volume) || !(volume > 0.0)) {
        return 0;
    }
    memset(rho_C_m3, 0, state->cell_count * sizeof(double));
    for (p = 0U; p < particle_count; ++p) {
        const double amount =
            particles[p].charge_C * particles[p].macro_weight / volume;
        if (!pic3d_particle_is_finite(&particles[p]) ||
            !isfinite(amount)) {
            return 0;
        }
        pic3d_deposit(
            state, rho_C_m3, particles[p].position_unwrapped_m,
            0.5, 0.5, 0.5, amount
        );
    }
    return pic3d_finite_array(rho_C_m3, state->cell_count);
}

static int pic3d_append_crossings(
    double start_m,
    double end_m,
    double spacing_m,
    double offset,
    double *points,
    size_t capacity,
    size_t *count
) {
    const double a = start_m / spacing_m - offset;
    const double b = end_m / spacing_m - offset;
    const double delta = b - a;
    const double low = fmin(a, b);
    const double high = fmax(a, b);
    long first;
    long last;
    long plane;
    if (points == NULL || count == NULL ||
        !isfinite(a) || !isfinite(b) ||
        !isfinite(spacing_m) || !(spacing_m > 0.0)) {
        return 0;
    }
    if (fabs(delta) <= DBL_MIN) {
        return 1;
    }
    if (low < (double)LONG_MIN + 2.0 ||
        high > (double)LONG_MAX - 2.0) {
        return 0;
    }
    first = (long)floor(low) + 1L;
    last = (long)ceil(high) - 1L;
    for (plane = first; plane <= last; ++plane) {
        const double t = ((double)plane - a) / delta;
        if (t > 64.0 * DBL_EPSILON &&
            t < 1.0 - 64.0 * DBL_EPSILON) {
            if (*count >= capacity) {
                return 0;
            }
            points[*count] = t;
            *count += 1U;
        }
    }
    return 1;
}

static int pic3d_compare_double(const void *left, const void *right) {
    const double a = *(const double *)left;
    const double b = *(const double *)right;
    return a < b ? -1 : (a > b ? 1 : 0);
}

static size_t pic3d_unique_breakpoints(
    double *points,
    size_t count
) {
    size_t read;
    size_t write;
    if (count == 0U) {
        return 0U;
    }
    qsort(points, count, sizeof(double), pic3d_compare_double);
    write = 1U;
    for (read = 1U; read < count; ++read) {
        if (fabs(points[read] - points[write - 1U]) >
            256.0 * DBL_EPSILON) {
            points[write] = points[read];
            write += 1U;
        }
    }
    return write;
}

static int pic3d_deposit_path_component(
    const swa_em_pic3d *state,
    swa_vec3 start_m,
    swa_vec3 end_m,
    double offset_x,
    double offset_y,
    double offset_z,
    double transported_charge_velocity_A_m,
    double *current_A_m2,
    size_t *segments
) {
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    const double volume = pic3d_cell_volume(state);
    const swa_vec3 displacement = swa_vsub(end_m, start_m);
    const size_t estimate_x =
        (size_t)ceil(fabs(displacement.x) / dx) + 3U;
    const size_t estimate_y =
        (size_t)ceil(fabs(displacement.y) / dy) + 3U;
    const size_t estimate_z =
        (size_t)ceil(fabs(displacement.z) / dz) + 3U;
    size_t capacity;
    double *points;
    size_t count = 0U;
    size_t interval;
    size_t used_segments;
    const double inverse_sqrt3 = 0.57735026918962576451;
    if (segments == NULL || current_A_m2 == NULL ||
        !isfinite(transported_charge_velocity_A_m) ||
        estimate_x > PIC3D_MAX_PATH_SEGMENTS ||
        estimate_y > PIC3D_MAX_PATH_SEGMENTS ||
        estimate_z > PIC3D_MAX_PATH_SEGMENTS ||
        estimate_x + estimate_y > PIC3D_MAX_PATH_SEGMENTS ||
        estimate_x + estimate_y + estimate_z + 2U >
            PIC3D_MAX_PATH_SEGMENTS ||
        !isfinite(volume) || !(volume > 0.0)) {
        return 0;
    }
    capacity = estimate_x + estimate_y + estimate_z + 2U;
    points = calloc(capacity, sizeof(double));
    if (points == NULL) {
        return 0;
    }
    points[count++] = 0.0;
    points[count++] = 1.0;
    if (!pic3d_append_crossings(
            start_m.x, end_m.x, dx, offset_x,
            points, capacity, &count
        ) ||
        !pic3d_append_crossings(
            start_m.y, end_m.y, dy, offset_y,
            points, capacity, &count
        ) ||
        !pic3d_append_crossings(
            start_m.z, end_m.z, dz, offset_z,
            points, capacity, &count
        )) {
        free(points);
        return 0;
    }
    count = pic3d_unique_breakpoints(points, count);
    if (count < 2U || count - 1U > PIC3D_MAX_PATH_SEGMENTS) {
        free(points);
        return 0;
    }
    used_segments = count - 1U;
    for (interval = 0U; interval + 1U < count; ++interval) {
        const double left = points[interval];
        const double right = points[interval + 1U];
        const double midpoint = 0.5 * (left + right);
        const double half = 0.5 * (right - left);
        const double nodes[2] = {
            midpoint - half * inverse_sqrt3,
            midpoint + half * inverse_sqrt3
        };
        size_t node;
        for (node = 0U; node < 2U; ++node) {
            const swa_vec3 position = swa_vadd(
                start_m, swa_vscale(displacement, nodes[node])
            );
            const double amount =
                transported_charge_velocity_A_m * half / volume;
            pic3d_deposit(
                state, current_A_m2, position,
                offset_x, offset_y, offset_z, amount
            );
        }
    }
    *segments = used_segments;
    free(points);
    return pic3d_finite_array(current_A_m2, state->cell_count);
}

int swa_em_pic3d_deposit_transport_current(
    const swa_em_pic3d *state,
    const swa_em_pic3d_particle *initial,
    const swa_em_pic3d_particle *final,
    size_t particle_count,
    double *jx,
    double *jy,
    double *jz,
    size_t *maximum_segments,
    swa_vec3 *transport_mean_current
) {
    const double domain_volume =
        state->length_x_m * state->length_y_m * state->length_z_m;
    swa_kahan transport_x = {0.0, 0.0};
    swa_kahan transport_y = {0.0, 0.0};
    swa_kahan transport_z = {0.0, 0.0};
    size_t p;
    if (!pic3d_valid_shape(state) ||
        initial == NULL || final == NULL ||
        particle_count > state->particle_count ||
        jx == NULL || jy == NULL || jz == NULL ||
        maximum_segments == NULL || transport_mean_current == NULL ||
        !isfinite(domain_volume) || !(domain_volume > 0.0)) {
        return 0;
    }
    memset(jx, 0, state->cell_count * sizeof(double));
    memset(jy, 0, state->cell_count * sizeof(double));
    memset(jz, 0, state->cell_count * sizeof(double));
    *maximum_segments = 0U;
    for (p = 0U; p < particle_count; ++p) {
        const swa_vec3 displacement = swa_vsub(
            final[p].position_unwrapped_m,
            initial[p].position_unwrapped_m
        );
        const double transported_charge =
            initial[p].charge_C * initial[p].macro_weight;
        const double factor = transported_charge / state->dt_s;
        size_t segments_x = 0U;
        size_t segments_y = 0U;
        size_t segments_z = 0U;
        if (!pic3d_particle_is_finite(&initial[p]) ||
            !pic3d_particle_is_finite(&final[p]) ||
            initial[p].charge_C != final[p].charge_C ||
            initial[p].mass_kg != final[p].mass_kg ||
            initial[p].macro_weight != final[p].macro_weight ||
            !pic3d_deposit_path_component(
                state,
                initial[p].position_unwrapped_m,
                final[p].position_unwrapped_m,
                0.0, 0.5, 0.5,
                factor * displacement.x, jx, &segments_x
            ) ||
            !pic3d_deposit_path_component(
                state,
                initial[p].position_unwrapped_m,
                final[p].position_unwrapped_m,
                0.5, 0.0, 0.5,
                factor * displacement.y, jy, &segments_y
            ) ||
            !pic3d_deposit_path_component(
                state,
                initial[p].position_unwrapped_m,
                final[p].position_unwrapped_m,
                0.5, 0.5, 0.0,
                factor * displacement.z, jz, &segments_z
            )) {
            return 0;
        }
        {
            size_t local_max = segments_x;
            if (segments_y > local_max) {
                local_max = segments_y;
            }
            if (segments_z > local_max) {
                local_max = segments_z;
            }
            if (local_max > *maximum_segments) {
                *maximum_segments = local_max;
            }
        }
        swa_kahan_add(
            &transport_x,
            transported_charge * displacement.x /
                (state->dt_s * domain_volume)
        );
        swa_kahan_add(
            &transport_y,
            transported_charge * displacement.y /
                (state->dt_s * domain_volume)
        );
        swa_kahan_add(
            &transport_z,
            transported_charge * displacement.z /
                (state->dt_s * domain_volume)
        );
    }
    *transport_mean_current = swa_v3(
        transport_x.sum, transport_y.sum, transport_z.sum
    );
    return pic3d_finite_array(jx, state->cell_count) &&
           pic3d_finite_array(jy, state->cell_count) &&
           pic3d_finite_array(jz, state->cell_count) &&
           pic3d_finite_vec3(*transport_mean_current);
}

static int pic3d_axis_move_x(
    const swa_em_pic3d *state,
    swa_vec3 start,
    swa_vec3 end,
    double transported_charge_C,
    double *jx
) {
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    const double volume = dx * dy * dz;
    size_t x0[2];
    size_t x1[2];
    size_t yi[2];
    size_t zi[2];
    double wx0[2];
    double wx1[2];
    double wy[2];
    double wz[2];
    double *delta;
    double *line;
    size_t a;
    size_t b;
    if (jx == NULL || !isfinite(transported_charge_C)) {
        return 0;
    }
    delta = calloc(2U * state->nx, sizeof(double));
    if (delta == NULL) {
        return 0;
    }
    line = delta + state->nx;
    pic3d_axis_weights(start.x, state->length_x_m, state->nx, 0.5,
                       &x0[0], &x0[1], &wx0[0], &wx0[1]);
    pic3d_axis_weights(end.x, state->length_x_m, state->nx, 0.5,
                       &x1[0], &x1[1], &wx1[0], &wx1[1]);
    pic3d_axis_weights(start.y, state->length_y_m, state->ny, 0.5,
                       &yi[0], &yi[1], &wy[0], &wy[1]);
    pic3d_axis_weights(start.z, state->length_z_m, state->nz, 0.5,
                       &zi[0], &zi[1], &wz[0], &wz[1]);
    for (b = 0U; b < 2U; ++b) {
        for (a = 0U; a < 2U; ++a) {
            const double transverse_weight = wy[a] * wz[b];
            const double density_scale =
                transported_charge_C * transverse_weight / volume;
            const double target_mean = transported_charge_C *
                transverse_weight * (end.x - start.x) /
                (state->dt_s * state->length_x_m * dy * dz);
            swa_kahan line_sum = {0.0, 0.0};
            double constant;
            size_t i;
            memset(delta, 0, state->nx * sizeof(double));
            memset(line, 0, state->nx * sizeof(double));
            delta[x0[0]] -= density_scale * wx0[0];
            delta[x0[1]] -= density_scale * wx0[1];
            delta[x1[0]] += density_scale * wx1[0];
            delta[x1[1]] += density_scale * wx1[1];
            for (i = 0U; i + 1U < state->nx; ++i) {
                line[i + 1U] = line[i] -
                    dx * delta[i] / state->dt_s;
            }
            for (i = 0U; i < state->nx; ++i) {
                swa_kahan_add(&line_sum, line[i]);
            }
            constant = target_mean -
                line_sum.sum / (double)state->nx;
            for (i = 0U; i < state->nx; ++i) {
                jx[pic3d_index(state, i, yi[a], zi[b])] +=
                    line[i] + constant;
            }
        }
    }
    free(delta);
    return 1;
}

static int pic3d_axis_move_y(
    const swa_em_pic3d *state,
    swa_vec3 start,
    swa_vec3 end,
    double transported_charge_C,
    double *jy
) {
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    const double volume = dx * dy * dz;
    size_t xi[2];
    size_t y0[2];
    size_t y1[2];
    size_t zi[2];
    double wx[2];
    double wy0[2];
    double wy1[2];
    double wz[2];
    double *delta;
    double *line;
    size_t a;
    size_t b;
    if (jy == NULL || !isfinite(transported_charge_C)) {
        return 0;
    }
    delta = calloc(2U * state->ny, sizeof(double));
    if (delta == NULL) {
        return 0;
    }
    line = delta + state->ny;
    pic3d_axis_weights(start.x, state->length_x_m, state->nx, 0.5,
                       &xi[0], &xi[1], &wx[0], &wx[1]);
    pic3d_axis_weights(start.y, state->length_y_m, state->ny, 0.5,
                       &y0[0], &y0[1], &wy0[0], &wy0[1]);
    pic3d_axis_weights(end.y, state->length_y_m, state->ny, 0.5,
                       &y1[0], &y1[1], &wy1[0], &wy1[1]);
    pic3d_axis_weights(start.z, state->length_z_m, state->nz, 0.5,
                       &zi[0], &zi[1], &wz[0], &wz[1]);
    for (b = 0U; b < 2U; ++b) {
        for (a = 0U; a < 2U; ++a) {
            const double transverse_weight = wx[a] * wz[b];
            const double density_scale =
                transported_charge_C * transverse_weight / volume;
            const double target_mean = transported_charge_C *
                transverse_weight * (end.y - start.y) /
                (state->dt_s * state->length_y_m * dx * dz);
            swa_kahan line_sum = {0.0, 0.0};
            double constant;
            size_t j;
            memset(delta, 0, state->ny * sizeof(double));
            memset(line, 0, state->ny * sizeof(double));
            delta[y0[0]] -= density_scale * wy0[0];
            delta[y0[1]] -= density_scale * wy0[1];
            delta[y1[0]] += density_scale * wy1[0];
            delta[y1[1]] += density_scale * wy1[1];
            for (j = 0U; j + 1U < state->ny; ++j) {
                line[j + 1U] = line[j] -
                    dy * delta[j] / state->dt_s;
            }
            for (j = 0U; j < state->ny; ++j) {
                swa_kahan_add(&line_sum, line[j]);
            }
            constant = target_mean -
                line_sum.sum / (double)state->ny;
            for (j = 0U; j < state->ny; ++j) {
                jy[pic3d_index(state, xi[a], j, zi[b])] +=
                    line[j] + constant;
            }
        }
    }
    free(delta);
    return 1;
}

static int pic3d_axis_move_z(
    const swa_em_pic3d *state,
    swa_vec3 start,
    swa_vec3 end,
    double transported_charge_C,
    double *jz
) {
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    const double volume = dx * dy * dz;
    size_t xi[2];
    size_t yi[2];
    size_t z0[2];
    size_t z1[2];
    double wx[2];
    double wy[2];
    double wz0[2];
    double wz1[2];
    double *delta;
    double *line;
    size_t a;
    size_t b;
    if (jz == NULL || !isfinite(transported_charge_C)) {
        return 0;
    }
    delta = calloc(2U * state->nz, sizeof(double));
    if (delta == NULL) {
        return 0;
    }
    line = delta + state->nz;
    pic3d_axis_weights(start.x, state->length_x_m, state->nx, 0.5,
                       &xi[0], &xi[1], &wx[0], &wx[1]);
    pic3d_axis_weights(start.y, state->length_y_m, state->ny, 0.5,
                       &yi[0], &yi[1], &wy[0], &wy[1]);
    pic3d_axis_weights(start.z, state->length_z_m, state->nz, 0.5,
                       &z0[0], &z0[1], &wz0[0], &wz0[1]);
    pic3d_axis_weights(end.z, state->length_z_m, state->nz, 0.5,
                       &z1[0], &z1[1], &wz1[0], &wz1[1]);
    for (b = 0U; b < 2U; ++b) {
        for (a = 0U; a < 2U; ++a) {
            const double transverse_weight = wx[a] * wy[b];
            const double density_scale =
                transported_charge_C * transverse_weight / volume;
            const double target_mean = transported_charge_C *
                transverse_weight * (end.z - start.z) /
                (state->dt_s * state->length_z_m * dx * dy);
            swa_kahan line_sum = {0.0, 0.0};
            double constant;
            size_t k;
            memset(delta, 0, state->nz * sizeof(double));
            memset(line, 0, state->nz * sizeof(double));
            delta[z0[0]] -= density_scale * wz0[0];
            delta[z0[1]] -= density_scale * wz0[1];
            delta[z1[0]] += density_scale * wz1[0];
            delta[z1[1]] += density_scale * wz1[1];
            for (k = 0U; k + 1U < state->nz; ++k) {
                line[k + 1U] = line[k] -
                    dz * delta[k] / state->dt_s;
            }
            for (k = 0U; k < state->nz; ++k) {
                swa_kahan_add(&line_sum, line[k]);
            }
            constant = target_mean -
                line_sum.sum / (double)state->nz;
            for (k = 0U; k < state->nz; ++k) {
                jz[pic3d_index(state, xi[a], yi[b], k)] +=
                    line[k] + constant;
            }
        }
    }
    free(delta);
    return 1;
}

int swa_em_pic3d_deposit_charge_conserving_current(
    const swa_em_pic3d *state,
    const swa_em_pic3d_particle *initial_particles,
    const swa_em_pic3d_particle *final_particles,
    size_t particle_count,
    double *jx_A_m2,
    double *jy_A_m2,
    double *jz_A_m2,
    swa_vec3 *transport_mean_current_A_m2
) {
    static const unsigned permutations[6][3] = {
        {0U, 1U, 2U}, {0U, 2U, 1U},
        {1U, 0U, 2U}, {1U, 2U, 0U},
        {2U, 0U, 1U}, {2U, 1U, 0U}
    };
    const double domain_volume =
        state != NULL ?
        state->length_x_m * state->length_y_m * state->length_z_m : NAN;
    swa_kahan mean_x = {0.0, 0.0};
    swa_kahan mean_y = {0.0, 0.0};
    swa_kahan mean_z = {0.0, 0.0};
    size_t p;
    if (!pic3d_valid_shape(state) ||
        initial_particles == NULL || final_particles == NULL ||
        particle_count > state->particle_count ||
        jx_A_m2 == NULL || jy_A_m2 == NULL || jz_A_m2 == NULL ||
        transport_mean_current_A_m2 == NULL ||
        !isfinite(domain_volume) || !(domain_volume > 0.0)) {
        return 0;
    }
    memset(jx_A_m2, 0, state->cell_count * sizeof(double));
    memset(jy_A_m2, 0, state->cell_count * sizeof(double));
    memset(jz_A_m2, 0, state->cell_count * sizeof(double));
    for (p = 0U; p < particle_count; ++p) {
        const double transported_charge =
            initial_particles[p].charge_C *
            initial_particles[p].macro_weight;
        const swa_vec3 displacement = swa_vsub(
            final_particles[p].position_unwrapped_m,
            initial_particles[p].position_unwrapped_m
        );
        size_t permutation;
        if (!pic3d_particle_is_finite(&initial_particles[p]) ||
            !pic3d_particle_is_finite(&final_particles[p]) ||
            initial_particles[p].charge_C != final_particles[p].charge_C ||
            initial_particles[p].mass_kg != final_particles[p].mass_kg ||
            initial_particles[p].macro_weight !=
                final_particles[p].macro_weight) {
            return 0;
        }
        for (permutation = 0U; permutation < 6U; ++permutation) {
            swa_vec3 point = initial_particles[p].position_unwrapped_m;
            size_t leg;
            for (leg = 0U; leg < 3U; ++leg) {
                const unsigned axis = permutations[permutation][leg];
                swa_vec3 next = point;
                if (axis == 0U) {
                    next.x = final_particles[p].position_unwrapped_m.x;
                    if (!pic3d_axis_move_x(
                            state, point, next,
                            transported_charge / 6.0, jx_A_m2)) {
                        return 0;
                    }
                } else if (axis == 1U) {
                    next.y = final_particles[p].position_unwrapped_m.y;
                    if (!pic3d_axis_move_y(
                            state, point, next,
                            transported_charge / 6.0, jy_A_m2)) {
                        return 0;
                    }
                } else {
                    next.z = final_particles[p].position_unwrapped_m.z;
                    if (!pic3d_axis_move_z(
                            state, point, next,
                            transported_charge / 6.0, jz_A_m2)) {
                        return 0;
                    }
                }
                point = next;
            }
        }
        swa_kahan_add(
            &mean_x,
            transported_charge * displacement.x /
                (state->dt_s * domain_volume)
        );
        swa_kahan_add(
            &mean_y,
            transported_charge * displacement.y /
                (state->dt_s * domain_volume)
        );
        swa_kahan_add(
            &mean_z,
            transported_charge * displacement.z /
                (state->dt_s * domain_volume)
        );
    }
    *transport_mean_current_A_m2 = swa_v3(
        mean_x.sum, mean_y.sum, mean_z.sum
    );
    return pic3d_finite_array(jx_A_m2, state->cell_count) &&
           pic3d_finite_array(jy_A_m2, state->cell_count) &&
           pic3d_finite_array(jz_A_m2, state->cell_count) &&
           pic3d_finite_vec3(*transport_mean_current_A_m2);
}
