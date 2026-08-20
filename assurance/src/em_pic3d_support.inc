
static double pic3d_max2(double a, double b) {
    return a > b ? a : b;
}

static double pic3d_max3(double a, double b, double c) {
    return pic3d_max2(pic3d_max2(a, b), c);
}

static double pic3d_max4(double a, double b, double c, double d) {
    return pic3d_max2(pic3d_max2(a, b), pic3d_max2(c, d));
}

static int pic3d_finite_vec3(swa_vec3 value) {
    return isfinite(value.x) && isfinite(value.y) && isfinite(value.z);
}

static int pic3d_finite_array(const double *values, size_t count) {
    size_t i;
    if (values == NULL) {
        return 0;
    }
    for (i = 0U; i < count; ++i) {
        if (!isfinite(values[i])) {
            return 0;
        }
    }
    return 1;
}

static int pic3d_checked_cell_count(
    size_t nx,
    size_t ny,
    size_t nz,
    size_t *out
) {
    size_t xy;
    if (out == NULL || nx == 0U || ny == 0U || nz == 0U ||
        nx > (size_t)LONG_MAX ||
        ny > (size_t)LONG_MAX ||
        nz > (size_t)LONG_MAX ||
        nx > SIZE_MAX / ny) {
        return 0;
    }
    xy = nx * ny;
    if (xy > SIZE_MAX / nz) {
        return 0;
    }
    *out = xy * nz;
    return *out <= SIZE_MAX / sizeof(double);
}

static size_t pic3d_index(
    const swa_em_pic3d *state,
    size_t i,
    size_t j,
    size_t k
) {
    return (k * state->ny + j) * state->nx + i;
}

static size_t pic3d_wrap_index(long index, size_t count) {
    const long modulus = (long)count;
    long wrapped = index % modulus;
    if (wrapped < 0L) {
        wrapped += modulus;
    }
    return (size_t)wrapped;
}

static double pic3d_wrap_coordinate(double value, double length) {
    double wrapped = fmod(value, length);
    if (wrapped < 0.0) {
        wrapped += length;
    }
    if (wrapped >= length) {
        wrapped = 0.0;
    }
    return wrapped;
}

static double pic3d_dx(const swa_em_pic3d *state) {
    return state->length_x_m / (double)state->nx;
}

static double pic3d_dy(const swa_em_pic3d *state) {
    return state->length_y_m / (double)state->ny;
}

static double pic3d_dz(const swa_em_pic3d *state) {
    return state->length_z_m / (double)state->nz;
}

static double pic3d_cell_volume(const swa_em_pic3d *state) {
    return pic3d_dx(state) * pic3d_dy(state) * pic3d_dz(state);
}

static int pic3d_valid_shape(const swa_em_pic3d *state) {
    size_t expected;
    return state != NULL &&
           pic3d_checked_cell_count(
               state->nx, state->ny, state->nz, &expected
           ) &&
           expected == state->cell_count &&
           isfinite(state->length_x_m) && state->length_x_m > 0.0 &&
           isfinite(state->length_y_m) && state->length_y_m > 0.0 &&
           isfinite(state->length_z_m) && state->length_z_m > 0.0 &&
           isfinite(state->dt_s) && state->dt_s > 0.0 &&
           state->ex_V_m != NULL &&
           state->ey_V_m != NULL &&
           state->ez_V_m != NULL &&
           state->bx_T != NULL &&
           state->by_T != NULL &&
           state->bz_T != NULL &&
           (state->particle_count == 0U || state->particles != NULL);
}

static int pic3d_particle_is_finite(
    const swa_em_pic3d_particle *particle
) {
    const double speed2 = particle != NULL ?
        swa_vdot(particle->velocity_mps, particle->velocity_mps) : NAN;
    return particle != NULL &&
           pic3d_finite_vec3(particle->position_unwrapped_m) &&
           pic3d_finite_vec3(particle->velocity_mps) &&
           isfinite(speed2) && speed2 >= 0.0 &&
           speed2 < SWA_C * SWA_C &&
           isfinite(particle->charge_C) &&
           isfinite(particle->mass_kg) && particle->mass_kg > 0.0 &&
           isfinite(particle->macro_weight) &&
           particle->macro_weight > 0.0;
}

static void pic3d_axis_weights(
    double position_m,
    double length_m,
    size_t count,
    double offset,
    size_t *index0,
    size_t *index1,
    double *weight0,
    double *weight1
) {
    const double spacing = length_m / (double)count;
    const double wrapped = pic3d_wrap_coordinate(position_m, length_m);
    const double coordinate = wrapped / spacing - offset;
    const double lower = floor(coordinate);
    const double fraction = coordinate - lower;
    if (count == 1U) {
        *index0 = 0U;
        *index1 = 0U;
        *weight0 = 1.0;
        *weight1 = 0.0;
        return;
    }
    *index0 = pic3d_wrap_index((long)lower, count);
    *index1 = (*index0 + 1U) % count;
    *weight0 = 1.0 - fraction;
    *weight1 = fraction;
}

static double pic3d_interpolate(
    const swa_em_pic3d *state,
    const double *values,
    swa_vec3 position_m,
    double offset_x,
    double offset_y,
    double offset_z
) {
    size_t ix[2];
    size_t iy[2];
    size_t iz[2];
    double wx[2];
    double wy[2];
    double wz[2];
    double value = 0.0;
    size_t a;
    size_t b;
    size_t c;
    pic3d_axis_weights(
        position_m.x, state->length_x_m, state->nx, offset_x,
        &ix[0], &ix[1], &wx[0], &wx[1]
    );
    pic3d_axis_weights(
        position_m.y, state->length_y_m, state->ny, offset_y,
        &iy[0], &iy[1], &wy[0], &wy[1]
    );
    pic3d_axis_weights(
        position_m.z, state->length_z_m, state->nz, offset_z,
        &iz[0], &iz[1], &wz[0], &wz[1]
    );
    for (c = 0U; c < 2U; ++c) {
        for (b = 0U; b < 2U; ++b) {
            for (a = 0U; a < 2U; ++a) {
                value += wx[a] * wy[b] * wz[c] *
                    values[pic3d_index(state, ix[a], iy[b], iz[c])];
            }
        }
    }
    return value;
}

static void pic3d_deposit(
    const swa_em_pic3d *state,
    double *values,
    swa_vec3 position_m,
    double offset_x,
    double offset_y,
    double offset_z,
    double amount
) {
    size_t ix[2];
    size_t iy[2];
    size_t iz[2];
    double wx[2];
    double wy[2];
    double wz[2];
    size_t a;
    size_t b;
    size_t c;
    pic3d_axis_weights(
        position_m.x, state->length_x_m, state->nx, offset_x,
        &ix[0], &ix[1], &wx[0], &wx[1]
    );
    pic3d_axis_weights(
        position_m.y, state->length_y_m, state->ny, offset_y,
        &iy[0], &iy[1], &wy[0], &wy[1]
    );
    pic3d_axis_weights(
        position_m.z, state->length_z_m, state->nz, offset_z,
        &iz[0], &iz[1], &wz[0], &wz[1]
    );
    for (c = 0U; c < 2U; ++c) {
        for (b = 0U; b < 2U; ++b) {
            for (a = 0U; a < 2U; ++a) {
                values[pic3d_index(state, ix[a], iy[b], iz[c])] +=
                    amount * wx[a] * wy[b] * wz[c];
            }
        }
    }
}

static swa_vec3 pic3d_interpolate_e(
    const swa_em_pic3d *state,
    const double *ex,
    const double *ey,
    const double *ez,
    swa_vec3 position_m
) {
    return swa_v3(
        pic3d_interpolate(state, ex, position_m, 0.0, 0.5, 0.5),
        pic3d_interpolate(state, ey, position_m, 0.5, 0.0, 0.5),
        pic3d_interpolate(state, ez, position_m, 0.5, 0.5, 0.0)
    );
}

static swa_vec3 pic3d_interpolate_b(
    const swa_em_pic3d *state,
    const double *bx,
    const double *by,
    const double *bz,
    swa_vec3 position_m
) {
    return swa_v3(
        pic3d_interpolate(state, bx, position_m, 0.5, 0.0, 0.0),
        pic3d_interpolate(state, by, position_m, 0.0, 0.5, 0.0),
        pic3d_interpolate(state, bz, position_m, 0.0, 0.0, 0.5)
    );
}

int swa_em_pic3d_init(
    swa_em_pic3d *state,
    size_t nx,
    size_t ny,
    size_t nz,
    size_t particle_count,
    double length_x_m,
    double length_y_m,
    double length_z_m,
    double dt_s
) {
    size_t cells;
    if (state == NULL ||
        !pic3d_checked_cell_count(nx, ny, nz, &cells) ||
        !isfinite(length_x_m) || !(length_x_m > 0.0) ||
        !isfinite(length_y_m) || !(length_y_m > 0.0) ||
        !isfinite(length_z_m) || !(length_z_m > 0.0) ||
        !isfinite(dt_s) || !(dt_s > 0.0) ||
        particle_count > SIZE_MAX / sizeof(swa_em_pic3d_particle)) {
        return 0;
    }
    memset(state, 0, sizeof(*state));
    state->nx = nx;
    state->ny = ny;
    state->nz = nz;
    state->cell_count = cells;
    state->particle_count = particle_count;
    state->length_x_m = length_x_m;
    state->length_y_m = length_y_m;
    state->length_z_m = length_z_m;
    state->dt_s = dt_s;
    if (particle_count > 0U) {
        state->particles = calloc(
            particle_count, sizeof(*state->particles)
        );
    }
    state->ex_V_m = calloc(cells, sizeof(double));
    state->ey_V_m = calloc(cells, sizeof(double));
    state->ez_V_m = calloc(cells, sizeof(double));
    state->bx_T = calloc(cells, sizeof(double));
    state->by_T = calloc(cells, sizeof(double));
    state->bz_T = calloc(cells, sizeof(double));
    if ((particle_count > 0U && state->particles == NULL) ||
        state->ex_V_m == NULL || state->ey_V_m == NULL ||
        state->ez_V_m == NULL || state->bx_T == NULL ||
        state->by_T == NULL || state->bz_T == NULL) {
        swa_em_pic3d_destroy(state);
        return 0;
    }
    return 1;
}

void swa_em_pic3d_destroy(swa_em_pic3d *state) {
    if (state == NULL) {
        return;
    }
    free(state->particles);
    free(state->ex_V_m);
    free(state->ey_V_m);
    free(state->ez_V_m);
    free(state->bx_T);
    free(state->by_T);
    free(state->bz_T);
    memset(state, 0, sizeof(*state));
}

int swa_em_pic3d_state_is_finite(const swa_em_pic3d *state) {
    size_t i;
    if (!pic3d_valid_shape(state) ||
        !pic3d_finite_array(state->ex_V_m, state->cell_count) ||
        !pic3d_finite_array(state->ey_V_m, state->cell_count) ||
        !pic3d_finite_array(state->ez_V_m, state->cell_count) ||
        !pic3d_finite_array(state->bx_T, state->cell_count) ||
        !pic3d_finite_array(state->by_T, state->cell_count) ||
        !pic3d_finite_array(state->bz_T, state->cell_count)) {
        return 0;
    }
    for (i = 0U; i < state->particle_count; ++i) {
        if (!pic3d_particle_is_finite(&state->particles[i])) {
            return 0;
        }
    }
    return 1;
}

int swa_em_pic3d_set_particle(
    swa_em_pic3d *state,
    size_t index,
    swa_vec3 position_unwrapped_m,
    swa_vec3 velocity_mps,
    double charge_C,
    double mass_kg,
    double macro_weight
) {
    swa_em_pic3d_particle particle;
    if (!pic3d_valid_shape(state) || index >= state->particle_count) {
        return 0;
    }
    particle.position_unwrapped_m = position_unwrapped_m;
    particle.velocity_mps = velocity_mps;
    particle.charge_C = charge_C;
    particle.mass_kg = mass_kg;
    particle.macro_weight = macro_weight;
    if (!pic3d_particle_is_finite(&particle)) {
        return 0;
    }
    state->particles[index] = particle;
    return 1;
}

int swa_em_pic3d_higuera_cary_push(
    swa_vec3 velocity_initial_mps,
    double charge_C,
    double mass_kg,
    swa_vec3 electric_field_V_m,
    swa_vec3 magnetic_field_T,
    double dt_s,
    swa_vec3 *velocity_final_mps
) {
    const double c2 = SWA_C * SWA_C;
    const double speed2 = swa_vdot(
        velocity_initial_mps, velocity_initial_mps
    );
    double gamma_initial;
    double qmdt2;
    swa_vec3 u_initial;
    swa_vec3 u_minus;
    swa_vec3 tau;
    double gamma_minus2;
    double tau2;
    double u_star;
    double sigma;
    double discriminant;
    double gamma_rotation;
    swa_vec3 t;
    swa_vec3 s;
    swa_vec3 u_prime;
    swa_vec3 u_plus;
    swa_vec3 u_final;
    double gamma_final;
    double final_speed2;
    if (velocity_final_mps == NULL ||
        !pic3d_finite_vec3(velocity_initial_mps) ||
        !isfinite(speed2) || speed2 < 0.0 || speed2 >= c2 ||
        !isfinite(charge_C) ||
        !isfinite(mass_kg) || !(mass_kg > 0.0) ||
        !pic3d_finite_vec3(electric_field_V_m) ||
        !pic3d_finite_vec3(magnetic_field_T) ||
        !isfinite(dt_s)) {
        return 0;
    }
    gamma_initial = 1.0 / sqrt(1.0 - speed2 / c2);
    u_initial = swa_vscale(velocity_initial_mps, gamma_initial);
    qmdt2 = charge_C * dt_s / (2.0 * mass_kg);
    u_minus = swa_vadd(
        u_initial, swa_vscale(electric_field_V_m, qmdt2)
    );
    tau = swa_vscale(magnetic_field_T, qmdt2);
    gamma_minus2 = 1.0 + swa_vdot(u_minus, u_minus) / c2;
    tau2 = swa_vdot(tau, tau);
    u_star = swa_vdot(u_minus, tau) / SWA_C;
    sigma = gamma_minus2 - tau2;
    discriminant = sigma * sigma + 4.0 * (tau2 + u_star * u_star);
    if (!isfinite(discriminant) || discriminant < 0.0) {
        return 0;
    }
    gamma_rotation = sqrt(
        0.5 * (sigma + sqrt(discriminant))
    );
    if (!isfinite(gamma_rotation) || !(gamma_rotation > 0.0)) {
        return 0;
    }
    t = swa_vscale(tau, 1.0 / gamma_rotation);
    s = swa_vscale(t, 2.0 / (1.0 + swa_vdot(t, t)));
    u_prime = swa_vadd(u_minus, swa_vcross(u_minus, t));
    u_plus = swa_vadd(u_minus, swa_vcross(u_prime, s));
    u_final = swa_vadd(
        u_plus, swa_vscale(electric_field_V_m, qmdt2)
    );
    gamma_final = sqrt(1.0 + swa_vdot(u_final, u_final) / c2);
    if (!isfinite(gamma_final) || !(gamma_final >= 1.0)) {
        return 0;
    }
    *velocity_final_mps = swa_vscale(u_final, 1.0 / gamma_final);
    final_speed2 = swa_vdot(*velocity_final_mps, *velocity_final_mps);
    return pic3d_finite_vec3(*velocity_final_mps) &&
           isfinite(final_speed2) && final_speed2 >= 0.0 &&
           final_speed2 < c2;
}
