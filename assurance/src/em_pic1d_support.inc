static double max2(double a, double b) { return a > b ? a : b; }
static double max3(double a, double b, double c) { return max2(max2(a, b), c); }
static size_t count_or_one(size_t count) { return count == 0U ? 1U : count; }

static int finite_vec(swa_vec3 value) {
    return isfinite(value.x) && isfinite(value.y) && isfinite(value.z);
}

static double wrap_periodic(double position_m, double length_m) {
    double wrapped = fmod(position_m, length_m);
    if (wrapped < 0.0) {
        wrapped += length_m;
    }
    if (wrapped >= length_m) {
        wrapped = 0.0;
    }
    return wrapped;
}

static size_t periodic_index(long long index, size_t count) {
    long long modulus;
    long long reduced;
    if (count == 0U || count > (size_t)LLONG_MAX) {
        return 0U;
    }
    modulus = (long long)count;
    reduced = index % modulus;
    if (reduced < 0) {
        reduced += modulus;
    }
    return (size_t)reduced;
}

static int finite_array(const double *values, size_t count) {
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

static int state_shape_valid(const swa_em_pic1d *state) {
    return state != NULL && state->cells >= 4U &&
           state->cells <= (size_t)LLONG_MAX &&
           isfinite(state->domain_length_m) &&
           state->domain_length_m > 0.0 &&
           isfinite(state->cross_section_area_m2) &&
           state->cross_section_area_m2 > 0.0 &&
           isfinite(state->dt_s) && state->dt_s > 0.0 &&
           isfinite(state->bx_uniform_T) &&
           state->particles != NULL &&
           state->ex_face_V_m != NULL &&
           state->ey_face_V_m != NULL &&
           state->ez_face_V_m != NULL &&
           state->by_cell_T != NULL &&
           state->bz_cell_T != NULL;
}

static int state_values_valid(const swa_em_pic1d *state) {
    size_t i;
    if (!state_shape_valid(state) ||
        !finite_array(state->ex_face_V_m, state->cells) ||
        !finite_array(state->ey_face_V_m, state->cells) ||
        !finite_array(state->ez_face_V_m, state->cells) ||
        !finite_array(state->by_cell_T, state->cells) ||
        !finite_array(state->bz_cell_T, state->cells)) {
        return 0;
    }
    for (i = 0U; i < state->particle_count; ++i) {
        const swa_em_pic1d_particle *particle = &state->particles[i];
        const double speed2 = swa_vdot(particle->velocity_mps,
                                       particle->velocity_mps);
        if (!isfinite(particle->position_unwrapped_m) ||
            !finite_vec(particle->velocity_mps) ||
            !isfinite(particle->charge_C) ||
            !isfinite(particle->mass_kg) ||
            !isfinite(particle->macro_weight) ||
            !(particle->mass_kg > 0.0) ||
            !(particle->macro_weight > 0.0) ||
            !(speed2 >= 0.0) || speed2 >= SWA_C * SWA_C) {
            return 0;
        }
    }
    return 1;
}

static int limits_valid(const swa_em_pic1d_limits *limits) {
    return limits != NULL &&
           isfinite(limits->continuity_relative_tolerance) &&
           isfinite(limits->gauss_relative_tolerance) &&
           isfinite(limits->magnetic_divergence_relative_tolerance) &&
           isfinite(limits->energy_relative_tolerance) &&
           isfinite(limits->momentum_relative_tolerance) &&
           isfinite(limits->maximum_courant) &&
           isfinite(limits->maximum_particle_cells_per_step) &&
           limits->continuity_relative_tolerance >= 0.0 &&
           limits->gauss_relative_tolerance >= 0.0 &&
           limits->magnetic_divergence_relative_tolerance >= 0.0 &&
           limits->energy_relative_tolerance >= 0.0 &&
           limits->momentum_relative_tolerance >= 0.0 &&
           limits->maximum_courant > 0.0 &&
           limits->maximum_courant <= 1.0 &&
           limits->maximum_particle_cells_per_step > 0.0;
}

int swa_em_pic1d_init(
    swa_em_pic1d *state,
    size_t cells,
    size_t particle_count,
    double domain_length_m,
    double cross_section_area_m2,
    double dt_s
) {
    if (state == NULL || cells < 4U || cells > (size_t)LLONG_MAX ||
        !isfinite(domain_length_m) || !(domain_length_m > 0.0) ||
        !isfinite(cross_section_area_m2) ||
        !(cross_section_area_m2 > 0.0) ||
        !isfinite(dt_s) || !(dt_s > 0.0) ||
        cells > SIZE_MAX / sizeof(double) ||
        particle_count > SIZE_MAX / sizeof(swa_em_pic1d_particle)) {
        return 0;
    }
    memset(state, 0, sizeof(*state));
    state->cells = cells;
    state->particle_count = particle_count;
    state->domain_length_m = domain_length_m;
    state->cross_section_area_m2 = cross_section_area_m2;
    state->dt_s = dt_s;
    state->particles = calloc(count_or_one(particle_count), sizeof(*state->particles));
    state->ex_face_V_m = calloc(cells, sizeof(double));
    state->ey_face_V_m = calloc(cells, sizeof(double));
    state->ez_face_V_m = calloc(cells, sizeof(double));
    state->by_cell_T = calloc(cells, sizeof(double));
    state->bz_cell_T = calloc(cells, sizeof(double));
    if (state->particles == NULL || state->ex_face_V_m == NULL ||
        state->ey_face_V_m == NULL || state->ez_face_V_m == NULL ||
        state->by_cell_T == NULL || state->bz_cell_T == NULL) {
        swa_em_pic1d_destroy(state);
        return 0;
    }
    return 1;
}

void swa_em_pic1d_destroy(swa_em_pic1d *state) {
    if (state == NULL) {
        return;
    }
    free(state->particles);
    free(state->ex_face_V_m);
    free(state->ey_face_V_m);
    free(state->ez_face_V_m);
    free(state->by_cell_T);
    free(state->bz_cell_T);
    memset(state, 0, sizeof(*state));
}

int swa_em_pic1d_set_particle(
    swa_em_pic1d *state,
    size_t index,
    double position_unwrapped_m,
    swa_vec3 velocity_mps,
    double charge_C,
    double mass_kg,
    double macro_weight
) {
    const double speed2 = swa_vdot(velocity_mps, velocity_mps);
    if (state == NULL || state->particles == NULL ||
        index >= state->particle_count ||
        !isfinite(position_unwrapped_m) || !finite_vec(velocity_mps) ||
        !isfinite(charge_C) || !isfinite(mass_kg) ||
        !isfinite(macro_weight) || !(mass_kg > 0.0) ||
        !(macro_weight > 0.0) || !(speed2 >= 0.0) ||
        speed2 >= SWA_C * SWA_C) {
        return 0;
    }
    state->particles[index].position_unwrapped_m = position_unwrapped_m;
    state->particles[index].velocity_mps = velocity_mps;
    state->particles[index].charge_C = charge_C;
    state->particles[index].mass_kg = mass_kg;
    state->particles[index].macro_weight = macro_weight;
    return 1;
}

