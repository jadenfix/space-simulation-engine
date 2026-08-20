
static double pic3d_particle_energy_J(
    const swa_em_pic3d_particle *particles,
    size_t particle_count
) {
    swa_kahan energy = {0.0, 0.0};
    size_t p;
    for (p = 0U; p < particle_count; ++p) {
        const double speed2 = swa_vdot(
            particles[p].velocity_mps, particles[p].velocity_mps
        );
        const double gamma = 1.0 / sqrt(
            1.0 - speed2 / (SWA_C * SWA_C)
        );
        swa_kahan_add(
            &energy,
            particles[p].macro_weight * particles[p].mass_kg *
            SWA_C * SWA_C * (gamma - 1.0)
        );
    }
    return energy.sum;
}

static swa_vec3 pic3d_particle_momentum_Ns(
    const swa_em_pic3d_particle *particles,
    size_t particle_count
) {
    swa_kahan px = {0.0, 0.0};
    swa_kahan py = {0.0, 0.0};
    swa_kahan pz = {0.0, 0.0};
    size_t p;
    for (p = 0U; p < particle_count; ++p) {
        const double speed2 = swa_vdot(
            particles[p].velocity_mps, particles[p].velocity_mps
        );
        const double gamma = 1.0 / sqrt(
            1.0 - speed2 / (SWA_C * SWA_C)
        );
        const double scale =
            particles[p].macro_weight * particles[p].mass_kg * gamma;
        swa_kahan_add(&px, scale * particles[p].velocity_mps.x);
        swa_kahan_add(&py, scale * particles[p].velocity_mps.y);
        swa_kahan_add(&pz, scale * particles[p].velocity_mps.z);
    }
    return swa_v3(px.sum, py.sum, pz.sum);
}

static double pic3d_field_array_energy_J(
    const swa_em_pic3d *state,
    const double *ex,
    const double *ey,
    const double *ez,
    const double *bx,
    const double *by,
    const double *bz
) {
    const double volume = pic3d_cell_volume(state);
    swa_kahan energy = {0.0, 0.0};
    size_t index;
    for (index = 0U; index < state->cell_count; ++index) {
        const double electric2 =
            ex[index] * ex[index] +
            ey[index] * ey[index] +
            ez[index] * ez[index];
        const double magnetic2 =
            bx[index] * bx[index] +
            by[index] * by[index] +
            bz[index] * bz[index];
        swa_kahan_add(
            &energy,
            0.5 * volume * (
                SWA_EPS0 * electric2 + magnetic2 / SWA_MU0
            )
        );
    }
    return energy.sum;
}

static swa_vec3 pic3d_field_array_momentum_Ns(
    const swa_em_pic3d *state,
    const double *ex,
    const double *ey,
    const double *ez,
    const double *bx,
    const double *by,
    const double *bz
) {
    const double volume = pic3d_cell_volume(state);
    swa_kahan px = {0.0, 0.0};
    swa_kahan py = {0.0, 0.0};
    swa_kahan pz = {0.0, 0.0};
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
                const swa_vec3 electric = swa_v3(
                    0.5 * (
                        ex[center] +
                        ex[pic3d_index(state, ip, j, k)]
                    ),
                    0.5 * (
                        ey[center] +
                        ey[pic3d_index(state, i, jp, k)]
                    ),
                    0.5 * (
                        ez[center] +
                        ez[pic3d_index(state, i, j, kp)]
                    )
                );
                const swa_vec3 magnetic = swa_v3(
                    0.25 * (
                        bx[center] +
                        bx[pic3d_index(state, i, jp, k)] +
                        bx[pic3d_index(state, i, j, kp)] +
                        bx[pic3d_index(state, i, jp, kp)]
                    ),
                    0.25 * (
                        by[center] +
                        by[pic3d_index(state, ip, j, k)] +
                        by[pic3d_index(state, i, j, kp)] +
                        by[pic3d_index(state, ip, j, kp)]
                    ),
                    0.25 * (
                        bz[center] +
                        bz[pic3d_index(state, ip, j, k)] +
                        bz[pic3d_index(state, i, jp, k)] +
                        bz[pic3d_index(state, ip, jp, k)]
                    )
                );
                const swa_vec3 density = swa_vscale(
                    swa_vcross(electric, magnetic), SWA_EPS0
                );
                swa_kahan_add(&px, density.x * volume);
                swa_kahan_add(&py, density.y * volume);
                swa_kahan_add(&pz, density.z * volume);
            }
        }
    }
    return swa_v3(px.sum, py.sum, pz.sum);
}

double swa_em_pic3d_field_energy_J(const swa_em_pic3d *state) {
    if (!swa_em_pic3d_state_is_finite(state)) {
        return NAN;
    }
    return pic3d_field_array_energy_J(
        state,
        state->ex_V_m, state->ey_V_m, state->ez_V_m,
        state->bx_T, state->by_T, state->bz_T
    );
}

swa_vec3 swa_em_pic3d_field_momentum_Ns(const swa_em_pic3d *state) {
    if (!swa_em_pic3d_state_is_finite(state)) {
        return swa_v3(NAN, NAN, NAN);
    }
    return pic3d_field_array_momentum_Ns(
        state,
        state->ex_V_m, state->ey_V_m, state->ez_V_m,
        state->bx_T, state->by_T, state->bz_T
    );
}

static void pic3d_electric_constraint_metrics(
    const swa_em_pic3d *state,
    const double *ex,
    const double *ey,
    const double *ez,
    const double *rho,
    double *maximum,
    double *rms,
    double *relative
) {
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    swa_kahan squared = {0.0, 0.0};
    double maximum_electric = 0.0;
    double maximum_charge = 0.0;
    double scale;
    size_t index;
    size_t i;
    size_t j;
    size_t k;
    *maximum = 0.0;
    for (index = 0U; index < state->cell_count; ++index) {
        maximum_electric = pic3d_max4(
            maximum_electric,
            fabs(ex[index]), fabs(ey[index]), fabs(ez[index])
        );
        maximum_charge = pic3d_max2(
            maximum_charge, fabs(rho[index])
        );
    }
    scale = pic3d_max3(
        maximum_charge,
        SWA_EPS0 * maximum_electric * sqrt(
            1.0 / (dx * dx) +
            1.0 / (dy * dy) +
            1.0 / (dz * dz)
        ),
        DBL_MIN
    );
    for (k = 0U; k < state->nz; ++k) {
        const size_t kp = (k + 1U) % state->nz;
        for (j = 0U; j < state->ny; ++j) {
            const size_t jp = (j + 1U) % state->ny;
            for (i = 0U; i < state->nx; ++i) {
                const size_t ip = (i + 1U) % state->nx;
                const size_t center = pic3d_index(state, i, j, k);
                const double divergence_charge = SWA_EPS0 * (
                    (ex[pic3d_index(state, ip, j, k)] - ex[center]) / dx +
                    (ey[pic3d_index(state, i, jp, k)] - ey[center]) / dy +
                    (ez[pic3d_index(state, i, j, kp)] - ez[center]) / dz
                );
                const double residual = divergence_charge - rho[center];
                scale = pic3d_max3(
                    scale, fabs(divergence_charge), fabs(rho[center])
                );
                *maximum = pic3d_max2(*maximum, fabs(residual));
                swa_kahan_add(&squared, residual * residual);
            }
        }
    }
    *rms = sqrt(squared.sum / (double)state->cell_count);
    *relative = *maximum / scale;
}

static void pic3d_magnetic_constraint_metrics(
    const swa_em_pic3d *state,
    const double *bx,
    const double *by,
    const double *bz,
    double *maximum,
    double *rms,
    double *relative
) {
    const double dx = pic3d_dx(state);
    const double dy = pic3d_dy(state);
    const double dz = pic3d_dz(state);
    swa_kahan squared = {0.0, 0.0};
    double maximum_magnetic = 0.0;
    double scale;
    size_t index;
    size_t i;
    size_t j;
    size_t k;
    *maximum = 0.0;
    for (index = 0U; index < state->cell_count; ++index) {
        maximum_magnetic = pic3d_max4(
            maximum_magnetic,
            fabs(bx[index]), fabs(by[index]), fabs(bz[index])
        );
    }
    scale = pic3d_max2(
        maximum_magnetic * sqrt(
            1.0 / (dx * dx) +
            1.0 / (dy * dy) +
            1.0 / (dz * dz)
        ),
        DBL_MIN
    );
    for (k = 0U; k < state->nz; ++k) {
        const size_t km = (k + state->nz - 1U) % state->nz;
        for (j = 0U; j < state->ny; ++j) {
            const size_t jm = (j + state->ny - 1U) % state->ny;
            for (i = 0U; i < state->nx; ++i) {
                const size_t im = (i + state->nx - 1U) % state->nx;
                const size_t center = pic3d_index(state, i, j, k);
                const double term_x =
                    (bx[center] -
                     bx[pic3d_index(state, im, j, k)]) / dx;
                const double term_y =
                    (by[center] -
                     by[pic3d_index(state, i, jm, k)]) / dy;
                const double term_z =
                    (bz[center] -
                     bz[pic3d_index(state, i, j, km)]) / dz;
                const double residual = term_x + term_y + term_z;
                scale = pic3d_max4(
                    scale, fabs(term_x), fabs(term_y), fabs(term_z)
                );
                *maximum = pic3d_max2(*maximum, fabs(residual));
                swa_kahan_add(&squared, residual * residual);
            }
        }
    }
    *rms = sqrt(squared.sum / (double)state->cell_count);
    *relative = *maximum / scale;
}

double swa_em_pic3d_electric_divergence_max_C_m3(
    const swa_em_pic3d *state,
    const double *rho_C_m3,
    double *relative
) {
    double maximum;
    double rms;
    if (!swa_em_pic3d_state_is_finite(state) ||
        rho_C_m3 == NULL || relative == NULL ||
        !pic3d_finite_array(rho_C_m3, state->cell_count)) {
        return NAN;
    }
    pic3d_electric_constraint_metrics(
        state,
        state->ex_V_m, state->ey_V_m, state->ez_V_m,
        rho_C_m3, &maximum, &rms, relative
    );
    return maximum;
}

double swa_em_pic3d_magnetic_divergence_max_T_m(
    const swa_em_pic3d *state,
    double *relative
) {
    double maximum;
    double rms;
    if (!swa_em_pic3d_state_is_finite(state) || relative == NULL) {
        return NAN;
    }
    pic3d_magnetic_constraint_metrics(
        state, state->bx_T, state->by_T, state->bz_T,
        &maximum, &rms, relative
    );
    return maximum;
}

static swa_vec3 pic3d_mean_current(
    const swa_em_pic3d *state,
    const double *jx,
    const double *jy,
    const double *jz
) {
    swa_kahan x = {0.0, 0.0};
    swa_kahan y = {0.0, 0.0};
    swa_kahan z = {0.0, 0.0};
    size_t index;
    for (index = 0U; index < state->cell_count; ++index) {
        swa_kahan_add(&x, jx[index]);
        swa_kahan_add(&y, jy[index]);
        swa_kahan_add(&z, jz[index]);
    }
    return swa_v3(
        x.sum / (double)state->cell_count,
        y.sum / (double)state->cell_count,
        z.sum / (double)state->cell_count
    );
}

static double pic3d_current_work_J(
    const swa_em_pic3d *state,
    const double *jx,
    const double *jy,
    const double *jz,
    const double *ex_mid,
    const double *ey_mid,
    const double *ez_mid
) {
    const double measure = state->dt_s * pic3d_cell_volume(state);
    swa_kahan work = {0.0, 0.0};
    size_t index;
    for (index = 0U; index < state->cell_count; ++index) {
        swa_kahan_add(
            &work,
            measure * (
                jx[index] * ex_mid[index] +
                jy[index] * ey_mid[index] +
                jz[index] * ez_mid[index]
            )
        );
    }
    return work.sum;
}
