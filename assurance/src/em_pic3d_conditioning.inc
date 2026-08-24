/*
 * Condition Gauss-law residuals against both the represented field/charge
 * amplitude and the absolute charge inventory.  The latter prevents a
 * physically neutral state from dividing a roundoff-sized residual by a
 * cancellation-sized net charge or a near-zero field at a phase crossing.
 */
static double pic3d_reported_scale(
    double maximum,
    double relative
) {
    double scale;
    if (!isfinite(maximum) || maximum < 0.0 ||
        !isfinite(relative) || !(relative > 0.0)) {
        return DBL_MIN;
    }
    scale = maximum / relative;
    return isfinite(scale) && scale > 0.0 ? scale : DBL_MIN;
}

static double pic3d_absolute_charge_density_C_m3(
    const swa_em_pic3d *state
) {
    const double volume =
        state->length_x_m * state->length_y_m * state->length_z_m;
    const double absolute_charge = pic3d_absolute_particle_charge_C(
        state->particles, state->particle_count
    );
    if (!isfinite(volume) || !(volume > 0.0) ||
        !isfinite(absolute_charge) || absolute_charge < 0.0) {
        return NAN;
    }
    return absolute_charge / volume;
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
    double unconditioned_relative;
    double represented_scale;
    double absolute_charge_density;
    double scale;

    pic3d_electric_constraint_metrics_unconditioned(
        state, ex, ey, ez, rho,
        maximum, rms, &unconditioned_relative
    );
    represented_scale = pic3d_reported_scale(
        *maximum, unconditioned_relative
    );
    absolute_charge_density =
        pic3d_absolute_charge_density_C_m3(state);
    scale = pic3d_max3(
        represented_scale,
        isfinite(absolute_charge_density) ?
            absolute_charge_density : DBL_MIN,
        DBL_MIN
    );
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
