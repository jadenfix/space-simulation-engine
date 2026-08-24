        const swa_vec3 position = swa_v3(
            ((double)i + 0.5) * dx,
            ((double)j + 0.5) * dy,
            ((double)k + 0.5) * dz
        );
        if (!swa_em_pic3d_set_particle(
                &state, 2U * cell, position,
                swa_v3(0.0, 0.0, 0.0),
                SWA_QE, SWA_ME, macro_weight
            ) ||
            !swa_em_pic3d_set_particle(
                &state, 2U * cell + 1U, position,
                swa_v3(0.0, 0.0, 0.0),
                SWA_QE, SWA_MP, macro_weight
            )) {
            configured = 0;
            break;
        }
    }
    gauss_initialized = configured ?
        swa_em_pic3d_initialize_gauss_field(
            &state, swa_v3(0.0, 0.0, 0.0), charge_tolerance
        ) : 0;
    swa_em_pic3d_destroy(&state);
    return configured && !gauss_initialized;
}

static void update_step_maxima(
    mode_result *result,
    const swa_em_pic3d_result *step
) {
    result->maximum_step_energy_relative = max2(
        result->maximum_step_energy_relative,
        step->total_energy_relative
    );
    result->maximum_particle_work_relative = max2(
        result->maximum_particle_work_relative,
        step->particle_work_relative
    );
    result->maximum_field_work_relative = max2(
        result->maximum_field_work_relative,
        step->field_work_relative
    );
    result->maximum_continuity_relative = max2(
        result->maximum_continuity_relative,
        step->continuity_relative_max
    );
    result->maximum_gauss_relative = max3(
        result->maximum_gauss_relative,
        step->initial_gauss_relative_max,
        step->final_gauss_relative_max
    );
    result->maximum_spectral_oracle_relative = max2(
        result->maximum_spectral_oracle_relative,
        step->spectral_oracle_difference_relative
    );
    result->maximum_momentum_relative = max2(
        result->maximum_momentum_relative,
        step->total_momentum_relative
    );
}

static int run_mode(
    double dt_s,
    double density_m3,
    double negative_mass_kg,
    double positive_mass_kg,
    double periods,
    mode_result *result
) {
    swa_em_pic3d state;
    swa_em_pic3d_limits gate = mode_limits();
    swa_em_pic3d_result step_result;
    const double length_x_m = 1.0;
    const double length_y_m = 1.0;
    const double length_z_m = 1.0;
    const double dx = length_x_m / (double)MODE_NX;
    const double dy = length_y_m / (double)MODE_NY;
    const double dz = length_z_m / (double)MODE_NZ;
    const double volume = dx * dy * dz;
    const double absolute_charge =
        2.0 * (double)(MODE_NX * MODE_NY * MODE_NZ) *
        density_m3 * volume * SWA_QE;
    const double charge_tolerance =
        max2(1.0e-13 * absolute_charge, 1.0e-30);
    double crossing_time_s[MODE_MAX_CROSSINGS];
    double previous_previous_e;
    double previous_e;
    double initial_energy;
    double requested_steps;
    size_t step;
    size_t period_index;
    swa_kahan period_sum = {0.0, 0.0};

    if (result == NULL ||
        !isfinite(dt_s) || !(dt_s > 0.0) ||
        !isfinite(density_m3) || !(density_m3 > 0.0) ||
        !isfinite(periods) || !(periods > 2.0)) {
        return 0;
    }
    memset(result, 0, sizeof(*result));
    memset(&state, 0, sizeof(state));
    memset(&step_result, 0, sizeof(step_result));

    result->density_m3 = density_m3;
    result->negative_mass_kg = negative_mass_kg;
    result->positive_mass_kg = positive_mass_kg;
    result->dt_s = dt_s;
    result->continuum_omega_rad_s = plasma_omega_rad_s(
        density_m3, negative_mass_kg, positive_mass_kg
    );
    result->midpoint_omega_rad_s = midpoint_omega_rad_s(
        result->continuum_omega_rad_s, dt_s
    );
