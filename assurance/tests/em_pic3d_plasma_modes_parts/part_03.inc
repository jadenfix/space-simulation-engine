    requested_steps =
        periods * 2.0 * SWA_PI /
        (result->continuum_omega_rad_s * dt_s);
    if (!isfinite(requested_steps) ||
        requested_steps < 1.0 || requested_steps > 100000.0) {
        return 0;
    }
    result->steps = (size_t)ceil(requested_steps) + 2U;

    if (!swa_em_pic3d_init(
            &state,
            MODE_NX, MODE_NY, MODE_NZ, MODE_PARTICLES,
            length_x_m, length_y_m, length_z_m, dt_s
        ) ||
        !configure_neutral_pair(
            &state, density_m3, negative_mass_kg, positive_mass_kg
        ) ||
        !swa_em_pic3d_initialize_gauss_field(
            &state,
            swa_v3(MODE_FIELD_AMPLITUDE_V_M, 0.0, 0.0),
            charge_tolerance
        )) {
        swa_em_pic3d_destroy(&state);
        return 0;
    }

    result->initial_mean_e_V_m = mean_ex_V_m(&state);
    initial_energy = total_energy_J(&state);
    previous_previous_e = result->initial_mean_e_V_m;
    previous_e = result->initial_mean_e_V_m;

    for (step = 1U; step <= result->steps; ++step) {
        double current_e;
        double current_energy;
        double energy_scale;
        double energy_drift;
        int executed;
        memset(&step_result, 0, sizeof(step_result));
        executed = swa_em_pic3d_step(&state, &gate, &step_result);
        if (!executed ||
            !step_result.finite || !step_result.passes ||
            !step_result.state_committed) {
            fprintf(
                stderr,
                "MODE step failure: step=%zu executed=%d finite=%d "
                "passes=%d committed=%d\n",
                step,
                executed,
                step_result.finite,
                step_result.passes,
                step_result.state_committed
            );
            fprintf(
                stderr,
                "MODE flags[cfl=%d pcfl=%d q=%d cont=%d g0=%d g1=%d "
                "d0=%d d1=%d proj=%d local=%d spectral=%d mean=%d "
                "nonlinear=%d sublum=%d energy=%d pwork=%d fwork=%d "
                "mom=%d]\n",
                step_result.courant_passes,
                step_result.particle_courant_passes,
                step_result.charge_closure_passes,
                step_result.continuity_passes,
                step_result.initial_gauss_passes,
                step_result.final_gauss_passes,
                step_result.initial_magnetic_divergence_passes,
                step_result.final_magnetic_divergence_passes,
                step_result.current_projection_passes,
                step_result.local_current_passes,
                step_result.spectral_oracle_passes,
                step_result.mean_current_passes,
                step_result.nonlinear_passes,
                step_result.subluminal_passes,
                step_result.energy_passes,
                step_result.particle_work_passes,
                step_result.field_work_passes,
                step_result.momentum_passes
            );
            fprintf(
                stderr,
                "MODE metrics[cont=%.17g g0=%.17g g1=%.17g local=%.17g "
                "spectral=%.17g nonlinear=%.17g energy=%.17g pwork=%.17g "
                "fwork=%.17g mom=%.17g disp=%.17g speed=%.17g]\n",
                step_result.continuity_relative_max,
                step_result.initial_gauss_relative_max,
                step_result.final_gauss_relative_max,
                step_result.local_current_difference_relative,
                step_result.spectral_oracle_difference_relative,
                step_result.nonlinear_relative_residual,
                step_result.total_energy_relative,
                step_result.particle_work_relative,
                step_result.field_work_relative,
                step_result.total_momentum_relative,
                step_result.maximum_particle_displacement_cells,
                step_result.maximum_speed_fraction_c
            );
            swa_em_pic3d_destroy(&state);
            return 0;
        }

        current_e = mean_ex_V_m(&state);
        current_energy = total_energy_J(&state);
        energy_scale = max3(
            fabs(initial_energy), fabs(current_energy), DBL_MIN
        );
        energy_drift = fabs(current_energy - initial_energy) / energy_scale;
        result->maximum_total_energy_relative_drift = max2(
            result->maximum_total_energy_relative_drift, energy_drift
        );
        update_step_maxima(result, &step_result);

        if (previous_e > 0.0 && current_e <= 0.0) {
            const double denominator = previous_e - current_e;
            if (!(denominator > 0.0) ||
                result->crossings >= MODE_MAX_CROSSINGS) {
                swa_em_pic3d_destroy(&state);
                return 0;
            }
            crossing_time_s[result->crossings] =
                ((double)(step - 1U) + previous_e / denominator) * dt_s;
            result->crossings += 1U;
        }

        if (step > 1U &&
            previous_e > previous_previous_e &&
            previous_e >= current_e &&
            previous_e > 0.0) {
            const double curvature =
                previous_previous_e - 2.0 * previous_e + current_e;
            const double difference = current_e - previous_previous_e;
            double interpolated_peak = previous_e;
            if (curvature < -DBL_MIN) {
                interpolated_peak -=
                    difference * difference / (8.0 * curvature);
            }
            if (result->positive_peaks == 0U) {
                result->first_positive_peak_V_m = interpolated_peak;
            }
            result->last_positive_peak_V_m = interpolated_peak;
            result->positive_peaks += 1U;
        }

        previous_previous_e = previous_e;
        previous_e = current_e;
    }

    if (result->crossings >= 2U) {
        for (period_index = 1U;
             period_index < result->crossings;
             ++period_index) {
            swa_kahan_add(
                &period_sum,
                crossing_time_s[period_index] -
                crossing_time_s[period_index - 1U]
