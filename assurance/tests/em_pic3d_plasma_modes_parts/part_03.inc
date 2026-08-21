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
