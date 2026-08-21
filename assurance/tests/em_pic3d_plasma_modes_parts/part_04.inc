            );
        }
        result->measured_omega_rad_s =
            2.0 * SWA_PI /
            (period_sum.sum / (double)(result->crossings - 1U));
    }
    if (result->positive_peaks >= 2U) {
        result->peak_relative_change = relative_scalar(
            result->first_positive_peak_V_m,
            result->last_positive_peak_V_m
        );
    }

    result->finite =
        isfinite(result->continuum_omega_rad_s) &&
        isfinite(result->midpoint_omega_rad_s) &&
        isfinite(result->measured_omega_rad_s) &&
        isfinite(result->initial_mean_e_V_m) &&
        isfinite(result->first_positive_peak_V_m) &&
        isfinite(result->last_positive_peak_V_m) &&
        isfinite(result->peak_relative_change) &&
        isfinite(result->maximum_total_energy_relative_drift) &&
        isfinite(result->maximum_step_energy_relative) &&
        isfinite(result->maximum_particle_work_relative) &&
        isfinite(result->maximum_field_work_relative) &&
        isfinite(result->maximum_continuity_relative) &&
        isfinite(result->maximum_gauss_relative) &&
        isfinite(result->maximum_momentum_relative);
    result->passes =
        result->finite &&
        result->crossings >= 4U &&
        result->positive_peaks >= 3U;

    swa_em_pic3d_destroy(&state);
    return 1;
}

static int write_case(
    FILE *fp,
    const char *name,
    const mode_result *value,
    int trailing_comma
) {
    return fprintf(
        fp,
        "    \"%s\": {"
        "\"density_m3\": %.17g, "
        "\"negative_mass_kg\": %.17g, "
        "\"positive_mass_kg\": %.17g, "
        "\"dt_s\": %.17g, "
        "\"continuum_omega_rad_s\": %.17g, "
        "\"midpoint_omega_rad_s\": %.17g, "
        "\"measured_omega_rad_s\": %.17g, "
        "\"crossings\": %zu, "
        "\"positive_peaks\": %zu, "
        "\"peak_relative_change\": %.17g, "
        "\"maximum_total_energy_relative_drift\": %.17g, "
        "\"maximum_step_energy_relative\": %.17g, "
        "\"maximum_particle_work_relative\": %.17g, "
        "\"maximum_field_work_relative\": %.17g, "
        "\"maximum_continuity_relative\": %.17g, "
        "\"maximum_gauss_relative\": %.17g, "
        "\"maximum_momentum_relative\": %.17g, "
        "\"passes\": %s}%s\n",
        name,
        value->density_m3,
        value->negative_mass_kg,
        value->positive_mass_kg,
        value->dt_s,
        value->continuum_omega_rad_s,
        value->midpoint_omega_rad_s,
        value->measured_omega_rad_s,
        value->crossings,
        value->positive_peaks,
        value->peak_relative_change,
        value->maximum_total_energy_relative_drift,
        value->maximum_step_energy_relative,
        value->maximum_particle_work_relative,
        value->maximum_field_work_relative,
        value->maximum_continuity_relative,
        value->maximum_gauss_relative,
        value->maximum_momentum_relative,
        value->passes ? "true" : "false",
        trailing_comma ? "," : ""
    ) >= 0;
}

static int write_receipt(
    const mode_result *coarse,
    const mode_result *medium,
    const mode_result *fine,
    const mode_result *high_density,
    const mode_result *heavy_electron,
    const mode_result *pair_plasma,
    double order_coarse_medium,
    double order_medium_fine,
    double density_ratio_error,
    double mass_ratio_error,
    double wrong_mass_model_error,
    int same_sign_rejected
) {
    FILE *fp = fopen("output/em_pic3d_plasma_modes_receipt.json", "w");
    int ok = fp != NULL;
    if (!ok) {
        return 0;
    }
    ok = ok && fprintf(
        fp,
        "{\n"
        "  \"schema\": \"spacewind.em_pic3d.plasma_modes.v1\",\n"
