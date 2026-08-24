        "  \"checks\": %zu,\n"
        "  \"failures\": %zu,\n"
        "  \"transcript_hash\": \"%016" PRIx64 "\",\n"
        "  \"cases\": {\n",
        checks, failures, transcript_hash
    ) >= 0;
    ok = ok && write_case(fp, "coarse", coarse, 1);
    ok = ok && write_case(fp, "medium", medium, 1);
    ok = ok && write_case(fp, "fine", fine, 1);
    ok = ok && write_case(fp, "high_density", high_density, 1);
    ok = ok && write_case(fp, "heavy_electron", heavy_electron, 1);
    ok = ok && write_case(fp, "pair_plasma", pair_plasma, 0);
    ok = ok && fprintf(
        fp,
        "  },\n"
        "  \"derived\": {"
        "\"order_coarse_medium\": %.17g, "
        "\"order_medium_fine\": %.17g, "
        "\"density_ratio_error\": %.17g, "
        "\"mass_ratio_error\": %.17g, "
        "\"wrong_mass_model_error\": %.17g},\n"
        "  \"adversarial\": {"
        "\"same_sign_charge_rejected\": %s},\n"
        "  \"passes\": %s\n"
        "}\n",
        order_coarse_medium,
        order_medium_fine,
        density_ratio_error,
        mass_ratio_error,
        wrong_mass_model_error,
        same_sign_rejected ? "true" : "false",
        failures == 0U ? "true" : "false"
    ) >= 0;
    if (fclose(fp) != 0) {
        ok = 0;
    }
    return ok;
}

int main(void) {
    const double length_x_m = 1.0;
    const double length_y_m = 1.0;
    const double length_z_m = 1.0;
    const double dx = length_x_m / (double)MODE_NX;
    const double dy = length_y_m / (double)MODE_NY;
    const double dz = length_z_m / (double)MODE_NZ;
    const double inverse_spacing2 =
        1.0 / (dx * dx) +
        1.0 / (dy * dy) +
        1.0 / (dz * dz);
    const double dt_coarse =
        0.12 / (SWA_C * sqrt(inverse_spacing2));
    const double target_omega =
        2.0 * SWA_PI / (64.0 * dt_coarse);
    const double base_density =
        target_omega * target_omega * SWA_EPS0 /
        (SWA_QE * SWA_QE * (1.0 / SWA_ME + 1.0 / SWA_MP));
    mode_result coarse;
    mode_result medium;
    mode_result fine;
    mode_result high_density;
    mode_result heavy_electron;
    mode_result pair_plasma;
    int ran_coarse;
    int ran_medium;
    int ran_fine;
    int ran_high_density;
    int ran_heavy_electron;
    int ran_pair_plasma;
    int same_sign_rejected;
    int all_ran;
    double coarse_continuum_error = 0.0;
    double medium_continuum_error = 0.0;
    double fine_continuum_error = 0.0;
    double coarse_midpoint_error = 0.0;
    double medium_midpoint_error = 0.0;
    double fine_midpoint_error = 0.0;
    double order_coarse_medium = 0.0;
    double order_medium_fine = 0.0;
    double density_ratio_error = 0.0;
    double mass_ratio_error = 0.0;
    double wrong_mass_model_error = 0.0;

    ran_coarse = run_mode(
        dt_coarse, base_density, SWA_ME, SWA_MP,
        MODE_PERIODS, &coarse
    );
    ran_medium = run_mode(
        0.5 * dt_coarse, base_density, SWA_ME, SWA_MP,
        MODE_PERIODS, &medium
    );
    ran_fine = run_mode(
        0.25 * dt_coarse, base_density, SWA_ME, SWA_MP,
        MODE_PERIODS, &fine
    );
    ran_high_density = run_mode(
        0.25 * dt_coarse, 4.0 * base_density, SWA_ME, SWA_MP,
        MODE_PERIODS, &high_density
    );
    ran_heavy_electron = run_mode(
        0.25 * dt_coarse, base_density, 4.0 * SWA_ME, SWA_MP,
        MODE_PERIODS, &heavy_electron
    );
    ran_pair_plasma = run_mode(
        0.25 * dt_coarse, base_density, SWA_ME, SWA_ME,
        MODE_PERIODS, &pair_plasma
    );
    same_sign_rejected = same_sign_charge_control_rejected(
        0.25 * dt_coarse, base_density
    );
