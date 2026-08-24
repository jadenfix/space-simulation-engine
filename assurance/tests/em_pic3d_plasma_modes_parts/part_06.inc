
    check_true(ran_coarse, "coarse plasma mode executes");
    check_true(ran_medium, "medium plasma mode executes");
    check_true(ran_fine, "fine plasma mode executes");
    check_true(ran_high_density, "high-density plasma mode executes");
    check_true(ran_heavy_electron, "heavy-electron plasma mode executes");
    check_true(ran_pair_plasma, "pair-plasma mode executes");
    check_true(
        same_sign_rejected,
        "same-sign charge control is rejected by periodic Gauss initialization"
    );

    all_ran =
        ran_coarse && ran_medium && ran_fine &&
        ran_high_density && ran_heavy_electron && ran_pair_plasma;
    if (all_ran) {
        const double measured_density_ratio =
            high_density.measured_omega_rad_s /
            fine.measured_omega_rad_s;
        const double predicted_density_ratio =
            high_density.midpoint_omega_rad_s /
            fine.midpoint_omega_rad_s;
        const double measured_mass_ratio =
            heavy_electron.measured_omega_rad_s /
            fine.measured_omega_rad_s;
        const double predicted_mass_ratio =
            heavy_electron.midpoint_omega_rad_s /
            fine.midpoint_omega_rad_s;
        const double wrong_pair_omega = midpoint_omega_rad_s(
            sqrt(
                base_density * SWA_QE * SWA_QE /
                (SWA_EPS0 * SWA_ME)
            ),
            pair_plasma.dt_s
        );

        check_true(coarse.passes, "coarse plasma mode has enough periods");
        check_true(medium.passes, "medium plasma mode has enough periods");
        check_true(fine.passes, "fine plasma mode has enough periods");
        check_true(
            high_density.passes,
            "high-density plasma mode has enough periods"
        );
        check_true(
            heavy_electron.passes,
            "heavy-electron plasma mode has enough periods"
        );
        check_true(
            pair_plasma.passes,
            "pair-plasma mode has enough periods"
        );

        coarse_continuum_error = relative_scalar(
            coarse.measured_omega_rad_s,
            coarse.continuum_omega_rad_s
        );
        medium_continuum_error = relative_scalar(
            medium.measured_omega_rad_s,
            medium.continuum_omega_rad_s
        );
        fine_continuum_error = relative_scalar(
            fine.measured_omega_rad_s,
            fine.continuum_omega_rad_s
        );
        coarse_midpoint_error = relative_scalar(
            coarse.measured_omega_rad_s,
            coarse.midpoint_omega_rad_s
        );
        medium_midpoint_error = relative_scalar(
            medium.measured_omega_rad_s,
            medium.midpoint_omega_rad_s
        );
        fine_midpoint_error = relative_scalar(
            fine.measured_omega_rad_s,
            fine.midpoint_omega_rad_s
        );
        if (coarse_continuum_error > DBL_MIN &&
            medium_continuum_error > DBL_MIN &&
            fine_continuum_error > DBL_MIN) {
            order_coarse_medium = log(
                coarse_continuum_error / medium_continuum_error
            ) / log(2.0);
            order_medium_fine = log(
                medium_continuum_error / fine_continuum_error
            ) / log(2.0);
        }
        density_ratio_error = relative_scalar(
            measured_density_ratio, predicted_density_ratio
        );
        mass_ratio_error = relative_scalar(
            measured_mass_ratio, predicted_mass_ratio
        );
        wrong_mass_model_error = relative_scalar(
            pair_plasma.measured_omega_rad_s, wrong_pair_omega
        );

        check_true(
            coarse_midpoint_error < 3.0e-4,
            "coarse frequency matches midpoint dispersion"
        );
        check_true(
            medium_midpoint_error < 1.0e-4,
            "medium frequency matches midpoint dispersion"
        );
        check_true(
            fine_midpoint_error < 5.0e-5,
            "fine frequency matches midpoint dispersion"
        );
        check_true(
            coarse_midpoint_error < coarse_continuum_error,
