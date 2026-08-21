            "coarse measurement prefers midpoint over continuum frequency"
        );
        check_true(
            medium_midpoint_error < medium_continuum_error,
            "medium measurement prefers midpoint over continuum frequency"
        );
        check_true(
            fine_midpoint_error < fine_continuum_error,
            "fine measurement prefers midpoint over continuum frequency"
        );
        check_true(
            coarse_continuum_error > medium_continuum_error &&
            medium_continuum_error > fine_continuum_error,
            "continuum phase error decreases under timestep refinement"
        );
        check_true(
            order_coarse_medium > 1.5 && order_coarse_medium < 2.5,
            "coarse-to-medium phase convergence is second order"
        );
        check_true(
            order_medium_fine > 1.5 && order_medium_fine < 2.5,
            "medium-to-fine phase convergence is second order"
        );
        check_true(
            coarse_continuum_error > 6.0e-4 &&
            fine_continuum_error < 8.0e-5,
            "refinement resolves and reduces midpoint phase dispersion"
        );

        check_true(
            density_ratio_error < 3.0e-4,
            "plasma frequency follows discrete square-root density scaling"
        );
        check_true(
            relative_scalar(measured_density_ratio, 2.0) < 3.0e-4,
            "fourfold density nearly doubles measured plasma frequency"
        );
        check_true(
            mass_ratio_error < 3.0e-4,
            "plasma frequency follows two-species charge-to-mass scaling"
        );
        check_true(
            measured_mass_ratio < 0.60 &&
            measured_mass_ratio > 0.45,
            "fourfold electron mass produces a resolvable frequency reduction"
        );
        check_true(
            relative_scalar(
                pair_plasma.measured_omega_rad_s,
                pair_plasma.midpoint_omega_rad_s
            ) < 2.0e-4,
            "pair-plasma frequency matches the two-mobile-species model"
        );
        check_true(
            wrong_mass_model_error > 0.20,
            "pair-plasma benchmark rejects an immobile-positive-species model"
        );

        check_true(
            coarse.maximum_continuity_relative < 2.0e-10 &&
            medium.maximum_continuity_relative < 2.0e-10 &&
            fine.maximum_continuity_relative < 2.0e-10 &&
            high_density.maximum_continuity_relative < 2.0e-10 &&
            heavy_electron.maximum_continuity_relative < 2.0e-10 &&
            pair_plasma.maximum_continuity_relative < 2.0e-10,
            "all plasma modes preserve local discrete continuity"
        );
        check_true(
            coarse.maximum_gauss_relative < 2.0e-10 &&
            medium.maximum_gauss_relative < 2.0e-10 &&
            fine.maximum_gauss_relative < 2.0e-10 &&
            high_density.maximum_gauss_relative < 2.0e-10 &&
            heavy_electron.maximum_gauss_relative < 2.0e-10 &&
            pair_plasma.maximum_gauss_relative < 2.0e-10,
            "all plasma modes preserve conditioned Gauss closure"
        );
        check_true(
            coarse.maximum_spectral_oracle_relative < 1.0e-8 &&
            medium.maximum_spectral_oracle_relative < 1.0e-8 &&
            fine.maximum_spectral_oracle_relative < 1.0e-8 &&
            high_density.maximum_spectral_oracle_relative < 1.0e-8 &&
            heavy_electron.maximum_spectral_oracle_relative < 1.0e-8 &&
            pair_plasma.maximum_spectral_oracle_relative < 1.0e-8,
            "all local currents remain within the spectral-oracle envelope"
        );

        check_true(
            coarse.peak_relative_change < 5.0e-4 &&
            medium.peak_relative_change < 5.0e-4 &&
            fine.peak_relative_change < 5.0e-4,
            "refinement cases preserve oscillation amplitude"
        );
        check_true(
            high_density.peak_relative_change < 5.0e-4 &&
            heavy_electron.peak_relative_change < 5.0e-4 &&
            pair_plasma.peak_relative_change < 5.0e-4,
            "scaling cases preserve oscillation amplitude"
        );
        check_true(
            coarse.maximum_total_energy_relative_drift < 2.0e-5 &&
            medium.maximum_total_energy_relative_drift < 2.0e-5 &&
            fine.maximum_total_energy_relative_drift < 2.0e-5,
            "refinement cases bound long-run total-energy drift"
        );
        check_true(
            high_density.maximum_total_energy_relative_drift < 2.0e-5 &&
            heavy_electron.maximum_total_energy_relative_drift < 2.0e-5 &&
            pair_plasma.maximum_total_energy_relative_drift < 2.0e-5,
            "scaling cases bound long-run total-energy drift"
        );
    }

    check_true(
        write_receipt(
            &coarse, &medium, &fine,
            &high_density, &heavy_electron, &pair_plasma,
            order_coarse_medium, order_medium_fine,
            density_ratio_error, mass_ratio_error,
            wrong_mass_model_error, same_sign_rejected
        ),
        "plasma-mode receipt writes"
    );

    fprintf(
        stderr,
        "spacewind 3D3V plasma modes: %zu checks, %zu failures, "
        "omega[c=%.17g m=%.17g f=%.17g] "
        "order[%.17g %.17g] hash=%016" PRIx64 "\n",
        checks, failures,
        coarse.measured_omega_rad_s,
        medium.measured_omega_rad_s,
        fine.measured_omega_rad_s,
        order_coarse_medium,
        order_medium_fine,
        transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
