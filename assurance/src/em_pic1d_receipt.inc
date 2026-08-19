int swa_em_pic1d_write_receipt(
    FILE *fp,
    const swa_em_pic1d *state,
    const swa_em_pic1d_limits *limits,
    const swa_em_pic1d_result *result
) {
    if (fp == NULL || state == NULL || limits == NULL || result == NULL) {
        return 0;
    }
    return fprintf(
        fp,
        "{\n"
        "  \"schema\": \"spacewind.em-pic1d/v1\",\n"
        "  \"cells\": %zu,\n"
        "  \"particles\": %zu,\n"
        "  \"dt_s\": %.17g,\n"
        "  \"courant_number\": %.17g,\n"
        "  \"maximum_particle_displacement_cells\": %.17g,\n"
        "  \"continuity_relative_max\": %.17g,\n"
        "  \"initial_gauss_relative_max\": %.17g,\n"
        "  \"final_gauss_relative_max\": %.17g,\n"
        "  \"magnetic_divergence_residual_T_m\": %.17g,\n"
        "  \"mean_jx_A_m2\": %.17g,\n"
        "  \"mean_jy_A_m2\": %.17g,\n"
        "  \"mean_jz_A_m2\": %.17g,\n"
        "  \"transverse_current_partition_relative_max\": %.17g,\n"
        "  \"harmonic_ampere_relative\": %.17g,\n"
        "  \"particle_energy_change_J\": %.17g,\n"
        "  \"field_energy_change_J\": %.17g,\n"
        "  \"midpoint_current_work_J\": %.17g,\n"
        "  \"particle_work_residual_J\": %.17g,\n"
        "  \"field_work_residual_J\": %.17g,\n"
        "  \"total_energy_relative\": %.17g,\n"
        "  \"total_momentum_relative\": %.17g,\n"
        "  \"maximum_speed_fraction_c\": %.17g,\n"
        "  \"continuity_passes\": %s,\n"
        "  \"gauss_passes\": %s,\n"
        "  \"energy_passes\": %s,\n"
        "  \"momentum_passes\": %s,\n"
        "  \"passes\": %s,\n"
        "  \"limits\": {\n"
        "    \"energy_relative\": %.17g,\n"
        "    \"momentum_relative\": %.17g,\n"
        "    \"maximum_courant\": %.17g\n"
        "  },\n"
        "  \"claim\": \"bounded one-dimensional three-velocity-component electromagnetic particle-field reference step\",\n"
        "  \"nonclaim\": \"this explicit periodic reference does not establish multidimensional kinetic accuracy, exact energy-momentum conservation, plasma-wing force, chamber similarity, or flight propulsion\"\n"
        "}\n",
        state->cells,
        state->particle_count,
        state->dt_s,
        result->courant_number,
        result->maximum_particle_displacement_cells,
        result->continuity_relative_max,
        result->initial_gauss_relative_max,
        result->final_gauss_relative_max,
        result->magnetic_divergence_residual_T_m,
        result->mean_jx_A_m2,
        result->mean_jy_A_m2,
        result->mean_jz_A_m2,
        result->transverse_current_partition_relative_max,
        result->harmonic_ampere_relative,
        result->particle_energy_change_J,
        result->field_energy_change_J,
        result->midpoint_current_work_J,
        result->particle_work_residual_J,
        result->field_work_residual_J,
        result->total_energy_relative,
        result->total_momentum_relative,
        result->maximum_speed_fraction_c,
        result->continuity_passes ? "true" : "false",
        (result->initial_gauss_passes && result->final_gauss_passes) ?
            "true" : "false",
        result->energy_passes ? "true" : "false",
        result->momentum_passes ? "true" : "false",
        result->passes ? "true" : "false",
        limits->energy_relative_tolerance,
        limits->momentum_relative_tolerance,
        limits->maximum_courant
    ) > 0;
}
