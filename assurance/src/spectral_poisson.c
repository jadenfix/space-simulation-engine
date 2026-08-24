#include "spacewind/spectral_poisson.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static double max2(double a, double b) { return a > b ? a : b; }
static double max3(double a, double b, double c) {
    return max2(max2(a, b), c);
}

static int finite_array(const double *values, size_t count) {
    size_t i;
    if (values == NULL) {
        return 0;
    }
    for (i = 0U; i < count; ++i) {
        if (!isfinite(values[i])) {
            return 0;
        }
    }
    return 1;
}

static int valid_nonnegative(double value) {
    return isfinite(value) && value >= 0.0;
}

static int residual_passes(
    double absolute_value,
    double relative_value,
    double absolute_tolerance,
    double relative_tolerance
) {
    return absolute_value <= absolute_tolerance ||
           relative_value <= relative_tolerance;
}

static void dft_coefficient(
    const double *values,
    size_t samples,
    size_t mode,
    double *real_out,
    double *imaginary_out
) {
    swa_kahan real_sum = {0.0, 0.0};
    swa_kahan imaginary_sum = {0.0, 0.0};
    size_t j;
    for (j = 0U; j < samples; ++j) {
        const double angle =
            -2.0 * SWA_PI * (double)mode * (double)j /
            (double)samples;
        swa_kahan_add(&real_sum, values[j] * cos(angle));
        swa_kahan_add(&imaginary_sum, values[j] * sin(angle));
    }
    *real_out = real_sum.sum;
    *imaginary_out = imaginary_sum.sum;
}

static double complex_norm(double real_value, double imaginary_value) {
    return hypot(real_value, imaginary_value);
}

void swa_default_spectral_poisson_tolerances(
    swa_spectral_poisson_tolerances *out
) {
    if (out == NULL) {
        return;
    }
    out->operator_relative_tolerance = 2.0e-11;
    out->charge_absolute_tolerance_C_m3 = 1.0e-18;
    out->electric_absolute_tolerance_V_m = 1.0e-9;
    out->energy_relative_tolerance = 2.0e-11;
    out->energy_absolute_tolerance_J_m2 = 1.0e-18;
    out->neutrality_relative_tolerance = 2.0e-12;
    out->neutrality_absolute_tolerance_C_m3 = 1.0e-18;
    out->harmonic_field_relative_tolerance = 2.0e-12;
    out->harmonic_field_absolute_tolerance_V_m = 1.0e-9;
    out->nyquist_charge_power_fraction_tolerance = 1.0e-18;
}

int swa_audit_periodic_spectral_poisson(
    const swa_spectral_poisson_input *input,
    const swa_spectral_poisson_tolerances *tolerances,
    swa_spectral_poisson_result *out
) {
    size_t mode;
    size_t sample;
    const int even_samples =
        input != NULL && (input->samples % 2U) == 0U;
    double charge_zero_real = 0.0;
    double charge_zero_imaginary = 0.0;
    double electric_zero_real = 0.0;
    double electric_zero_imaginary = 0.0;
    double nyquist_charge_power = 0.0;
    double nonzero_charge_power = 0.0;
    double poisson_residual_power = 0.0;
    double poisson_charge_scale_power = 0.0;
    double poisson_laplacian_scale_power = 0.0;
    double poisson_residual_maximum = 0.0;
    double poisson_scale_maximum = 0.0;
    double gradient_residual_power = 0.0;
    double gradient_electric_scale_power = 0.0;
    double gradient_potential_scale_power = 0.0;
    double gradient_residual_maximum = 0.0;
    double gradient_scale_maximum = 0.0;
    double gauss_residual_power = 0.0;
    double gauss_charge_scale_power = 0.0;
    double gauss_divergence_scale_power = 0.0;
    double gauss_residual_maximum = 0.0;
    double gauss_scale_maximum = 0.0;
    double electric_spectral_power = 0.0;
    swa_kahan charge_square = {0.0, 0.0};
    swa_kahan electric_square = {0.0, 0.0};
    swa_kahan charge_potential = {0.0, 0.0};
    double cell_width_m;
    double charge_scale;
    double electric_scale;
    double energy_scale;

    if (input == NULL || tolerances == NULL || out == NULL ||
        input->samples < 3U ||
        input->samples > (size_t)INT64_MAX ||
        !(input->length_m > 0.0) ||
        !(input->permittivity_F_m > 0.0) ||
        !isfinite(input->length_m) ||
        !isfinite(input->permittivity_F_m) ||
        !finite_array(input->charge_density_C_m3, input->samples) ||
        !finite_array(input->potential_V, input->samples) ||
        !finite_array(input->electric_field_V_m, input->samples) ||
        !valid_nonnegative(tolerances->operator_relative_tolerance) ||
        !valid_nonnegative(
            tolerances->charge_absolute_tolerance_C_m3) ||
        !valid_nonnegative(
            tolerances->electric_absolute_tolerance_V_m) ||
        !valid_nonnegative(tolerances->energy_relative_tolerance) ||
        !valid_nonnegative(
            tolerances->energy_absolute_tolerance_J_m2) ||
        !valid_nonnegative(
            tolerances->neutrality_relative_tolerance) ||
        !valid_nonnegative(
            tolerances->neutrality_absolute_tolerance_C_m3) ||
        !valid_nonnegative(
            tolerances->harmonic_field_relative_tolerance) ||
        !valid_nonnegative(
            tolerances->harmonic_field_absolute_tolerance_V_m) ||
        !valid_nonnegative(
            tolerances->nyquist_charge_power_fraction_tolerance) ||
        tolerances->nyquist_charge_power_fraction_tolerance > 1.0) {
        return 0;
    }

    memset(out, 0, sizeof(*out));
    out->samples = input->samples;
    out->nyquist_present = even_samples;
    cell_width_m = input->length_m / (double)input->samples;
    if (!(cell_width_m > 0.0) || !isfinite(cell_width_m)) {
        return 0;
    }

    for (sample = 0U; sample < input->samples; ++sample) {
        const double charge = input->charge_density_C_m3[sample];
        const double electric = input->electric_field_V_m[sample];
        swa_kahan_add(&charge_square, charge * charge);
        swa_kahan_add(&electric_square, electric * electric);
        swa_kahan_add(
            &charge_potential,
            charge * input->potential_V[sample]
        );
    }
    out->charge_rms_C_m3 =
        sqrt(charge_square.sum / (double)input->samples);
    out->electric_rms_V_m =
        sqrt(electric_square.sum / (double)input->samples);

    for (mode = 0U; mode < input->samples; ++mode) {
        double charge_real;
        double charge_imaginary;
        double potential_real;
        double potential_imaginary;
        double electric_real;
        double electric_imaginary;
        double electric_power;
        int64_t signed_mode;
        double wave_number;
        double poisson_real;
        double poisson_imaginary;
        double poisson_norm;
        double poisson_charge_norm;
        double poisson_laplacian_norm;
        int is_nyquist;

        dft_coefficient(
            input->charge_density_C_m3,
            input->samples,
            mode,
            &charge_real,
            &charge_imaginary
        );
        dft_coefficient(
            input->potential_V,
            input->samples,
            mode,
            &potential_real,
            &potential_imaginary
        );
        dft_coefficient(
            input->electric_field_V_m,
            input->samples,
            mode,
            &electric_real,
            &electric_imaginary
        );
        electric_power =
            electric_real * electric_real +
            electric_imaginary * electric_imaginary;
        electric_spectral_power += electric_power;

        if (mode == 0U) {
            charge_zero_real = charge_real;
            charge_zero_imaginary = charge_imaginary;
            electric_zero_real = electric_real;
            electric_zero_imaginary = electric_imaginary;
            continue;
        }

        signed_mode = mode <= input->samples / 2U ?
            (int64_t)mode :
            (int64_t)mode - (int64_t)input->samples;
        wave_number =
            2.0 * SWA_PI * (double)signed_mode /
            input->length_m;
        is_nyquist =
            even_samples && mode == input->samples / 2U;

        poisson_real =
            input->permittivity_F_m *
                wave_number * wave_number * potential_real -
            charge_real;
        poisson_imaginary =
            input->permittivity_F_m *
                wave_number * wave_number * potential_imaginary -
            charge_imaginary;
        poisson_norm = complex_norm(
            poisson_real, poisson_imaginary
        );
        poisson_charge_norm = complex_norm(
            charge_real, charge_imaginary
        );
        poisson_laplacian_norm = complex_norm(
            input->permittivity_F_m *
                wave_number * wave_number * potential_real,
            input->permittivity_F_m *
                wave_number * wave_number * potential_imaginary
        );
        poisson_residual_power += poisson_norm * poisson_norm;
        poisson_charge_scale_power +=
            poisson_charge_norm * poisson_charge_norm;
        poisson_laplacian_scale_power +=
            poisson_laplacian_norm * poisson_laplacian_norm;
        poisson_residual_maximum = max2(
            poisson_residual_maximum, poisson_norm
        );
        poisson_scale_maximum = max3(
            poisson_scale_maximum,
            poisson_charge_norm,
            poisson_laplacian_norm
        );
        nonzero_charge_power +=
            poisson_charge_norm * poisson_charge_norm;

        if (is_nyquist) {
            const double gradient_real =
                electric_real -
                wave_number * potential_imaginary;
            const double gradient_imaginary =
                electric_imaginary +
                wave_number * potential_real;
            const double gauss_real =
                -input->permittivity_F_m *
                    wave_number * electric_imaginary -
                charge_real;
            const double gauss_imaginary =
                input->permittivity_F_m *
                    wave_number * electric_real -
                charge_imaginary;
            nyquist_charge_power =
                poisson_charge_norm * poisson_charge_norm;
            out->nyquist_charge_rms_C_m3 =
                poisson_charge_norm / (double)input->samples;
            out->nyquist_gradient_residual_V_m =
                complex_norm(
                    gradient_real, gradient_imaginary
                ) / (double)input->samples;
            out->nyquist_gauss_residual_C_m3 =
                complex_norm(
                    gauss_real, gauss_imaginary
                ) / (double)input->samples;
            continue;
        }

        {
            const double potential_gradient_real =
                -wave_number * potential_imaginary;
            const double potential_gradient_imaginary =
                wave_number * potential_real;
            const double gradient_real =
                electric_real + potential_gradient_real;
            const double gradient_imaginary =
                electric_imaginary +
                potential_gradient_imaginary;
            const double gradient_norm = complex_norm(
                gradient_real, gradient_imaginary
            );
            const double electric_norm = complex_norm(
                electric_real, electric_imaginary
            );
            const double potential_gradient_norm = complex_norm(
                potential_gradient_real,
                potential_gradient_imaginary
            );
            const double gauss_divergence_real =
                -input->permittivity_F_m *
                    wave_number * electric_imaginary;
            const double gauss_divergence_imaginary =
                input->permittivity_F_m *
                    wave_number * electric_real;
            const double gauss_real =
                gauss_divergence_real - charge_real;
            const double gauss_imaginary =
                gauss_divergence_imaginary -
                charge_imaginary;
            const double gauss_norm = complex_norm(
                gauss_real, gauss_imaginary
            );
            const double gauss_divergence_norm = complex_norm(
                gauss_divergence_real,
                gauss_divergence_imaginary
            );

            gradient_residual_power +=
                gradient_norm * gradient_norm;
            gradient_electric_scale_power +=
                electric_norm * electric_norm;
            gradient_potential_scale_power +=
                potential_gradient_norm *
                potential_gradient_norm;
            gradient_residual_maximum = max2(
                gradient_residual_maximum,
                gradient_norm
            );
            gradient_scale_maximum = max3(
                gradient_scale_maximum,
                electric_norm,
                potential_gradient_norm
            );

            gauss_residual_power +=
                gauss_norm * gauss_norm;
            gauss_charge_scale_power +=
                poisson_charge_norm * poisson_charge_norm;
            gauss_divergence_scale_power +=
                gauss_divergence_norm *
                gauss_divergence_norm;
            gauss_residual_maximum = max2(
                gauss_residual_maximum, gauss_norm
            );
            gauss_scale_maximum = max3(
                gauss_scale_maximum,
                poisson_charge_norm,
                gauss_divergence_norm
            );
            ++out->resolved_nonzero_modes;
        }
    }

    out->charge_mean_C_m3 =
        charge_zero_real / (double)input->samples;
    out->electric_mean_V_m =
        electric_zero_real / (double)input->samples;
    out->charge_zero_mode_relative =
        complex_norm(
            charge_zero_real, charge_zero_imaginary
        ) /
        ((double)input->samples *
         max2(out->charge_rms_C_m3, DBL_MIN));
    out->harmonic_field_relative =
        complex_norm(
            electric_zero_real, electric_zero_imaginary
        ) /
        ((double)input->samples *
         max2(out->electric_rms_V_m, DBL_MIN));
    out->nyquist_charge_power_fraction =
        nonzero_charge_power > 0.0 ?
        nyquist_charge_power / nonzero_charge_power : 0.0;

    out->poisson_residual_rms_C_m3 =
        sqrt(poisson_residual_power) /
        (double)input->samples;
    out->poisson_residual_max_spectral_C_m3 =
        poisson_residual_maximum /
        (double)input->samples;
    charge_scale = max2(
        sqrt(poisson_charge_scale_power) /
            (double)input->samples,
        sqrt(poisson_laplacian_scale_power) /
            (double)input->samples
    );
    out->poisson_relative_rms =
        out->poisson_residual_rms_C_m3 /
        max2(charge_scale, DBL_MIN);
    out->poisson_relative_maximum =
        out->poisson_residual_max_spectral_C_m3 /
        max2(poisson_scale_maximum /
                 (double)input->samples,
             DBL_MIN);

    out->gradient_residual_rms_V_m =
        sqrt(gradient_residual_power) /
        (double)input->samples;
    out->gradient_residual_max_spectral_V_m =
        gradient_residual_maximum /
        (double)input->samples;
    electric_scale = max2(
        sqrt(gradient_electric_scale_power) /
            (double)input->samples,
        sqrt(gradient_potential_scale_power) /
            (double)input->samples
    );
    out->gradient_relative_rms =
        out->gradient_residual_rms_V_m /
        max2(electric_scale, DBL_MIN);
    out->gradient_relative_maximum =
        out->gradient_residual_max_spectral_V_m /
        max2(gradient_scale_maximum /
                 (double)input->samples,
             DBL_MIN);

    out->gauss_residual_rms_C_m3 =
        sqrt(gauss_residual_power) /
        (double)input->samples;
    out->gauss_residual_max_spectral_C_m3 =
        gauss_residual_maximum /
        (double)input->samples;
    charge_scale = max2(
        sqrt(gauss_charge_scale_power) /
            (double)input->samples,
        sqrt(gauss_divergence_scale_power) /
            (double)input->samples
    );
    out->gauss_relative_rms =
        out->gauss_residual_rms_C_m3 /
        max2(charge_scale, DBL_MIN);
    out->gauss_relative_maximum =
        out->gauss_residual_max_spectral_C_m3 /
        max2(gauss_scale_maximum /
                 (double)input->samples,
             DBL_MIN);

    out->field_energy_sample_J_m2 =
        0.5 * input->permittivity_F_m *
        cell_width_m * electric_square.sum;
    out->field_energy_parseval_J_m2 =
        0.5 * input->permittivity_F_m *
        cell_width_m * electric_spectral_power /
        (double)input->samples;
    out->charge_potential_energy_J_m2 =
        0.5 * cell_width_m * charge_potential.sum;
    energy_scale = max3(
        fabs(out->field_energy_sample_J_m2),
        fabs(out->field_energy_parseval_J_m2),
        DBL_MIN
    );
    out->parseval_energy_relative_error =
        fabs(out->field_energy_sample_J_m2 -
             out->field_energy_parseval_J_m2) /
        energy_scale;
    energy_scale = max3(
        fabs(out->field_energy_sample_J_m2),
        fabs(out->charge_potential_energy_J_m2),
        DBL_MIN
    );
    out->electrostatic_energy_relative_error =
        fabs(out->field_energy_sample_J_m2 -
             out->charge_potential_energy_J_m2) /
        energy_scale;

    out->finite =
        isfinite(out->charge_mean_C_m3) &&
        isfinite(out->charge_rms_C_m3) &&
        isfinite(out->electric_mean_V_m) &&
        isfinite(out->electric_rms_V_m) &&
        isfinite(out->charge_zero_mode_relative) &&
        isfinite(out->harmonic_field_relative) &&
        isfinite(out->nyquist_charge_rms_C_m3) &&
        isfinite(out->nyquist_charge_power_fraction) &&
        isfinite(out->nyquist_gradient_residual_V_m) &&
        isfinite(out->nyquist_gauss_residual_C_m3) &&
        isfinite(out->poisson_residual_rms_C_m3) &&
        isfinite(out->poisson_residual_max_spectral_C_m3) &&
        isfinite(out->poisson_relative_rms) &&
        isfinite(out->poisson_relative_maximum) &&
        isfinite(out->gradient_residual_rms_V_m) &&
        isfinite(out->gradient_residual_max_spectral_V_m) &&
        isfinite(out->gradient_relative_rms) &&
        isfinite(out->gradient_relative_maximum) &&
        isfinite(out->gauss_residual_rms_C_m3) &&
        isfinite(out->gauss_residual_max_spectral_C_m3) &&
        isfinite(out->gauss_relative_rms) &&
        isfinite(out->gauss_relative_maximum) &&
        isfinite(out->field_energy_sample_J_m2) &&
        isfinite(out->field_energy_parseval_J_m2) &&
        isfinite(out->charge_potential_energy_J_m2) &&
        isfinite(out->parseval_energy_relative_error) &&
        isfinite(out->electrostatic_energy_relative_error);

    out->neutrality_passes = residual_passes(
        fabs(out->charge_mean_C_m3),
        out->charge_zero_mode_relative,
        tolerances->neutrality_absolute_tolerance_C_m3,
        tolerances->neutrality_relative_tolerance
    );
    out->harmonic_field_passes = residual_passes(
        fabs(out->electric_mean_V_m),
        out->harmonic_field_relative,
        tolerances->harmonic_field_absolute_tolerance_V_m,
        tolerances->harmonic_field_relative_tolerance
    );
    out->nyquist_passes =
        !out->nyquist_present ||
        out->nyquist_charge_power_fraction <=
            tolerances->nyquist_charge_power_fraction_tolerance;
    out->poisson_passes =
        residual_passes(
            out->poisson_residual_rms_C_m3,
            out->poisson_relative_rms,
            tolerances->charge_absolute_tolerance_C_m3,
            tolerances->operator_relative_tolerance
        ) &&
        residual_passes(
            out->poisson_residual_max_spectral_C_m3,
            out->poisson_relative_maximum,
            tolerances->charge_absolute_tolerance_C_m3,
            tolerances->operator_relative_tolerance
        );
    out->gradient_passes =
        residual_passes(
            out->gradient_residual_rms_V_m,
            out->gradient_relative_rms,
            tolerances->electric_absolute_tolerance_V_m,
            tolerances->operator_relative_tolerance
        ) &&
        residual_passes(
            out->gradient_residual_max_spectral_V_m,
            out->gradient_relative_maximum,
            tolerances->electric_absolute_tolerance_V_m,
            tolerances->operator_relative_tolerance
        );
    out->gauss_passes =
        residual_passes(
            out->gauss_residual_rms_C_m3,
            out->gauss_relative_rms,
            tolerances->charge_absolute_tolerance_C_m3,
            tolerances->operator_relative_tolerance
        ) &&
        residual_passes(
            out->gauss_residual_max_spectral_C_m3,
            out->gauss_relative_maximum,
            tolerances->charge_absolute_tolerance_C_m3,
            tolerances->operator_relative_tolerance
        );
    out->parseval_passes = residual_passes(
        fabs(out->field_energy_sample_J_m2 -
             out->field_energy_parseval_J_m2),
        out->parseval_energy_relative_error,
        tolerances->energy_absolute_tolerance_J_m2,
        tolerances->energy_relative_tolerance
    );
    out->electrostatic_energy_passes = residual_passes(
        fabs(out->field_energy_sample_J_m2 -
             out->charge_potential_energy_J_m2),
        out->electrostatic_energy_relative_error,
        tolerances->energy_absolute_tolerance_J_m2,
        tolerances->energy_relative_tolerance
    );
    out->passes = out->finite &&
                  out->neutrality_passes &&
                  out->harmonic_field_passes &&
                  out->nyquist_passes &&
                  out->poisson_passes &&
                  out->gradient_passes &&
                  out->gauss_passes &&
                  out->parseval_passes &&
                  out->electrostatic_energy_passes;
    return 1;
}

int swa_write_spectral_poisson_receipt(
    FILE *fp,
    const swa_spectral_poisson_input *input,
    const swa_spectral_poisson_tolerances *tolerances,
    const swa_spectral_poisson_result *result
) {
    if (fp == NULL || input == NULL ||
        tolerances == NULL || result == NULL) {
        return 0;
    }
    return fprintf(
        fp,
        "{\n"
        "  \"schema\": \"spacewind.spectral-poisson-audit/v1\",\n"
        "  \"samples\": %zu,\n"
        "  \"length_m\": %.17g,\n"
        "  \"permittivity_F_m\": %.17g,\n"
        "  \"resolved_nonzero_modes\": %zu,\n"
        "  \"charge_mean_C_m3\": %.17g,\n"
        "  \"charge_zero_mode_relative\": %.17g,\n"
        "  \"electric_mean_V_m\": %.17g,\n"
        "  \"harmonic_field_relative\": %.17g,\n"
        "  \"nyquist_present\": %s,\n"
        "  \"nyquist_charge_power_fraction\": %.17g,\n"
        "  \"nyquist_gradient_residual_V_m\": %.17g,\n"
        "  \"poisson_relative_rms\": %.17g,\n"
        "  \"gradient_relative_rms\": %.17g,\n"
        "  \"gauss_relative_rms\": %.17g,\n"
        "  \"field_energy_sample_J_m2\": %.17g,\n"
        "  \"field_energy_parseval_J_m2\": %.17g,\n"
        "  \"charge_potential_energy_J_m2\": %.17g,\n"
        "  \"parseval_energy_relative_error\": %.17g,\n"
        "  \"electrostatic_energy_relative_error\": %.17g,\n"
        "  \"operator_relative_tolerance\": %.17g,\n"
        "  \"nyquist_charge_power_fraction_tolerance\": %.17g,\n"
        "  \"neutrality_passes\": %s,\n"
        "  \"harmonic_field_passes\": %s,\n"
        "  \"nyquist_passes\": %s,\n"
        "  \"poisson_passes\": %s,\n"
        "  \"gradient_passes\": %s,\n"
        "  \"gauss_passes\": %s,\n"
        "  \"parseval_passes\": %s,\n"
        "  \"electrostatic_energy_passes\": %s,\n"
        "  \"passes\": %s,\n"
        "  \"nonclaim\": \"spectral residual closure verifies the declared periodic discrete operator; unresolved mean or Nyquist content, plasma-model fidelity, and propulsion performance remain separate questions\"\n"
        "}\n",
        input->samples,
        input->length_m,
        input->permittivity_F_m,
        result->resolved_nonzero_modes,
        result->charge_mean_C_m3,
        result->charge_zero_mode_relative,
        result->electric_mean_V_m,
        result->harmonic_field_relative,
        result->nyquist_present ? "true" : "false",
        result->nyquist_charge_power_fraction,
        result->nyquist_gradient_residual_V_m,
        result->poisson_relative_rms,
        result->gradient_relative_rms,
        result->gauss_relative_rms,
        result->field_energy_sample_J_m2,
        result->field_energy_parseval_J_m2,
        result->charge_potential_energy_J_m2,
        result->parseval_energy_relative_error,
        result->electrostatic_energy_relative_error,
        tolerances->operator_relative_tolerance,
        tolerances->nyquist_charge_power_fraction_tolerance,
        result->neutrality_passes ? "true" : "false",
        result->harmonic_field_passes ? "true" : "false",
        result->nyquist_passes ? "true" : "false",
        result->poisson_passes ? "true" : "false",
        result->gradient_passes ? "true" : "false",
        result->gauss_passes ? "true" : "false",
        result->parseval_passes ? "true" : "false",
        result->electrostatic_energy_passes ? "true" : "false",
        result->passes ? "true" : "false"
    ) > 0;
}
