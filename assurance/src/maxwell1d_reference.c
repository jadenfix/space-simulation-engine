#include "spacewind/maxwell1d_reference.h"

#include <complex.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
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

static int valid_state_shape(const swa_maxwell1d_state *state) {
    return state != NULL && state->samples >= 5U &&
           (state->samples & 1U) == 1U &&
           isfinite(state->length_m) && state->length_m > 0.0 &&
           isfinite(state->cross_section_area_m2) &&
           state->cross_section_area_m2 > 0.0 &&
           state->electric_y_V_m != NULL &&
           state->electric_z_V_m != NULL &&
           state->magnetic_y_T != NULL &&
           state->magnetic_z_T != NULL;
}

static int pointers_alias_state(
    const swa_maxwell1d_state *state,
    const double *values
) {
    return values != NULL &&
           (values == state->electric_y_V_m ||
            values == state->electric_z_V_m ||
            values == state->magnetic_y_T ||
            values == state->magnetic_z_T);
}

int swa_maxwell1d_init(
    swa_maxwell1d_state *state,
    size_t samples,
    double length_m,
    double cross_section_area_m2
) {
    if (state == NULL || samples < 5U || (samples & 1U) == 0U ||
        !isfinite(length_m) || !(length_m > 0.0) ||
        !isfinite(cross_section_area_m2) ||
        !(cross_section_area_m2 > 0.0) ||
        samples > SIZE_MAX / sizeof(double)) {
        return 0;
    }
    memset(state, 0, sizeof(*state));
    state->samples = samples;
    state->length_m = length_m;
    state->cross_section_area_m2 = cross_section_area_m2;
    state->electric_y_V_m = calloc(samples, sizeof(double));
    state->electric_z_V_m = calloc(samples, sizeof(double));
    state->magnetic_y_T = calloc(samples, sizeof(double));
    state->magnetic_z_T = calloc(samples, sizeof(double));
    if (state->electric_y_V_m == NULL ||
        state->electric_z_V_m == NULL ||
        state->magnetic_y_T == NULL ||
        state->magnetic_z_T == NULL) {
        swa_maxwell1d_destroy(state);
        return 0;
    }
    return 1;
}

void swa_maxwell1d_destroy(swa_maxwell1d_state *state) {
    if (state == NULL) {
        return;
    }
    free(state->electric_y_V_m);
    free(state->electric_z_V_m);
    free(state->magnetic_y_T);
    free(state->magnetic_z_T);
    memset(state, 0, sizeof(*state));
}

int swa_maxwell1d_copy(
    const swa_maxwell1d_state *source,
    swa_maxwell1d_state *destination
) {
    size_t bytes;
    if (!swa_maxwell1d_state_is_finite(source) ||
        destination == NULL || destination == source ||
        !swa_maxwell1d_init(destination, source->samples,
                            source->length_m,
                            source->cross_section_area_m2)) {
        return 0;
    }
    bytes = source->samples * sizeof(double);
    memcpy(destination->electric_y_V_m, source->electric_y_V_m, bytes);
    memcpy(destination->electric_z_V_m, source->electric_z_V_m, bytes);
    memcpy(destination->magnetic_y_T, source->magnetic_y_T, bytes);
    memcpy(destination->magnetic_z_T, source->magnetic_z_T, bytes);
    return 1;
}

int swa_maxwell1d_state_is_finite(const swa_maxwell1d_state *state) {
    return valid_state_shape(state) &&
           finite_array(state->electric_y_V_m, state->samples) &&
           finite_array(state->electric_z_V_m, state->samples) &&
           finite_array(state->magnetic_y_T, state->samples) &&
           finite_array(state->magnetic_z_T, state->samples);
}

static int diagnostics(
    const swa_maxwell1d_state *state,
    double *energy_J,
    double *momentum_Ns,
    double *right_energy_J,
    double *left_energy_J
) {
    const double dx = state->length_m / (double)state->samples;
    const double volume = state->cross_section_area_m2 * dx;
    swa_kahan energy = {0.0, 0.0};
    swa_kahan momentum = {0.0, 0.0};
    swa_kahan right = {0.0, 0.0};
    swa_kahan left = {0.0, 0.0};
    size_t i;
    if (!swa_maxwell1d_state_is_finite(state) ||
        energy_J == NULL || momentum_Ns == NULL ||
        right_energy_J == NULL || left_energy_J == NULL ||
        !isfinite(dx) || !(dx > 0.0) ||
        !isfinite(volume) || !(volume > 0.0)) {
        return 0;
    }
    for (i = 0U; i < state->samples; ++i) {
        const double ey = state->electric_y_V_m[i];
        const double ez = state->electric_z_V_m[i];
        const double by = state->magnetic_y_T[i];
        const double bz = state->magnetic_z_T[i];
        const double cby = SWA_C * by;
        const double cbz = SWA_C * bz;
        const double density =
            0.5 * SWA_EPS0 * (ey * ey + ez * ez) +
            0.5 / SWA_MU0 * (by * by + bz * bz);
        const double poynting_x = (ey * bz - ez * by) / SWA_MU0;
        const double right_density = 0.25 * SWA_EPS0 *
            ((ey + cbz) * (ey + cbz) +
             (ez - cby) * (ez - cby));
        const double left_density = 0.25 * SWA_EPS0 *
            ((ey - cbz) * (ey - cbz) +
             (ez + cby) * (ez + cby));
        swa_kahan_add(&energy, density * volume);
        swa_kahan_add(&momentum,
                      poynting_x * volume / (SWA_C * SWA_C));
        swa_kahan_add(&right, right_density * volume);
        swa_kahan_add(&left, left_density * volume);
    }
    *energy_J = energy.sum;
    *momentum_Ns = momentum.sum;
    *right_energy_J = right.sum;
    *left_energy_J = left.sum;
    return isfinite(*energy_J) && isfinite(*momentum_Ns) &&
           isfinite(*right_energy_J) && isfinite(*left_energy_J);
}

double swa_maxwell1d_field_energy_J(const swa_maxwell1d_state *state) {
    double energy = NAN;
    double momentum;
    double right;
    double left;
    if (!diagnostics(state, &energy, &momentum, &right, &left)) {
        return NAN;
    }
    return energy;
}

double swa_maxwell1d_field_momentum_Ns(const swa_maxwell1d_state *state) {
    double energy;
    double momentum = NAN;
    double right;
    double left;
    if (!diagnostics(state, &energy, &momentum, &right, &left)) {
        return NAN;
    }
    return momentum;
}

static int spectral_shift_real(
    const double *input,
    size_t samples,
    double shift_fraction_of_domain,
    double *output,
    double *maximum_imaginary
) {
    double complex *coefficient;
    size_t k;
    size_t j;
    double imaginary_maximum = 0.0;
    if (input == NULL || output == NULL || maximum_imaginary == NULL ||
        samples < 5U || (samples & 1U) == 0U ||
        !isfinite(shift_fraction_of_domain) ||
        !finite_array(input, samples)) {
        return 0;
    }
    coefficient = calloc(samples, sizeof(*coefficient));
    if (coefficient == NULL) {
        return 0;
    }
    for (k = 0U; k < samples; ++k) {
        double complex sum = 0.0 + 0.0 * I;
        const long mode = k <= samples / 2U ?
            (long)k : (long)k - (long)samples;
        const double phase = -2.0 * SWA_PI *
                             (double)mode *
                             shift_fraction_of_domain;
        const double complex translation = cos(phase) + I * sin(phase);
        for (j = 0U; j < samples; ++j) {
            const double angle = -2.0 * SWA_PI *
                                 (double)(k * j) /
                                 (double)samples;
            sum += input[j] * (cos(angle) + I * sin(angle));
        }
        coefficient[k] = sum * translation;
    }
    for (j = 0U; j < samples; ++j) {
        double complex sum = 0.0 + 0.0 * I;
        for (k = 0U; k < samples; ++k) {
            const double angle = 2.0 * SWA_PI *
                                 (double)(k * j) /
                                 (double)samples;
            sum += coefficient[k] * (cos(angle) + I * sin(angle));
        }
        sum /= (double)samples;
        output[j] = creal(sum);
        imaginary_maximum = max2(imaginary_maximum, fabs(cimag(sum)));
    }
    free(coefficient);
    *maximum_imaginary = imaginary_maximum;
    return finite_array(output, samples) && isfinite(imaginary_maximum);
}

int swa_maxwell1d_propagate_vacuum(
    swa_maxwell1d_state *state,
    double dt_s,
    double *spectral_imaginary_maximum
) {
    const size_t n = state != NULL ? state->samples : 0U;
    double *right_y;
    double *left_y;
    double *right_z;
    double *left_z;
    double *shifted_right_y;
    double *shifted_left_y;
    double *shifted_right_z;
    double *shifted_left_z;
    double max_imaginary = 0.0;
    double local_imaginary;
    double shift_fraction;
    size_t i;
    if (!swa_maxwell1d_state_is_finite(state) ||
        spectral_imaginary_maximum == NULL || !isfinite(dt_s) ||
        n > SIZE_MAX / (8U * sizeof(double))) {
        return 0;
    }
    right_y = calloc(8U * n, sizeof(double));
    if (right_y == NULL) {
        return 0;
    }
    left_y = right_y + n;
    right_z = left_y + n;
    left_z = right_z + n;
    shifted_right_y = left_z + n;
    shifted_left_y = shifted_right_y + n;
    shifted_right_z = shifted_left_y + n;
    shifted_left_z = shifted_right_z + n;

    for (i = 0U; i < n; ++i) {
        right_y[i] = state->electric_y_V_m[i] +
                     SWA_C * state->magnetic_z_T[i];
        left_y[i] = state->electric_y_V_m[i] -
                    SWA_C * state->magnetic_z_T[i];
        right_z[i] = state->electric_z_V_m[i] -
                     SWA_C * state->magnetic_y_T[i];
        left_z[i] = state->electric_z_V_m[i] +
                    SWA_C * state->magnetic_y_T[i];
    }
    shift_fraction = SWA_C * dt_s / state->length_m;
    if (!isfinite(shift_fraction) ||
        !spectral_shift_real(right_y, n, shift_fraction,
                             shifted_right_y, &local_imaginary)) {
        free(right_y);
        return 0;
    }
    max_imaginary = max2(max_imaginary, local_imaginary);
    if (!spectral_shift_real(left_y, n, -shift_fraction,
                             shifted_left_y, &local_imaginary)) {
        free(right_y);
        return 0;
    }
    max_imaginary = max2(max_imaginary, local_imaginary);
    if (!spectral_shift_real(right_z, n, shift_fraction,
                             shifted_right_z, &local_imaginary)) {
        free(right_y);
        return 0;
    }
    max_imaginary = max2(max_imaginary, local_imaginary);
    if (!spectral_shift_real(left_z, n, -shift_fraction,
                             shifted_left_z, &local_imaginary)) {
        free(right_y);
        return 0;
    }
    max_imaginary = max2(max_imaginary, local_imaginary);

    for (i = 0U; i < n; ++i) {
        state->electric_y_V_m[i] =
            0.5 * (shifted_right_y[i] + shifted_left_y[i]);
        state->magnetic_z_T[i] =
            0.5 * (shifted_right_y[i] - shifted_left_y[i]) / SWA_C;
        state->electric_z_V_m[i] =
            0.5 * (shifted_right_z[i] + shifted_left_z[i]);
        state->magnetic_y_T[i] =
            0.5 * (shifted_left_z[i] - shifted_right_z[i]) / SWA_C;
    }
    free(right_y);
    *spectral_imaginary_maximum = max_imaginary;
    return swa_maxwell1d_state_is_finite(state);
}

static int source_half_step(
    swa_maxwell1d_state *state,
    const double *current_y_A_m2,
    const double *current_z_A_m2,
    double half_dt_s,
    double *current_work_J
) {
    const double dx = state->length_m / (double)state->samples;
    const double volume = state->cross_section_area_m2 * dx;
    swa_kahan work = {0.0, 0.0};
    size_t i;
    if (current_work_J == NULL || !isfinite(half_dt_s) ||
        !isfinite(volume) || !(volume > 0.0)) {
        return 0;
    }
    for (i = 0U; i < state->samples; ++i) {
        const double jy = current_y_A_m2 == NULL ? 0.0 : current_y_A_m2[i];
        const double jz = current_z_A_m2 == NULL ? 0.0 : current_z_A_m2[i];
        const double ey_initial = state->electric_y_V_m[i];
        const double ez_initial = state->electric_z_V_m[i];
        double ey_final;
        double ez_final;
        if (!isfinite(jy) || !isfinite(jz)) {
            return 0;
        }
        ey_final = ey_initial - half_dt_s * jy / SWA_EPS0;
        ez_final = ez_initial - half_dt_s * jz / SWA_EPS0;
        swa_kahan_add(
            &work,
            half_dt_s * volume *
            (jy * 0.5 * (ey_initial + ey_final) +
             jz * 0.5 * (ez_initial + ez_final))
        );
        state->electric_y_V_m[i] = ey_final;
        state->electric_z_V_m[i] = ez_final;
    }
    *current_work_J = work.sum;
    return swa_maxwell1d_state_is_finite(state) && isfinite(work.sum);
}

int swa_maxwell1d_audit_energy_balance(
    double field_energy_initial_J,
    double field_energy_final_J,
    double current_work_J,
    double relative_tolerance,
    double absolute_tolerance_J,
    double *residual_J,
    double *relative_residual,
    int *passes
) {
    double scale;
    if (!isfinite(field_energy_initial_J) ||
        !isfinite(field_energy_final_J) ||
        !isfinite(current_work_J) ||
        !isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
        !isfinite(absolute_tolerance_J) || absolute_tolerance_J < 0.0 ||
        residual_J == NULL || relative_residual == NULL || passes == NULL) {
        return 0;
    }
    *residual_J = field_energy_final_J -
                  field_energy_initial_J + current_work_J;
    scale = max3(fabs(field_energy_initial_J),
                 fabs(field_energy_final_J),
                 max2(fabs(current_work_J), DBL_MIN));
    *relative_residual = fabs(*residual_J) / scale;
    *passes = fabs(*residual_J) <= absolute_tolerance_J ||
              *relative_residual <= relative_tolerance;
    return isfinite(*residual_J) && isfinite(*relative_residual);
}

int swa_maxwell1d_advance(
    swa_maxwell1d_state *state,
    const double *current_y_A_m2,
    const double *current_z_A_m2,
    double dt_s,
    double relative_tolerance,
    double absolute_energy_tolerance_J,
    double absolute_momentum_tolerance_Ns,
    swa_maxwell1d_step_result *out
) {
    double first_work;
    double second_work;
    double momentum_scale;
    if (!swa_maxwell1d_state_is_finite(state) || out == NULL ||
        !isfinite(dt_s) ||
        !isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
        !isfinite(absolute_energy_tolerance_J) ||
        absolute_energy_tolerance_J < 0.0 ||
        !isfinite(absolute_momentum_tolerance_Ns) ||
        absolute_momentum_tolerance_Ns < 0.0 ||
        pointers_alias_state(state, current_y_A_m2) ||
        pointers_alias_state(state, current_z_A_m2) ||
        (current_y_A_m2 != NULL &&
         !finite_array(current_y_A_m2, state->samples)) ||
        (current_z_A_m2 != NULL &&
         !finite_array(current_z_A_m2, state->samples))) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    out->vacuum_step = current_y_A_m2 == NULL && current_z_A_m2 == NULL;
    if (!diagnostics(state,
                     &out->field_energy_initial_J,
                     &out->field_momentum_initial_Ns,
                     &out->right_moving_energy_initial_J,
                     &out->left_moving_energy_initial_J) ||
        !source_half_step(state, current_y_A_m2, current_z_A_m2,
                          0.5 * dt_s, &first_work) ||
        !swa_maxwell1d_propagate_vacuum(
            state, dt_s, &out->spectral_imaginary_maximum
        ) ||
        !source_half_step(state, current_y_A_m2, current_z_A_m2,
                          0.5 * dt_s, &second_work) ||
        !diagnostics(state,
                     &out->field_energy_final_J,
                     &out->field_momentum_final_Ns,
                     &out->right_moving_energy_final_J,
                     &out->left_moving_energy_final_J)) {
        return 0;
    }
    out->current_work_J = first_work + second_work;
    if (!swa_maxwell1d_audit_energy_balance(
            out->field_energy_initial_J,
            out->field_energy_final_J,
            out->current_work_J,
            relative_tolerance,
            absolute_energy_tolerance_J,
            &out->energy_residual_J,
            &out->energy_relative_residual,
            &out->energy_closes)) {
        return 0;
    }
    out->vacuum_momentum_residual_Ns =
        out->field_momentum_final_Ns -
        out->field_momentum_initial_Ns;
    momentum_scale = max3(
        fabs(out->field_momentum_initial_Ns),
        fabs(out->field_momentum_final_Ns),
        DBL_MIN
    );
    out->vacuum_momentum_relative_residual =
        fabs(out->vacuum_momentum_residual_Ns) / momentum_scale;
    out->vacuum_momentum_closes =
        !out->vacuum_step ||
        fabs(out->vacuum_momentum_residual_Ns) <=
            absolute_momentum_tolerance_Ns ||
        out->vacuum_momentum_relative_residual <= relative_tolerance;
    out->finite =
        swa_maxwell1d_state_is_finite(state) &&
        isfinite(out->field_energy_initial_J) &&
        isfinite(out->field_energy_final_J) &&
        isfinite(out->field_momentum_initial_Ns) &&
        isfinite(out->field_momentum_final_Ns) &&
        isfinite(out->current_work_J) &&
        isfinite(out->energy_relative_residual) &&
        isfinite(out->vacuum_momentum_relative_residual) &&
        isfinite(out->spectral_imaginary_maximum);
    out->passes = out->finite && out->energy_closes &&
                  out->vacuum_momentum_closes;
    return 1;
}

int swa_write_maxwell1d_reference_receipt(
    FILE *fp,
    const swa_maxwell1d_state *state,
    const swa_maxwell1d_step_result *result
) {
    if (fp == NULL || state == NULL || result == NULL) {
        return 0;
    }
    return fprintf(
        fp,
        "{\n"
        "  \"schema\": \"spacewind.maxwell1d-reference/v1\",\n"
        "  \"samples\": %zu,\n"
        "  \"length_m\": %.17g,\n"
        "  \"cross_section_area_m2\": %.17g,\n"
        "  \"vacuum_step\": %s,\n"
        "  \"field_energy_initial_J\": %.17g,\n"
        "  \"field_energy_final_J\": %.17g,\n"
        "  \"right_moving_energy_initial_J\": %.17g,\n"
        "  \"right_moving_energy_final_J\": %.17g,\n"
        "  \"left_moving_energy_initial_J\": %.17g,\n"
        "  \"left_moving_energy_final_J\": %.17g,\n"
        "  \"field_momentum_initial_Ns\": %.17g,\n"
        "  \"field_momentum_final_Ns\": %.17g,\n"
        "  \"current_work_J\": %.17g,\n"
        "  \"energy_residual_J\": %.17g,\n"
        "  \"energy_relative_residual\": %.17g,\n"
        "  \"vacuum_momentum_residual_Ns\": %.17g,\n"
        "  \"vacuum_momentum_relative_residual\": %.17g,\n"
        "  \"spectral_imaginary_maximum\": %.17g,\n"
        "  \"energy_closes\": %s,\n"
        "  \"vacuum_momentum_closes\": %s,\n"
        "  \"passes\": %s,\n"
        "  \"nonclaim\": \"this direct-DFT periodic transverse Maxwell reference is an independent numerical oracle; prescribed current is not a self-consistent particle current and the result is not plasma-wing thrust evidence\"\n"
        "}\n",
        state->samples,
        state->length_m,
        state->cross_section_area_m2,
        result->vacuum_step ? "true" : "false",
        result->field_energy_initial_J,
        result->field_energy_final_J,
        result->right_moving_energy_initial_J,
        result->right_moving_energy_final_J,
        result->left_moving_energy_initial_J,
        result->left_moving_energy_final_J,
        result->field_momentum_initial_Ns,
        result->field_momentum_final_Ns,
        result->current_work_J,
        result->energy_residual_J,
        result->energy_relative_residual,
        result->vacuum_momentum_residual_Ns,
        result->vacuum_momentum_relative_residual,
        result->spectral_imaginary_maximum,
        result->energy_closes ? "true" : "false",
        result->vacuum_momentum_closes ? "true" : "false",
        result->passes ? "true" : "false"
    ) > 0;
}
