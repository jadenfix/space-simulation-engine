#include "spacewind/charge_deposition.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static double max2(double a, double b) { return a > b ? a : b; }

static int valid_nonnegative(double value) {
    return isfinite(value) && value >= 0.0;
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

static double wrap_periodic(double position_m, double length_m) {
    double wrapped = fmod(position_m, length_m);
    if (wrapped < 0.0) {
        wrapped += length_m;
    }
    if (wrapped >= length_m) {
        wrapped = 0.0;
    }
    return wrapped;
}

static void compensated_add(double *sum, double *correction, double value) {
    const double y = value - *correction;
    const double t = *sum + y;
    *correction = (t - *sum) - y;
    *sum = t;
}

static int deposit_cic(
    double position_unwrapped_m,
    double charge_C,
    size_t cells,
    double domain_length_m,
    double cell_width_m,
    double inverse_cell_volume_m3,
    double *density_C_m3,
    double *compensation_C_m3
) {
    const double wrapped = wrap_periodic(position_unwrapped_m,
                                         domain_length_m);
    const double coordinate = wrapped / cell_width_m - 0.5;
    const double lower_real = floor(coordinate);
    const double upper_weight = coordinate - lower_real;
    const double lower_weight = 1.0 - upper_weight;
    int64_t lower_signed;
    size_t lower;
    size_t upper;
    double scale;
    if (!isfinite(wrapped) || !isfinite(coordinate) ||
        lower_real < -1.0 ||
        lower_real > (double)(cells - 1U)) {
        return 0;
    }
    lower_signed = (int64_t)lower_real;
    lower = lower_signed < 0 ? cells - 1U : (size_t)lower_signed;
    upper = lower + 1U;
    if (upper == cells) {
        upper = 0U;
    }
    scale = charge_C * inverse_cell_volume_m3;
    compensated_add(&density_C_m3[lower],
                    &compensation_C_m3[lower],
                    scale * lower_weight);
    compensated_add(&density_C_m3[upper],
                    &compensation_C_m3[upper],
                    scale * upper_weight);
    return 1;
}

int swa_deposit_periodic_cic_transport(
    const swa_periodic_cic_transport_input *input,
    const swa_periodic_cic_transport_buffers *buffers,
    swa_periodic_cic_transport_result *out
) {
    size_t i;
    double cell_width_m;
    double cell_volume_m3;
    double inverse_cell_volume_m3;
    double *initial_compensation;
    double *final_compensation;
    double *jy_zero;
    swa_kahan particle_charge = {0.0, 0.0};
    swa_kahan transported_charge_turns = {0.0, 0.0};
    swa_kahan deposited_initial = {0.0, 0.0};
    swa_kahan deposited_final = {0.0, 0.0};
    swa_kahan raw_current_mean = {0.0, 0.0};
    swa_kahan deposited_current_mean = {0.0, 0.0};
    double current_shift;
    double current_scale;
    double charge_scale;
    swa_charge_continuity_input continuity_input;

    if (input == NULL || buffers == NULL || out == NULL ||
        input->cells < 2U || input->particles == 0U ||
        input->cells > (size_t)INT64_MAX ||
        !(input->domain_length_m > 0.0) ||
        !(input->cross_section_area_m2 > 0.0) ||
        !(input->dt_s > 0.0) ||
        !isfinite(input->domain_length_m) ||
        !isfinite(input->cross_section_area_m2) ||
        !isfinite(input->dt_s) ||
        !valid_nonnegative(input->continuity_relative_tolerance) ||
        !valid_nonnegative(input->continuity_absolute_tolerance_A_m3) ||
        !valid_nonnegative(
            input->continuity_global_absolute_tolerance_C_m) ||
        !valid_nonnegative(input->charge_absolute_tolerance_C) ||
        !valid_nonnegative(input->mean_current_relative_tolerance) ||
        !valid_nonnegative(
            input->mean_current_absolute_tolerance_A_m2) ||
        buffers->rho_initial_C_m3 == NULL ||
        buffers->rho_final_C_m3 == NULL ||
        buffers->current_face_A_m2 == NULL ||
        buffers->rho_capacity < input->cells ||
        buffers->current_face_capacity < input->cells + 1U ||
        !finite_array(input->initial_position_unwrapped_m,
                      input->particles) ||
        !finite_array(input->final_position_unwrapped_m,
                      input->particles) ||
        !finite_array(input->particle_charge_C,
                      input->particles) ||
        input->cells > SIZE_MAX / 2U) {
        return 0;
    }

    cell_width_m = input->domain_length_m / (double)input->cells;
    cell_volume_m3 = cell_width_m * input->cross_section_area_m2;
    if (!(cell_width_m > 0.0) || !(cell_volume_m3 > 0.0) ||
        !isfinite(cell_width_m) || !isfinite(cell_volume_m3)) {
        return 0;
    }
    inverse_cell_volume_m3 = 1.0 / cell_volume_m3;
    if (!isfinite(inverse_cell_volume_m3)) {
        return 0;
    }

    initial_compensation = calloc(input->cells, sizeof(double));
    final_compensation = calloc(input->cells, sizeof(double));
    jy_zero = calloc(2U * input->cells, sizeof(double));
    if (initial_compensation == NULL ||
        final_compensation == NULL ||
        jy_zero == NULL) {
        free(initial_compensation);
        free(final_compensation);
        free(jy_zero);
        return 0;
    }

    memset(out, 0, sizeof(*out));
    memset(buffers->rho_initial_C_m3, 0,
           input->cells * sizeof(double));
    memset(buffers->rho_final_C_m3, 0,
           input->cells * sizeof(double));
    memset(buffers->current_face_A_m2, 0,
           (input->cells + 1U) * sizeof(double));

    for (i = 0U; i < input->particles; ++i) {
        const double displacement_m =
            input->final_position_unwrapped_m[i] -
            input->initial_position_unwrapped_m[i];
        const double charge_C = input->particle_charge_C[i];
        if (!isfinite(displacement_m) ||
            !deposit_cic(input->initial_position_unwrapped_m[i],
                         charge_C, input->cells,
                         input->domain_length_m, cell_width_m,
                         inverse_cell_volume_m3,
                         buffers->rho_initial_C_m3,
                         initial_compensation) ||
            !deposit_cic(input->final_position_unwrapped_m[i],
                         charge_C, input->cells,
                         input->domain_length_m, cell_width_m,
                         inverse_cell_volume_m3,
                         buffers->rho_final_C_m3,
                         final_compensation)) {
            free(initial_compensation);
            free(final_compensation);
            free(jy_zero);
            return 0;
        }
        swa_kahan_add(&particle_charge, charge_C);
        swa_kahan_add(
            &transported_charge_turns,
            charge_C * displacement_m / input->domain_length_m
        );
        out->maximum_particle_displacement_cells = max2(
            out->maximum_particle_displacement_cells,
            fabs(displacement_m) / cell_width_m
        );
    }

    for (i = 0U; i < input->cells; ++i) {
        swa_kahan_add(
            &deposited_initial,
            buffers->rho_initial_C_m3[i] * cell_volume_m3
        );
        swa_kahan_add(
            &deposited_final,
            buffers->rho_final_C_m3[i] * cell_volume_m3
        );
    }

    buffers->current_face_A_m2[0] = 0.0;
    for (i = 0U; i + 1U < input->cells; ++i) {
        const double density_change =
            buffers->rho_final_C_m3[i] -
            buffers->rho_initial_C_m3[i];
        buffers->current_face_A_m2[i + 1U] =
            buffers->current_face_A_m2[i] -
            cell_width_m * density_change / input->dt_s;
    }
    for (i = 0U; i < input->cells; ++i) {
        swa_kahan_add(&raw_current_mean,
                      buffers->current_face_A_m2[i]);
    }

    out->particle_charge_sum_C = particle_charge.sum;
    out->deposited_initial_charge_C = deposited_initial.sum;
    out->deposited_final_charge_C = deposited_final.sum;
    out->initial_charge_error_C =
        out->deposited_initial_charge_C -
        out->particle_charge_sum_C;
    out->final_charge_error_C =
        out->deposited_final_charge_C -
        out->particle_charge_sum_C;
    out->transported_charge_turns_C =
        transported_charge_turns.sum;
    out->target_mean_current_A_m2 =
        out->transported_charge_turns_C /
        (input->cross_section_area_m2 * input->dt_s);

    current_shift =
        out->target_mean_current_A_m2 -
        raw_current_mean.sum / (double)input->cells;
    for (i = 0U; i < input->cells; ++i) {
        buffers->current_face_A_m2[i] += current_shift;
        swa_kahan_add(&deposited_current_mean,
                      buffers->current_face_A_m2[i]);
    }
    buffers->current_face_A_m2[input->cells] =
        buffers->current_face_A_m2[0];
    out->deposited_mean_current_A_m2 =
        deposited_current_mean.sum / (double)input->cells;
    out->mean_current_error_A_m2 =
        out->deposited_mean_current_A_m2 -
        out->target_mean_current_A_m2;

    continuity_input.nx = input->cells;
    continuity_input.ny = 1U;
    continuity_input.dx_m = cell_width_m;
    continuity_input.dy_m = 1.0;
    continuity_input.dt_s = input->dt_s;
    continuity_input.rho_initial_C_m3 =
        buffers->rho_initial_C_m3;
    continuity_input.rho_final_C_m3 =
        buffers->rho_final_C_m3;
    continuity_input.jx_face_A_m2 =
        buffers->current_face_A_m2;
    continuity_input.jy_face_A_m2 = jy_zero;

    if (!swa_audit_charge_continuity(
            &continuity_input,
            input->continuity_relative_tolerance,
            input->continuity_absolute_tolerance_A_m3,
            input->continuity_global_absolute_tolerance_C_m,
            &out->continuity)) {
        free(initial_compensation);
        free(final_compensation);
        free(jy_zero);
        return 0;
    }

    charge_scale = max2(
        max2(fabs(out->particle_charge_sum_C),
             fabs(out->deposited_initial_charge_C)),
        max2(fabs(out->deposited_final_charge_C), DBL_MIN)
    );
    out->charge_closure_passes =
        fabs(out->initial_charge_error_C) <=
            input->charge_absolute_tolerance_C ||
        fabs(out->initial_charge_error_C) /
                charge_scale <=
            input->continuity_relative_tolerance;
    out->charge_closure_passes =
        out->charge_closure_passes &&
        (fabs(out->final_charge_error_C) <=
             input->charge_absolute_tolerance_C ||
         fabs(out->final_charge_error_C) /
                 charge_scale <=
             input->continuity_relative_tolerance);

    current_scale = max2(
        max2(fabs(out->target_mean_current_A_m2),
             fabs(out->deposited_mean_current_A_m2)),
        DBL_MIN
    );
    out->mean_current_passes =
        fabs(out->mean_current_error_A_m2) <=
            input->mean_current_absolute_tolerance_A_m2 ||
        fabs(out->mean_current_error_A_m2) /
                current_scale <=
            input->mean_current_relative_tolerance;
    out->finite =
        isfinite(out->particle_charge_sum_C) &&
        isfinite(out->deposited_initial_charge_C) &&
        isfinite(out->deposited_final_charge_C) &&
        isfinite(out->initial_charge_error_C) &&
        isfinite(out->final_charge_error_C) &&
        isfinite(out->transported_charge_turns_C) &&
        isfinite(out->target_mean_current_A_m2) &&
        isfinite(out->deposited_mean_current_A_m2) &&
        isfinite(out->mean_current_error_A_m2) &&
        isfinite(out->maximum_particle_displacement_cells) &&
        out->continuity.finite;
    out->passes = out->finite &&
                  out->charge_closure_passes &&
                  out->mean_current_passes &&
                  out->continuity.passes;

    free(initial_compensation);
    free(final_compensation);
    free(jy_zero);
    return 1;
}

int swa_write_periodic_cic_transport_receipt(
    FILE *fp,
    const swa_periodic_cic_transport_input *input,
    const swa_periodic_cic_transport_result *result
) {
    if (fp == NULL || input == NULL || result == NULL) {
        return 0;
    }
    return fprintf(
        fp,
        "{\n"
        "  \"schema\": \"spacewind.periodic-cic-transport/v1\",\n"
        "  \"cells\": %zu,\n"
        "  \"particles\": %zu,\n"
        "  \"domain_length_m\": %.17g,\n"
        "  \"cross_section_area_m2\": %.17g,\n"
        "  \"dt_s\": %.17g,\n"
        "  \"particle_charge_sum_C\": %.17g,\n"
        "  \"deposited_initial_charge_C\": %.17g,\n"
        "  \"deposited_final_charge_C\": %.17g,\n"
        "  \"transported_charge_turns_C\": %.17g,\n"
        "  \"target_mean_current_A_m2\": %.17g,\n"
        "  \"deposited_mean_current_A_m2\": %.17g,\n"
        "  \"mean_current_error_A_m2\": %.17g,\n"
        "  \"maximum_particle_displacement_cells\": %.17g,\n"
        "  \"continuity_rms_A_m3\": %.17g,\n"
        "  \"continuity_max_A_m3\": %.17g,\n"
        "  \"continuity_global_residual_C_m\": %.17g,\n"
        "  \"charge_closure_passes\": %s,\n"
        "  \"mean_current_passes\": %s,\n"
        "  \"continuity_passes\": %s,\n"
        "  \"passes\": %s,\n"
        "  \"nonclaim\": \"the reference reconstruction enforces a declared discrete continuity contract; it is not yet a local production electromagnetic PIC current deposition scheme\"\n"
        "}\n",
        input->cells,
        input->particles,
        input->domain_length_m,
        input->cross_section_area_m2,
        input->dt_s,
        result->particle_charge_sum_C,
        result->deposited_initial_charge_C,
        result->deposited_final_charge_C,
        result->transported_charge_turns_C,
        result->target_mean_current_A_m2,
        result->deposited_mean_current_A_m2,
        result->mean_current_error_A_m2,
        result->maximum_particle_displacement_cells,
        result->continuity.residual_rms_A_m3,
        result->continuity.residual_maximum_A_m3,
        result->continuity.global_charge_residual_C_m,
        result->charge_closure_passes ? "true" : "false",
        result->mean_current_passes ? "true" : "false",
        result->continuity.passes ? "true" : "false",
        result->passes ? "true" : "false"
    ) > 0;
}
