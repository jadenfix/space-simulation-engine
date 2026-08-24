#include "spacewind/field_ledger.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static double max2(double a, double b) { return a > b ? a : b; }
static double max5(double a, double b, double c, double d, double e) {
    return max2(max2(max2(a, b), max2(c, d)), e);
}
static int finite_v3(swa_vec3 v) {
    return isfinite(v.x) && isfinite(v.y) && isfinite(v.z);
}
static size_t cell_index(size_t i, size_t j, size_t nx) {
    return j * nx + i;
}
static size_t xface_index(size_t i, size_t j, size_t nx) {
    return j * (nx + 1U) + i;
}
static size_t yface_index(size_t i, size_t j, size_t nx) {
    return j * nx + i;
}
static int valid_grid(size_t nx, size_t ny) {
    if (nx == 0U || ny == 0U || nx == SIZE_MAX || ny == SIZE_MAX) {
        return 0;
    }
    if (ny > SIZE_MAX / nx) {
        return 0;
    }
    if (ny > SIZE_MAX / (nx + 1U)) {
        return 0;
    }
    if ((ny + 1U) > SIZE_MAX / nx) {
        return 0;
    }
    return 1;
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
static int valid_tolerances(double relative, double absolute) {
    return isfinite(relative) && isfinite(absolute) &&
           relative >= 0.0 && absolute >= 0.0;
}
static int residual_within(double absolute_value, double relative_value,
                           double absolute_tolerance,
                           double relative_tolerance) {
    return absolute_value <= absolute_tolerance ||
           relative_value <= relative_tolerance;
}

int swa_audit_charge_continuity(
    const swa_charge_continuity_input *input,
    double relative_tolerance,
    double local_absolute_tolerance_A_m3,
    double global_absolute_tolerance_C_m,
    swa_charge_continuity_result *out
) {
    size_t cells;
    size_t xfaces;
    size_t yfaces;
    size_t i;
    size_t j;
    double cell_area;
    double source_rms;
    double source_maximum = 0.0;
    double residual_maximum = 0.0;
    double boundary_current;
    double global_scale;
    swa_kahan initial_charge = {0.0, 0.0};
    swa_kahan final_charge = {0.0, 0.0};
    swa_kahan residual_abs = {0.0, 0.0};
    swa_kahan residual_square = {0.0, 0.0};
    swa_kahan source_square = {0.0, 0.0};
    swa_kahan residual_integral = {0.0, 0.0};
    swa_kahan boundary = {0.0, 0.0};

    if (input == NULL || out == NULL ||
        !valid_grid(input->nx, input->ny) ||
        !(input->dx_m > 0.0) || !(input->dy_m > 0.0) ||
        !(input->dt_s > 0.0) ||
        !isfinite(input->dx_m) || !isfinite(input->dy_m) ||
        !isfinite(input->dt_s) ||
        !valid_tolerances(relative_tolerance,
                          local_absolute_tolerance_A_m3) ||
        !isfinite(global_absolute_tolerance_C_m) ||
        global_absolute_tolerance_C_m < 0.0) {
        return 0;
    }
    cells = input->nx * input->ny;
    xfaces = (input->nx + 1U) * input->ny;
    yfaces = input->nx * (input->ny + 1U);
    if (!finite_array(input->rho_initial_C_m3, cells) ||
        !finite_array(input->rho_final_C_m3, cells) ||
        !finite_array(input->jx_face_A_m2, xfaces) ||
        !finite_array(input->jy_face_A_m2, yfaces)) {
        return 0;
    }

    memset(out, 0, sizeof(*out));
    cell_area = input->dx_m * input->dy_m;
    if (!(cell_area > 0.0) || !isfinite(cell_area)) {
        return 0;
    }

    for (j = 0U; j < input->ny; ++j) {
        for (i = 0U; i < input->nx; ++i) {
            const size_t c = cell_index(i, j, input->nx);
            const double temporal =
                (input->rho_final_C_m3[c] -
                 input->rho_initial_C_m3[c]) / input->dt_s;
            const double divergence =
                (input->jx_face_A_m2[xface_index(i + 1U, j, input->nx)] -
                 input->jx_face_A_m2[xface_index(i, j, input->nx)]) /
                    input->dx_m +
                (input->jy_face_A_m2[yface_index(i, j + 1U, input->nx)] -
                 input->jy_face_A_m2[yface_index(i, j, input->nx)]) /
                    input->dy_m;
            const double residual = temporal + divergence;
            const double source = max2(fabs(temporal), fabs(divergence));
            swa_kahan_add(&initial_charge,
                          input->rho_initial_C_m3[c] * cell_area);
            swa_kahan_add(&final_charge,
                          input->rho_final_C_m3[c] * cell_area);
            swa_kahan_add(&residual_abs, fabs(residual));
            swa_kahan_add(&residual_square, residual * residual);
            swa_kahan_add(&source_square, source * source);
            swa_kahan_add(&residual_integral, residual * cell_area);
            residual_maximum = max2(residual_maximum, fabs(residual));
            source_maximum = max2(source_maximum, source);
        }
    }

    for (j = 0U; j < input->ny; ++j) {
        swa_kahan_add(
            &boundary,
            (input->jx_face_A_m2[xface_index(input->nx, j, input->nx)] -
             input->jx_face_A_m2[xface_index(0U, j, input->nx)]) *
                input->dy_m
        );
    }
    for (i = 0U; i < input->nx; ++i) {
        swa_kahan_add(
            &boundary,
            (input->jy_face_A_m2[yface_index(i, input->ny, input->nx)] -
             input->jy_face_A_m2[yface_index(i, 0U, input->nx)]) *
                input->dx_m
        );
    }

    out->initial_charge_C_m = initial_charge.sum;
    out->final_charge_C_m = final_charge.sum;
    out->outward_boundary_current_A_m = boundary.sum;
    out->residual_mean_absolute_A_m3 =
        residual_abs.sum / (double)cells;
    out->residual_rms_A_m3 =
        sqrt(residual_square.sum / (double)cells);
    out->residual_maximum_A_m3 = residual_maximum;
    source_rms = sqrt(source_square.sum / (double)cells);
    out->residual_relative_rms =
        out->residual_rms_A_m3 / max2(source_rms, DBL_MIN);
    out->residual_relative_maximum =
        residual_maximum / max2(source_maximum, DBL_MIN);
    out->residual_integral_A_m = residual_integral.sum;
    boundary_current = out->outward_boundary_current_A_m;
    out->global_charge_residual_C_m =
        out->final_charge_C_m - out->initial_charge_C_m +
        input->dt_s * boundary_current;
    global_scale = max5(
        fabs(out->initial_charge_C_m),
        fabs(out->final_charge_C_m),
        fabs(input->dt_s * boundary_current),
        fabs(input->dt_s * out->residual_integral_A_m),
        DBL_MIN
    );
    out->global_charge_relative =
        fabs(out->global_charge_residual_C_m) / global_scale;
    out->divergence_theorem_error_C_m =
        out->global_charge_residual_C_m -
        input->dt_s * out->residual_integral_A_m;
    out->finite =
        isfinite(out->initial_charge_C_m) &&
        isfinite(out->final_charge_C_m) &&
        isfinite(out->outward_boundary_current_A_m) &&
        isfinite(out->residual_mean_absolute_A_m3) &&
        isfinite(out->residual_rms_A_m3) &&
        isfinite(out->residual_maximum_A_m3) &&
        isfinite(out->residual_relative_rms) &&
        isfinite(out->residual_relative_maximum) &&
        isfinite(out->residual_integral_A_m) &&
        isfinite(out->global_charge_residual_C_m) &&
        isfinite(out->global_charge_relative) &&
        isfinite(out->divergence_theorem_error_C_m);
    out->local_passes =
        residual_within(out->residual_rms_A_m3,
                        out->residual_relative_rms,
                        local_absolute_tolerance_A_m3,
                        relative_tolerance) &&
        residual_within(out->residual_maximum_A_m3,
                        out->residual_relative_maximum,
                        local_absolute_tolerance_A_m3,
                        relative_tolerance);
    out->global_passes =
        residual_within(fabs(out->global_charge_residual_C_m),
                        out->global_charge_relative,
                        global_absolute_tolerance_C_m,
                        relative_tolerance);
    out->passes = out->finite && out->local_passes &&
                  out->global_passes;
    return 1;
}

static int audit_staggered_constraint(
    const swa_staggered_scalar_constraint_input *input,
    double multiplier,
    int subtract_source,
    double relative_tolerance,
    double absolute_tolerance,
    swa_field_constraint_result *out
) {
    size_t cells;
    size_t xfaces;
    size_t yfaces;
    size_t i;
    size_t j;
    double residual_maximum = 0.0;
    double source_maximum = 0.0;
    double source_rms;
    double area;
    swa_kahan residual_abs = {0.0, 0.0};
    swa_kahan residual_square = {0.0, 0.0};
    swa_kahan source_square = {0.0, 0.0};
    swa_kahan signed_integral = {0.0, 0.0};

    if (input == NULL || out == NULL ||
        !valid_grid(input->nx, input->ny) ||
        !(input->dx_m > 0.0) || !(input->dy_m > 0.0) ||
        !isfinite(input->dx_m) || !isfinite(input->dy_m) ||
        !isfinite(multiplier) ||
        !valid_tolerances(relative_tolerance, absolute_tolerance)) {
        return 0;
    }
    cells = input->nx * input->ny;
    xfaces = (input->nx + 1U) * input->ny;
    yfaces = input->nx * (input->ny + 1U);
    if (!finite_array(input->x_face_field, xfaces) ||
        !finite_array(input->y_face_field, yfaces) ||
        (subtract_source && !finite_array(input->cell_source, cells))) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    area = input->dx_m * input->dy_m;
    if (!(area > 0.0) || !isfinite(area)) {
        return 0;
    }

    for (j = 0U; j < input->ny; ++j) {
        for (i = 0U; i < input->nx; ++i) {
            const size_t c = cell_index(i, j, input->nx);
            const double divergence =
                (input->x_face_field[xface_index(i + 1U, j, input->nx)] -
                 input->x_face_field[xface_index(i, j, input->nx)]) /
                    input->dx_m +
                (input->y_face_field[yface_index(i, j + 1U, input->nx)] -
                 input->y_face_field[yface_index(i, j, input->nx)]) /
                    input->dy_m;
            const double lhs = multiplier * divergence;
            const double rhs = subtract_source ? input->cell_source[c] : 0.0;
            const double residual = lhs - rhs;
            const double scale = max2(fabs(lhs), fabs(rhs));
            swa_kahan_add(&residual_abs, fabs(residual));
            swa_kahan_add(&residual_square, residual * residual);
            swa_kahan_add(&source_square, scale * scale);
            swa_kahan_add(&signed_integral, residual * area);
            residual_maximum = max2(residual_maximum, fabs(residual));
            source_maximum = max2(source_maximum, scale);
        }
    }

    out->residual_mean_absolute =
        residual_abs.sum / (double)cells;
    out->residual_rms = sqrt(residual_square.sum / (double)cells);
    out->residual_maximum = residual_maximum;
    source_rms = sqrt(source_square.sum / (double)cells);
    out->residual_relative_rms =
        out->residual_rms / max2(source_rms, DBL_MIN);
    out->residual_relative_maximum =
        residual_maximum / max2(source_maximum, DBL_MIN);
    out->residual_signed_integral = signed_integral.sum;
    out->finite =
        isfinite(out->residual_mean_absolute) &&
        isfinite(out->residual_rms) &&
        isfinite(out->residual_maximum) &&
        isfinite(out->residual_relative_rms) &&
        isfinite(out->residual_relative_maximum) &&
        isfinite(out->residual_signed_integral);
    out->passes = out->finite &&
        residual_within(out->residual_rms,
                        out->residual_relative_rms,
                        absolute_tolerance,
                        relative_tolerance) &&
        residual_within(out->residual_maximum,
                        out->residual_relative_maximum,
                        absolute_tolerance,
                        relative_tolerance);
    return 1;
}

int swa_audit_gauss_law(
    const swa_staggered_scalar_constraint_input *input,
    double relative_tolerance,
    double absolute_tolerance_C_m3,
    swa_field_constraint_result *out
) {
    return audit_staggered_constraint(
        input,
        SWA_EPS0,
        1,
        relative_tolerance,
        absolute_tolerance_C_m3,
        out
    );
}

int swa_audit_magnetic_divergence(
    const swa_staggered_scalar_constraint_input *input,
    double relative_tolerance,
    double absolute_tolerance_T_m,
    swa_field_constraint_result *out
) {
    return audit_staggered_constraint(
        input,
        1.0,
        0,
        relative_tolerance,
        absolute_tolerance_T_m,
        out
    );
}

int swa_audit_poynting_ledger(
    const swa_poynting_ledger *ledger,
    double relative_tolerance,
    double absolute_tolerance_J,
    swa_poynting_ledger_result *out
) {
    double scale;
    if (ledger == NULL || out == NULL ||
        !valid_tolerances(relative_tolerance, absolute_tolerance_J) ||
        !isfinite(ledger->field_energy_initial_J) ||
        !isfinite(ledger->field_energy_final_J) ||
        !isfinite(ledger->particle_work_J) ||
        !isfinite(ledger->outward_poynting_energy_J) ||
        !isfinite(ledger->impressed_source_energy_J)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    out->residual_J =
        ledger->field_energy_final_J -
        ledger->field_energy_initial_J +
        ledger->particle_work_J +
        ledger->outward_poynting_energy_J -
        ledger->impressed_source_energy_J;
    scale = max5(
        fabs(ledger->field_energy_initial_J),
        fabs(ledger->field_energy_final_J),
        fabs(ledger->particle_work_J),
        fabs(ledger->outward_poynting_energy_J),
        max2(fabs(ledger->impressed_source_energy_J), DBL_MIN)
    );
    out->relative_residual = fabs(out->residual_J) / scale;
    out->finite = isfinite(out->residual_J) &&
                  isfinite(out->relative_residual);
    out->passes = out->finite &&
        residual_within(fabs(out->residual_J),
                        out->relative_residual,
                        absolute_tolerance_J,
                        relative_tolerance);
    return 1;
}

int swa_audit_field_momentum(
    const swa_field_momentum_ledger *ledger,
    double relative_tolerance,
    double absolute_tolerance_Ns,
    swa_field_momentum_result *out
) {
    double scale;
    swa_vec3 field_change;
    if (ledger == NULL || out == NULL ||
        !valid_tolerances(relative_tolerance, absolute_tolerance_Ns) ||
        !finite_v3(ledger->field_momentum_initial_Ns) ||
        !finite_v3(ledger->field_momentum_final_Ns) ||
        !finite_v3(ledger->mechanical_impulse_Ns) ||
        !finite_v3(ledger->outward_maxwell_impulse_Ns) ||
        !finite_v3(ledger->impressed_external_impulse_Ns)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    field_change = swa_vsub(
        ledger->field_momentum_final_Ns,
        ledger->field_momentum_initial_Ns
    );
    out->residual_Ns = swa_vsub(
        swa_vadd(
            swa_vadd(field_change, ledger->mechanical_impulse_Ns),
            ledger->outward_maxwell_impulse_Ns
        ),
        ledger->impressed_external_impulse_Ns
    );
    out->residual_norm_Ns = swa_vnorm(out->residual_Ns);
    scale = max5(
        swa_vnorm(ledger->field_momentum_initial_Ns),
        swa_vnorm(ledger->field_momentum_final_Ns),
        swa_vnorm(ledger->mechanical_impulse_Ns),
        swa_vnorm(ledger->outward_maxwell_impulse_Ns),
        max2(swa_vnorm(ledger->impressed_external_impulse_Ns), DBL_MIN)
    );
    out->relative_residual = out->residual_norm_Ns / scale;
    out->finite = finite_v3(out->residual_Ns) &&
                  isfinite(out->residual_norm_Ns) &&
                  isfinite(out->relative_residual);
    out->passes = out->finite &&
        residual_within(out->residual_norm_Ns,
                        out->relative_residual,
                        absolute_tolerance_Ns,
                        relative_tolerance);
    return 1;
}

int swa_local_electromagnetic_state(
    swa_vec3 electric_field_V_m,
    swa_vec3 magnetic_field_T,
    swa_local_em_result *out
) {
    double electric_square;
    double magnetic_square;
    double energy_scale;
    if (out == NULL || !finite_v3(electric_field_V_m) ||
        !finite_v3(magnetic_field_T)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    electric_square = swa_vdot(electric_field_V_m, electric_field_V_m);
    magnetic_square = swa_vdot(magnetic_field_T, magnetic_field_T);
    out->energy_density_J_m3 =
        0.5 * (SWA_EPS0 * electric_square +
               magnetic_square / SWA_MU0);
    out->poynting_flux_W_m2 =
        swa_vscale(
            swa_vcross(electric_field_V_m, magnetic_field_T),
            1.0 / SWA_MU0
        );
    out->momentum_density_Ns_m3 =
        swa_vscale(
            out->poynting_flux_W_m2,
            1.0 / (SWA_C * SWA_C)
        );
    out->invariant_B2_minus_E2_over_c2_T2 =
        magnetic_square -
        electric_square / (SWA_C * SWA_C);
    out->invariant_E_dot_B_over_c =
        swa_vdot(electric_field_V_m, magnetic_field_T) / SWA_C;
    energy_scale = SWA_C * out->energy_density_J_m3;
    out->dominant_energy_ratio =
        swa_vnorm(out->poynting_flux_W_m2) /
        max2(energy_scale, DBL_MIN);
    out->finite =
        isfinite(out->energy_density_J_m3) &&
        finite_v3(out->poynting_flux_W_m2) &&
        finite_v3(out->momentum_density_Ns_m3) &&
        isfinite(out->invariant_B2_minus_E2_over_c2_T2) &&
        isfinite(out->invariant_E_dot_B_over_c) &&
        isfinite(out->dominant_energy_ratio);
    return 1;
}

int swa_maxwell_traction(
    swa_vec3 electric_field_V_m,
    swa_vec3 magnetic_field_T,
    swa_vec3 surface_normal,
    swa_maxwell_traction_result *out
) {
    double normal_norm;
    double electric_square;
    double magnetic_square;
    swa_vec3 unit;
    swa_vec3 electric_traction;
    swa_vec3 magnetic_traction;
    if (out == NULL || !finite_v3(electric_field_V_m) ||
        !finite_v3(magnetic_field_T) || !finite_v3(surface_normal)) {
        return 0;
    }
    normal_norm = swa_vnorm(surface_normal);
    if (!(normal_norm > 0.0) || !isfinite(normal_norm)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    unit = swa_vscale(surface_normal, 1.0 / normal_norm);
    electric_square = swa_vdot(electric_field_V_m,
                               electric_field_V_m);
    magnetic_square = swa_vdot(magnetic_field_T,
                               magnetic_field_T);
    electric_traction = swa_vscale(
        swa_vsub(
            swa_vscale(electric_field_V_m,
                       swa_vdot(electric_field_V_m, unit)),
            swa_vscale(unit, 0.5 * electric_square)
        ),
        SWA_EPS0
    );
    magnetic_traction = swa_vscale(
        swa_vsub(
            swa_vscale(magnetic_field_T,
                       swa_vdot(magnetic_field_T, unit)),
            swa_vscale(unit, 0.5 * magnetic_square)
        ),
        1.0 / SWA_MU0
    );
    out->unit_normal = unit;
    out->traction_N_m2 =
        swa_vadd(electric_traction, magnetic_traction);
    out->normal_pressure_Pa =
        swa_vdot(out->traction_N_m2, unit);
    out->finite = finite_v3(out->unit_normal) &&
                  finite_v3(out->traction_N_m2) &&
                  isfinite(out->normal_pressure_Pa);
    return 1;
}

int swa_write_field_ledger_receipt(
    FILE *fp,
    const swa_charge_continuity_result *continuity,
    const swa_field_constraint_result *gauss,
    const swa_field_constraint_result *magnetic_divergence,
    const swa_poynting_ledger_result *energy,
    const swa_field_momentum_result *momentum
) {
    if (fp == NULL || continuity == NULL || gauss == NULL ||
        magnetic_divergence == NULL || energy == NULL ||
        momentum == NULL) {
        return 0;
    }
    return fprintf(
        fp,
        "{\n"
        "  \"schema\": \"spacewind.field-particle-ledger/v1\",\n"
        "  \"charge_continuity\": {\n"
        "    \"local_rms_A_m3\": %.17g,\n"
        "    \"local_max_A_m3\": %.17g,\n"
        "    \"global_residual_C_m\": %.17g,\n"
        "    \"divergence_theorem_error_C_m\": %.17g,\n"
        "    \"passes\": %s\n"
        "  },\n"
        "  \"gauss_law\": {\n"
        "    \"rms_C_m3\": %.17g,\n"
        "    \"max_C_m3\": %.17g,\n"
        "    \"passes\": %s\n"
        "  },\n"
        "  \"magnetic_divergence\": {\n"
        "    \"rms_T_m\": %.17g,\n"
        "    \"max_T_m\": %.17g,\n"
        "    \"passes\": %s\n"
        "  },\n"
        "  \"poynting_energy\": {\n"
        "    \"residual_J\": %.17g,\n"
        "    \"relative_residual\": %.17g,\n"
        "    \"passes\": %s\n"
        "  },\n"
        "  \"field_momentum\": {\n"
        "    \"residual_Ns\": [%.17g, %.17g, %.17g],\n"
        "    \"relative_residual\": %.17g,\n"
        "    \"passes\": %s\n"
        "  },\n"
        "  \"nonclaim\": \"ledger closure verifies declared discrete bookkeeping only; it does not establish a realizable plasma-wing force\"\n"
        "}\n",
        continuity->residual_rms_A_m3,
        continuity->residual_maximum_A_m3,
        continuity->global_charge_residual_C_m,
        continuity->divergence_theorem_error_C_m,
        continuity->passes ? "true" : "false",
        gauss->residual_rms,
        gauss->residual_maximum,
        gauss->passes ? "true" : "false",
        magnetic_divergence->residual_rms,
        magnetic_divergence->residual_maximum,
        magnetic_divergence->passes ? "true" : "false",
        energy->residual_J,
        energy->relative_residual,
        energy->passes ? "true" : "false",
        momentum->residual_Ns.x,
        momentum->residual_Ns.y,
        momentum->residual_Ns.z,
        momentum->relative_residual,
        momentum->passes ? "true" : "false"
    ) > 0;
}
