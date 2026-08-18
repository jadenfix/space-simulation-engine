#include "spacewind/covariant_em.h"

#include <float.h>
#include <math.h>
#include <string.h>

static double max2(double a, double b) { return a > b ? a : b; }
static double max3(double a, double b, double c) {
    return max2(max2(a, b), c);
}

static int finite_vec3(swa_vec3 value) {
    return isfinite(value.x) && isfinite(value.y) && isfinite(value.z);
}

static int finite_fields(const swa_em_fields *fields) {
    return fields != NULL &&
           finite_vec3(fields->electric_field_V_m) &&
           finite_vec3(fields->magnetic_field_T);
}

static int finite_tensor(const swa_tensor4 *tensor) {
    size_t mu;
    size_t nu;
    if (tensor == NULL) {
        return 0;
    }
    for (mu = 0U; mu < 4U; ++mu) {
        for (nu = 0U; nu < 4U; ++nu) {
            if (!isfinite(tensor->component[mu][nu])) {
                return 0;
            }
        }
    }
    return 1;
}

swa_four swa_four_make(double time_component, double x, double y, double z) {
    swa_four value = {{time_component, x, y, z}};
    return value;
}

swa_four swa_four_add(swa_four a, swa_four b) {
    size_t i;
    swa_four out;
    for (i = 0U; i < 4U; ++i) {
        out.component[i] = a.component[i] + b.component[i];
    }
    return out;
}

swa_four swa_four_sub(swa_four a, swa_four b) {
    size_t i;
    swa_four out;
    for (i = 0U; i < 4U; ++i) {
        out.component[i] = a.component[i] - b.component[i];
    }
    return out;
}

swa_four swa_four_scale(swa_four a, double scale) {
    size_t i;
    swa_four out;
    for (i = 0U; i < 4U; ++i) {
        out.component[i] = a.component[i] * scale;
    }
    return out;
}

double swa_four_minkowski_dot(swa_four a, swa_four b) {
    return -a.component[0] * b.component[0] +
            a.component[1] * b.component[1] +
            a.component[2] * b.component[2] +
            a.component[3] * b.component[3];
}

double swa_four_euclidean_norm(swa_four a) {
    return hypot(hypot(a.component[0], a.component[1]),
                 hypot(a.component[2], a.component[3]));
}

int swa_four_is_finite(swa_four value) {
    size_t i;
    for (i = 0U; i < 4U; ++i) {
        if (!isfinite(value.component[i])) {
            return 0;
        }
    }
    return 1;
}

int swa_lorentz_boost_matrix(
    swa_vec3 frame_velocity_mps,
    swa_tensor4 *out
) {
    double beta[3];
    double beta2;
    double gamma;
    double spatial_factor;
    size_t i;
    size_t j;
    if (out == NULL || !finite_vec3(frame_velocity_mps)) {
        return 0;
    }
    beta[0] = frame_velocity_mps.x / SWA_C;
    beta[1] = frame_velocity_mps.y / SWA_C;
    beta[2] = frame_velocity_mps.z / SWA_C;
    beta2 = beta[0] * beta[0] + beta[1] * beta[1] +
            beta[2] * beta[2];
    if (!isfinite(beta2) || beta2 < 0.0 || beta2 >= 1.0) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    gamma = 1.0 / sqrt(1.0 - beta2);
    out->component[0][0] = gamma;
    for (i = 0U; i < 3U; ++i) {
        out->component[0][i + 1U] = -gamma * beta[i];
        out->component[i + 1U][0] = -gamma * beta[i];
    }
    spatial_factor = beta2 > 0.0 ? (gamma - 1.0) / beta2 : 0.0;
    for (i = 0U; i < 3U; ++i) {
        for (j = 0U; j < 3U; ++j) {
            out->component[i + 1U][j + 1U] =
                (i == j ? 1.0 : 0.0) +
                spatial_factor * beta[i] * beta[j];
        }
    }
    return finite_tensor(out);
}

int swa_lorentz_transform_four(
    swa_vec3 frame_velocity_mps,
    swa_four input,
    swa_four *out
) {
    swa_tensor4 matrix;
    swa_four transformed = {{0.0, 0.0, 0.0, 0.0}};
    size_t mu;
    size_t nu;
    if (out == NULL || !swa_four_is_finite(input) ||
        !swa_lorentz_boost_matrix(frame_velocity_mps, &matrix)) {
        return 0;
    }
    for (mu = 0U; mu < 4U; ++mu) {
        for (nu = 0U; nu < 4U; ++nu) {
            transformed.component[mu] +=
                matrix.component[mu][nu] * input.component[nu];
        }
    }
    *out = transformed;
    return swa_four_is_finite(*out);
}

int swa_lorentz_transform_tensor(
    swa_vec3 frame_velocity_mps,
    const swa_tensor4 *input,
    swa_tensor4 *out
) {
    swa_tensor4 matrix;
    swa_tensor4 transformed;
    size_t mu;
    size_t nu;
    size_t alpha;
    size_t beta;
    if (input == NULL || out == NULL || !finite_tensor(input) ||
        !swa_lorentz_boost_matrix(frame_velocity_mps, &matrix)) {
        return 0;
    }
    memset(&transformed, 0, sizeof(transformed));
    for (mu = 0U; mu < 4U; ++mu) {
        for (nu = 0U; nu < 4U; ++nu) {
            for (alpha = 0U; alpha < 4U; ++alpha) {
                for (beta = 0U; beta < 4U; ++beta) {
                    transformed.component[mu][nu] +=
                        matrix.component[mu][alpha] *
                        matrix.component[nu][beta] *
                        input->component[alpha][beta];
                }
            }
        }
    }
    *out = transformed;
    return finite_tensor(out);
}

int swa_em_field_tensor(
    const swa_em_fields *fields,
    swa_tensor4 *out
) {
    const double ex = fields != NULL ? fields->electric_field_V_m.x : 0.0;
    const double ey = fields != NULL ? fields->electric_field_V_m.y : 0.0;
    const double ez = fields != NULL ? fields->electric_field_V_m.z : 0.0;
    const double bx = fields != NULL ? fields->magnetic_field_T.x : 0.0;
    const double by = fields != NULL ? fields->magnetic_field_T.y : 0.0;
    const double bz = fields != NULL ? fields->magnetic_field_T.z : 0.0;
    if (!finite_fields(fields) || out == NULL) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    out->component[0][1] = ex / SWA_C;
    out->component[0][2] = ey / SWA_C;
    out->component[0][3] = ez / SWA_C;
    out->component[1][0] = -ex / SWA_C;
    out->component[2][0] = -ey / SWA_C;
    out->component[3][0] = -ez / SWA_C;
    out->component[1][2] = bz;
    out->component[2][1] = -bz;
    out->component[1][3] = -by;
    out->component[3][1] = by;
    out->component[2][3] = bx;
    out->component[3][2] = -bx;
    return finite_tensor(out);
}

int swa_em_stress_energy_tensor(
    const swa_em_fields *fields,
    swa_tensor4 *out,
    swa_covariant_em_diagnostics *diagnostics
) {
    static const double metric[4] = {-1.0, 1.0, 1.0, 1.0};
    swa_tensor4 field_tensor;
    double invariant = 0.0;
    double trace = 0.0;
    size_t mu;
    size_t nu;
    size_t alpha;
    size_t beta;
    if (!finite_fields(fields) || out == NULL || diagnostics == NULL ||
        !swa_em_field_tensor(fields, &field_tensor)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    memset(diagnostics, 0, sizeof(*diagnostics));
    for (alpha = 0U; alpha < 4U; ++alpha) {
        for (beta = 0U; beta < 4U; ++beta) {
            invariant += field_tensor.component[alpha][beta] *
                         field_tensor.component[alpha][beta] *
                         metric[alpha] * metric[beta];
        }
    }
    for (mu = 0U; mu < 4U; ++mu) {
        for (nu = 0U; nu < 4U; ++nu) {
            double contraction = 0.0;
            for (alpha = 0U; alpha < 4U; ++alpha) {
                contraction += field_tensor.component[mu][alpha] *
                               field_tensor.component[nu][alpha] *
                               metric[alpha];
            }
            out->component[mu][nu] =
                (contraction -
                 (mu == nu ? 0.25 * metric[mu] * invariant : 0.0)) /
                SWA_MU0;
        }
    }
    diagnostics->invariant_B2_minus_E2_over_c2_T2 =
        swa_vdot(fields->magnetic_field_T, fields->magnetic_field_T) -
        swa_vdot(fields->electric_field_V_m,
                 fields->electric_field_V_m) /
            (SWA_C * SWA_C);
    diagnostics->invariant_E_dot_B_over_c =
        swa_vdot(fields->electric_field_V_m,
                 fields->magnetic_field_T) /
        SWA_C;
    diagnostics->energy_density_J_m3 = out->component[0][0];
    diagnostics->poynting_flux_W_m2 = swa_vscale(
        swa_vcross(fields->electric_field_V_m,
                   fields->magnetic_field_T),
        1.0 / SWA_MU0
    );
    diagnostics->dominant_energy_ratio =
        swa_vnorm(diagnostics->poynting_flux_W_m2) /
        max2(SWA_C * diagnostics->energy_density_J_m3, DBL_MIN);
    for (mu = 0U; mu < 4U; ++mu) {
        trace += metric[mu] * out->component[mu][mu];
    }
    diagnostics->tensor_trace_J_m3 = trace;
    diagnostics->finite = finite_tensor(out) &&
        isfinite(diagnostics->invariant_B2_minus_E2_over_c2_T2) &&
        isfinite(diagnostics->invariant_E_dot_B_over_c) &&
        isfinite(diagnostics->energy_density_J_m3) &&
        finite_vec3(diagnostics->poynting_flux_W_m2) &&
        isfinite(diagnostics->dominant_energy_ratio) &&
        isfinite(diagnostics->tensor_trace_J_m3);
    diagnostics->dominant_energy_condition =
        diagnostics->finite &&
        diagnostics->energy_density_J_m3 >= 0.0 &&
        diagnostics->dominant_energy_ratio <= 1.0 + 64.0 * DBL_EPSILON;
    return diagnostics->finite;
}

int swa_lorentz_transform_em_fields(
    swa_vec3 frame_velocity_mps,
    const swa_em_fields *input,
    swa_em_fields *out
) {
    const double speed2 = swa_vdot(frame_velocity_mps,
                                   frame_velocity_mps);
    double gamma;
    double coefficient;
    swa_vec3 transformed_electric;
    swa_vec3 transformed_magnetic;
    if (!finite_fields(input) || out == NULL ||
        !finite_vec3(frame_velocity_mps) ||
        !isfinite(speed2) || speed2 < 0.0 ||
        speed2 >= SWA_C * SWA_C) {
        return 0;
    }
    gamma = 1.0 / sqrt(1.0 - speed2 / (SWA_C * SWA_C));
    coefficient = gamma * gamma / (gamma + 1.0) /
                  (SWA_C * SWA_C);
    transformed_electric = swa_vsub(
        swa_vscale(
            swa_vadd(input->electric_field_V_m,
                     swa_vcross(frame_velocity_mps,
                                input->magnetic_field_T)),
            gamma
        ),
        swa_vscale(frame_velocity_mps,
                   coefficient *
                   swa_vdot(frame_velocity_mps,
                            input->electric_field_V_m))
    );
    transformed_magnetic = swa_vsub(
        swa_vscale(
            swa_vsub(
                input->magnetic_field_T,
                swa_vscale(
                    swa_vcross(frame_velocity_mps,
                               input->electric_field_V_m),
                    1.0 / (SWA_C * SWA_C)
                )
            ),
            gamma
        ),
        swa_vscale(frame_velocity_mps,
                   coefficient *
                   swa_vdot(frame_velocity_mps,
                            input->magnetic_field_T))
    );
    out->electric_field_V_m = transformed_electric;
    out->magnetic_field_T = transformed_magnetic;
    return finite_fields(out);
}

int swa_four_current(
    double charge_density_C_m3,
    swa_vec3 current_density_A_m2,
    swa_four *out
) {
    if (out == NULL || !isfinite(charge_density_C_m3) ||
        !finite_vec3(current_density_A_m2)) {
        return 0;
    }
    *out = swa_four_make(
        SWA_C * charge_density_C_m3,
        current_density_A_m2.x,
        current_density_A_m2.y,
        current_density_A_m2.z
    );
    return swa_four_is_finite(*out);
}

int swa_lorentz_force_density(
    const swa_em_fields *fields,
    double charge_density_C_m3,
    swa_vec3 current_density_A_m2,
    swa_four *out_force_density_N_m3
) {
    static const double metric[4] = {-1.0, 1.0, 1.0, 1.0};
    swa_tensor4 field_tensor;
    swa_four current;
    swa_four force = {{0.0, 0.0, 0.0, 0.0}};
    size_t mu;
    size_t nu;
    if (out_force_density_N_m3 == NULL ||
        !swa_em_field_tensor(fields, &field_tensor) ||
        !swa_four_current(charge_density_C_m3,
                          current_density_A_m2, &current)) {
        return 0;
    }
    for (mu = 0U; mu < 4U; ++mu) {
        for (nu = 0U; nu < 4U; ++nu) {
            force.component[mu] +=
                field_tensor.component[mu][nu] *
                metric[nu] * current.component[nu];
        }
    }
    *out_force_density_N_m3 = force;
    return swa_four_is_finite(*out_force_density_N_m3);
}

int swa_particle_four_momentum(
    double rest_mass_kg,
    swa_vec3 velocity_mps,
    swa_four *out_momentum_Ns
) {
    const double speed2 = swa_vdot(velocity_mps, velocity_mps);
    double gamma;
    if (out_momentum_Ns == NULL || !isfinite(rest_mass_kg) ||
        !(rest_mass_kg > 0.0) || !finite_vec3(velocity_mps) ||
        !isfinite(speed2) || speed2 < 0.0 ||
        speed2 >= SWA_C * SWA_C) {
        return 0;
    }
    gamma = 1.0 / sqrt(1.0 - speed2 / (SWA_C * SWA_C));
    *out_momentum_Ns = swa_four_make(
        gamma * rest_mass_kg * SWA_C,
        gamma * rest_mass_kg * velocity_mps.x,
        gamma * rest_mass_kg * velocity_mps.y,
        gamma * rest_mass_kg * velocity_mps.z
    );
    return swa_four_is_finite(*out_momentum_Ns);
}

int swa_audit_particle_mass_shell(
    double rest_mass_kg,
    swa_four momentum_Ns,
    double relative_tolerance,
    double absolute_tolerance_kg2_m2_s2,
    double *out_residual_kg2_m2_s2,
    double *out_relative_residual,
    int *out_passes
) {
    const double expected = -rest_mass_kg * rest_mass_kg *
                            SWA_C * SWA_C;
    double observed;
    double scale;
    if (!isfinite(rest_mass_kg) || !(rest_mass_kg > 0.0) ||
        !swa_four_is_finite(momentum_Ns) ||
        !isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
        !isfinite(absolute_tolerance_kg2_m2_s2) ||
        absolute_tolerance_kg2_m2_s2 < 0.0 ||
        out_residual_kg2_m2_s2 == NULL ||
        out_relative_residual == NULL || out_passes == NULL) {
        return 0;
    }
    observed = swa_four_minkowski_dot(momentum_Ns, momentum_Ns);
    *out_residual_kg2_m2_s2 = observed - expected;
    scale = max3(fabs(observed), fabs(expected), DBL_MIN);
    *out_relative_residual =
        fabs(*out_residual_kg2_m2_s2) / scale;
    *out_passes =
        fabs(*out_residual_kg2_m2_s2) <=
            absolute_tolerance_kg2_m2_s2 ||
        *out_relative_residual <= relative_tolerance;
    return isfinite(*out_residual_kg2_m2_s2) &&
           isfinite(*out_relative_residual);
}

int swa_audit_four_momentum_ledger(
    const swa_four_momentum_ledger *ledger,
    double relative_tolerance,
    double absolute_tolerance_Ns,
    swa_four_momentum_ledger_result *out
) {
    swa_four expected;
    swa_four initial_total;
    swa_four final_total;
    double scale;
    if (ledger == NULL || out == NULL ||
        !isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
        !isfinite(absolute_tolerance_Ns) ||
        absolute_tolerance_Ns < 0.0 ||
        !swa_four_is_finite(ledger->field_initial_Ns) ||
        !swa_four_is_finite(ledger->field_final_Ns) ||
        !swa_four_is_finite(ledger->matter_initial_Ns) ||
        !swa_four_is_finite(ledger->matter_final_Ns) ||
        !swa_four_is_finite(ledger->external_impulse_Ns) ||
        !swa_four_is_finite(ledger->momentum_in_Ns) ||
        !swa_four_is_finite(ledger->momentum_out_Ns)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    initial_total = swa_four_add(
        ledger->field_initial_Ns,
        ledger->matter_initial_Ns
    );
    final_total = swa_four_add(
        ledger->field_final_Ns,
        ledger->matter_final_Ns
    );
    expected = swa_four_add(
        initial_total,
        swa_four_add(
            ledger->external_impulse_Ns,
            swa_four_sub(ledger->momentum_in_Ns,
                         ledger->momentum_out_Ns)
        )
    );
    out->residual_Ns = swa_four_sub(final_total, expected);
    out->residual_norm_Ns = swa_four_euclidean_norm(out->residual_Ns);
    scale = max3(swa_four_euclidean_norm(initial_total),
                 swa_four_euclidean_norm(final_total),
                 max2(swa_four_euclidean_norm(expected), DBL_MIN));
    out->relative_residual = out->residual_norm_Ns / scale;
    out->finite = swa_four_is_finite(out->residual_Ns) &&
                  isfinite(out->residual_norm_Ns) &&
                  isfinite(out->relative_residual);
    out->passes = out->finite &&
        (out->residual_norm_Ns <= absolute_tolerance_Ns ||
         out->relative_residual <= relative_tolerance);
    return 1;
}

int swa_write_covariant_em_receipt(
    FILE *fp,
    size_t checks,
    size_t failures,
    uint64_t deterministic_hash,
    double worst_field_invariant_relative,
    double worst_tensor_covariance_relative,
    double worst_four_force_covariance_relative,
    double worst_mass_shell_relative
) {
    if (fp == NULL || !isfinite(worst_field_invariant_relative) ||
        !isfinite(worst_tensor_covariance_relative) ||
        !isfinite(worst_four_force_covariance_relative) ||
        !isfinite(worst_mass_shell_relative)) {
        return 0;
    }
    return fprintf(
        fp,
        "{\n"
        "  \"schema\": \"spacewind.covariant-em-assurance/v1\",\n"
        "  \"checks\": %zu,\n"
        "  \"failures\": %zu,\n"
        "  \"deterministic_hash\": \"%016llx\",\n"
        "  \"worst_field_invariant_relative\": %.17g,\n"
        "  \"worst_tensor_covariance_relative\": %.17g,\n"
        "  \"worst_four_force_covariance_relative\": %.17g,\n"
        "  \"worst_mass_shell_relative\": %.17g,\n"
        "  \"passes\": %s,\n"
        "  \"nonclaim\": \"finite Lorentz-covariance and stress-energy checks do not establish a complete relativistic plasma solver, net propulsion, or chamber validity\"\n"
        "}\n",
        checks,
        failures,
        (unsigned long long)deterministic_hash,
        worst_field_invariant_relative,
        worst_tensor_covariance_relative,
        worst_four_force_covariance_relative,
        worst_mass_shell_relative,
        failures == 0U ? "true" : "false"
    ) > 0;
}
