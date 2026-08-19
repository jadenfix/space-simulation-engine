#include "spacewind/open_boundary.h"

#include <float.h>
#include <math.h>
#include <string.h>

static double max2(double a, double b) { return a > b ? a : b; }

static int finite_vec3(swa_vec3 value) {
    return isfinite(value.x) && isfinite(value.y) && isfinite(value.z);
}

static int normalize(swa_vec3 input, swa_vec3 *out) {
    const double norm = swa_vnorm(input);
    if (out == NULL || !finite_vec3(input) ||
        !isfinite(norm) || !(norm > 0.0)) {
        return 0;
    }
    *out = swa_vscale(input, 1.0 / norm);
    return finite_vec3(*out);
}

int swa_em_boundary_flux_density(
    const swa_em_fields *fields,
    swa_vec3 outward_normal,
    swa_em_boundary_flux *out
) {
    swa_tensor4 stress_energy;
    swa_covariant_em_diagnostics diagnostics;
    swa_vec3 unit;
    swa_four flux = {{0.0, 0.0, 0.0, 0.0}};
    size_t nu;
    size_t i;
    if (out == NULL || !normalize(outward_normal, &unit) ||
        !swa_em_stress_energy_tensor(fields, &stress_energy,
                                     &diagnostics)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    for (nu = 0U; nu < 4U; ++nu) {
        for (i = 0U; i < 3U; ++i) {
            const double normal_component =
                i == 0U ? unit.x : (i == 1U ? unit.y : unit.z);
            flux.component[nu] += normal_component *
                                  stress_energy.component[i + 1U][nu];
        }
    }
    out->unit_normal = unit;
    out->outward_four_momentum_flux_N_m2 = flux;
    out->outward_power_W_m2 = SWA_C * flux.component[0];
    out->outward_momentum_flux_N_m2 = swa_v3(
        flux.component[1], flux.component[2], flux.component[3]
    );
    out->finite = swa_four_is_finite(flux) &&
                  finite_vec3(out->unit_normal) &&
                  isfinite(out->outward_power_W_m2) &&
                  finite_vec3(out->outward_momentum_flux_N_m2);
    return out->finite;
}

int swa_integrate_em_boundary_flux(
    const swa_em_boundary_flux *flux,
    double area_m2,
    double duration_s,
    swa_four *out_four_impulse_Ns
) {
    double measure;
    if (flux == NULL || out_four_impulse_Ns == NULL || !flux->finite ||
        !swa_four_is_finite(flux->outward_four_momentum_flux_N_m2) ||
        !isfinite(area_m2) || !(area_m2 >= 0.0) ||
        !isfinite(duration_s) || !(duration_s >= 0.0)) {
        return 0;
    }
    measure = area_m2 * duration_s;
    if (!isfinite(measure)) {
        return 0;
    }
    *out_four_impulse_Ns = swa_four_scale(
        flux->outward_four_momentum_flux_N_m2,
        measure
    );
    return swa_four_is_finite(*out_four_impulse_Ns);
}

int swa_particle_boundary_crossing(
    double rest_mass_kg,
    swa_vec3 velocity_mps,
    double macro_weight,
    swa_vec3 outward_normal,
    swa_particle_boundary_transport *out
) {
    swa_vec3 unit;
    swa_four momentum;
    double normal_velocity;
    if (out == NULL || !normalize(outward_normal, &unit) ||
        !isfinite(rest_mass_kg) || !(rest_mass_kg > 0.0) ||
        !isfinite(macro_weight) || !(macro_weight > 0.0) ||
        !finite_vec3(velocity_mps) ||
        !swa_particle_four_momentum(rest_mass_kg,
                                    velocity_mps, &momentum)) {
        return 0;
    }
    memset(out, 0, sizeof(*out));
    momentum = swa_four_scale(momentum, macro_weight);
    normal_velocity = swa_vdot(velocity_mps, unit);
    if (!isfinite(normal_velocity)) {
        return 0;
    }
    out->normal_velocity_mps = normal_velocity;
    if (normal_velocity > 0.0) {
        out->outward_Ns = momentum;
        out->is_outward = 1;
    } else if (normal_velocity < 0.0) {
        out->inward_Ns = momentum;
        out->is_inward = 1;
    }
    out->finite = swa_four_is_finite(out->outward_Ns) &&
                  swa_four_is_finite(out->inward_Ns) &&
                  isfinite(out->normal_velocity_mps);
    return out->finite;
}

int swa_audit_open_boundary_four_momentum(
    swa_four field_initial_Ns,
    swa_four field_final_Ns,
    swa_four matter_initial_Ns,
    swa_four matter_final_Ns,
    swa_four external_impulse_Ns,
    swa_four field_outward_Ns,
    swa_four particle_outward_Ns,
    swa_four particle_inward_Ns,
    double relative_tolerance,
    double absolute_tolerance_Ns,
    swa_four_momentum_ledger_result *out
) {
    swa_four_momentum_ledger ledger;
    if (out == NULL || !swa_four_is_finite(field_initial_Ns) ||
        !swa_four_is_finite(field_final_Ns) ||
        !swa_four_is_finite(matter_initial_Ns) ||
        !swa_four_is_finite(matter_final_Ns) ||
        !swa_four_is_finite(external_impulse_Ns) ||
        !swa_four_is_finite(field_outward_Ns) ||
        !swa_four_is_finite(particle_outward_Ns) ||
        !swa_four_is_finite(particle_inward_Ns) ||
        !isfinite(relative_tolerance) || relative_tolerance < 0.0 ||
        !isfinite(absolute_tolerance_Ns) || absolute_tolerance_Ns < 0.0) {
        return 0;
    }
    ledger.field_initial_Ns = field_initial_Ns;
    ledger.field_final_Ns = field_final_Ns;
    ledger.matter_initial_Ns = matter_initial_Ns;
    ledger.matter_final_Ns = matter_final_Ns;
    ledger.external_impulse_Ns = external_impulse_Ns;
    ledger.momentum_in_Ns = particle_inward_Ns;
    ledger.momentum_out_Ns = swa_four_add(
        field_outward_Ns, particle_outward_Ns
    );
    return swa_audit_four_momentum_ledger(
        &ledger, relative_tolerance, absolute_tolerance_Ns, out
    );
}

int swa_write_open_boundary_receipt(
    FILE *fp,
    const swa_em_boundary_flux *right_wave,
    const swa_em_boundary_flux *left_wave,
    double worst_traction_relative,
    double worst_control_volume_relative,
    size_t checks,
    size_t failures,
    uint64_t deterministic_hash
) {
    if (fp == NULL || right_wave == NULL || left_wave == NULL ||
        !right_wave->finite || !left_wave->finite ||
        !isfinite(worst_traction_relative) ||
        !isfinite(worst_control_volume_relative)) {
        return 0;
    }
    return fprintf(
        fp,
        "{\n"
        "  \"schema\": \"spacewind.open-boundary-four-momentum/v1\",\n"
        "  \"checks\": %zu,\n"
        "  \"failures\": %zu,\n"
        "  \"deterministic_hash\": \"%016llx\",\n"
        "  \"right_wave_power_W_m2\": %.17g,\n"
        "  \"right_wave_momentum_flux_N_m2\": [%.17g, %.17g, %.17g],\n"
        "  \"left_wave_power_W_m2\": %.17g,\n"
        "  \"left_wave_momentum_flux_N_m2\": [%.17g, %.17g, %.17g],\n"
        "  \"worst_traction_relative\": %.17g,\n"
        "  \"worst_control_volume_relative\": %.17g,\n"
        "  \"passes\": %s,\n"
        "  \"nonclaim\": \"open-boundary stress-energy bookkeeping is a finite assurance layer; it does not establish a multidimensional plasma-wing force or propulsion\"\n"
        "}\n",
        checks,
        failures,
        (unsigned long long)deterministic_hash,
        right_wave->outward_power_W_m2,
        right_wave->outward_momentum_flux_N_m2.x,
        right_wave->outward_momentum_flux_N_m2.y,
        right_wave->outward_momentum_flux_N_m2.z,
        left_wave->outward_power_W_m2,
        left_wave->outward_momentum_flux_N_m2.x,
        left_wave->outward_momentum_flux_N_m2.y,
        left_wave->outward_momentum_flux_N_m2.z,
        worst_traction_relative,
        worst_control_volume_relative,
        failures == 0U ? "true" : "false"
    ) > 0;
}
