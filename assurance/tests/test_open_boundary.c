#include "spacewind/open_boundary.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);
static double worst_traction_relative = 0.0;
static double worst_control_volume_relative = 0.0;
static swa_em_boundary_flux right_receipt;
static swa_em_boundary_flux left_receipt;

static double max2(double a, double b) { return a > b ? a : b; }
static double max3(double a, double b, double c) { return max2(max2(a, b), c); }

static void check_true(int condition, const char *name) {
    ++checks;
    transcript_hash ^= swa_fnv1a64(name, strlen(name));
    transcript_hash *= UINT64_C(1099511628211);
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}

static swa_four zero_four(void) {
    return swa_four_make(0.0, 0.0, 0.0, 0.0);
}

static double four_relative(swa_four value, double scale) {
    return swa_four_euclidean_norm(value) / max2(scale, DBL_MIN);
}

static void test_local_flux_tensor_consistency(void) {
    size_t i;
    for (i = 1U; i <= 3000U; ++i) {
        swa_em_fields fields;
        swa_em_boundary_flux flux;
        swa_maxwell_traction_result traction;
        swa_local_em_result local;
        const swa_vec3 normal = swa_v3(
            0.1 + swa_halton(i, 17U),
            -0.8 + 1.6 * swa_halton(i, 19U),
            -0.7 + 1.4 * swa_halton(i, 23U)
        );
        double expected_power;
        double traction_error;
        double traction_scale;

        fields.electric_field_V_m = swa_v3(
            4.0e5 * (swa_halton(i, 2U) - 0.5),
            4.0e5 * (swa_halton(i, 3U) - 0.5),
            4.0e5 * (swa_halton(i, 5U) - 0.5)
        );
        fields.magnetic_field_T = swa_v3(
            4.0e-6 * (swa_halton(i, 7U) - 0.5),
            4.0e-6 * (swa_halton(i, 11U) - 0.5),
            4.0e-6 * (swa_halton(i, 13U) - 0.5)
        );
        check_true(
            swa_em_boundary_flux_density(&fields, normal, &flux),
            "boundary stress-energy flux executes"
        );
        check_true(
            swa_maxwell_traction(fields.electric_field_V_m,
                                 fields.magnetic_field_T,
                                 normal, &traction),
            "independent Maxwell traction executes"
        );
        check_true(
            swa_local_electromagnetic_state(fields.electric_field_V_m,
                                            fields.magnetic_field_T,
                                            &local),
            "independent local electromagnetic state executes"
        );
        expected_power = swa_vdot(local.poynting_flux_W_m2,
                                  flux.unit_normal);
        check_true(
            fabs(flux.outward_power_W_m2 - expected_power) <=
                256.0 * DBL_EPSILON *
                max3(fabs(flux.outward_power_W_m2),
                     fabs(expected_power),
                     swa_vnorm(local.poynting_flux_W_m2)),
            "stress-energy boundary power equals normal Poynting flux"
        );
        traction_error = swa_vnorm(swa_vadd(
            flux.outward_momentum_flux_N_m2,
            traction.traction_N_m2
        ));
        traction_scale = max3(
            swa_vnorm(flux.outward_momentum_flux_N_m2),
            swa_vnorm(traction.traction_N_m2),
            max2(local.energy_density_J_m3, DBL_MIN)
        );
        worst_traction_relative = max2(
            worst_traction_relative,
            traction_error / traction_scale
        );
        check_true(
            traction_error <= 512.0 * DBL_EPSILON * traction_scale,
            "stress-energy spatial flux is negative conventional Maxwell traction"
        );
    }

    {
        swa_em_fields fields;
        swa_em_boundary_flux flux;
        fields.electric_field_V_m = swa_v3(1.0, 0.0, 0.0);
        fields.magnetic_field_T = swa_v3(0.0, 0.0, 0.0);
        check_true(
            !swa_em_boundary_flux_density(
                &fields, swa_v3(0.0, 0.0, 0.0), &flux
            ),
            "zero boundary normal is rejected"
        );
    }
}

static void test_plane_wave_escape(int direction) {
    const double amplitude = 2.5e5;
    const double area = 3.0;
    const double duration = 2.0e-8;
    const double slab_length = SWA_C * duration;
    swa_em_fields fields;
    swa_em_boundary_flux flux;
    swa_tensor4 tensor;
    swa_covariant_em_diagnostics diagnostics;
    swa_four outward;
    swa_four field_initial;
    swa_four_momentum_ledger_result result;
    swa_vec3 normal;
    double volume;
    double poynting_x;
    double energy;
    double momentum_x;
    double physical_scale;

    fields.electric_field_V_m = swa_v3(0.0, amplitude, 0.0);
    fields.magnetic_field_T = swa_v3(
        0.0, 0.0, (double)direction * amplitude / SWA_C
    );
    normal = swa_v3((double)direction, 0.0, 0.0);
    check_true(
        swa_em_boundary_flux_density(&fields, normal, &flux),
        "plane-wave boundary flux executes"
    );
    check_true(
        swa_em_stress_energy_tensor(&fields, &tensor, &diagnostics),
        "plane-wave stress-energy tensor executes"
    );
    check_true(
        swa_integrate_em_boundary_flux(&flux, area, duration, &outward),
        "integrate plane-wave boundary four-flux"
    );
    volume = area * slab_length;
    energy = diagnostics.energy_density_J_m3 * volume;
    poynting_x = diagnostics.poynting_flux_W_m2.x;
    momentum_x = poynting_x * volume / (SWA_C * SWA_C);
    field_initial = swa_four_make(
        energy / SWA_C, momentum_x, 0.0, 0.0
    );
    physical_scale = max2(energy / SWA_C, DBL_MIN);
    check_true(
        swa_audit_open_boundary_four_momentum(
            field_initial, zero_four(), zero_four(), zero_four(),
            zero_four(), outward, zero_four(), zero_four(),
            2.0e-13, 1.0e-24, &result
        ),
        "plane-wave open control-volume audit executes"
    );
    worst_control_volume_relative = max2(
        worst_control_volume_relative,
        four_relative(result.residual_Ns, physical_scale)
    );
    check_true(result.passes,
               "plane-wave field four-momentum exits through boundary");
    check_true(flux.outward_power_W_m2 > 0.0,
               "outgoing plane wave has positive outward power");
    check_true(
        (double)direction * flux.outward_momentum_flux_N_m2.x > 0.0,
        "outgoing plane-wave momentum points out of domain"
    );

    if (direction > 0) {
        right_receipt = flux;
    } else {
        left_receipt = flux;
    }
}

static void test_particle_boundary_transport(void) {
    const double mass = 3.0e-20;
    const double weight = 7.0;
    const swa_vec3 outward_normal = swa_v3(1.0, 0.0, 0.0);
    const swa_vec3 outward_velocity = swa_v3(
        0.12 * SWA_C, 0.03 * SWA_C, -0.01 * SWA_C
    );
    const swa_vec3 inward_velocity = swa_v3(
        -0.08 * SWA_C, 0.02 * SWA_C, 0.01 * SWA_C
    );
    swa_particle_boundary_transport outward;
    swa_particle_boundary_transport inward;
    swa_four_momentum_ledger_result result;

    check_true(
        swa_particle_boundary_crossing(
            mass, outward_velocity, weight, outward_normal, &outward
        ),
        "outward particle crossing executes"
    );
    check_true(outward.is_outward && !outward.is_inward,
               "outward particle classified correctly");
    check_true(
        swa_audit_open_boundary_four_momentum(
            zero_four(), zero_four(),
            outward.outward_Ns, zero_four(), zero_four(),
            zero_four(), outward.outward_Ns, zero_four(),
            1.0e-14, 1.0e-24, &result
        ),
        "outward particle control-volume audit executes"
    );
    check_true(result.passes,
               "outward particle four-momentum leaves control volume");

    check_true(
        swa_particle_boundary_crossing(
            mass, inward_velocity, weight, outward_normal, &inward
        ),
        "inward particle crossing executes"
    );
    check_true(inward.is_inward && !inward.is_outward,
               "inward particle classified correctly");
    check_true(
        swa_audit_open_boundary_four_momentum(
            zero_four(), zero_four(),
            zero_four(), inward.inward_Ns, zero_four(),
            zero_four(), zero_four(), inward.inward_Ns,
            1.0e-14, 1.0e-24, &result
        ),
        "inward particle control-volume audit executes"
    );
    check_true(result.passes,
               "inward particle four-momentum enters control volume");

    check_true(
        !swa_particle_boundary_crossing(
            mass, swa_v3(SWA_C, 0.0, 0.0), weight,
            outward_normal, &outward
        ),
        "luminal particle boundary crossing rejected"
    );
}

static void test_combined_and_corrupted_control_volume(void) {
    const double amplitude = 1.8e5;
    const double area = 2.5;
    const double duration = 1.5e-8;
    const double volume = area * SWA_C * duration;
    const double mass = 2.0e-20;
    const double weight = 4.0;
    swa_em_fields fields;
    swa_em_boundary_flux flux;
    swa_tensor4 tensor;
    swa_covariant_em_diagnostics diagnostics;
    swa_particle_boundary_transport particle;
    swa_four field_out;
    swa_four field_initial;
    swa_four matter_initial;
    swa_four_momentum_ledger_result result;
    swa_four corrupted_field_out;
    double energy;
    double momentum_x;
    double scale;

    fields.electric_field_V_m = swa_v3(0.0, amplitude, 0.0);
    fields.magnetic_field_T = swa_v3(0.0, 0.0, amplitude / SWA_C);
    check_true(
        swa_em_boundary_flux_density(
            &fields, swa_v3(1.0, 0.0, 0.0), &flux
        ),
        "combined field boundary flux executes"
    );
    check_true(
        swa_em_stress_energy_tensor(&fields, &tensor, &diagnostics),
        "combined field stress-energy executes"
    );
    check_true(
        swa_integrate_em_boundary_flux(&flux, area, duration, &field_out),
        "combined field boundary impulse integrates"
    );
    energy = diagnostics.energy_density_J_m3 * volume;
    momentum_x = diagnostics.poynting_flux_W_m2.x * volume /
                 (SWA_C * SWA_C);
    field_initial = swa_four_make(energy / SWA_C, momentum_x, 0.0, 0.0);

    check_true(
        swa_particle_boundary_crossing(
            mass,
            swa_v3(0.10 * SWA_C, -0.02 * SWA_C, 0.01 * SWA_C),
            weight,
            swa_v3(1.0, 0.0, 0.0),
            &particle
        ),
        "combined outward particle crossing executes"
    );
    matter_initial = particle.outward_Ns;
    check_true(
        swa_audit_open_boundary_four_momentum(
            field_initial, zero_four(), matter_initial, zero_four(),
            zero_four(), field_out, particle.outward_Ns, zero_four(),
            2.0e-13, 1.0e-24, &result
        ),
        "combined field-particle boundary audit executes"
    );
    check_true(result.passes,
               "combined field and particle boundary transport closes");

    scale = max3(
        swa_four_euclidean_norm(field_initial),
        swa_four_euclidean_norm(matter_initial),
        1.0e-24
    );
    corrupted_field_out = field_out;
    corrupted_field_out.component[1] += 0.05 * scale;
    check_true(
        swa_audit_open_boundary_four_momentum(
            field_initial, zero_four(), matter_initial, zero_four(),
            zero_four(), corrupted_field_out,
            particle.outward_Ns, zero_four(),
            1.0e-6, 1.0e-24, &result
        ),
        "corrupted field-particle boundary audit executes"
    );
    check_true(!result.passes,
               "untracked Maxwell-stress boundary impulse is rejected");
}

int main(void) {
    FILE *fp;
    memset(&right_receipt, 0, sizeof(right_receipt));
    memset(&left_receipt, 0, sizeof(left_receipt));
    test_local_flux_tensor_consistency();
    test_plane_wave_escape(1);
    test_plane_wave_escape(-1);
    test_particle_boundary_transport();
    test_combined_and_corrupted_control_volume();

    fp = fopen("output/open_boundary_receipt.json", "w");
    check_true(fp != NULL, "open boundary receipt");
    if (fp != NULL) {
        check_true(
            swa_write_open_boundary_receipt(
                fp, &right_receipt, &left_receipt,
                worst_traction_relative,
                worst_control_volume_relative,
                checks, failures, transcript_hash
            ),
            "write open boundary receipt"
        );
        (void)fclose(fp);
    }

    printf(
        "spacewind open-boundary assurance: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
