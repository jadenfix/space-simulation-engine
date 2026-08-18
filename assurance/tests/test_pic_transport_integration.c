#include "spacewind/charge_deposition.h"
#include "spacewind/constants.h"
#include "spacewind/plasma.h"

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);

static double max2(double a, double b) { return a > b ? a : b; }

static void check_true(int condition, const char *name) {
    ++checks;
    transcript_hash ^= swa_fnv1a64(name, strlen(name));
    transcript_hash *= UINT64_C(1099511628211);
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}

static double wrap_periodic(double x, double length) {
    double wrapped = fmod(x, length);
    if (wrapped < 0.0) {
        wrapped += length;
    }
    return wrapped;
}

static double nearest_periodic_displacement(
    double old_unwrapped,
    double new_wrapped,
    double length
) {
    const double old_wrapped = wrap_periodic(old_unwrapped, length);
    double displacement = new_wrapped - old_wrapped;
    if (displacement > 0.5 * length) {
        displacement -= length;
    } else if (displacement < -0.5 * length) {
        displacement += length;
    }
    return displacement;
}

static swa_periodic_cic_transport_input transport_input(
    size_t cells,
    size_t particles,
    double length,
    double dt,
    const double *initial,
    const double *final,
    const double *charge
) {
    swa_periodic_cic_transport_input input;
    memset(&input, 0, sizeof(input));
    input.cells = cells;
    input.particles = particles;
    input.domain_length_m = length;
    input.cross_section_area_m2 = 1.0;
    input.dt_s = dt;
    input.initial_position_unwrapped_m = initial;
    input.final_position_unwrapped_m = final;
    input.particle_charge_C = charge;
    input.continuity_relative_tolerance = 5.0e-12;
    input.continuity_absolute_tolerance_A_m3 = 1.0e-18;
    input.continuity_global_absolute_tolerance_C_m = 1.0e-22;
    input.charge_absolute_tolerance_C = 1.0e-22;
    input.mean_current_relative_tolerance = 5.0e-12;
    input.mean_current_absolute_tolerance_A_m2 = 1.0e-20;
    return input;
}

static void test_grid_staggering_contract(void) {
    enum { CELLS = 8, PARTICLES = 3 };
    sw_pic1d pic;
    double shifted_initial[PARTICLES];
    double shifted_final[PARTICLES];
    double unshifted_initial[PARTICLES];
    double unshifted_final[PARTICLES];
    double effective_charge[PARTICLES];
    double rho_initial[CELLS];
    double rho_final[CELLS];
    double current[CELLS + 1];
    double wrong_rho_initial[CELLS];
    double wrong_rho_final[CELLS];
    double wrong_current[CELLS + 1];
    swa_periodic_cic_transport_input input;
    swa_periodic_cic_transport_buffers buffers;
    swa_periodic_cic_transport_buffers wrong_buffers;
    swa_periodic_cic_transport_result result;
    swa_periodic_cic_transport_result wrong_result;
    double maximum_aligned_error = 0.0;
    double maximum_wrong_error = 0.0;
    size_t i;

    if (!sw_pic1d_init(&pic, CELLS, PARTICLES, 8.0, 0.1)) {
        check_true(0, "initialize staggering-contract PIC");
        return;
    }
    check_true(1, "initialize staggering-contract PIC");
    pic.neutralizing_background_c_m3 = 0.0;
    pic.particles[0].position_m = sw_v3(0.1, 0.0, 0.0);
    pic.particles[1].position_m = sw_v3(2.7, 0.0, 0.0);
    pic.particles[2].position_m = sw_v3(7.8, 0.0, 0.0);
    pic.particles[0].charge_c = 1.0e-9;
    pic.particles[1].charge_c = -0.7e-9;
    pic.particles[2].charge_c = 0.4e-9;
    pic.particles[0].macro_weight = 1.0;
    pic.particles[1].macro_weight = 1.3;
    pic.particles[2].macro_weight = 0.8;
    sw_pic1d_deposit_charge(&pic);

    for (i = 0U; i < PARTICLES; ++i) {
        const double position = pic.particles[i].position_m.x;
        shifted_initial[i] = position + 0.5 * pic.dx_m;
        shifted_final[i] = shifted_initial[i];
        unshifted_initial[i] = position;
        unshifted_final[i] = position;
        effective_charge[i] = pic.particles[i].charge_c *
                              pic.particles[i].macro_weight;
    }
    input = transport_input(CELLS, PARTICLES, pic.length_m,
                            pic.step_s, shifted_initial,
                            shifted_final, effective_charge);
    buffers.rho_initial_C_m3 = rho_initial;
    buffers.rho_final_C_m3 = rho_final;
    buffers.rho_capacity = CELLS;
    buffers.current_face_A_m2 = current;
    buffers.current_face_capacity = CELLS + 1U;
    check_true(
        swa_deposit_periodic_cic_transport(
            &input, &buffers, &result
        ),
        "aligned root/reference deposition executes"
    );
    check_true(result.passes,
               "aligned root/reference deposition closes");
    for (i = 0U; i < CELLS; ++i) {
        maximum_aligned_error = max2(
            maximum_aligned_error,
            fabs(pic.charge_density_c_m3[i] - rho_initial[i])
        );
    }
    check_true(maximum_aligned_error < 1.0e-23,
               "half-cell topology transform matches root nodal CIC");

    input = transport_input(CELLS, PARTICLES, pic.length_m,
                            pic.step_s, unshifted_initial,
                            unshifted_final, effective_charge);
    wrong_buffers.rho_initial_C_m3 = wrong_rho_initial;
    wrong_buffers.rho_final_C_m3 = wrong_rho_final;
    wrong_buffers.rho_capacity = CELLS;
    wrong_buffers.current_face_A_m2 = wrong_current;
    wrong_buffers.current_face_capacity = CELLS + 1U;
    check_true(
        swa_deposit_periodic_cic_transport(
            &input, &wrong_buffers, &wrong_result
        ),
        "misaligned reference deposition executes"
    );
    for (i = 0U; i < CELLS; ++i) {
        maximum_wrong_error = max2(
            maximum_wrong_error,
            fabs(pic.charge_density_c_m3[i] -
                 wrong_rho_initial[i])
        );
    }
    check_true(maximum_wrong_error > 1.0e-11,
               "omitting half-cell transform is detected");
    sw_pic1d_destroy(&pic);
}

static void test_root_pic_transport_sequence(void) {
    enum { CELLS = 32, PARTICLES = 512, STEPS = 100 };
    sw_pic1d pic;
    double unwrapped[PARTICLES];
    double reference_initial_position[PARTICLES];
    double reference_final_position[PARTICLES];
    double effective_charge[PARTICLES];
    double root_initial_density[CELLS];
    double reference_initial_density[CELLS];
    double reference_final_density[CELLS];
    double reference_current[CELLS + 1];
    swa_periodic_cic_transport_buffers buffers;
    swa_periodic_cic_transport_result result;
    double worst_initial_density_error = 0.0;
    double worst_final_density_error = 0.0;
    double worst_density_relative = 0.0;
    double worst_continuity_relative = 0.0;
    double worst_mean_current_error = 0.0;
    double maximum_unwrapped_displacement_cells = 0.0;
    double initial_energy;
    double final_energy;
    double relative_energy_change;
    size_t step;
    size_t p;
    size_t cell;
    int sequence_passes = 1;
    FILE *fp;

    if (!sw_pic1d_init(&pic, CELLS, PARTICLES,
                       1000.0, 1.0e-7)) {
        check_true(0, "initialize root PIC integration sequence");
        return;
    }
    check_true(1, "initialize root PIC integration sequence");
    sw_pic1d_quiet_start(&pic, -SW_QE, SW_ME,
                         1.0e6, 1.0e5, 0.01, 1U);
    sw_pic1d_deposit_charge(&pic);
    sw_pic1d_solve_poisson_spectral(&pic);
    initial_energy = sw_pic1d_kinetic_energy_j(&pic) +
                     sw_pic1d_field_energy_j(&pic);
    for (p = 0U; p < PARTICLES; ++p) {
        unwrapped[p] = pic.particles[p].position_m.x;
        effective_charge[p] = pic.particles[p].charge_c *
                              pic.particles[p].macro_weight;
    }
    buffers.rho_initial_C_m3 = reference_initial_density;
    buffers.rho_final_C_m3 = reference_final_density;
    buffers.rho_capacity = CELLS;
    buffers.current_face_A_m2 = reference_current;
    buffers.current_face_capacity = CELLS + 1U;

    for (step = 0U; step < STEPS; ++step) {
        swa_periodic_cic_transport_input input;
        memcpy(root_initial_density,
               pic.charge_density_c_m3,
               sizeof(root_initial_density));
        for (p = 0U; p < PARTICLES; ++p) {
            reference_initial_position[p] =
                unwrapped[p] + 0.5 * pic.dx_m;
        }

        sw_pic1d_step(&pic);

        for (p = 0U; p < PARTICLES; ++p) {
            const double displacement =
                nearest_periodic_displacement(
                    unwrapped[p],
                    pic.particles[p].position_m.x,
                    pic.length_m
                );
            unwrapped[p] += displacement;
            reference_final_position[p] =
                unwrapped[p] + 0.5 * pic.dx_m;
        }
        input = transport_input(
            CELLS, PARTICLES, pic.length_m, pic.step_s,
            reference_initial_position,
            reference_final_position,
            effective_charge
        );
        if (!swa_deposit_periodic_cic_transport(
                &input, &buffers, &result)) {
            sequence_passes = 0;
            check_true(0, "root PIC reference transport executes");
            break;
        }
        check_true(1, "root PIC reference transport executes");
        check_true(result.passes,
                   "root PIC reference transport closes");
        if (!result.passes) {
            sequence_passes = 0;
        }
        maximum_unwrapped_displacement_cells = max2(
            maximum_unwrapped_displacement_cells,
            result.maximum_particle_displacement_cells
        );
        worst_continuity_relative = max2(
            worst_continuity_relative,
            result.continuity.residual_relative_maximum
        );
        worst_mean_current_error = max2(
            worst_mean_current_error,
            fabs(result.mean_current_error_A_m2)
        );
        for (cell = 0U; cell < CELLS; ++cell) {
            const double root_initial_particle_density =
                root_initial_density[cell] -
                pic.neutralizing_background_c_m3;
            const double root_final_particle_density =
                pic.charge_density_c_m3[cell] -
                pic.neutralizing_background_c_m3;
            const double initial_error = fabs(
                root_initial_particle_density -
                reference_initial_density[cell]
            );
            const double final_error = fabs(
                root_final_particle_density -
                reference_final_density[cell]
            );
            const double density_scale = max2(
                fabs(pic.neutralizing_background_c_m3),
                DBL_MIN
            );
            worst_initial_density_error = max2(
                worst_initial_density_error,
                initial_error
            );
            worst_final_density_error = max2(
                worst_final_density_error,
                final_error
            );
            worst_density_relative = max2(
                worst_density_relative,
                max2(initial_error, final_error) / density_scale
            );
        }
        check_true(worst_density_relative < 2.0e-12,
                   "root and independent CIC densities agree");
        if (!(worst_density_relative < 2.0e-12)) {
            sequence_passes = 0;
        }
    }

    final_energy = sw_pic1d_kinetic_energy_j(&pic) +
                   sw_pic1d_field_energy_j(&pic);
    relative_energy_change =
        (final_energy - initial_energy) /
        max2(fabs(initial_energy), DBL_MIN);
    check_true(isfinite(initial_energy) && isfinite(final_energy) &&
               isfinite(relative_energy_change),
               "root PIC sequence energy diagnostic is finite");
    check_true(maximum_unwrapped_displacement_cells < 0.01,
               "nearest-image unwrapping remains unambiguous");
    check_true(sequence_passes,
               "complete root PIC transport sequence passes");

    fp = fopen("output/pic_transport_integration_receipt.json", "w");
    check_true(fp != NULL, "open PIC integration receipt");
    if (fp != NULL) {
        check_true(
            fprintf(
                fp,
                "{\n"
                "  \"schema\": \"spacewind.pic-transport-integration/v1\",\n"
                "  \"steps\": %u,\n"
                "  \"cells\": %u,\n"
                "  \"particles\": %u,\n"
                "  \"worst_initial_density_error_C_m3\": %.17g,\n"
                "  \"worst_final_density_error_C_m3\": %.17g,\n"
                "  \"worst_density_relative_to_background\": %.17g,\n"
                "  \"worst_continuity_relative\": %.17g,\n"
                "  \"worst_mean_current_error_A_m2\": %.17g,\n"
                "  \"maximum_unwrapped_displacement_cells\": %.17g,\n"
                "  \"relative_energy_change_diagnostic\": %.17g,\n"
                "  \"passes\": %s,\n"
                "  \"nonclaim\": \"this cross-layer check validates charge topology and continuity bookkeeping for the bounded electrostatic PIC sequence; it does not establish electromagnetic field-energy closure or plasma-wing thrust\"\n"
                "}\n",
                (unsigned)STEPS,
                (unsigned)CELLS,
                (unsigned)PARTICLES,
                worst_initial_density_error,
                worst_final_density_error,
                worst_density_relative,
                worst_continuity_relative,
                worst_mean_current_error,
                maximum_unwrapped_displacement_cells,
                relative_energy_change,
                sequence_passes ? "true" : "false"
            ) > 0,
            "write PIC integration receipt"
        );
        (void)fclose(fp);
    }
    sw_pic1d_destroy(&pic);
}

int main(void) {
    test_grid_staggering_contract();
    test_root_pic_transport_sequence();
    printf(
        "spacewind PIC transport integration: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
