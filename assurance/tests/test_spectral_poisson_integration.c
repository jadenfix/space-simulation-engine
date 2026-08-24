#include "spacewind/constants.h"
#include "spacewind/plasma.h"
#include "spacewind/spectral_poisson.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t transcript_hash = UINT64_C(14695981039346656037);

static void check_true(int condition, const char *name) {
    ++checks;
    transcript_hash ^= swa_fnv1a64(name, strlen(name));
    transcript_hash *= UINT64_C(1099511628211);
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}

static swa_spectral_poisson_input root_input(const sw_pic1d *pic) {
    swa_spectral_poisson_input input;
    input.samples = pic->grid_count;
    input.length_m = pic->length_m;
    input.permittivity_F_m = SW_EPS0;
    input.charge_density_C_m3 = pic->charge_density_c_m3;
    input.potential_V = pic->potential_v;
    input.electric_field_V_m = pic->electric_field_v_m;
    return input;
}

static void set_manufactured_charge(sw_pic1d *pic, int add_mean) {
    size_t j;
    const double wave_a = 2.0 * SW_PI * 3.0 / pic->length_m;
    const double wave_b = 2.0 * SW_PI * 11.0 / pic->length_m;
    for (j = 0U; j < pic->grid_count; ++j) {
        const double x = pic->length_m * (double)j /
                         (double)pic->grid_count;
        pic->charge_density_c_m3[j] =
            3.0e-7 * cos(wave_a * x) -
            1.5e-7 * sin(wave_b * x) +
            (add_mean ? 2.0e-7 : 0.0);
    }
}

static void test_root_manufactured_modes(
    swa_spectral_poisson_result *receipt_result,
    swa_spectral_poisson_input *receipt_input,
    swa_spectral_poisson_tolerances *receipt_tolerances,
    double *receipt_charge,
    double *receipt_potential,
    double *receipt_electric
) {
    enum { SAMPLES = 64 };
    sw_pic1d pic;
    swa_spectral_poisson_input input;
    swa_spectral_poisson_tolerances tolerances;
    swa_spectral_poisson_result result;
    swa_spectral_poisson_result corrupted;
    size_t j;

    if (!sw_pic1d_init(&pic, SAMPLES, 1U, 12.8, 1.0e-7)) {
        check_true(0, "initialize root spectral integration grid");
        return;
    }
    check_true(1, "initialize root spectral integration grid");
    swa_default_spectral_poisson_tolerances(&tolerances);
    set_manufactured_charge(&pic, 0);
    sw_pic1d_solve_poisson_spectral(&pic);
    input = root_input(&pic);
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &result
        ),
        "root manufactured spectral audit executes"
    );
    check_true(result.passes,
               "root manufactured spectral solve passes");
    check_true(result.poisson_passes,
               "root potential closes spectral Poisson equation");
    check_true(result.gradient_passes,
               "root electric field is negative spectral gradient");
    check_true(result.gauss_passes,
               "root electric field closes spectral Gauss law");
    check_true(result.electrostatic_energy_passes,
               "root rho-phi and field energy agree");

    pic.electric_field_v_m[9] += 100.0;
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &corrupted
        ),
        "root electric corruption audit executes"
    );
    check_true(!corrupted.gradient_passes,
               "root electric corruption is detected");
    check_true(!corrupted.gauss_passes,
               "root electric corruption breaks Gauss law");
    pic.electric_field_v_m[9] -= 100.0;

    set_manufactured_charge(&pic, 1);
    sw_pic1d_solve_poisson_spectral(&pic);
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &corrupted
        ),
        "root nonneutral spectral audit executes"
    );
    check_true(!corrupted.neutrality_passes,
               "root skipped zero mode is reported explicitly");
    check_true(corrupted.poisson_passes,
               "resolved nonzero root modes remain correct");
    check_true(!corrupted.passes,
               "nonneutral root solution is not promoted");

    set_manufactured_charge(&pic, 0);
    sw_pic1d_solve_poisson_spectral(&pic);
    input = root_input(&pic);
    *receipt_result = result;
    *receipt_input = input;
    *receipt_tolerances = tolerances;
    for (j = 0U; j < SAMPLES; ++j) {
        receipt_charge[j] = pic.charge_density_c_m3[j];
        receipt_potential[j] = pic.potential_v[j];
        receipt_electric[j] = pic.electric_field_v_m[j];
    }
    receipt_input->charge_density_C_m3 = receipt_charge;
    receipt_input->potential_V = receipt_potential;
    receipt_input->electric_field_V_m = receipt_electric;
    sw_pic1d_destroy(&pic);
}

static void test_root_nyquist_boundary(void) {
    enum { SAMPLES = 64 };
    sw_pic1d pic;
    swa_spectral_poisson_input input;
    swa_spectral_poisson_tolerances tolerances;
    swa_spectral_poisson_result result;
    double maximum_electric = 0.0;
    size_t j;

    if (!sw_pic1d_init(&pic, SAMPLES, 1U, 6.4, 1.0e-7)) {
        check_true(0, "initialize root Nyquist grid");
        return;
    }
    check_true(1, "initialize root Nyquist grid");
    for (j = 0U; j < SAMPLES; ++j) {
        pic.charge_density_c_m3[j] =
            (j & 1U) == 0U ? 2.0e-7 : -2.0e-7;
    }
    sw_pic1d_solve_poisson_spectral(&pic);
    for (j = 0U; j < SAMPLES; ++j) {
        maximum_electric = fmax(
            maximum_electric,
            fabs(pic.electric_field_v_m[j])
        );
    }
    input = root_input(&pic);
    swa_default_spectral_poisson_tolerances(&tolerances);
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &result
        ),
        "root Nyquist spectral audit executes"
    );
    check_true(result.poisson_passes,
               "root Nyquist potential closes Poisson equation");
    check_true(maximum_electric < 1.0e-9,
               "collocated root Nyquist derivative vanishes");
    check_true(!result.nyquist_passes,
               "root Nyquist charge is marked unresolved");
    check_true(!result.electrostatic_energy_passes,
               "root Nyquist energy mismatch is exposed");
    check_true(!result.passes,
               "root Nyquist state cannot pass the audit");
    sw_pic1d_destroy(&pic);
}

static void test_actual_root_pic_state(void) {
    enum { SAMPLES = 64, PARTICLES = 1024, STEPS = 25 };
    sw_pic1d pic;
    swa_spectral_poisson_input input;
    swa_spectral_poisson_tolerances tolerances;
    swa_spectral_poisson_result result;
    size_t step;

    if (!sw_pic1d_init(&pic, SAMPLES, PARTICLES,
                       1000.0, 1.0e-7)) {
        check_true(0, "initialize actual root PIC spectral audit");
        return;
    }
    check_true(1, "initialize actual root PIC spectral audit");
    sw_pic1d_quiet_start(&pic, -SW_QE, SW_ME,
                         1.0e6, 1.0e5, 0.01, 1U);
    for (step = 0U; step < STEPS; ++step) {
        sw_pic1d_step(&pic);
    }
    input = root_input(&pic);
    swa_default_spectral_poisson_tolerances(&tolerances);
    tolerances.nyquist_charge_power_fraction_tolerance = 1.0e-8;
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &result
        ),
        "actual root PIC spectral audit executes"
    );
    check_true(result.neutrality_passes,
               "actual root PIC remains neutral");
    check_true(result.poisson_passes,
               "actual root PIC Poisson solve closes");
    check_true(result.gradient_passes,
               "actual root PIC electric gradient closes");
    check_true(result.gauss_passes,
               "actual root PIC Gauss law closes");
    check_true(result.parseval_passes,
               "actual root PIC Parseval identity closes");
    check_true(result.passes,
               "actual bounded root PIC state passes spectral audit");
    sw_pic1d_destroy(&pic);
}

int main(void) {
    enum { RECEIPT_SAMPLES = 64 };
    double receipt_charge[RECEIPT_SAMPLES];
    double receipt_potential[RECEIPT_SAMPLES];
    double receipt_electric[RECEIPT_SAMPLES];
    swa_spectral_poisson_input receipt_input;
    swa_spectral_poisson_tolerances receipt_tolerances;
    swa_spectral_poisson_result receipt_result;
    FILE *fp;

    test_root_manufactured_modes(
        &receipt_result,
        &receipt_input,
        &receipt_tolerances,
        receipt_charge,
        receipt_potential,
        receipt_electric
    );
    test_root_nyquist_boundary();
    test_actual_root_pic_state();

    fp = fopen("output/root_spectral_poisson_receipt.json", "w");
    check_true(fp != NULL, "open root spectral Poisson receipt");
    if (fp != NULL) {
        check_true(
            swa_write_spectral_poisson_receipt(
                fp,
                &receipt_input,
                &receipt_tolerances,
                &receipt_result
            ),
            "write root spectral Poisson receipt"
        );
        (void)fclose(fp);
    }

    printf(
        "spacewind root spectral Poisson integration: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
