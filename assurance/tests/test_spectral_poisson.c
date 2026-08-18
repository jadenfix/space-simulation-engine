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

static void manufacture_superposition(
    size_t samples,
    double length_m,
    double scale,
    double *charge,
    double *potential,
    double *electric
) {
    size_t j;
    const double amplitude_cos = 2.5e-7 * scale;
    const double amplitude_sin = -1.7e-7 * scale;
    const double wave_cos = 2.0 * SWA_PI * 3.0 / length_m;
    const double wave_sin = 2.0 * SWA_PI * 7.0 / length_m;
    for (j = 0U; j < samples; ++j) {
        const double x = length_m * (double)j /
                         (double)samples;
        charge[j] =
            amplitude_cos * cos(wave_cos * x) +
            amplitude_sin * sin(wave_sin * x);
        potential[j] =
            amplitude_cos * cos(wave_cos * x) /
                (SWA_EPS0 * wave_cos * wave_cos) +
            amplitude_sin * sin(wave_sin * x) /
                (SWA_EPS0 * wave_sin * wave_sin);
        electric[j] =
            amplitude_cos * sin(wave_cos * x) /
                (SWA_EPS0 * wave_cos) -
            amplitude_sin * cos(wave_sin * x) /
                (SWA_EPS0 * wave_sin);
    }
}

static void test_manufactured_and_corruptions(
    swa_spectral_poisson_result *receipt_result,
    swa_spectral_poisson_input *receipt_input,
    swa_spectral_poisson_tolerances *receipt_tolerances,
    double *receipt_charge,
    double *receipt_potential,
    double *receipt_electric
) {
    enum { SAMPLES = 64 };
    const double length = 7.3;
    double charge[SAMPLES];
    double potential[SAMPLES];
    double electric[SAMPLES];
    swa_spectral_poisson_input input;
    swa_spectral_poisson_tolerances tolerances;
    swa_spectral_poisson_result result;
    swa_spectral_poisson_result gauge_result;
    swa_spectral_poisson_result corrupted;
    size_t j;

    manufacture_superposition(
        SAMPLES, length, 1.0,
        charge, potential, electric
    );
    input.samples = SAMPLES;
    input.length_m = length;
    input.permittivity_F_m = SWA_EPS0;
    input.charge_density_C_m3 = charge;
    input.potential_V = potential;
    input.electric_field_V_m = electric;
    swa_default_spectral_poisson_tolerances(&tolerances);

    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &result
        ),
        "manufactured spectral audit executes"
    );
    check_true(result.finite,
               "manufactured spectral audit finite");
    check_true(result.passes,
               "manufactured spectral Poisson system passes");
    check_true(result.neutrality_passes,
               "manufactured charge is neutral");
    check_true(result.nyquist_passes,
               "manufactured charge avoids Nyquist mode");
    check_true(result.poisson_passes,
               "manufactured Poisson equation closes");
    check_true(result.gradient_passes,
               "manufactured electric gradient closes");
    check_true(result.gauss_passes,
               "manufactured spectral Gauss law closes");
    check_true(result.parseval_passes,
               "manufactured Parseval energy closes");
    check_true(result.electrostatic_energy_passes,
               "manufactured rho-phi energy closes");

    for (j = 0U; j < SAMPLES; ++j) {
        potential[j] += 1234.5;
    }
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &gauge_result
        ),
        "gauge-shifted spectral audit executes"
    );
    check_true(gauge_result.passes,
               "constant potential gauge is ignored");
    for (j = 0U; j < SAMPLES; ++j) {
        potential[j] -= 1234.5;
    }

    electric[11] += 50.0;
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &corrupted
        ),
        "electric-corruption audit executes"
    );
    check_true(!corrupted.gradient_passes,
               "electric corruption breaks gradient identity");
    check_true(!corrupted.gauss_passes,
               "electric corruption breaks Gauss identity");
    check_true(!corrupted.passes,
               "electric corruption is rejected");
    electric[11] -= 50.0;

    potential[23] += 10.0;
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &corrupted
        ),
        "potential-corruption audit executes"
    );
    check_true(!corrupted.poisson_passes,
               "potential corruption breaks Poisson identity");
    check_true(!corrupted.passes,
               "potential corruption is rejected");
    potential[23] -= 10.0;

    charge[3] += 1.0e-5;
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &corrupted
        ),
        "charge-corruption audit executes"
    );
    check_true(!corrupted.poisson_passes,
               "charge corruption breaks Poisson identity");
    check_true(!corrupted.gauss_passes,
               "charge corruption breaks Gauss identity");
    check_true(!corrupted.passes,
               "charge corruption is rejected");
    charge[3] -= 1.0e-5;

    *receipt_result = result;
    *receipt_input = input;
    *receipt_tolerances = tolerances;
    memcpy(receipt_charge, charge, sizeof(charge));
    memcpy(receipt_potential, potential, sizeof(potential));
    memcpy(receipt_electric, electric, sizeof(electric));
    receipt_input->charge_density_C_m3 = receipt_charge;
    receipt_input->potential_V = receipt_potential;
    receipt_input->electric_field_V_m = receipt_electric;
}

static void test_scaling(void) {
    enum { SAMPLES = 63 };
    const double length = 9.1;
    const double scales[] = {1.0e-6, 1.0, 1.0e6};
    double charge[SAMPLES];
    double potential[SAMPLES];
    double electric[SAMPLES];
    swa_spectral_poisson_input input;
    swa_spectral_poisson_tolerances tolerances;
    swa_spectral_poisson_result result;
    size_t scale_index;

    swa_default_spectral_poisson_tolerances(&tolerances);
    input.samples = SAMPLES;
    input.length_m = length;
    input.permittivity_F_m = SWA_EPS0;
    input.charge_density_C_m3 = charge;
    input.potential_V = potential;
    input.electric_field_V_m = electric;
    for (scale_index = 0U;
         scale_index < sizeof(scales) / sizeof(scales[0]);
         ++scale_index) {
        manufacture_superposition(
            SAMPLES, length, scales[scale_index],
            charge, potential, electric
        );
        check_true(
            swa_audit_periodic_spectral_poisson(
                &input, &tolerances, &result
            ),
            "scaled spectral audit executes"
        );
        check_true(result.passes,
                   "scaled spectral identity remains closed");
        check_true(!result.nyquist_present,
                   "odd grid has no Nyquist singleton");
    }
}

static void test_mean_and_harmonic_fail_closed(void) {
    enum { SAMPLES = 64 };
    const double length = 8.0;
    double charge[SAMPLES];
    double potential[SAMPLES];
    double electric[SAMPLES];
    swa_spectral_poisson_input input;
    swa_spectral_poisson_tolerances tolerances;
    swa_spectral_poisson_result result;
    size_t j;

    manufacture_superposition(
        SAMPLES, length, 1.0,
        charge, potential, electric
    );
    input.samples = SAMPLES;
    input.length_m = length;
    input.permittivity_F_m = SWA_EPS0;
    input.charge_density_C_m3 = charge;
    input.potential_V = potential;
    input.electric_field_V_m = electric;
    swa_default_spectral_poisson_tolerances(&tolerances);

    for (j = 0U; j < SAMPLES; ++j) {
        charge[j] += 1.0e-6;
    }
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &result
        ),
        "nonneutral spectral audit executes"
    );
    check_true(!result.neutrality_passes,
               "unresolved net charge is rejected");
    check_true(!result.passes,
               "nonneutral periodic Poisson result is not promoted");
    for (j = 0U; j < SAMPLES; ++j) {
        charge[j] -= 1.0e-6;
        electric[j] += 100.0;
    }
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &result
        ),
        "harmonic-field spectral audit executes"
    );
    check_true(!result.harmonic_field_passes,
               "undeclared uniform electric field is rejected");
    check_true(!result.passes,
               "harmonic field blocks Poisson promotion");
}

static void test_nyquist_boundary(void) {
    enum { SAMPLES = 64 };
    const double length = 6.4;
    const double amplitude = 3.0e-7;
    const double wave_number =
        SWA_PI * (double)SAMPLES / length;
    double charge[SAMPLES];
    double potential[SAMPLES];
    double electric[SAMPLES];
    swa_spectral_poisson_input input;
    swa_spectral_poisson_tolerances tolerances;
    swa_spectral_poisson_result result;
    size_t j;

    for (j = 0U; j < SAMPLES; ++j) {
        charge[j] = (j & 1U) == 0U ? amplitude : -amplitude;
        potential[j] = charge[j] /
            (SWA_EPS0 * wave_number * wave_number);
        electric[j] = 0.0;
    }
    input.samples = SAMPLES;
    input.length_m = length;
    input.permittivity_F_m = SWA_EPS0;
    input.charge_density_C_m3 = charge;
    input.potential_V = potential;
    input.electric_field_V_m = electric;
    swa_default_spectral_poisson_tolerances(&tolerances);
    check_true(
        swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &result
        ),
        "Nyquist spectral audit executes"
    );
    check_true(result.nyquist_present,
               "even grid reports Nyquist singleton");
    check_true(result.poisson_passes,
               "Nyquist potential still closes Poisson equation");
    check_true(result.nyquist_charge_power_fraction > 0.999999,
               "pure Nyquist charge is identified");
    check_true(!result.nyquist_passes,
               "unresolved collocated Nyquist charge is rejected");
    check_true(!result.electrostatic_energy_passes,
               "missing Nyquist electric field breaks energy identity");
    check_true(!result.passes,
               "Nyquist boundary cannot be promoted");
}

static void test_invalid_inputs(void) {
    double charge[4] = {0.0, 0.0, 0.0, 0.0};
    double potential[4] = {0.0, 0.0, 0.0, 0.0};
    double electric[4] = {0.0, 0.0, 0.0, 0.0};
    swa_spectral_poisson_input input;
    swa_spectral_poisson_tolerances tolerances;
    swa_spectral_poisson_result result;

    input.samples = 4U;
    input.length_m = 1.0;
    input.permittivity_F_m = SWA_EPS0;
    input.charge_density_C_m3 = charge;
    input.potential_V = potential;
    input.electric_field_V_m = electric;
    swa_default_spectral_poisson_tolerances(&tolerances);

    input.length_m = 0.0;
    check_true(
        !swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &result
        ),
        "zero spectral domain rejected"
    );
    input.length_m = 1.0;
    tolerances.nyquist_charge_power_fraction_tolerance = 2.0;
    check_true(
        !swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &result
        ),
        "invalid Nyquist tolerance rejected"
    );
    tolerances.nyquist_charge_power_fraction_tolerance = 0.0;
    electric[0] = NAN;
    check_true(
        !swa_audit_periodic_spectral_poisson(
            &input, &tolerances, &result
        ),
        "nonfinite spectral field rejected"
    );
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

    test_manufactured_and_corruptions(
        &receipt_result,
        &receipt_input,
        &receipt_tolerances,
        receipt_charge,
        receipt_potential,
        receipt_electric
    );
    test_scaling();
    test_mean_and_harmonic_fail_closed();
    test_nyquist_boundary();
    test_invalid_inputs();

    fp = fopen("output/spectral_poisson_receipt.json", "w");
    check_true(fp != NULL, "open spectral Poisson receipt");
    if (fp != NULL) {
        check_true(
            swa_write_spectral_poisson_receipt(
                fp,
                &receipt_input,
                &receipt_tolerances,
                &receipt_result
            ),
            "write spectral Poisson receipt"
        );
        (void)fclose(fp);
    }

    printf(
        "spacewind spectral Poisson assurance: %zu checks, %zu failures, hash=%016llx\n",
        checks,
        failures,
        (unsigned long long)transcript_hash
    );
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
