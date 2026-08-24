#include "spacewind/invariants.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static size_t checks = 0U;
static size_t failures = 0U;
static uint64_t receipt_hash = UINT64_C(14695981039346656037);

static void mix_bytes(const void *data, size_t size) {
    receipt_hash ^= swa_fnv1a64(data, size);
    receipt_hash *= UINT64_C(1099511628211);
}

static void check_true(int condition, const char *name) {
    ++checks;
    mix_bytes(name, strlen(name));
    if (!condition) {
        ++failures;
        fprintf(stderr, "FAIL: %s\n", name);
    }
}

static int near(double a, double b, double relative, double absolute) {
    const double scale = fmax(fabs(a), fabs(b));
    return fabs(a - b) <= fmax(absolute, relative * scale);
}

static swa_vec3 random_vector(size_t i, double scale) {
    return swa_v3(scale * (-1.0 + 2.0 * swa_halton(i, 2U)),
                  scale * (-1.0 + 2.0 * swa_halton(i, 3U)),
                  scale * (-1.0 + 2.0 * swa_halton(i, 5U)));
}

static void test_four_velocity(void) {
    size_t i;
    for (i = 1U; i <= 2000U; ++i) {
        swa_vec3 direction = random_vector(i, 1.0);
        swa_four_vector u;
        const double norm = swa_vnorm(direction);
        const double speed = 0.95 * SWA_C * swa_halton(i, 7U);
        if (norm == 0.0) {
            direction = swa_v3(1.0, 0.0, 0.0);
        } else {
            direction = swa_vscale(direction, 1.0 / norm);
        }
        check_true(swa_four_velocity_from_three_velocity(swa_vscale(direction, speed), &u),
                   "construct subluminal four velocity");
        check_true(swa_audit_four_velocity(u, 2e-13),
                   "four velocity preserves Minkowski norm");
        mix_bytes(&u, sizeof(u));
    }
    {
        swa_four_vector u;
        check_true(!swa_four_velocity_from_three_velocity(swa_v3(SWA_C, 0.0, 0.0), &u),
                   "reject luminal three velocity");
        check_true(!swa_four_velocity_from_three_velocity(swa_v3(1.01 * SWA_C, 0.0, 0.0), &u),
                   "reject superluminal three velocity");
    }
}

static void test_electromagnetic_invariants(void) {
    size_t i;
    for (i = 1U; i <= 2000U; ++i) {
        swa_em_state s;
        swa_em_state rotated;
        swa_em_audit a;
        swa_em_audit b;
        const double angle = 2.0 * SWA_PI * swa_halton(i, 17U);
        s.electric_V_m = random_vector(i + 11U, 0.02);
        s.magnetic_T = random_vector(i + 23U, 15e-9);
        s.velocity_mps = random_vector(i + 37U, 7e5);
        s.charge_C = (i % 2U == 0U) ? SWA_QE : -SWA_QE;
        check_true(swa_audit_em_state(&s, 2e-14, &a), "audit electromagnetic state");
        check_true(a.finite && a.magnetic_work_free, "magnetic Lorentz force performs no work");
        check_true(near(a.total_particle_power_W, a.electric_particle_power_W, 2e-13, 1e-30),
                   "Lorentz power equals electric power");
        check_true(a.field_energy_density_J_m3 >= 0.0, "field energy density nonnegative");
        rotated = s;
        rotated.electric_V_m = swa_rotate_z(s.electric_V_m, angle);
        rotated.magnetic_T = swa_rotate_z(s.magnetic_T, angle);
        rotated.velocity_mps = swa_rotate_z(s.velocity_mps, angle);
        check_true(swa_audit_em_state(&rotated, 2e-14, &b), "audit rotated electromagnetic state");
        check_true(near(a.field_energy_density_J_m3, b.field_energy_density_J_m3, 2e-13, 1e-30),
                   "field energy rotation invariant");
        check_true(near(a.invariant_B2_minus_E2_over_c2_T2,
                        b.invariant_B2_minus_E2_over_c2_T2, 3e-13, 1e-32),
                   "first electromagnetic invariant rotation invariant");
        check_true(near(a.invariant_E_dot_B_over_c,
                        b.invariant_E_dot_B_over_c, 3e-13, 1e-32),
                   "second electromagnetic invariant rotation invariant");
        check_true(swa_vnorm(swa_vsub(swa_rotate_z(a.poynting_W_m2, angle), b.poynting_W_m2)) <=
                       3e-13 * fmax(swa_vnorm(a.poynting_W_m2), 1e-30),
                   "Poynting vector rotates covariantly");
        mix_bytes(&a, sizeof(a));
    }
}

static swa_parker_sample parker_sample(double radius_m, double density_1au,
                                       double radial_speed_mps,
                                       double radial_field_1au_T) {
    swa_parker_sample s;
    const double ratio = SWA_AU / radius_m;
    const double solar_rotation_rad_s = 2.86533e-6;
    const double br = radial_field_1au_T * ratio * ratio;
    const double bphi = -br * solar_rotation_rad_s * radius_m / radial_speed_mps;
    memset(&s, 0, sizeof(s));
    s.radius_m = radius_m;
    s.number_density_m3 = density_1au * ratio * ratio;
    s.radial_speed_mps = radial_speed_mps;
    s.radial_magnetic_field_T = br;
    s.bulk_velocity_mps = swa_v3(radial_speed_mps, 0.0, 0.0);
    s.magnetic_field_T = swa_v3(br, bphi, 0.0);
    s.electric_field_V_m = swa_vscale(swa_vcross(s.bulk_velocity_mps, s.magnetic_field_T), -1.0);
    return s;
}

static void test_parker_contract(void) {
    size_t i;
    for (i = 1U; i <= 1000U; ++i) {
        const double r1 = (0.2 + 4.0 * swa_halton(i, 2U)) * SWA_AU;
        const double r2 = (0.2 + 4.0 * swa_halton(i, 3U)) * SWA_AU;
        const double speed = 3e5 + 5e5 * swa_halton(i, 5U);
        swa_parker_sample a = parker_sample(r1, 5e6, speed, 3e-9);
        swa_parker_sample b = parker_sample(r2, 5e6, speed, 3e-9);
        swa_parker_audit result;
        check_true(swa_audit_parker_pair(&a, &b, 5e-13, &result),
                   "audit Parker sample pair");
        check_true(result.passes, "Parker mass, magnetic, and motional-field invariants close");
        b.number_density_m3 *= 1.01;
        check_true(swa_audit_parker_pair(&a, &b, 1e-4, &result) && !result.passes,
                   "Parker mass-flux mismatch detected");
        mix_bytes(&result, sizeof(result));
    }
}

static void test_force_bounds(void) {
    size_t i;
    for (i = 1U; i <= 1000U; ++i) {
        swa_force_bound bound;
        const double r = (0.3 + 8.0 * swa_halton(i, 2U)) * SWA_AU;
        const double irradiance = SWA_SOLAR_IRRADIANCE_1AU * (SWA_AU / r) * (SWA_AU / r);
        const double area = 1.0 + 10000.0 * swa_halton(i, 3U);
        const double reflectivity = swa_halton(i, 5U);
        const double photon_limit = (1.0 + reflectivity) * irradiance * area / SWA_C;
        check_true(swa_audit_photon_force(irradiance, area, reflectivity,
                                          0.999 * photon_limit, 1e-12, &bound) && bound.passes,
                   "photon momentum force below bound");
        check_true(swa_audit_photon_force(irradiance, area, reflectivity,
                                          1.01 * photon_limit, 1e-12, &bound) && !bound.passes,
                   "photon momentum excess rejected");
        {
            const double density = 1e-22 + 1e-19 * swa_halton(i, 7U);
            const double speed = 2e5 + 8e5 * swa_halton(i, 11U);
            const double interaction_area = 10.0 + 1e8 * swa_halton(i, 13U);
            const double multiplier = 2.0;
            const double plasma_limit = multiplier * density * speed * speed * interaction_area;
            check_true(swa_audit_plasma_momentum_force(density, speed, interaction_area,
                                                        multiplier, 0.99 * plasma_limit,
                                                        1e-12, &bound) && bound.passes,
                       "plasma momentum force below declared bound");
            check_true(swa_audit_plasma_momentum_force(density, speed, interaction_area,
                                                        multiplier, 1.02 * plasma_limit,
                                                        1e-12, &bound) && !bound.passes,
                       "plasma momentum excess rejected");
        }
    }
}

static void test_orbit_invariants(void) {
    size_t i;
    const double mu = 1.32712440018e20;
    for (i = 1U; i <= 1000U; ++i) {
        const double radius = (0.2 + 5.0 * swa_halton(i, 2U)) * SWA_AU;
        const double angle = 2.0 * SWA_PI * swa_halton(i, 3U);
        const swa_vec3 r = swa_rotate_z(swa_v3(radius, 0.0, 0.0), angle);
        const swa_vec3 v = swa_rotate_z(swa_v3(0.0, sqrt(mu / radius), 0.0), angle);
        swa_orbit_invariants a;
        swa_orbit_invariants b;
        const double extra_rotation = 2.0 * SWA_PI * swa_halton(i, 5U);
        check_true(swa_compute_orbit_invariants(mu, r, v, &a) && a.finite,
                   "compute circular-orbit invariants");
        check_true(a.bound && a.eccentricity < 2e-14,
                   "circular orbit is bound and nearly zero eccentricity");
        check_true(near(a.semimajor_axis_m, radius, 3e-14, 1e-4),
                   "circular semimajor axis equals radius");
        check_true(swa_compute_orbit_invariants(mu, swa_rotate_z(r, extra_rotation),
                                                swa_rotate_z(v, extra_rotation), &b),
                   "compute rotated orbit invariants");
        check_true(near(a.specific_energy_J_kg, b.specific_energy_J_kg, 3e-14, 1e-8),
                   "orbital energy rotation invariant");
        check_true(near(a.eccentricity, b.eccentricity, 1e-8, 2e-14),
                   "orbital eccentricity rotation invariant");
        mix_bytes(&a, sizeof(a));
    }
}

static void test_nbody_invariants(void) {
    size_t sample;
    for (sample = 1U; sample <= 500U; ++sample) {
        enum { N = 5 };
        double mass[N];
        swa_vec3 position[N];
        swa_vec3 velocity[N];
        swa_vec3 shifted_position[N];
        swa_vec3 boosted_velocity[N];
        double mass_reversed[N];
        swa_vec3 position_reversed[N];
        swa_vec3 velocity_reversed[N];
        swa_nbody_invariants base;
        swa_nbody_invariants shifted;
        swa_nbody_invariants boosted;
        swa_nbody_invariants reversed;
        const swa_vec3 translation = random_vector(sample + 4000U, 1e9);
        const swa_vec3 boost = random_vector(sample + 5000U, 1e3);
        size_t i;
        for (i = 0U; i < N; ++i) {
            mass[i] = 1e20 * (1.0 + 100.0 * swa_halton(sample + i * 11U, 2U));
            position[i] = random_vector(sample + i * 17U, 1e9);
            position[i].x += 3e9 * (double)i;
            velocity[i] = random_vector(sample + i * 23U, 1e4);
            shifted_position[i] = swa_vadd(position[i], translation);
            boosted_velocity[i] = swa_vadd(velocity[i], boost);
            mass_reversed[N - 1U - i] = mass[i];
            position_reversed[N - 1U - i] = position[i];
            velocity_reversed[N - 1U - i] = velocity[i];
        }
        check_true(swa_compute_nbody_invariants(N, mass, position, velocity,
                                                SWA_G, 100.0, &base) && base.finite,
                   "compute N-body invariants");
        check_true(swa_compute_nbody_invariants(N, mass, shifted_position, velocity,
                                                SWA_G, 100.0, &shifted),
                   "compute translated N-body invariants");
        check_true(near(base.total_energy_J, shifted.total_energy_J, 3e-13, 1e-3),
                   "N-body energy translation invariant");
        check_true(swa_vnorm(swa_vsub(shifted.center_of_mass_m,
                                      swa_vadd(base.center_of_mass_m, translation))) < 1e-5,
                   "center of mass translates covariantly");
        check_true(swa_vnorm(swa_vsub(base.total_momentum_Ns,
                                      shifted.total_momentum_Ns)) <=
                       2e-13 * fmax(swa_vnorm(base.total_momentum_Ns), 1.0),
                   "total momentum translation invariant");
        check_true(swa_compute_nbody_invariants(N, mass, position, boosted_velocity,
                                                SWA_G, 100.0, &boosted),
                   "compute Galilean-boosted N-body invariants");
        check_true(swa_vnorm(swa_vsub(boosted.total_momentum_Ns,
                                      swa_vadd(base.total_momentum_Ns,
                                               swa_vscale(boost, base.total_mass_kg)))) <=
                       4e-13 * fmax(swa_vnorm(boosted.total_momentum_Ns), 1.0),
                   "N-body momentum follows Galilean boost law");
        check_true(near(boosted.kinetic_energy_J,
                        base.kinetic_energy_J + swa_vdot(boost, base.total_momentum_Ns) +
                            0.5 * base.total_mass_kg * swa_vdot(boost, boost),
                        5e-13, 1e-3),
                   "N-body kinetic energy follows Galilean boost law");
        check_true(swa_compute_nbody_invariants(N, mass_reversed, position_reversed,
                                                velocity_reversed, SWA_G, 100.0, &reversed),
                   "compute permuted N-body invariants");
        check_true(near(base.total_energy_J, reversed.total_energy_J, 5e-13, 1e-3),
                   "N-body energy permutation invariant");
        check_true(swa_vnorm(swa_vsub(base.total_momentum_Ns,
                                      reversed.total_momentum_Ns)) <=
                       5e-13 * fmax(swa_vnorm(base.total_momentum_Ns), 1.0),
                   "N-body momentum permutation invariant");
        mix_bytes(&base, sizeof(base));
    }
}

int main(void) {
    FILE *fp;
    test_four_velocity();
    test_electromagnetic_invariants();
    test_parker_contract();
    test_force_bounds();
    test_orbit_invariants();
    test_nbody_invariants();
    fp = fopen("output/invariants_receipt.json", "w");
    if (fp != NULL) {
        (void)swa_write_assurance_receipt(fp, checks, failures, receipt_hash);
        (void)fclose(fp);
    } else {
        ++failures;
    }
    printf("spacewind invariants: %zu checks, %zu failures, hash=%016llx\n",
           checks, failures, (unsigned long long)receipt_hash);
    return failures == 0U ? EXIT_SUCCESS : EXIT_FAILURE;
}
