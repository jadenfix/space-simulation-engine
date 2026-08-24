#include "spacewind/invariants.h"

#include <float.h>
#include <math.h>
#include <string.h>

static double max2(double a, double b) { return a > b ? a : b; }
static double max3(double a, double b, double c) { return max2(max2(a, b), c); }
static int finite_v3(swa_vec3 v) {
    return isfinite(v.x) && isfinite(v.y) && isfinite(v.z);
}
static double relative_error(double observed, double expected) {
    return fabs(observed - expected) / max3(fabs(observed), fabs(expected), DBL_MIN);
}

double swa_minkowski_norm(swa_four_vector v) {
    return -v.t * v.t + v.x * v.x + v.y * v.y + v.z * v.z;
}

int swa_four_velocity_from_three_velocity(swa_vec3 velocity_mps, swa_four_vector *out) {
    const double speed2 = swa_vdot(velocity_mps, velocity_mps);
    double gamma;
    if (out == NULL || !finite_v3(velocity_mps) || !(speed2 >= 0.0) ||
        speed2 >= SWA_C * SWA_C) {
        return 0;
    }
    gamma = 1.0 / sqrt(1.0 - speed2 / (SWA_C * SWA_C));
    out->t = gamma * SWA_C;
    out->x = gamma * velocity_mps.x;
    out->y = gamma * velocity_mps.y;
    out->z = gamma * velocity_mps.z;
    return isfinite(out->t) && isfinite(out->x) && isfinite(out->y) && isfinite(out->z);
}

int swa_audit_four_velocity(swa_four_vector u, double relative_tolerance) {
    const double norm = swa_minkowski_norm(u);
    if (!(relative_tolerance >= 0.0) || !isfinite(norm)) {
        return 0;
    }
    return relative_error(norm, -SWA_C * SWA_C) <= relative_tolerance;
}

int swa_audit_em_state(const swa_em_state *s, double relative_tolerance,
                       swa_em_audit *o) {
    double e2, b2, magnetic_power_scale;
    if (s == NULL || o == NULL || !(relative_tolerance >= 0.0) ||
        !finite_v3(s->electric_V_m) || !finite_v3(s->magnetic_T) ||
        !finite_v3(s->velocity_mps) || !isfinite(s->charge_C)) {
        return 0;
    }
    memset(o, 0, sizeof(*o));
    o->electric_force_N = swa_vscale(s->electric_V_m, s->charge_C);
    o->magnetic_force_N = swa_vscale(swa_vcross(s->velocity_mps, s->magnetic_T), s->charge_C);
    o->lorentz_force_N = swa_vadd(o->electric_force_N, o->magnetic_force_N);
    o->poynting_W_m2 = swa_vscale(swa_vcross(s->electric_V_m, s->magnetic_T), 1.0 / SWA_MU0);
    e2 = swa_vdot(s->electric_V_m, s->electric_V_m);
    b2 = swa_vdot(s->magnetic_T, s->magnetic_T);
    o->field_energy_density_J_m3 = 0.5 * (SWA_EPS0 * e2 + b2 / SWA_MU0);
    o->invariant_B2_minus_E2_over_c2_T2 = b2 - e2 / (SWA_C * SWA_C);
    o->invariant_E_dot_B_over_c = swa_vdot(s->electric_V_m, s->magnetic_T) / SWA_C;
    o->electric_particle_power_W = swa_vdot(o->electric_force_N, s->velocity_mps);
    o->magnetic_particle_power_W = swa_vdot(o->magnetic_force_N, s->velocity_mps);
    o->total_particle_power_W = swa_vdot(o->lorentz_force_N, s->velocity_mps);
    magnetic_power_scale = max3(
        fabs(s->charge_C) * swa_vnorm(s->velocity_mps) *
            swa_vnorm(swa_vcross(s->velocity_mps, s->magnetic_T)),
        fabs(o->total_particle_power_W),
        1.0);
    o->magnetic_power_relative = fabs(o->magnetic_particle_power_W) / magnetic_power_scale;
    o->finite = finite_v3(o->lorentz_force_N) && finite_v3(o->poynting_W_m2) &&
                isfinite(o->field_energy_density_J_m3) &&
                isfinite(o->invariant_B2_minus_E2_over_c2_T2) &&
                isfinite(o->invariant_E_dot_B_over_c) &&
                isfinite(o->magnetic_power_relative);
    o->magnetic_work_free = o->finite && o->magnetic_power_relative <= relative_tolerance;
    return 1;
}

static void parker_sample_residuals(const swa_parker_sample *s,
                                    double *motional_relative,
                                    double *e_dot_b_relative,
                                    double *e_dot_v_relative) {
    const swa_vec3 expected_e = swa_vscale(swa_vcross(s->bulk_velocity_mps, s->magnetic_field_T), -1.0);
    const double e_scale = max3(swa_vnorm(s->electric_field_V_m), swa_vnorm(expected_e), DBL_MIN);
    const double eb_scale = max3(swa_vnorm(s->electric_field_V_m) * swa_vnorm(s->magnetic_field_T),
                                 fabs(swa_vdot(s->electric_field_V_m, s->magnetic_field_T)), DBL_MIN);
    const double ev_scale = max3(swa_vnorm(s->electric_field_V_m) * swa_vnorm(s->bulk_velocity_mps),
                                 fabs(swa_vdot(s->electric_field_V_m, s->bulk_velocity_mps)), DBL_MIN);
    *motional_relative = swa_vnorm(swa_vsub(s->electric_field_V_m, expected_e)) / e_scale;
    *e_dot_b_relative = fabs(swa_vdot(s->electric_field_V_m, s->magnetic_field_T)) / eb_scale;
    *e_dot_v_relative = fabs(swa_vdot(s->electric_field_V_m, s->bulk_velocity_mps)) / ev_scale;
}

int swa_audit_parker_pair(const swa_parker_sample *a, const swa_parker_sample *b,
                          double relative_tolerance, swa_parker_audit *o) {
    double mass_flux_a, mass_flux_b, magnetic_flux_a, magnetic_flux_b;
    if (a == NULL || b == NULL || o == NULL || !(relative_tolerance >= 0.0) ||
        !(a->radius_m > 0.0) || !(b->radius_m > 0.0) ||
        !(a->number_density_m3 > 0.0) || !(b->number_density_m3 > 0.0) ||
        !finite_v3(a->bulk_velocity_mps) || !finite_v3(b->bulk_velocity_mps) ||
        !finite_v3(a->magnetic_field_T) || !finite_v3(b->magnetic_field_T) ||
        !finite_v3(a->electric_field_V_m) || !finite_v3(b->electric_field_V_m)) {
        return 0;
    }
    memset(o, 0, sizeof(*o));
    mass_flux_a = a->number_density_m3 * a->radial_speed_mps * a->radius_m * a->radius_m;
    mass_flux_b = b->number_density_m3 * b->radial_speed_mps * b->radius_m * b->radius_m;
    magnetic_flux_a = a->radial_magnetic_field_T * a->radius_m * a->radius_m;
    magnetic_flux_b = b->radial_magnetic_field_T * b->radius_m * b->radius_m;
    o->mass_flux_relative_error = relative_error(mass_flux_a, mass_flux_b);
    o->magnetic_flux_relative_error = relative_error(magnetic_flux_a, magnetic_flux_b);
    parker_sample_residuals(a, &o->motional_field_relative_error_a,
                            &o->electric_dot_magnetic_relative_a,
                            &o->electric_dot_velocity_relative_a);
    parker_sample_residuals(b, &o->motional_field_relative_error_b,
                            &o->electric_dot_magnetic_relative_b,
                            &o->electric_dot_velocity_relative_b);
    o->finite = isfinite(o->mass_flux_relative_error) &&
                isfinite(o->magnetic_flux_relative_error) &&
                isfinite(o->motional_field_relative_error_a) &&
                isfinite(o->motional_field_relative_error_b) &&
                isfinite(o->electric_dot_magnetic_relative_a) &&
                isfinite(o->electric_dot_magnetic_relative_b) &&
                isfinite(o->electric_dot_velocity_relative_a) &&
                isfinite(o->electric_dot_velocity_relative_b);
    o->passes = o->finite &&
                o->mass_flux_relative_error <= relative_tolerance &&
                o->magnetic_flux_relative_error <= relative_tolerance &&
                o->motional_field_relative_error_a <= relative_tolerance &&
                o->motional_field_relative_error_b <= relative_tolerance &&
                o->electric_dot_magnetic_relative_a <= relative_tolerance &&
                o->electric_dot_magnetic_relative_b <= relative_tolerance &&
                o->electric_dot_velocity_relative_a <= relative_tolerance &&
                o->electric_dot_velocity_relative_b <= relative_tolerance;
    return 1;
}

int swa_audit_photon_force(double irradiance_W_m2, double area_m2,
                           double reflectivity, double requested_force_N,
                           double relative_tolerance, swa_force_bound *o) {
    double momentum_coefficient;
    if (o == NULL || !(irradiance_W_m2 >= 0.0) || !(area_m2 >= 0.0) ||
        !(reflectivity >= 0.0 && reflectivity <= 1.0) ||
        !(requested_force_N >= 0.0) || !(relative_tolerance >= 0.0)) {
        return 0;
    }
    momentum_coefficient = 1.0 + reflectivity;
    o->requested_force_N = requested_force_N;
    o->upper_bound_force_N = momentum_coefficient * irradiance_W_m2 * area_m2 / SWA_C;
    o->utilization = requested_force_N / max2(o->upper_bound_force_N, DBL_MIN);
    o->finite = isfinite(o->upper_bound_force_N) && isfinite(o->utilization);
    o->passes = o->finite && requested_force_N <=
                o->upper_bound_force_N * (1.0 + relative_tolerance);
    return 1;
}

int swa_audit_plasma_momentum_force(double mass_density_kg_m3,
                                    double relative_speed_mps,
                                    double interaction_area_m2,
                                    double momentum_multiplier,
                                    double requested_force_N,
                                    double relative_tolerance,
                                    swa_force_bound *o) {
    if (o == NULL || !(mass_density_kg_m3 >= 0.0) || !(relative_speed_mps >= 0.0) ||
        !(interaction_area_m2 >= 0.0) || !(momentum_multiplier >= 0.0) ||
        !(requested_force_N >= 0.0) || !(relative_tolerance >= 0.0)) {
        return 0;
    }
    o->requested_force_N = requested_force_N;
    o->upper_bound_force_N = momentum_multiplier * mass_density_kg_m3 *
                             relative_speed_mps * relative_speed_mps *
                             interaction_area_m2;
    o->utilization = requested_force_N / max2(o->upper_bound_force_N, DBL_MIN);
    o->finite = isfinite(o->upper_bound_force_N) && isfinite(o->utilization);
    o->passes = o->finite && requested_force_N <=
                o->upper_bound_force_N * (1.0 + relative_tolerance);
    return 1;
}

int swa_compute_orbit_invariants(double mu, swa_vec3 r, swa_vec3 v,
                                 swa_orbit_invariants *o) {
    double radius, speed2, rv, energy;
    swa_vec3 h, eccentricity;
    if (o == NULL || !(mu > 0.0) || !finite_v3(r) || !finite_v3(v)) {
        return 0;
    }
    radius = swa_vnorm(r);
    if (!(radius > 0.0)) {
        return 0;
    }
    memset(o, 0, sizeof(*o));
    speed2 = swa_vdot(v, v);
    rv = swa_vdot(r, v);
    energy = 0.5 * speed2 - mu / radius;
    h = swa_vcross(r, v);
    eccentricity = swa_vsub(
        swa_vscale(r, (speed2 - mu / radius) / mu),
        swa_vscale(v, rv / mu));
    o->specific_energy_J_kg = energy;
    o->specific_angular_momentum_m2_s = h;
    o->eccentricity_vector = eccentricity;
    o->eccentricity = swa_vnorm(eccentricity);
    o->bound = energy < 0.0;
    if (energy != 0.0) {
        o->semimajor_axis_m = -mu / (2.0 * energy);
    } else {
        o->semimajor_axis_m = INFINITY;
    }
    if (o->bound) {
        o->periapsis_m = o->semimajor_axis_m * (1.0 - o->eccentricity);
        o->apoapsis_m = o->semimajor_axis_m * (1.0 + o->eccentricity);
    } else {
        o->periapsis_m = (swa_vdot(h, h) / mu) / (1.0 + o->eccentricity);
        o->apoapsis_m = INFINITY;
    }
    o->finite = isfinite(o->specific_energy_J_kg) && finite_v3(h) &&
                finite_v3(eccentricity) && isfinite(o->eccentricity) &&
                isfinite(o->periapsis_m) &&
                (isfinite(o->apoapsis_m) || isinf(o->apoapsis_m));
    return 1;
}

int swa_compute_nbody_invariants(size_t count, const double *mass_kg,
                                 const swa_vec3 *position_m,
                                 const swa_vec3 *velocity_mps,
                                 double gravitational_constant,
                                 double softening_m,
                                 swa_nbody_invariants *o) {
    size_t i, j;
    double total_mass = 0.0, kinetic = 0.0, potential = 0.0;
    swa_vec3 weighted_position = {0.0, 0.0, 0.0};
    swa_vec3 momentum = {0.0, 0.0, 0.0};
    swa_vec3 angular = {0.0, 0.0, 0.0};
    if (o == NULL || count == 0U || mass_kg == NULL || position_m == NULL ||
        velocity_mps == NULL || !(gravitational_constant > 0.0) ||
        !(softening_m >= 0.0)) {
        return 0;
    }
    memset(o, 0, sizeof(*o));
    for (i = 0; i < count; ++i) {
        const double m = mass_kg[i];
        if (!(m > 0.0) || !finite_v3(position_m[i]) || !finite_v3(velocity_mps[i])) {
            return 0;
        }
        total_mass += m;
        weighted_position = swa_vadd(weighted_position, swa_vscale(position_m[i], m));
        momentum = swa_vadd(momentum, swa_vscale(velocity_mps[i], m));
        angular = swa_vadd(angular, swa_vscale(swa_vcross(position_m[i], velocity_mps[i]), m));
        kinetic += 0.5 * m * swa_vdot(velocity_mps[i], velocity_mps[i]);
    }
    for (i = 0; i < count; ++i) {
        for (j = i + 1U; j < count; ++j) {
            const swa_vec3 dr = swa_vsub(position_m[j], position_m[i]);
            const double distance = sqrt(swa_vdot(dr, dr) + softening_m * softening_m);
            if (!(distance > 0.0)) {
                return 0;
            }
            potential -= gravitational_constant * mass_kg[i] * mass_kg[j] / distance;
        }
    }
    o->total_mass_kg = total_mass;
    o->center_of_mass_m = swa_vscale(weighted_position, 1.0 / total_mass);
    o->total_momentum_Ns = momentum;
    o->angular_momentum_kg_m2_s = angular;
    o->kinetic_energy_J = kinetic;
    o->potential_energy_J = potential;
    o->total_energy_J = kinetic + potential;
    o->finite = isfinite(total_mass) && finite_v3(o->center_of_mass_m) &&
                finite_v3(momentum) && finite_v3(angular) &&
                isfinite(kinetic) && isfinite(potential) && isfinite(o->total_energy_J);
    return 1;
}
