#include "spacewind/fields.h"
#include "spacewind/constants.h"
#include "spacewind/rng.h"

#include <float.h>
#include <math.h>
#include <stddef.h>
#include <string.h>

static double sw_clamp(double x, double lo, double hi) {
    return fmax(lo, fmin(hi, x));
}

void sw_solar_wind_default(sw_solar_wind_model *model) {
    memset(model, 0, sizeof(*model));
    model->central_mu = SW_SOLAR_MU;
    model->proton_temperature_k = 8.0e5;
    model->electron_temperature_k = 8.0e5;
    model->reference_radius_m = SW_AU;
    model->reference_number_density_m3 = 5.0e6;
    model->reference_radial_field_t = 3.0e-9;
    model->rotation_rate_rad_s = SW_SOLAR_ROTATION_RATE;
    model->source_surface_radius_m = 2.5 * SW_SOLAR_RADIUS;
    model->stream_fraction = 0.25;
    model->stream_arms = 2U;
    model->stream_phase_rad = 0.0;
    model->stream_latitude_width_rad = 0.8;
    model->shear_enabled = false;
    model->shear_normal = sw_v3(0.0, 1.0, 0.0);
    model->shear_width_m = 2.0e8;
    model->shear_flow_direction = sw_v3(1.0, 0.0, 0.0);
    model->cme_enabled = false;
    model->cme_launch_radius_m = 20.0 * SW_SOLAR_RADIUS;
    model->cme_speed_m_s = 900000.0;
    model->cme_width_m = 0.03 * SW_AU;
    model->cme_density_multiplier = 4.0;
    model->cme_field_multiplier = 3.0;
    model->cme_speed_increment_m_s = 300000.0;
}

void sw_solar_wind_init_turbulence(
    sw_solar_wind_model *model,
    uint64_t seed,
    size_t mode_count,
    double rms_magnetic_field_t,
    double min_wavelength_m,
    double max_wavelength_m
) {
    sw_rng rng;
    size_t i;
    double amplitude_norm = 0.0;

    mode_count = mode_count > SW_MAX_TURBULENCE_MODES ? SW_MAX_TURBULENCE_MODES : mode_count;
    if (mode_count == 0U || min_wavelength_m <= 0.0 || max_wavelength_m < min_wavelength_m) {
        model->turbulence_mode_count = 0U;
        return;
    }

    sw_rng_seed(&rng, seed, 17U);
    model->turbulence_mode_count = mode_count;
    for (i = 0U; i < mode_count; ++i) {
        const double u = (mode_count == 1U) ? 0.5 : (double)i / (double)(mode_count - 1U);
        const double wavelength = min_wavelength_m
            * pow(max_wavelength_m / min_wavelength_m, u);
        const double cos_theta = 2.0 * sw_rng_uniform(&rng) - 1.0;
        const double sin_theta = sqrt(fmax(0.0, 1.0 - cos_theta * cos_theta));
        const double phi = SW_TWO_PI * sw_rng_uniform(&rng);
        const sw_vec3 khat = sw_v3(
            sin_theta * cos(phi),
            sin_theta * sin(phi),
            cos_theta
        );
        sw_vec3 trial = fabs(khat.z) < 0.8 ? sw_v3(0.0, 0.0, 1.0) : sw_v3(0.0, 1.0, 0.0);
        sw_vec3 polarization = sw_v3_normalize(sw_v3_cross(khat, trial));
        const double raw_amplitude = pow(wavelength / min_wavelength_m, 1.0 / 6.0);

        model->turbulence_modes[i].wave_vector = sw_v3_scale(khat, SW_TWO_PI / wavelength);
        model->turbulence_modes[i].polarization = polarization;
        model->turbulence_modes[i].angular_frequency = 400000.0 * SW_TWO_PI / wavelength;
        model->turbulence_modes[i].phase = SW_TWO_PI * sw_rng_uniform(&rng);
        model->turbulence_modes[i].magnetic_amplitude_t = raw_amplitude;
        model->turbulence_modes[i].velocity_sign = sw_rng_uniform(&rng) < 0.5 ? -1.0 : 1.0;
        amplitude_norm += raw_amplitude * raw_amplitude;
    }

    amplitude_norm = sqrt(amplitude_norm);
    if (amplitude_norm > 0.0) {
        for (i = 0U; i < mode_count; ++i) {
            model->turbulence_modes[i].magnetic_amplitude_t *= rms_magnetic_field_t / amplitude_norm;
        }
    }
}

double sw_parker_sound_speed(const sw_solar_wind_model *model) {
    const double total_temperature = model->proton_temperature_k + model->electron_temperature_k;
    return sqrt(SW_KB * total_temperature / SW_MP);
}

double sw_parker_critical_radius(const sw_solar_wind_model *model) {
    const double cs = sw_parker_sound_speed(model);
    return model->central_mu / (2.0 * cs * cs);
}

static double sw_parker_equation(double y, double rhs) {
    return y - log(y) - rhs;
}

double sw_parker_speed(const sw_solar_wind_model *model, double radius_m) {
    const double cs = sw_parker_sound_speed(model);
    const double rc = sw_parker_critical_radius(model);
    double rhs;
    double lo;
    double hi;
    double f_lo;
    double f_hi;
    unsigned iter;

    if (radius_m <= 0.0 || cs <= 0.0 || rc <= 0.0) {
        return 0.0;
    }
    if (fabs(radius_m - rc) <= 1.0e-12 * rc) {
        return cs;
    }

    rhs = 4.0 * log(radius_m / rc) + 4.0 * rc / radius_m - 3.0;
    rhs = fmax(rhs, 1.0);

    if (radius_m < rc) {
        lo = DBL_MIN * 16.0;
        hi = 1.0;
    } else {
        lo = 1.0;
        hi = fmax(2.0, rhs + log(rhs + 1.0) + 2.0);
        while (sw_parker_equation(hi, rhs) <= 0.0 && hi < 1.0e12) {
            hi *= 2.0;
        }
    }

    f_lo = sw_parker_equation(lo, rhs);
    f_hi = sw_parker_equation(hi, rhs);
    if (f_lo == 0.0) {
        return cs * sqrt(lo);
    }
    if (f_hi == 0.0) {
        return cs * sqrt(hi);
    }

    for (iter = 0U; iter < 120U; ++iter) {
        const double mid = 0.5 * (lo + hi);
        const double f_mid = sw_parker_equation(mid, rhs);
        if (fabs(f_mid) < 1.0e-13) {
            return cs * sqrt(mid);
        }
        if ((f_lo > 0.0 && f_mid > 0.0) || (f_lo < 0.0 && f_mid < 0.0)) {
            lo = mid;
            f_lo = f_mid;
        } else {
            hi = mid;
            f_hi = f_mid;
        }
    }
    (void)f_hi;
    return cs * sqrt(0.5 * (lo + hi));
}

static void sw_spherical_basis(sw_vec3 position, sw_vec3 *er, sw_vec3 *etheta, sw_vec3 *ephi, double *theta, double *phi) {
    const double r = sw_v3_norm(position);
    const double rho_xy = hypot(position.x, position.y);
    *er = r > 0.0 ? sw_v3_scale(position, 1.0 / r) : sw_v3(1.0, 0.0, 0.0);
    *theta = atan2(rho_xy, position.z);
    *phi = atan2(position.y, position.x);
    *etheta = sw_v3(cos(*theta) * cos(*phi), cos(*theta) * sin(*phi), -sin(*theta));
    *ephi = sw_v3(-sin(*phi), cos(*phi), 0.0);
}

bool sw_environment_evaluate(
    const sw_solar_wind_model *model,
    double time_s,
    sw_vec3 position_m,
    sw_environment_sample *sample
) {
    sw_vec3 er;
    sw_vec3 etheta;
    sw_vec3 ephi;
    double theta;
    double phi;
    double base_speed;
    double ref_speed;
    double stream_factor = 1.0;
    double speed;
    double density;
    double br;
    double bphi;
    sw_vec3 velocity;
    sw_vec3 magnetic;
    size_t i;

    memset(sample, 0, sizeof(*sample));
    sample->radius_m = sw_v3_norm(position_m);
    if (sample->radius_m <= model->source_surface_radius_m || sample->radius_m <= 0.0) {
        return false;
    }

    sw_spherical_basis(position_m, &er, &etheta, &ephi, &theta, &phi);
    (void)etheta;
    base_speed = sw_parker_speed(model, sample->radius_m);
    ref_speed = sw_parker_speed(model, model->reference_radius_m);

    if (model->stream_arms > 0U && model->stream_fraction != 0.0) {
        const double latitude = 0.5 * SW_PI - theta;
        const double latitude_envelope = exp(
            -0.5 * (latitude / fmax(model->stream_latitude_width_rad, 1.0e-6))
            * (latitude / fmax(model->stream_latitude_width_rad, 1.0e-6))
        );
        const double phase = (double)model->stream_arms
            * (phi - model->rotation_rate_rad_s * time_s) + model->stream_phase_rad;
        stream_factor += model->stream_fraction * latitude_envelope * sin(phase);
    }
    stream_factor = sw_clamp(stream_factor, 0.2, 3.0);
    speed = base_speed * stream_factor;
    velocity = sw_v3_scale(er, speed);

    sample->shear_normal = sw_v3_normalize(model->shear_normal);
    sample->shear_coordinate_m = sw_v3_dot(position_m, sample->shear_normal) - model->shear_offset_m;
    if (model->shear_enabled) {
        const double blend = 0.5 * (
            1.0 + tanh(sample->shear_coordinate_m / fmax(model->shear_width_m, 1.0))
        );
        const sw_vec3 direction = sw_v3_normalize(model->shear_flow_direction);
        velocity = sw_v3_add(velocity, sw_v3_scale(direction, model->shear_delta_speed_m_s * blend));
    }

    density = model->reference_number_density_m3
        * pow(model->reference_radius_m / sample->radius_m, 2.0)
        * ref_speed / fmax(base_speed, 1.0);
    density /= stream_factor;

    br = model->reference_radial_field_t
        * pow(model->reference_radius_m / sample->radius_m, 2.0);
    bphi = -br * model->rotation_rate_rad_s
        * (sample->radius_m - model->source_surface_radius_m)
        * sin(theta) / fmax(base_speed, 1.0);
    magnetic = sw_v3_add(sw_v3_scale(er, br), sw_v3_scale(ephi, bphi));

    if (model->cme_enabled && time_s >= model->cme_launch_time_s) {
        const double shell_radius = model->cme_launch_radius_m
            + model->cme_speed_m_s * (time_s - model->cme_launch_time_s);
        const double z = (sample->radius_m - shell_radius) / fmax(model->cme_width_m, 1.0);
        const double envelope = exp(-0.5 * z * z);
        density *= 1.0 + (model->cme_density_multiplier - 1.0) * envelope;
        magnetic = sw_v3_scale(
            magnetic,
            1.0 + (model->cme_field_multiplier - 1.0) * envelope
        );
        velocity = sw_v3_add(
            velocity,
            sw_v3_scale(er, model->cme_speed_increment_m_s * envelope)
        );
    }

    sample->mass_density_kg_m3 = density * SW_MP;
    for (i = 0U; i < model->turbulence_mode_count; ++i) {
        const sw_turbulence_mode *mode = &model->turbulence_modes[i];
        const double phase = sw_v3_dot(mode->wave_vector, position_m)
            - mode->angular_frequency * time_s + mode->phase;
        const sw_vec3 db = sw_v3_scale(
            mode->polarization,
            mode->magnetic_amplitude_t * cos(phase)
        );
        const double alfven_denominator = sqrt(
            fmax(SW_MU0 * sample->mass_density_kg_m3, DBL_MIN)
        );
        magnetic = sw_v3_add(magnetic, db);
        velocity = sw_v3_add(
            velocity,
            sw_v3_scale(db, mode->velocity_sign / alfven_denominator)
        );
    }

    sample->proton_number_density_m3 = density;
    sample->proton_temperature_k = model->proton_temperature_k;
    sample->electron_temperature_k = model->electron_temperature_k;
    sample->wind_velocity_m_s = velocity;
    sample->magnetic_field_t = magnetic;
    sample->electric_field_v_m = sw_v3_scale(sw_v3_cross(velocity, magnetic), -1.0);
    sample->photon_flux_w_m2 = SW_SOLAR_CONSTANT_1AU
        * pow(SW_AU / sample->radius_m, 2.0);
    sample->thermal_pressure_pa = density * SW_KB
        * (model->proton_temperature_k + model->electron_temperature_k);
    sample->dynamic_pressure_pa = sample->mass_density_kg_m3 * sw_v3_norm2(velocity);
    sample->magnetic_pressure_pa = sw_v3_norm2(magnetic) / (2.0 * SW_MU0);
    sample->debye_length_m = sqrt(
        SW_EPS0 * SW_KB * fmax(model->electron_temperature_k, 1.0)
        / fmax(density * SW_QE * SW_QE, DBL_MIN)
    );
    sample->alfven_speed_m_s = sw_v3_norm(magnetic)
        / sqrt(fmax(SW_MU0 * sample->mass_density_kg_m3, DBL_MIN));
    sample->sound_speed_m_s = sw_parker_sound_speed(model);
    sample->fast_magnetosonic_speed_m_s = sqrt(
        sample->sound_speed_m_s * sample->sound_speed_m_s
        + sample->alfven_speed_m_s * sample->alfven_speed_m_s
    );
    return sw_v3_isfinite(sample->wind_velocity_m_s)
        && sw_v3_isfinite(sample->magnetic_field_t)
        && isfinite(sample->proton_number_density_m3);
}

sw_vec3 sw_gravity_acceleration(
    const sw_gravity_system *system,
    sw_vec3 position_m
) {
    sw_vec3 acceleration = sw_v3(0.0, 0.0, 0.0);
    size_t i;
    for (i = 0U; i < system->count; ++i) {
        const sw_body *body = &system->bodies[i];
        const sw_vec3 delta = sw_v3_sub(body->position_m, position_m);
        const double r2 = sw_v3_norm2(delta);
        const double r = sqrt(r2);
        if (r > body->radius_m * 1.0e-9 && r2 > 0.0) {
            acceleration = sw_v3_add(
                acceleration,
                sw_v3_scale(delta, body->gravitational_mu / (r2 * r))
            );
        }
    }
    return acceleration;
}

sw_vec3 sw_gravity_acceleration_with_j2(
    const sw_gravity_system *system,
    sw_vec3 position_m
) {
    sw_vec3 acceleration = sw_gravity_acceleration(system, position_m);
    size_t i;
    for (i = 0U; i < system->count; ++i) {
        const sw_body *body = &system->bodies[i];
        if (body->j2 == 0.0 || body->radius_m <= 0.0) {
            continue;
        }
        {
            const sw_vec3 rel = sw_v3_sub(position_m, body->position_m);
            const double r = sw_v3_norm(rel);
            const sw_vec3 axis = sw_v3_normalize(body->spin_axis);
            if (r <= body->radius_m * 1.0e-9) {
                continue;
            }
            {
                const double z = sw_v3_dot(rel, axis);
                const double z2_over_r2 = (z * z) / (r * r);
                const double scale = 1.5 * body->j2 * body->gravitational_mu
                    * body->radius_m * body->radius_m / pow(r, 5.0);
                const sw_vec3 term = sw_v3_sub(
                    sw_v3_scale(rel, 5.0 * z2_over_r2 - 1.0),
                    sw_v3_scale(axis, 2.0 * z)
                );
                acceleration = sw_v3_add(acceleration, sw_v3_scale(term, scale));
            }
        }
    }
    return acceleration;
}

double sw_specific_orbital_energy(double central_mu, sw_vec3 position_m, sw_vec3 velocity_m_s) {
    const double r = sw_v3_norm(position_m);
    if (r <= 0.0) {
        return INFINITY;
    }
    return 0.5 * sw_v3_norm2(velocity_m_s) - central_mu / r;
}

sw_vec3 sw_specific_angular_momentum(sw_vec3 position_m, sw_vec3 velocity_m_s) {
    return sw_v3_cross(position_m, velocity_m_s);
}
