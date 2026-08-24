#include "spacewind/nbody.h"
#include "spacewind/constants.h"

#include <math.h>
#include <stddef.h>

void sw_nbody_accelerations(sw_nbody_system *system) {
    size_t i;
    size_t j;
    double eps2;
    if (system == NULL || system->count > SW_NBODY_MAX) {
        return;
    }
    eps2 = system->softening_m * system->softening_m;
    for (i = 0U; i < system->count; ++i) {
        system->objects[i].acceleration_m_s2 = sw_v3(0.0, 0.0, 0.0);
    }
    for (i = 0U; i < system->count; ++i) {
        for (j = i + 1U; j < system->count; ++j) {
            sw_nbody_object *a = &system->objects[i];
            sw_nbody_object *b = &system->objects[j];
            const sw_vec3 delta = sw_v3_sub(b->position_m, a->position_m);
            const double distance2 = sw_v3_norm2(delta) + eps2;
            const double inv_distance = 1.0 / sqrt(distance2);
            const double inv_distance3 = inv_distance * inv_distance * inv_distance;
            const sw_vec3 direction_factor = sw_v3_scale(delta, SW_G * inv_distance3);
            a->acceleration_m_s2 = sw_v3_add(
                a->acceleration_m_s2,
                sw_v3_scale(direction_factor, b->mass_kg)
            );
            b->acceleration_m_s2 = sw_v3_sub(
                b->acceleration_m_s2,
                sw_v3_scale(direction_factor, a->mass_kg)
            );
        }
    }
}

bool sw_nbody_leapfrog_step(sw_nbody_system *system, double step_s) {
    size_t i;
    if (system == NULL || system->count == 0U || system->count > SW_NBODY_MAX
        || step_s <= 0.0 || !isfinite(step_s)) {
        return false;
    }
    sw_nbody_accelerations(system);
    for (i = 0U; i < system->count; ++i) {
        sw_nbody_object *object = &system->objects[i];
        object->velocity_m_s = sw_v3_add(
            object->velocity_m_s,
            sw_v3_scale(object->acceleration_m_s2, 0.5 * step_s)
        );
        object->position_m = sw_v3_add(
            object->position_m,
            sw_v3_scale(object->velocity_m_s, step_s)
        );
    }
    sw_nbody_accelerations(system);
    for (i = 0U; i < system->count; ++i) {
        sw_nbody_object *object = &system->objects[i];
        object->velocity_m_s = sw_v3_add(
            object->velocity_m_s,
            sw_v3_scale(object->acceleration_m_s2, 0.5 * step_s)
        );
        if (!sw_v3_isfinite(object->position_m) || !sw_v3_isfinite(object->velocity_m_s)) {
            return false;
        }
    }
    return true;
}

double sw_nbody_total_energy_j(const sw_nbody_system *system) {
    double energy = 0.0;
    size_t i;
    size_t j;
    double eps2;
    if (system == NULL || system->count > SW_NBODY_MAX) {
        return NAN;
    }
    eps2 = system->softening_m * system->softening_m;
    for (i = 0U; i < system->count; ++i) {
        const sw_nbody_object *object = &system->objects[i];
        energy += 0.5 * object->mass_kg * sw_v3_norm2(object->velocity_m_s);
    }
    for (i = 0U; i < system->count; ++i) {
        for (j = i + 1U; j < system->count; ++j) {
            const sw_vec3 delta = sw_v3_sub(
                system->objects[j].position_m,
                system->objects[i].position_m
            );
            energy -= SW_G * system->objects[i].mass_kg * system->objects[j].mass_kg
                / sqrt(sw_v3_norm2(delta) + eps2);
        }
    }
    return energy;
}

sw_vec3 sw_nbody_total_momentum_kg_m_s(const sw_nbody_system *system) {
    sw_vec3 momentum = sw_v3(0.0, 0.0, 0.0);
    size_t i;
    if (system == NULL || system->count > SW_NBODY_MAX) {
        return momentum;
    }
    for (i = 0U; i < system->count; ++i) {
        momentum = sw_v3_add(
            momentum,
            sw_v3_scale(system->objects[i].velocity_m_s, system->objects[i].mass_kg)
        );
    }
    return momentum;
}

sw_vec3 sw_nbody_center_of_mass_m(const sw_nbody_system *system) {
    sw_vec3 weighted = sw_v3(0.0, 0.0, 0.0);
    double total_mass = 0.0;
    size_t i;
    if (system == NULL || system->count > SW_NBODY_MAX) {
        return weighted;
    }
    for (i = 0U; i < system->count; ++i) {
        weighted = sw_v3_add(
            weighted,
            sw_v3_scale(system->objects[i].position_m, system->objects[i].mass_kg)
        );
        total_mass += system->objects[i].mass_kg;
    }
    return total_mass > 0.0 ? sw_v3_scale(weighted, 1.0 / total_mass) : sw_v3(0.0, 0.0, 0.0);
}

size_t sw_nbody_merge_overlaps(sw_nbody_system *system) {
    size_t merges = 0U;
    size_t i = 0U;
    if (system == NULL || system->count > SW_NBODY_MAX) {
        return 0U;
    }
    while (i < system->count) {
        size_t j = i + 1U;
        while (j < system->count) {
            sw_nbody_object *a = &system->objects[i];
            sw_nbody_object *b = &system->objects[j];
            const double separation = sw_v3_norm(sw_v3_sub(b->position_m, a->position_m));
            if (separation <= a->radius_m + b->radius_m) {
                const double total_mass = a->mass_kg + b->mass_kg;
                const sw_vec3 momentum = sw_v3_add(
                    sw_v3_scale(a->velocity_m_s, a->mass_kg),
                    sw_v3_scale(b->velocity_m_s, b->mass_kg)
                );
                const sw_vec3 center = sw_v3_scale(
                    sw_v3_add(
                        sw_v3_scale(a->position_m, a->mass_kg),
                        sw_v3_scale(b->position_m, b->mass_kg)
                    ),
                    1.0 / total_mass
                );
                a->position_m = center;
                a->velocity_m_s = sw_v3_scale(momentum, 1.0 / total_mass);
                a->mass_kg = total_mass;
                a->radius_m = cbrt(a->radius_m * a->radius_m * a->radius_m
                    + b->radius_m * b->radius_m * b->radius_m);
                system->objects[j] = system->objects[system->count - 1U];
                --system->count;
                ++merges;
            } else {
                ++j;
            }
        }
        ++i;
    }
    return merges;
}
