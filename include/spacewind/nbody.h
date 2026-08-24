#ifndef SPACEWIND_NBODY_H
#define SPACEWIND_NBODY_H

#include <stdbool.h>
#include <stddef.h>

#include "spacewind/math3.h"

#define SW_NBODY_MAX 128U

typedef struct {
    double mass_kg;
    double radius_m;
    sw_vec3 position_m;
    sw_vec3 velocity_m_s;
    sw_vec3 acceleration_m_s2;
    unsigned id;
} sw_nbody_object;

typedef struct {
    size_t count;
    double softening_m;
    sw_nbody_object objects[SW_NBODY_MAX];
} sw_nbody_system;

void sw_nbody_accelerations(sw_nbody_system *system);
bool sw_nbody_leapfrog_step(sw_nbody_system *system, double step_s);
double sw_nbody_total_energy_j(const sw_nbody_system *system);
sw_vec3 sw_nbody_total_momentum_kg_m_s(const sw_nbody_system *system);
sw_vec3 sw_nbody_center_of_mass_m(const sw_nbody_system *system);
size_t sw_nbody_merge_overlaps(sw_nbody_system *system);

#endif
