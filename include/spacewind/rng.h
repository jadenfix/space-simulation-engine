#ifndef SPACEWIND_RNG_H
#define SPACEWIND_RNG_H

#include <stdint.h>

typedef struct {
    uint64_t state;
    uint64_t inc;
    int has_spare;
    double spare;
} sw_rng;

void sw_rng_seed(sw_rng *rng, uint64_t seed, uint64_t stream);
uint32_t sw_rng_u32(sw_rng *rng);
double sw_rng_uniform(sw_rng *rng);
double sw_rng_normal(sw_rng *rng);

#endif
