#include "spacewind/rng.h"
#include "spacewind/constants.h"

#include <math.h>

static uint32_t sw_rotr32(uint32_t x, uint32_t r) {
    return (x >> r) | (x << ((0U - r) & 31U));
}

void sw_rng_seed(sw_rng *rng, uint64_t seed, uint64_t stream) {
    rng->state = 0U;
    rng->inc = (stream << 1U) | 1U;
    rng->has_spare = 0;
    rng->spare = 0.0;
    (void)sw_rng_u32(rng);
    rng->state += seed;
    (void)sw_rng_u32(rng);
}

uint32_t sw_rng_u32(sw_rng *rng) {
    const uint64_t oldstate = rng->state;
    const uint32_t xorshifted = (uint32_t)(((oldstate >> 18U) ^ oldstate) >> 27U);
    const uint32_t rot = (uint32_t)(oldstate >> 59U);
    rng->state = oldstate * 6364136223846793005ULL + rng->inc;
    return sw_rotr32(xorshifted, rot);
}

double sw_rng_uniform(sw_rng *rng) {
    return ((double)sw_rng_u32(rng) + 0.5) / 4294967296.0;
}

double sw_rng_normal(sw_rng *rng) {
    if (rng->has_spare != 0) {
        rng->has_spare = 0;
        return rng->spare;
    }
    {
        const double u1 = sw_rng_uniform(rng);
        const double u2 = sw_rng_uniform(rng);
        const double radius = sqrt(-2.0 * log(u1));
        const double theta = SW_TWO_PI * u2;
        rng->spare = radius * sin(theta);
        rng->has_spare = 1;
        return radius * cos(theta);
    }
}
