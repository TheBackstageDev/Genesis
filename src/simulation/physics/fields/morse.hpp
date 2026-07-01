#pragma once

#include "common.hpp"

namespace sim
{
    // De = bond dissociation energy
    // re = equilibrium bond length
    // a  = width parameter (controls steepness)
    struct MorseParams
    {
        float De;
        float re;
        float a;
        char order;
    };

    glm::vec3 computeMorseForce(uint32_t i, uint32_t j,
                                const float* __restrict x,
                                const float* __restrict y,
                                const float* __restrict z,
                                const float* __restrict De,
                                const float* __restrict re,
                                const float* __restrict a);

    float computeMorsePotential(uint32_t i, uint32_t j,
                            const float* __restrict x,
                            const float* __restrict y,
                            const float* __restrict z,
                            const float* __restrict De,
                            const float* __restrict re,
                            const float* __restrict a);
} // namespace sim
