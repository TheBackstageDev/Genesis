#pragma once

#include "common.hpp"

namespace sim
{
    struct LJParams
    {
        float sigma;
        float epsilon;
    };

    glm::vec3 computeLJforce(uint32_t i, uint32_t j,
                             const float* __restrict x,
                             const float* __restrict y,
                             const float* __restrict z,
                             const float* __restrict ljParams);

    float computeLJpotential(uint32_t i, uint32_t j,
                             const float* __restrict x,
                             const float* __restrict y,
                             const float* __restrict z,
                             const float* __restrict ljParams);
} // namespace sim
