#pragma once

#include "common.hpp"

namespace sim
{
    struct TersoffParams
    {
        // Pairwise terms
        float A, B;             // repulsive and attractive prefactors
        float lambda1, lambda2; // exponential decay constants
        float lambda3;          // angular decay constant

        // Cutoff radii
        float R, D;             // cutoff radius and smoothing distance

        // Bond order / angular dependence
        float beta;             // bond order scaling
        float n;                // bond order exponent
        float c, d;             // angular parameters
        float h;                // preferred cos(theta0)
        float gamma;            // angular weight

        float m;
    };

    constexpr size_t PARAM_SIZE = 13;

    glm::vec3 computeTersoffForce(
        uint32_t i, uint32_t j,
        const float* __restrict x,
        const float* __restrict y,
        const float* __restrict z,
        const float* __restrict params,
        const std::vector<uint32_t>& neighbors_i);
} // namespace sim
