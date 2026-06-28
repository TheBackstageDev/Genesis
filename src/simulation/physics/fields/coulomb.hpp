#pragma once

#include "common.hpp"

namespace sim
{
    constexpr float COULOMB_CUTOFF = 12.f;
    constexpr float COULOMB_K = 1389.354576f; // kJ·mol⁻¹· Å ·e⁻²

    glm::vec3 computeCoulombForce(
        uint32_t i, uint32_t j,
        const float* __restrict x,
        const float* __restrict y,
        const float* __restrict z,
        const float* __restrict q);
} // namespace sim
