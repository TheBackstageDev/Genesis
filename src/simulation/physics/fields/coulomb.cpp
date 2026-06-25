#include "coulomb.hpp"

namespace sim
{
    glm::vec3 sim::computeCoulombForce(
        uint32_t i, uint32_t j,
        const float* __restrict x,
        const float* __restrict y,
        const float* __restrict z,
        const float* __restrict q)
    {
        float dx = x[j] - x[i];
        float dy = y[j] - y[i];
        float dz = z[j] - z[i];

        float dr2 = dx*dx + dy*dy + dz*dz;

        if (dr2 < std::numeric_limits<float>::epsilon() || dr2 > COULOMB_CUTOFF * COULOMB_CUTOFF)
            return glm::vec3{0.f};

        float qq = q[i] * q[j];
        if (qq == 0.f)
            return glm::vec3{0.f};

        float inv_r = 1.0f / sqrtf(dr2);
        float forceMag = COULOMB_K * qq * inv_r * inv_r * inv_r;

        return glm::vec3{-forceMag * dx, -forceMag * dy, -forceMag * dz};
    }
} // namespace sim
