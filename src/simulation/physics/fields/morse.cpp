#include "morse.hpp"
#include <cmath>

namespace sim
{
    glm::vec3 sim::computeMorseForce(uint32_t i, uint32_t j,
                                const float* __restrict x,
                                const float* __restrict y,
                                const float* __restrict z,
                                const float* __restrict De,
                                const float* __restrict re,
                                const float* __restrict a)
    {
        float dx = x[i] - x[j];
        float dy = y[i] - y[j];
        float dz = z[i] - z[j];
        float r2 = dx*dx + dy*dy + dz*dz;
        float r  = std::sqrt(r2);

        if (r < 1e-8f) return glm::vec3(0.0f);

        float De_ = De[i];
        float re_ = re[i];
        float a_  = a[i];

        float expTerm = std::exp(-a_ * (r - re_));
        float Fmag = 2.0f * a_ * De_ * (1.0f - expTerm) * expTerm;

        float invr = 1.0f / r;
        return glm::vec3(dx * invr * Fmag,
                         dy * invr * Fmag,
                         dz * invr * Fmag);
    }

    float sim::computeMorsePotential(uint32_t i, uint32_t j,
                                const float* __restrict x,
                                const float* __restrict y,
                                const float* __restrict z,
                                const float* __restrict De,
                                const float* __restrict re,
                                const float* __restrict a)
    {
        float dx = x[i] - x[j];
        float dy = y[i] - y[j];
        float dz = z[i] - z[j];
        float r  = std::sqrt(dx*dx + dy*dy + dz*dz);

        float De_ = De[i];
        float re_ = re[i];
        float a_  = a[i];

        float term = 1.0f - std::exp(-a_ * (r - re_));
        return De_ * term * term;
    }
}
