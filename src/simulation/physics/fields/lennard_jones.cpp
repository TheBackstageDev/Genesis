#include "lennard_jones.hpp"

namespace sim
{
    glm::vec3 sim::computeLJforce(
        uint32_t i, uint32_t j,
        const float& dx,
        const float& dy,
        const float& dz,
        const float* __restrict ljParams)
    {
        float dr2 = dx*dx + dy*dy + dz*dz;

        float sigma_i   = ljParams[i*2];
        float epsilon_i = ljParams[i*2 + 1];
        float sigma_j   = ljParams[j*2];
        float epsilon_j = ljParams[j*2 + 1];

        float sigma = 0.5f * (sigma_i + sigma_j);
        float epsilon = sqrtf(epsilon_i * epsilon_j);

        if (dr2 < std::numeric_limits<float>::epsilon())
            return glm::vec3{0.f};

        float cutoff = 2.5f * sigma;
        if (dr2 > cutoff * cutoff)
            return glm::vec3{0.f};

        float inv_r2 = 1.f / dr2;
        float inv_r6 = inv_r2 * inv_r2 * inv_r2;
        float inv_r12 = inv_r6 * inv_r6;

        float sig2 = sigma * sigma;
        float sig6 = sig2 * sig2 * sig2;
        float sig12 = sig6 * sig6;

        float force_mag = 24.f * epsilon * inv_r2 * (2.f * sig12 * inv_r12 - sig6 * inv_r6);
        
        return glm::vec3{-force_mag * dx, -force_mag * dy, -force_mag * dz};
    }

    float sim::computeLJpotential(uint32_t i, uint32_t j,
                             const float* __restrict x,
                             const float* __restrict y,
                             const float* __restrict z,
                             const float* __restrict ljParams)
    {
        float dx = x[i] - x[j];
        float dy = y[i] - y[j];
        float dz = z[i] - z[j];
        float r2 = dx*dx + dy*dy + dz*dz;
        float r  = std::sqrt(r2);
        
        float sigma_i   = ljParams[i*2];
        float epsilon_i = ljParams[i*2 + 1];
        float sigma_j   = ljParams[j*2];
        float epsilon_j = ljParams[j*2 + 1];
        
        float sigma = 0.5f * (sigma_i + sigma_j);
        float epsilon = sqrtf(epsilon_i * epsilon_j);
        
        if (r > sigma * 2.5f) return 0.0f;

        float sr6 = std::pow(sigma / r, 6);
        float sr12 = sr6 * sr6;

        return 4.0f * epsilon * (sr12 - sr6);
    }
} // namespace sim
