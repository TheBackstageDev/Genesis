#include "tersoff.hpp"

#include "stdexcept"
#include <glm/gtc/constants.hpp>

namespace sim
{
    inline float cutoff(float r, float R, float D) {
        if (r < R - D) return 1.0f;
        else if (r > R + D) return 0.0f;
        else return 0.5f + 0.5f * std::cos(glm::pi<float>() * (r - (R - D)) / (2.0f * D));
    }

    glm::vec3 computeTersoffForce(
        uint32_t i, uint32_t j,
        const float* __restrict x,
        const float* __restrict y,
        const float* __restrict z,
        const float* __restrict params,
        const std::vector<uint32_t>& neighbors_i)
    {
        float A       = params[0];
        float B       = params[1];
        float lambda1 = params[2];
        float lambda2 = params[3];
        float lambda3 = params[4];
        float R       = params[5];
        float D       = params[6];
        float beta    = params[7];
        float n       = params[8];
        float c       = params[9];
        float d       = params[10];
        float h       = params[11];
        float gamma   = params[12];
        float m   = params[13];

        float dx = x[j] - x[i];
        float dy = y[j] - y[i];
        float dz = z[j] - z[i];
        float r2 = dx*dx + dy*dy + dz*dz;
        float r  = std::sqrt(r2);

        float fc = cutoff(r, R, D);

        float VR = A * std::exp(-lambda1 * r);
        float VA = -B * std::exp(-lambda2 * r);

        float zeta = 0.0f;
        for (uint32_t k : neighbors_i) 
        {
            if (k == j) continue;

            float dx2 = x[k] - x[i];
            float dy2 = y[k] - y[i];
            float dz2 = z[k] - z[i];
            float r2k = dx2*dx2 + dy2*dy2 + dz2*dz2;
            float rk  = std::sqrt(r2k);

            float fcik = cutoff(rk, R, D);
            if (fcik == 0.0f) continue;

            float cosTheta = (dx*dx2 + dy*dy2 + dz*dz2) / (r * rk);

            float gTheta = 1.0f + (c*c)/(d*d)
                - (c*c)/(d*d + (h - cosTheta)*(h - cosTheta));

            float expTerm = std::exp(lambda3 * (r - rk)*(r - rk)*(r - rk));

            zeta += fcik * gTheta * expTerm;
        }

        float bij = std::pow(1.0f + std::pow(beta * zeta, n), -1.0f/(2.0f*m));

        float dVRdr = -lambda1 * VR;
        float dVAdr = -lambda2 * VA;
        float dEdr  = fc * (dVRdr + bij * dVAdr);

        float F = -dEdr / r;
        return glm::vec3(F * dx, F * dy, F * dz);
    }
} // namespace sim
