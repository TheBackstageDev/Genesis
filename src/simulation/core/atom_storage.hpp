#pragma once

#include <vector>
#include <glm/glm.hpp>

namespace sim
{
    struct atomView
    {
        size_t i = 0;
        float px, py, pz;
        float vx, vy, vz;
        float fx, fy, fz;
        float invMass;
        float q;
        float sigma;
        float epsilon;

        // Morse parameters
        float De = 0.f;
        float re = 0.f;
        float a = 0.f;
    };
    
    class atomStorage 
    {
    public:
        atomStorage(size_t n);
        atomStorage() = default;

        size_t mobileCount() const { return N; }

        void clear();
        void clearForces();
        void resize(size_t n);
        void setAtom(const atomView& aView);

        float* xData() { return x.data(); }
        float* yData() { return y.data(); }
        float* zData() { return z.data(); }

        float* vxData() { return vx.data(); }
        float* vyData() { return vy.data(); }
        float* vzData() { return vz.data(); }

        float* fxData() { return fx.data(); }
        float* fyData() { return fy.data(); }
        float* fzData() { return fz.data(); }

        float* qData() { return q.data(); }
        float* ljParamsData() { return ljParams.data(); }
        float* invMassData() { return invMass.data(); }

        float* DeData() { return De.data(); }
        float* reData() { return re.data(); }
        float* aData()  { return a.data(); }
        
        glm::vec3 position(size_t i) { return glm::vec3(x[i], y[i], z[i]); }
        glm::vec3 force(size_t i) { return glm::vec3(fx[i], fy[i], fz[i]); }
        glm::vec3 velocity(size_t i) { return glm::vec3(vx[i], vy[i], vz[i]); }
        float& charge(size_t i) { return q[i]; }
    private:
        size_t N = 0;

        std::vector<float> x, y, z;
        std::vector<float> vx, vy, vz;
        std::vector<float> fx, fy, fz;
        std::vector<float> invMass;
        std::vector<float> q, ljParams;

        std::vector<float> De;
        std::vector<float> re;
        std::vector<float> a;
    };
} // namespace sim
