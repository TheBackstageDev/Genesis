#include "barostats.hpp"

namespace sim
{
    float computePressure(fun::universe& u, float virial)
    {
        auto& data = u.getData();
        auto& atomData = u.getAtomData();
        
        if (atomData.atoms.empty())
            return 0.0f;

        float kinetic = u.calculateKineticEnergy();
        glm::vec3 box = u.boxSizes();
        float volume = box.x * box.y * box.z;
        
        if (volume <= 0.0f) return 0.0f;

        size_t N = atomData.atoms.size();
        float ndof = 3.0f * N - 3.0f; 
        if (ndof < 1.0f) ndof = 1.0f;

        float temperature = (2.0f * kinetic) / (ndof * KB);

        float P_ideal = (N * KB * temperature) / volume;
        float P_virial = -virial / (3.0f * volume);

        const float conv_factor = 1660.53906660f;
        float pressure_kPa = (P_ideal + P_virial) * conv_factor;

        return -pressure_kPa;
    }
} // namespace sim