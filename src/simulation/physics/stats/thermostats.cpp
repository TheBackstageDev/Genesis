#include "thermostats.hpp"

namespace sim
{
    // Bussi–Donadio–Parrinello
    float computeTemperature(fun::universe &u, float kelvin)
    {
        auto& atomData = u.getAtomData();
        auto& data = u.getData();
        size_t N = atomData.atoms.size();
        if (N == 0) return 0.f;

        float ndof = 3.0f * static_cast<float>(N) - 3.0f;
        if (ndof <= 0.0f) ndof = 1.0f;

        float KE = u.calculateKineticEnergy();
        return (2.0f * KE) / (ndof * KB);
    }
} // namespace sim