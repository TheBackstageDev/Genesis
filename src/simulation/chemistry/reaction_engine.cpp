#include "reaction_engine.hpp"
#include <iostream>
#include <numeric>

namespace sim
{
    reaction_engine::reaction_engine()
    {
        std::filesystem::path folder = "resource";
    }

    void reaction_engine::update(fun::universe &u, core::verlet_list &verlet_list)
    {
        auto &atomData = u.getAtomData();
        auto &data = u.getData();

        float* __restrict x = data.xData();
        float* __restrict y = data.yData();
        float* __restrict z = data.zData();

        float* __restrict fx = data.fxData();
        float* __restrict fy = data.fyData();
        float* __restrict fz = data.fzData();

        float* __restrict De = data.DeData();
        float* __restrict re = data.reData();
        float* __restrict a  = data.aData();

        for (int32_t i = 0; i < atomData.atoms.size(); ++i)
        {
            for (auto& j : verlet_list.verlet[i])
            {
                glm::vec3 force = computeMorseForce(i, j, x, y, z, De, re, a);

                fx[i] += force.x;
                fy[i] += force.y;
                fz[i] += force.z;

                fx[j] -= force.x;
                fy[j] -= force.y;
                fz[j] -= force.z;
            }
        }
    }
} // namespace sim
