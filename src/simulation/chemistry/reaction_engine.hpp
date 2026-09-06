#pragma once

#include "simulation/fundamental_structures.hpp"
#include "simulation/physics/fields/parameter_table.hpp"
#include "simulation/physics/fields/morse.hpp"
#include "simulation/core/universe.hpp"

#include "core/utils/verletlist.hpp"
#include "core/utils/spatialgrid.hpp"

#include <vector>
#include <functional>

namespace sim
{
    class reaction_engine
    {
    public:
        reaction_engine();
        ~reaction_engine() = default;

        void update(fun::universe& u, core::verlet_list& verlet_list);
    private:
        std::vector<uint32_t> bondOrders;
        bool isRadical(fun::universe& u, int32_t atom)
        {
            return bondOrders[atom] < constants::getUsualBonds(u.getAtomData().atoms[atom].ZIndex);
        }

        bool inReactiveRange(fun::universe& u, int32_t a, int32_t b)
        {
            auto& data = u.getData();
            auto& atomData = u.getAtomData();

            float* __restrict x = data.xData();
            float* __restrict y = data.yData();
            float* __restrict z = data.zData();

            float r = glm::length(glm::vec3(x[a] - x[b], y[a] - y[b], z[a] - z[b]));
            float range = (constants::VDW_RADII[atomData.atoms[a].ZIndex] + constants::VDW_RADII[atomData.atoms[b].ZIndex]);
            
            return r < range;
        }

        float computeUniverseEnergy(fun::universe& u);
        void findLowestEnergyStateAtom(fun::universe& u, int32_t atom, std::vector<int32_t> atomNeighbours);
        void propagateValence(fun::universe& u, int32_t atom, std::vector<int32_t> visited);
    
        void createBond(fun::universe& u, int32_t a, int32_t b);
        void breakBond(fun::universe& u, int32_t a, int32_t b);
    };
} // namespace sim
