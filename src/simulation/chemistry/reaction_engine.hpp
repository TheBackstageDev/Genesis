#pragma once

#include "simulation/fundamental_structures.hpp"
#include "simulation/physics/fields/parameter_table.hpp"
#include "simulation/physics/fields/morse.hpp"
#include "simulation/core/universe.hpp"

#include "core/utils/verletlist.hpp"
#include "core/utils/spatialgrid.hpp"

#include <vector>
#include <unordered_set>
#include <functional>

namespace sim
{
    using atomId_t = uint32_t;

    struct bondEdge
    {
        atomId_t a, b;
        int32_t order;
    };

    struct functionalGroup
    {
        struct PatternLevel
        {
            uint32_t Z;                     // expected ZIndex
            int32_t count;                  // at this level, how many should be
        };

        std::vector<std::vector<PatternLevel>> pattern;
    };

    struct moleculeGraph
    {
        std::vector<atomId_t> nodes;
        std::vector<uint32_t> depthOffset;
        std::vector<uint16_t> groups;

        std::unordered_set<atomId_t> bonds;
    };

    struct reactionRule
    {
        std::vector<atomId_t> requiredGroups;
        void execute(atomId_t molA, atomId_t molB);
    };

    class reaction_engine
    {
    public:
        reaction_engine();
        ~reaction_engine() = default;

        void update(fun::universe& u, core::verlet_list& verlet_list);
        void updateGraph() { m_updateGraph = true; }
        void resetSortedFlag() { m_sortedIndices = true; }
    private:
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

        std::vector<moleculeGraph> m_moleculeGraphs{};
        std::vector<functionalGroup> m_functionalGroups{};      // Functional groups which are 'molds' that are used to find on the graph the groups
        std::vector<functionalGroup> m_functionalGroupsGraph{}; // Functional groups that are present on the molecules graph
        std::vector<bondEdge> m_bondGraph{};

        bool m_updateGraph = true;
        bool m_sortedIndices = true;

        float computeUniverseEnergy(fun::universe& u);

        void createMoleculeGraph(fun::universe& u);
        void computeMoleculeGraph(const fun::universe& u, const atomId_t current, moleculeGraph& graph);
        void findFunctionalGroups(fun::universe& u, moleculeGraph& graph);

        void createBond(fun::universe& u, int32_t a, int32_t b);
        void breakBond(fun::universe& u, int32_t a, int32_t b);

        void setupFunctionalGroups();
        void setupReactionRules();
    };
} // namespace sim
