#include "reaction_engine.hpp"
#include <iostream>
#include <numeric>
#include <queue>

namespace sim
{
    reaction_engine::reaction_engine()
    {
        std::filesystem::path folder = "resource";

        setupFunctionalGroups();
        setupReactionRules();
    }

    void getAtomCount(const std::vector<uint32_t> ZIndices, std::unordered_map<uint32_t, int32_t>& counts)
    {
        for (uint32_t i : ZIndices)
            counts[i]++;
    }

    void reaction_engine::setupFunctionalGroups()
    {
        
    }
    
    void reaction_engine::setupReactionRules()
    {
    
    }

    void reaction_engine::update(fun::universe &u, core::verlet_list &verlet_list)
    {
        auto &atomData = u.getAtomData();
        auto &data = u.getData();

        if (m_updateGraph)
        {
            createMoleculeGraph(u);
            m_updateGraph = false;
            m_sortedIndices = false; // cannot assume anymore the molecules will contain linear (1, 2, 3, 4) indices
        }

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
                float dx = x[i] - x[j];
                float dy = y[i] - y[j];
                float dz = z[i] - z[j];

                glm::vec3 force = computeMorseForce(dx, dy, dz, De[i], re[i], a[i]);

                fx[i] += force.x;
                fy[i] += force.y;
                fz[i] += force.z;

                fx[j] -= force.x;
                fy[j] -= force.y;
                fz[j] -= force.z;
            }
        }
    }

    void reaction_engine::computeMoleculeGraph(const fun::universe& u, const atomId_t current, moleculeGraph& graph)
    {
        graph.nodes.clear();
        graph.depthOffset.clear();
        graph.bonds.clear();
        graph.groups.clear();

        const size_t atomCount = u.numAtoms();

        std::vector<bool> visited(atomCount, false);
        visited[current] = true;
        
        std::queue<std::pair<atomId_t, int32_t>> q;
        q.push(std::make_pair(current, 0));
        
        int32_t currentDepth = -1;

        while (!q.empty())
        {
            auto [current, depth] = q.front();
            q.pop();

            if (depth != currentDepth)
            {
                graph.depthOffset.push_back(static_cast<uint32_t>(graph.nodes.size()));
                currentDepth = depth;
            }

            for (size_t i = 0; i < m_bondGraph.size(); ++i)
            {
                const auto& bond = m_bondGraph[i];
                atomId_t neighbor = UINT32_MAX;

                if (bond.a == current)
                    neighbor = bond.b;
                else if (bond.b == current)
                    neighbor = bond.a;
                else
                    continue;

                graph.bonds.emplace(static_cast<uint32_t>(i));

                if (!visited[neighbor])
                {
                    visited[neighbor] = true;
                    q.push({neighbor, depth + 1});
                }
            }
        }
    }

    bool functionalGroupPossible(const std::unordered_map<uint32_t, int32_t>& atomCounts, const std::unordered_map<uint32_t, int32_t>& moldCounts)
    {
        for (auto& [ZIndex, count] : atomCounts)
            if (!moldCounts.contains(ZIndex) || moldCounts.at(ZIndex) < count) return false;

        return true;
    }

    void reaction_engine::findFunctionalGroups(fun::universe& u, moleculeGraph& graph)
    {
        std::unordered_map<uint32_t, int32_t> atomCounts{};

        auto& data = u.getAtomData();

        for (int32_t i = 0; i < graph.nodes.size(); ++i)
            atomCounts[data.atoms[i].ZIndex]++;

        for (int32_t m = 0; m < m_functionalGroups.size(); ++m)
        {
            const auto& mold = m_functionalGroups[m];
            
            if (graph.depthOffset.size() < mold.pattern.size() + 1) continue;

            for (size_t depth = 0; depth < mold.pattern.size(); ++depth)
            {
                const uint32_t start = graph.depthOffset[depth];
                const uint32_t end   = graph.depthOffset[depth + 1];
            }
        }
    }

    void reaction_engine::createMoleculeGraph(fun::universe& u)
    {
        auto &atomData = u.getAtomData();
        
        for (int32_t i = 0; i < atomData.bonds.size(); ++i)
        {
            bondEdge edge{};
            edge.a = atomData.bonds[i].centralAtom;
            edge.b = atomData.bonds[i].bondedAtom;
            edge.order = static_cast<int32_t>(atomData.bonds[i].type);

            m_bondGraph.emplace_back(std::move(edge));
        }

        std::vector<atomId_t> excludeList{}; // atoms alreadly known to be in molecules

        for (int32_t i = 0; i < atomData.atoms.size(); ++i)
        {
            if (std::find(excludeList.begin(), excludeList.end(), i) != excludeList.end()) continue;
            
            moleculeGraph m{};

            computeMoleculeGraph(u, i, m);
            findFunctionalGroups(u, m);

            if (m_sortedIndices)
                i += static_cast<int32_t>(m.nodes.size());
            else
                excludeList.insert(excludeList.end(), m.nodes.begin(), m.nodes.end());

            m_moleculeGraphs.emplace_back(std::move(m));
        }
    }
} // namespace sim
