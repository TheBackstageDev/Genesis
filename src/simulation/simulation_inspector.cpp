#include "simulation_inspector.hpp"

#include <stack>

namespace sim
{
    simulation_inspector::simulation_inspector()
    {
    
    }

    std::vector<point> simulation_inspector::getRDFgraph(fun::universe& u, uint32_t Zi, uint32_t Zj)
    {
        return {};
    }

    std::string buildFormula(const std::unordered_map<uint8_t, uint32_t>& counts)
    {
        std::vector<std::pair<std::string, uint32_t>> elems;
        for (auto& kv : counts)
            elems.emplace_back(constants::getAtomLetter(kv.first), kv.second);

        std::sort(elems.begin(), elems.end(), [&](auto& a, auto& b) {
            if (a.first == "C") return true;
            if (b.first == "C") return false;
            if (a.first == "H") return true;
            if (b.first == "H") return false;
            if (a.first == "S") return true;
            if (b.first == "S") return false;
            return a.first < b.first;
        });

        std::ostringstream formula;
        for (auto& kv : elems)
        {
            formula << kv.first;
            if (kv.second > 1) formula << kv.second;
        }
        return formula.str();
    }

    std::vector<std::pair<std::string, uint32_t>> simulation_inspector::getMoleculesPresent(fun::universe& u)
    {
        const auto& atomData = u.getAtomData();
        const auto& bonds = atomData.bonds;

        std::unordered_map<int32_t, std::vector<int32_t>> bondTree{};
        for (int32_t i = 0; i < bonds.size(); ++i)
        {
            const auto& bond = bonds[i];
            bondTree[bond.centralAtom].emplace_back(bond.bondedAtom);
            bondTree[bond.bondedAtom].emplace_back(bond.centralAtom);
        }

        std::vector<bool> visited(atomData.atoms.size(), false);
        std::vector<std::pair<std::string, uint32_t>> molecules;

        for (int32_t i = 0; i < atomData.atoms.size(); ++i)
        {
            if (!visited[i])
            {
                std::vector<int32_t> molAtoms;
                std::stack<int32_t> stack;
                stack.push(i);

                while (!stack.empty())
                {
                    int32_t current = stack.top();
                    stack.pop();

                    if (visited[current]) continue;
                    visited[current] = true;
                    molAtoms.push_back(current);

                    for (int32_t neighbor : bondTree[current])
                    {
                        if (!visited[neighbor]) stack.push(neighbor);
                    }
                }

                std::unordered_map<uint8_t, uint32_t> counts;
                for (int32_t idx : molAtoms)
                    counts[atomData.atoms[idx].ZIndex]++;

                const std::string formula = buildFormula(counts);

                auto it = std::find_if(molecules.begin(), molecules.end(),
                    [&](const std::pair<std::string, uint32_t>& pair)
                    {
                        return pair.first == formula;
                    });

                if (it != molecules.end())
                    it->second++;
                else
                    molecules.emplace_back(formula, 1);
            }
        }

        return molecules;
    }

    float simulation_inspector::computeDensity(fun::universe& u, uint32_t Z, uint32_t frame)
    {
        return 0.f;
    }
} // namespace sim
