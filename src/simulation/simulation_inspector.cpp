#include "simulation_inspector.hpp"

#include <stack>

namespace sim
{
    const float simulation_inspector::r_max_gr = 10.0f;
    const int32_t simulation_inspector::nbins_gr = 200;
    const float simulation_inspector::dr_gr = r_max_gr / nbins_gr;

    simulation_inspector::simulation_inspector()
    {
    
    }

    std::vector<point> simulation_inspector::getRDFgraph(fun::universe& u, uint32_t Zi, uint32_t Zj)
    {
        constexpr float r_max = 10.0f;
        constexpr int32_t nbins = 200;
        constexpr float dr = r_max / nbins;

        std::vector<int> histogram(nbins, 0);

        auto& data = u.getData();
        auto& atoms = u.getAtomData().atoms;
        std::vector<glm::vec3>& positions = data.positions;
        size_t n_atoms = positions.size();

        std::vector<int> indices_i, indices_j;
        indices_i.reserve(n_atoms);
        indices_j.reserve(n_atoms);

        for (int i = 0; i < n_atoms; ++i) {
            uint32_t z = atoms[i].ZIndex;
            if (z == Zi) indices_i.push_back(i);
            if (z == Zj) indices_j.push_back(i);
        }

        glm::vec3 boxSizes = u.boxSizes();
        float volume = boxSizes.x * boxSizes.y * boxSizes.z;

        for (int idx_i : indices_i) {
            for (int idx_j : indices_j) {
                if (Zi == Zj && idx_i == idx_j) continue;

                float dist = glm::length(u.minImageVec(positions[idx_i] - positions[idx_j]));
                if (dist < r_max && dist > 0.001f) {
                    int bin = static_cast<int>(dist / dr);
                    if (bin < nbins) histogram[bin]++;
                }
            }
        }

        std::vector<point> rdf;
        rdf.reserve(nbins);

        size_t Ni = indices_i.size();
        size_t Nj = indices_j.size();
        float rho_j = static_cast<float>(Nj) / volume;

        for (int bin = 0; bin < nbins; ++bin) 
        {
            float r = (bin + 0.5f) * dr;
            float shell_volume = 4.0f * M_PI * r * r * dr;

            float count = static_cast<float>(histogram[bin]);
            if (Zi == Zj) count *= 0.5f;

            // Normalized RDF
            float g = count / (Ni * shell_volume * rho_j);

            rdf.push_back({r, g});
        }

        return rdf;
    }

    void simulation_inspector::beginRDFaccumulation(uint32_t Zi, uint32_t Zj, uint32_t frames, fun::universe& u) 
    {
        m_Zi = Zi;
        m_Zj = Zj;
        m_targetFrames = frames;
        m_accumulatedFrames = 0;
        m_rdfHistogram.assign(nbins_gr, 0.0);
        m_accumulating = true;

        auto& atoms = u.getAtomData().atoms;
        m_Ni = std::count_if(atoms.begin(), atoms.end(), [&](auto& a){ return a.ZIndex == Zi; });
        m_Nj = std::count_if(atoms.begin(), atoms.end(), [&](auto& a){ return a.ZIndex == Zj; });

        m_rdfInProgress = true;
    }

    void simulation_inspector::accumulateRDFFrame(fun::universe& u) {
        if (!m_accumulating) return;

        std::vector<int> frame_hist(nbins_gr, 0);
        auto& data = u.getData();
        auto& atoms = u.getAtomData().atoms;
        auto& positions = data.positions;

        std::vector<int> indices_i, indices_j;
        for (int i = 0; i < atoms.size(); ++i) 
        {
            if (atoms[i].ZIndex == m_Zi) indices_i.push_back(i);
            if (atoms[i].ZIndex == m_Zj) indices_j.push_back(i);
        }

        for (int idx_i : indices_i) 
        {
            for (int idx_j : indices_j) 
            {
                if (m_Zi == m_Zj && idx_i == idx_j) continue;
                float dist = glm::length(u.minImageVec(positions[idx_i] - positions[idx_j]));
                if (dist < r_max_gr && dist > 0.001f) {
                    int bin = static_cast<int>(dist / dr_gr);
                    if (bin < nbins_gr) frame_hist[bin]++;
                }
            }
        }

        for (int b = 0; b < nbins_gr; ++b)
            m_rdfHistogram[b] += frame_hist[b];

        m_accumulatedFrames++;
        if (m_accumulatedFrames >= m_targetFrames)
            m_accumulating = false;
    }

    std::vector<point> simulation_inspector::finalizeRDF(fun::universe& u) 
    {
        std::vector<point> rdf;
        if (m_accumulatedFrames == 0) return rdf;

        for (float& h : m_rdfHistogram)
            h /= static_cast<float>(m_accumulatedFrames);

        glm::vec3 boxSizes = u.boxSizes();
        float volume = boxSizes.x * boxSizes.y * boxSizes.z;
        float rho_j = static_cast<float>(m_Nj) / volume;

        for (int bin = 0; bin < nbins_gr; ++bin) 
        {
            float r = (bin + 0.5f) * dr_gr;
            float shell_volume = 4.0f * M_PI * r * r * dr_gr;

            float count = m_rdfHistogram[bin];
            if (m_Zi == m_Zj) count *= 0.5f;

            float g = count / (m_Ni * shell_volume * rho_j);
            rdf.push_back({r, g});
        }

        m_rdfInProgress = false;

        return rdf;
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

    float simulation_inspector::calculateKineticEnergy(fun::universe& u)
    {
        return u.calculateKineticEnergy();
    }

    float simulation_inspector::calculatePotentialEnergy(fun::universe& u)
    {
        float totalPotential = 0.f;

        return totalPotential;
    }

    float simulation_inspector::meanParticleSpeed(fun::universe& u)
    {
        const auto& velocities = u.getData().velocities;
        float totalVelocity = 0.f;

        for (const glm::vec3& velocity : velocities)
            totalVelocity += glm::length(velocity);

        return totalVelocity / static_cast<float>(velocities.size()) * 100.f; // 100.f to convert from A/Ps to m/s
    }

    // Energies from Forces

    float simulation_inspector::ljPotential(fun::universe& u, int32_t i, int32_t j)
    {
        const auto& data = u.getData();
        return 0.f;
    }

    float simulation_inspector::coulombPotential(fun::universe& u, int32_t i, int32_t j)
    {
        const auto& data = u.getData();

        float r = glm::length(data.positions[i] - data.positions[j]);

        return COULOMB_K * data.q[i] * data.q[j] / r;
    }
} // namespace sim
