#include "simulation_inspector.hpp"

#include <stack>

namespace sim
{
    const float r_max_gr = 13.0f;
    const int32_t nbins_gr = 200;
    const float dr_gr = r_max_gr / nbins_gr;

    simulation_inspector::simulation_inspector()
    {
    
    }

    void simulation_inspector::clear()
    {
        m_rdfData.clear();
        m_activeRDFPairs.clear();
        m_totalFrames = 0;
    }

    void histogram::push(const std::vector<float>& frame_hist)
    {
        if (capacity == 0) capacity = 50;

        if (data.size() != capacity * nbins_gr) 
        {
            data.assign(capacity * nbins_gr, 0.0f);
            write_idx = 0;
            total_frames = 0;
        }

        std::copy(frame_hist.begin(), frame_hist.begin() + std::min(frame_hist.size(), (size_t)nbins_gr),
                  data.begin() + write_idx * nbins_gr);

        write_idx = (write_idx + 1) % capacity;
        total_frames++;
    }

    std::vector<float> histogram::get_average() const
    {
        if (total_frames == 0) return std::vector<float>(nbins_gr, 0.0f);

        std::vector<float> avg(nbins_gr, 0.0f);
        uint32_t frames_in_window = std::min(static_cast<uint32_t>(total_frames), capacity);

        for (uint32_t i = 0; i < capacity; ++i) 
        {
            const float* frame = data.data() + i * nbins_gr;
            for (int b = 0; b < nbins_gr; ++b)
                avg[b] += frame[b];
        }

        for (auto& v : avg) v /= frames_in_window;
        return avg;
    }

    void simulation_inspector::startRDF(uint32_t Zi, uint32_t Zj)
    {
        uint64_t key = makeKey(Zi, Zj);
        if (m_rdfData.find(key) != m_rdfData.end()) return;

        RDFPair rp;
        rp.Zi = Zi;
        rp.Zj = Zj;
        rp.m_histogram.capacity = m_analysisWindow;
        m_rdfData[key] = std::move(rp);

        auto pair = std::make_pair(Zi, Zj);
        if (std::find(m_activeRDFPairs.begin(), m_activeRDFPairs.end(), pair) == m_activeRDFPairs.end())
            m_activeRDFPairs.push_back(pair);
    }

    void simulation_inspector::update(fun::universe& u)
    {
        m_totalFrames++;
        m_totalFrames %= (m_analysisWindow + 1);

        for (const auto& p : m_activeRDFPairs)
        {
            uint64_t key = makeKey(p.first, p.second);
            auto it = m_rdfData.find(key);
            if (it == m_rdfData.end()) continue;

            RDFPair& rdfp = it->second;
            std::vector<float> frame = computeFrameRDF(u, p.first, p.second);
            rdfp.m_histogram.push(frame);

            if (rdfp.Ni == 0 || rdfp.Nj == 0) 
            {
                const auto& atoms = u.getAtomData().atoms;
                rdfp.Ni = std::count_if(atoms.begin(), atoms.end(), [&](const auto& a){ return a.ZIndex == p.first; });
                rdfp.Nj = std::count_if(atoms.begin(), atoms.end(), [&](const auto& a){ return a.ZIndex == p.second; });
            }
        }
    }

    std::vector<float> simulation_inspector::computeFrameRDF(fun::universe& u, uint32_t Zi, uint32_t Zj) const
    {
        std::vector<float> hist(nbins_gr, 0.0f);
        const auto& atoms = u.getAtomData().atoms;
        const auto& positions = u.getData().positions;

        std::vector<int> idx_i, idx_j;
        for (size_t i = 0; i < atoms.size(); ++i) {
            if (atoms[i].ZIndex == Zi) idx_i.push_back(i);
            if (atoms[i].ZIndex == Zj) idx_j.push_back(i);
        }

        for (int i : idx_i) 
        {
            for (int j : idx_j) 
            {
                if (Zi == Zj && i == j) continue;
                float dist = glm::length(u.minImageVec(positions[i] - positions[j]));
                if (dist < r_max_gr && dist > 0.001f) 
                {
                    int bin = static_cast<int>(dist / dr_gr);
                    if (bin >= 0 && bin < nbins_gr)
                        hist[bin] += 1.0f;
                }
            }
        }
        return hist;
    }

    std::vector<std::string> simulation_inspector::getActiveRDFs() const
    {
        std::vector<std::string> pairs;
        pairs.reserve(m_activeRDFPairs.size());

        for (const auto& p : m_activeRDFPairs)
        {
            uint32_t Zi = p.first;
            uint32_t Zj = p.second;

            std::string symbolA = constants::getAtomLetter(Zi);
            std::string symbolB = constants::getAtomLetter(Zj);

            std::string label = symbolA + "-" + symbolB;

            if (symbolA > symbolB)
                label = symbolB + "-" + symbolA;

            pairs.push_back(label);
        }

        return pairs;
    }

    std::vector<std::vector<point>> simulation_inspector::getRDFs(fun::universe& u) const
    {
        std::vector<std::vector<point>> result{};

        for (const auto& p : m_activeRDFPairs)
            result.push_back(getRDF(p.first, p.second, u));
            
        return result;
    }

    std::vector<point> simulation_inspector::getRDF(uint32_t Zi, uint32_t Zj, fun::universe& u) const
    {
        uint64_t key = makeKey(Zi, Zj);
        auto it = m_rdfData.find(key);
        if (it == m_rdfData.end()) return {};

        const RDFPair& rp = it->second;
        std::vector<float> avg = rp.m_histogram.get_average();

        glm::vec3 box = u.boxSizes();
        float volume = box.x * box.y * box.z;
        float rho_j = rp.Nj > 0 ? static_cast<float>(rp.Nj) / volume : 0.0f;

        std::vector<point> rdf;
        rdf.reserve(nbins_gr);

        for (int bin = 0; bin < nbins_gr; ++bin) 
        {
            float r = (bin + 0.5f) * dr_gr;
            float shell_vol = 4.0f * M_PI * r * r * dr_gr;
            float count = avg[bin];
            if (Zi == Zj) count *= 0.5f;

            float g = (rp.Ni > 0 && shell_vol > 0.0f) ? count / (rp.Ni * shell_vol * rho_j) : 0.0f;
            rdf.push_back({r, g});
        }

        float sum_tail = 0.0f;
        int32_t tail_bins = 0;
        
        for (int b = nbins_gr * 0.7; b < nbins_gr; ++b)
        {
            sum_tail += rdf[b].y;
            tail_bins++;
        }

        if (tail_bins > 0)
        {
            float avg_tail = sum_tail / tail_bins;
            if (avg_tail > 0.01f)
            {
                for (auto& p : rdf)
                    p.y /= avg_tail;
            }
        }

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

    std::vector<std::pair<std::string, uint32_t>> simulation_inspector::getMoleculesPresent(fun::universe& u) const
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

    float simulation_inspector::calculateKineticEnergy(fun::universe& u) const
    {
        return u.calculateKineticEnergy();
    }

    float simulation_inspector::calculatePotentialEnergy(fun::universe& u) const
    {
        float totalPotential = 0.f;

        return totalPotential;
    }

    float simulation_inspector::meanParticleSpeed(fun::universe& u) const
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
