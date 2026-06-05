#pragma once

#include "universe.hpp"

#include <vector>

namespace sim
{
    struct point
    {
        float x, y;
    };

    class simulation_inspector
    {
    public:
        simulation_inspector();

        // between what atoms
        void beginRDFaccumulation(uint32_t Zi, uint32_t Zj, uint32_t frames, fun::universe& u);
        void accumulateRDFFrame(fun::universe& u);
        std::vector<point> finalizeRDF(fun::universe& u);
        std::vector<point> getRDFgraph(fun::universe& u, uint32_t Zi, uint32_t Zj);
        std::vector<std::pair<std::string, uint32_t>> getMoleculesPresent(fun::universe& u);
        
        float calculateKineticEnergy(fun::universe& u);
        float calculatePotentialEnergy(fun::universe& u);
        float meanParticleSpeed(fun::universe& u);

        bool accumulating() { return m_accumulating; } 
        bool rdfInProgress() { return m_rdfInProgress; }

        static const float r_max_gr;
        static const int32_t nbins_gr;
        static const float dr_gr;

    private:
        float ljPotential(fun::universe& u, int32_t i, int32_t j);
        float coulombPotential(fun::universe& u, int32_t i, int32_t j);

        std::vector<float> m_rdfHistogram;
        uint32_t m_accumulatedFrames = 0;
        uint32_t m_targetFrames = 1;

        uint32_t m_Zi = 1, m_Zj = 1;
        uint32_t m_Ni = 0, m_Nj = 0;

        bool m_accumulating = false;
        bool m_rdfInProgress = false;
    }; 
} // namespace sim
