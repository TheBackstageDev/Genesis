#pragma once

#include "universe.hpp"

#include <filesystem>
#include <fstream>
#include <vector>

namespace sim
{
    struct point
    {
        float x, y;
    };

    struct histogram
    {
        std::vector<float> data;
        uint32_t capacity = 0;
        uint32_t write_idx = 0;
        uint32_t total_frames = 0;

        void push(const std::vector<float>& frame_hist);
        std::vector<float> get_average() const;
        void clear()
        {
            data.clear();
            total_frames = 0;
            write_idx = 0;
        }
    };

    class simulation_inspector
    {
    public:
        simulation_inspector();

        void set_analysis_window(uint32_t steps) { m_analysisWindow = steps; }
        uint32_t get_analysis_window() { return m_analysisWindow; }

        void startRDF(uint32_t Zi, uint32_t Zj);
        void update(fun::universe& u);

        void clear();

        // Graphs
        std::vector<point> getRDF(uint32_t Zi, uint32_t Zj, fun::universe& u) const;
        std::vector<std::vector<point>> getRDFs(fun::universe& u) const;
        std::vector<std::string> getActiveRDFs() const;

        std::vector<std::pair<std::string, uint32_t>> getMoleculesPresent(fun::universe& u) const;
        
        float calculateKineticEnergy(fun::universe& u) const;
        float calculatePotentialEnergy(fun::universe& u) const;
        float meanParticleSpeed(fun::universe& u) const;

        bool isRDFready() const { return m_activeRDFPairs.empty() && m_totalFrames == m_analysisWindow; }

        void exportData(const std::filesystem::path where);
    private:
        float ljPotential(fun::universe& u, int32_t i, int32_t j);
        float coulombPotential(fun::universe& u, int32_t i, int32_t j);

        uint32_t m_analysisWindow = 50;
        uint32_t m_totalFrames    = 0;

        struct RDFPair
        {
            uint32_t Zi, Zj;
            histogram m_histogram;
            uint32_t Ni = 0, Nj = 0;
        };

        std::unordered_map<uint64_t, RDFPair> m_rdfData;
        std::vector<std::pair<uint32_t, uint32_t>> m_activeRDFPairs;

        std::vector<float> computeFrameRDF(fun::universe& u, uint32_t Zi, uint32_t Zj) const;

        static uint64_t makeKey(uint32_t a, uint32_t b)
        {
            return (static_cast<uint64_t>(a) << 32) | b;
        }
    }; 
} // namespace sim
