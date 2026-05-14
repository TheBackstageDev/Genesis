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
        std::vector<point> getRDFgraph(fun::universe& u, uint32_t Zi, uint32_t Zj);
        std::vector<std::pair<std::string, uint32_t>> getMoleculesPresent(fun::universe& u);
        //std::vector<std::pair<std::string, uint32_t>> getMoleculesPresent(fun::universe& u, bool reactive = true);
        
    private:
        float computeDensity(fun::universe& u, uint32_t Z, uint32_t frame = UINT32_MAX);
    }; 
} // namespace sim
