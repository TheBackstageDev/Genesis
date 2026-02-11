#pragma once

#include "universe.hpp"

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
        
    private:
        float computeDensity(fun::universe& u, uint32_t Z, uint32_t frame = UINT32_MAX);
    }; 
} // namespace sim
