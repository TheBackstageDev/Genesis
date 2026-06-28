#pragma once

#include "spatialgrid.hpp"
#include "simulation/core/universe.hpp"

#include <vector>
#include <glm/glm.hpp>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>

namespace core
{
    struct verlet_list
    {
        std::vector<float> old_x, old_y, old_z;

        std::vector<std::vector<uint32_t>> verlet;
        float cutoff = 10.0f;
        float skin = 1.f;
        
        void construct(const core::SpatialGrid& grid, sim::fun::universe &u);
        bool needsRebuild(sim::fun::universe& u);
    };
} // namespace core
