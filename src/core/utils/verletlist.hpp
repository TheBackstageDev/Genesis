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
        std::vector<std::vector<uint32_t>> verlet;
        float cutoff = 0.0f;
        float skin = 1.0f;
        
        void construct(const core::SpatialGrid& grid, sim::fun::universe &u);
        bool needsRebuild(const std::vector<glm::vec3>& old_positions,
                          const std::vector<glm::vec3>& new_positions);
    };
} // namespace core
