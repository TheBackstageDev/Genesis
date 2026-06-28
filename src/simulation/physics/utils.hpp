#pragma once

#include "common.hpp"

namespace sim
{
    inline glm::vec3 minImageVec(glm::vec3 dr, glm::vec3 box, bool wall_col, bool rf_col)
    {
        if (wall_col && rf_col) return dr;

        if (!rf_col)
            dr.z -= box.z * std::floor(dr.z / box.z + 0.5f);

        if (!wall_col)
        {
            dr.x -= box.x * std::floor(dr.x / box.x + 0.5f);
            dr.y -= box.y * std::floor(dr.y / box.y + 0.5f);
        }

        return dr;
    }
} // namespace sim
