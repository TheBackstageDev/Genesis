#pragma once

#include <glm/glm.hpp>
#include "simulation/constants.hpp"

namespace core
{
    struct SpatialGrid
    {
        glm::ivec3 gridDimensions;

        std::vector<uint32_t> particleIndices; // first number is always the ammount of particles in the cell
        std::vector<uint32_t> cellOffsets;

        void rebuild(const std::vector<glm::vec3>& positions, const glm::vec3& box, float cutoff = CELL_CUTOFF)
        {
            gridDimensions = glm::ivec3(
                std::ceil(box.x / cutoff),
                std::ceil(box.y / cutoff),
                std::ceil(box.z / cutoff)
            );

            size_t numCells = static_cast<size_t>(gridDimensions.x) *
                            gridDimensions.y * gridDimensions.z;

            std::vector<std::vector<uint32_t>> tempLists(numCells);

            for (size_t i = 0; i < positions.size(); ++i)
            {
                glm::vec3 p = positions[i];

                p.x -= cutoff * std::floor(p.x / cutoff);
                p.y -= cutoff * std::floor(p.y / cutoff);
                p.z -= cutoff * std::floor(p.z / cutoff);

                int32_t ix = static_cast<int32_t>(p.x / cutoff);
                int32_t iy = static_cast<int32_t>(p.y / cutoff);
                int32_t iz = static_cast<int32_t>(p.z / cutoff);

                ix = std::clamp(ix, 0, gridDimensions.x - 1);
                iy = std::clamp(iy, 0, gridDimensions.y - 1);
                iz = std::clamp(iz, 0, gridDimensions.z - 1);

                size_t id = static_cast<size_t>(ix + gridDimensions.x * (iy + gridDimensions.y * iz));

                tempLists[id].push_back(static_cast<uint32_t>(i));
            }

            particleIndices.clear();
            cellOffsets.clear();

            size_t totalSize = numCells + positions.size();
            particleIndices.reserve(totalSize);
            cellOffsets.reserve(numCells);

            uint32_t offset = 0;
            for (size_t c = 0; c < numCells; ++c)
            {
                const auto& lst = tempLists[c];
                uint32_t count = static_cast<uint32_t>(lst.size());

                cellOffsets.push_back(offset);
                offset += count + 1;

                particleIndices.push_back(count);
                particleIndices.insert(particleIndices.end(), lst.begin(), lst.end());
            }
        }

        size_t cellToIndex(const uint32_t& x, const uint32_t& y, const uint32_t& z) const {
            return static_cast<size_t>(z * gridDimensions.x * gridDimensions.y +
                                       y * gridDimensions.x + x);
        }

        glm::ivec3 indexToCell(uint32_t index) const
        {
            if (gridDimensions.x == 0 || gridDimensions.y == 0 || gridDimensions.z == 0) {
                return glm::ivec3(0);
            }

            uint32_t nxny = gridDimensions.x * gridDimensions.y;

            uint32_t z = index / nxny;
            uint32_t remainder = index % nxny;
            uint32_t y = remainder / gridDimensions.x;
            uint32_t x = remainder % gridDimensions.x;

            return glm::ivec3(static_cast<int32_t>(x),
                            static_cast<int32_t>(y),
                            static_cast<int32_t>(z));
        }

        template<typename Func>
        void foreach(size_t cellId, Func&& func,
                    uint32_t target_start = UINT32_MAX,
                    uint32_t target_end   = UINT32_MAX) const
        {
            if (cellId + 1 >= cellOffsets.size()) return;

            uint32_t cellStartOffset = cellOffsets[cellId];
            uint32_t cellCount       = particleIndices[cellStartOffset];

            uint32_t start = (target_start == UINT32_MAX) ? cellStartOffset + 1 : target_start;
            uint32_t end   = (target_end   == UINT32_MAX) ? cellStartOffset + 1 + cellCount : target_end + 1;

            start = std::max(start, cellStartOffset + 1);
            end   = std::min(end,   cellStartOffset + 1 + cellCount);

            if (start >= end) return;

            for (uint32_t idx = start; idx < end; ++idx)
            {
                uint32_t atomId = particleIndices[idx];
                func(atomId);
            }
        }
    };
} // namespace core
