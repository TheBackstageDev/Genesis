#pragma once

#include <glm/glm.hpp>

#include <algorithm>
#include "simulation/constants.hpp"
#include "simulation/physics/utils.hpp"

namespace core
{
    struct SpatialGrid
    {
        glm::ivec3 gridDimensions;

        std::vector<uint32_t> particleIndices; // first number is always the ammount of particles in the cell
        std::vector<uint32_t> cellOffsets;

        float cell_size = CELL_CUTOFF;

        void rebuild(const float* __restrict x, const float* __restrict y, const float* __restrict z, size_t N, const glm::vec3& box, float cutoff = CELL_CUTOFF)
        {
            float one_over_cutoff = 1.f / cutoff;
            gridDimensions = glm::ivec3(
                std::ceil(box.x * one_over_cutoff),
                std::ceil(box.y * one_over_cutoff),
                std::ceil(box.z * one_over_cutoff)
            );

            cell_size = cutoff;

            size_t numCells = static_cast<size_t>(gridDimensions.x) *
                            gridDimensions.y * gridDimensions.z;

            std::vector<std::vector<uint32_t>> tempLists(numCells);

            for (size_t i = 0; i < N; ++i)
            {
                glm::vec3 p = glm::vec3(x[i], y[i], z[i]);

                p.x -= cutoff * int32_t(p.x * one_over_cutoff);
                p.y -= cutoff * int32_t(p.y * one_over_cutoff);
                p.z -= cutoff * int32_t(p.z * one_over_cutoff);

                int32_t ix = static_cast<int32_t>(p.x * one_over_cutoff);
                int32_t iy = static_cast<int32_t>(p.y * one_over_cutoff);
                int32_t iz = static_cast<int32_t>(p.z * one_over_cutoff);

                ix = std::clamp(ix, 0, gridDimensions.x - 1);
                iy = std::clamp(iy, 0, gridDimensions.y - 1);
                iz = std::clamp(iz, 0, gridDimensions.z - 1);

                size_t id = static_cast<size_t>(ix + gridDimensions.x * (iy + gridDimensions.y * iz));

                tempLists[id].push_back(static_cast<uint32_t>(i));
            }

            particleIndices.clear();
            cellOffsets.clear();

            size_t totalSize = numCells + N;
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
                    uint32_t target_end = UINT32_MAX) const
        {
            if (cellId >= cellOffsets.size()) return;

            uint32_t cellStart = cellOffsets[cellId];
            uint32_t count = particleIndices[cellStart];
            if (count == 0) return;

            uint32_t start = (target_start == UINT32_MAX) ? cellStart + 1 : target_start;
            uint32_t end   = (target_end == UINT32_MAX) ? cellStart + 1 + count : target_end + 1;

            start = std::max(start, cellStart + 1);
            end   = std::min(end, cellStart + 1 + count);

            for (uint32_t idx = start; idx < end; ++idx)
            {
                func(particleIndices[idx]);
            }
        }

        glm::ivec3 positionToCell(const glm::vec3& pos) const
        {
            float one_over_cellsize = 1.f / cell_size;

            return glm::ivec3(
                static_cast<int32_t>(std::floor(pos.x * one_over_cellsize)),
                static_cast<int32_t>(std::floor(pos.y * one_over_cellsize)),
                static_cast<int32_t>(std::floor(pos.z * one_over_cellsize))
            );
        }
    };
} // namespace core
