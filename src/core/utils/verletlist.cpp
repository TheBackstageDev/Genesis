#include "verletlist.hpp"

namespace core
{
    void verlet_list::construct(const core::SpatialGrid& grid, sim::fun::universe &u)
    {
        verlet.assign(u.numAtoms(), {});
        auto& data = u.getData();
        const float cutoff2 = cutoff * cutoff;

        #pragma omp parallel
        {
            std::vector<std::vector<uint32_t>> local_lists(u.numAtoms());

            #pragma omp for schedule(dynamic, 8)
            for (int32_t cell = 0; cell < grid.cellOffsets.size(); ++cell)
            {
                glm::ivec3 cellDims = grid.indexToCell(cell);

                grid.foreach(cell, [&](const uint32_t &i)
                {
                    const glm::vec3 pi = data.position(i);

                    for (int dx = -1; dx <= 1; ++dx)
                    for (int dy = -1; dy <= 1; ++dy)
                    for (int dz = -1; dz <= 1; ++dz)
                    {
                        int32_t n_ix = (cellDims.x + dx + grid.gridDimensions.x) % grid.gridDimensions.x;
                        int32_t n_iy = (cellDims.y + dy + grid.gridDimensions.y) % grid.gridDimensions.y;
                        int32_t n_iz = (cellDims.z + dz + grid.gridDimensions.z) % grid.gridDimensions.z;

                        int32_t neighbor_id = grid.cellToIndex(n_ix, n_iy, n_iz);
                        bool same_cell = (dx == 0 && dy == 0 && dz == 0);

                        grid.foreach(neighbor_id, [&](const uint32_t& j)
                        {
                            if (j <= i) return;

                            glm::vec3 dr = u.minImageVec(data.position(j) - pi);
                            float dr2 = glm::dot(dr, dr);

                            if (dr2 <= cutoff2)
                                local_lists[i].push_back(j);
                        }, same_cell ? i + 1 : 0, -1);
                    }
                });
            }

            #pragma omp critical
            {
                for (size_t i = 0; i < local_lists.size(); ++i)
                {
                    verlet[i].insert(verlet[i].end(),
                                     local_lists[i].begin(),
                                     local_lists[i].end());
                }
            }
        }
    }

    bool verlet_list::needsRebuild(const std::vector<glm::vec3>& old_positions,
                                   const std::vector<glm::vec3>& new_positions)
    {
        return false;
    }
}; // core