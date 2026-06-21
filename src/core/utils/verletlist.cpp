#include "verletlist.hpp"

namespace core
{
    void verlet_list::construct(const core::SpatialGrid& grid, sim::fun::universe &u)
    {
        verlet.assign(u.numAtoms(), {});
        auto& data = u.getData();

        for (int32_t cell = 0; cell < grid.cellOffsets.size(); ++cell)
        {
            glm::ivec3 cellDims = grid.indexToCell(cell);

            grid.foreach(cell, [&](const uint32_t &i)
            {
                glm::vec3 pi = data.positions[i];

                for (int32_t dx = -1; dx <= 1; ++dx)
                for (int32_t dy = -1; dy <= 1; ++dy)
                for (int32_t dz = -1; dz <= 1; ++dz)
                {
                    int32_t n_ix = cellDims.x + dx;
                    int32_t n_iy = cellDims.y + dy;
                    int32_t n_iz = cellDims.z + dz;

                    n_ix = (n_ix % grid.gridDimensions.x + grid.gridDimensions.x) % grid.gridDimensions.x;
                    n_iy = (n_iy % grid.gridDimensions.y + grid.gridDimensions.y) % grid.gridDimensions.y;
                    n_iz = (n_iz % grid.gridDimensions.z + grid.gridDimensions.z) % grid.gridDimensions.z;

                    int32_t neighbor_id = grid.cellToIndex(n_ix, n_iy, n_iz);

                    if (!u.wallcollision() && !u.rooffloorcollision() && neighbor_id < -1) 
                        continue;

                    bool same_cell = (dx == 0 && dy == 0 && dz == 0);

                    grid.foreach(neighbor_id, [&](const uint32_t& j)
                    {
                        if (j <= i) return;

                        glm::vec3 dr = u.minImageVec(data.positions[j] - pi);
                        float dr2 = glm::length2(dr);

                        if (dr2 <= cutoff * cutoff)
                            verlet[i].push_back(j);
                    }, same_cell ? i + 1 : 0, -1);
                }
            });
        }

        for (auto& lst : verlet)
            std::sort(lst.begin(), lst.end());
    }

    bool verlet_list::needsRebuild(const std::vector<glm::vec3>& old_positions,
                                   const std::vector<glm::vec3>& new_positions)
    {
        return false;
    }
}; // core