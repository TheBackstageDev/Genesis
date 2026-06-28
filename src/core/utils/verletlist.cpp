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

    bool verlet_list::needsRebuild(sim::fun::universe& u)
    {
        auto& storage = u.getData();
        
        size_t N = storage.mobileCount();

        if (old_x.size() != N)
        {
            old_x.resize(N);
            old_y.resize(N);
            old_z.resize(N);
            
            return true;
        }

        const float* __restrict x = storage.xData();
        const float* __restrict y = storage.yData();
        const float* __restrict z = storage.zData();

        const float rebuild_threshold = (skin > 0.0f) ? (skin * 0.5f) : (cutoff * 0.5f);
        const float thresh2 = rebuild_threshold * rebuild_threshold;

        for (size_t i = 0; i < N; ++i)
        {
            glm::vec3 d = u.minImageVec(glm::vec3(x[i], y[i], z[i]) - glm::vec3(old_x[i], old_y[i], old_z[i]));
            if (glm::dot(d, d) > thresh2)
                return true;
        }

        return false;
    }
}; // core