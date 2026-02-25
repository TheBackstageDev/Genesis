#include "dynamics.hpp"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>

namespace sim
{
    sim_dynamics::sim_dynamics(fun::universe &u)
        : m_universe(u)
    {
    }

    void sim_dynamics::createComputeShaders()
    {
        if (unbonded_program.id() != 0 && bonded_program.id() != 0 && integrate_program.id() != 0)
            return;

        const std::filesystem::path shaders_path = "resource/shaders";

        unbonded_program = core::glProgram{
            core::glShader{GL_COMPUTE_SHADER, shaders_path / "unbonded.comp"},
        };

        bonded_program = core::glProgram{
            core::glShader{GL_COMPUTE_SHADER, shaders_path / "bonded.comp"},
        };

        integrate_program = core::glProgram{
            core::glShader{GL_COMPUTE_SHADER, shaders_path / "integrate_program.comp"},
        };
    }

    void sim_dynamics::updateSSBOs()
    {
        auto &data = m_universe.getData();
        auto &atoms = m_universe.getAtoms();
        size_t N = m_universe.numAtoms();

        if (N == 0)
            return;

        std::vector<glm::vec4> pos_data(N);
        for (size_t i = 0; i < N; ++i)
            pos_data[i] = glm::vec4(data.positions[i], atoms[i].mass);

        if (!ssbo_pos.id() ||
            ssbo_pos.getSize() != static_cast<GLsizeiptr>(N * sizeof(glm::vec4)))
        {
            ssbo_pos = core::glBuffer(
                GL_SHADER_STORAGE_BUFFER,
                pos_data.data(),
                pos_data.size() * sizeof(glm::vec4),
                GL_DYNAMIC_DRAW);
        }
        else
        {
            ssbo_pos.update(pos_data.data(), pos_data.size() * sizeof(glm::vec4));
        }

        if (!ssbo_lj.id() ||
            ssbo_lj.getSize() != static_cast<GLsizeiptr>(data.lj_params.size() * sizeof(float)))
        {
            ssbo_lj = core::glBuffer(
                GL_SHADER_STORAGE_BUFFER,
                data.lj_params.data(),
                data.lj_params.size() * sizeof(float),
                GL_STATIC_DRAW);
        }

        if (!ssbo_q.id() ||
            ssbo_q.getSize() != static_cast<GLsizeiptr>(data.q.size() * sizeof(float)))
        {
            ssbo_q = core::glBuffer(
                GL_SHADER_STORAGE_BUFFER,
                data.q.data(),
                data.q.size() * sizeof(float),
                GL_DYNAMIC_DRAW);
        }

        const GLsizeiptr force_size = N * sizeof(glm::vec4);

        bool need_recreate_force = !ssbo_force.id() || ssbo_force.getSize() != force_size;

        if (need_recreate_force)
        {
            ssbo_force = core::glBuffer(GL_SHADER_STORAGE_BUFFER, nullptr, force_size, GL_DYNAMIC_DRAW);

            glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo_force.id());
            glClearBufferData(GL_SHADER_STORAGE_BUFFER, GL_RGBA32F, GL_RGBA, GL_FLOAT, nullptr);
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
        }
    }

    void sim_dynamics::computeUnbondedGPU()
    {
        auto &data = m_universe.getData();

        const uint32_t N = static_cast<uint32_t>(m_universe.numAtoms());
        const uint32_t local_size = 256;
        const uint32_t num_groups = (N + local_size - 1) / local_size;

        unbonded_program.use();
        unbonded_program.setUniform("numParticles", N);
        unbonded_program.setUniform("box", m_universe.boxSizes());
        unbonded_program.setUniform("r_cut2", CELL_CUTOFF * CELL_CUTOFF);
        unbonded_program.setUniform("collide_floor_roof", m_universe.rooffloorcollision());
        unbonded_program.setUniform("collide_walls", m_universe.wallcollision());

        ssbo_pos.bindBase(0);
        ssbo_lj.bindBase(1);
        ssbo_q.bindBase(2);
        ssbo_force.bindBase(3);

        glDispatchCompute(num_groups, 1, 1);

        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        core::glBuffer::unbind(GL_SHADER_STORAGE_BUFFER);

        ssbo_force.bind();
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, N * sizeof(glm::vec4), data.forces.data());
        ssbo_force.unbind();
    }

    void sim_dynamics::computeBondedGPU()
    {
    }

    std::vector<glm::vec3> sim_dynamics::processCellUnbonded(int32_t ix, int32_t iy, int32_t iz, int32_t atom_start, int32_t atom_end)
    {
        float local_virial = 0.f;

        auto &data = m_universe.getData();
        auto &atomData = m_universe.getAtomData();

        int32_t cellID = universe_grid.cellToIndex(ix, iy, iz);

        thread_local std::vector<glm::vec3> local_forces(atomData.atoms.size(), {0, 0, 0});

        if (local_forces.size() != atomData.atoms.size())
            local_forces.assign(atomData.atoms.size(), {0.0f, 0.0f, 0.0f});

        std::fill(local_forces.begin(), local_forces.end(), glm::vec3{0, 0, 0});

        universe_grid.foreach (cellID, [&](const uint32_t &i)
                               {
                universe_grid.foreach(cellID, [&](const uint32_t& j)
                {
                    if (j <= i)
                        return;

                    glm::vec3 dr = m_universe.minImageVec(data.positions[j] - data.positions[i]);

                    if (m_universe.areBonded(i, j))
                        return;

                    glm::vec3 cForce = computeCoulombForce(i, j, dr);
                    glm::vec3 lForce = computeLJforce(i, j, dr);

                    glm::vec3 total_force = cForce + lForce;
                    local_forces[i] += total_force;
                    local_forces[j] -= total_force;

                    local_virial += dr.x * total_force.x +
                                    dr.y * total_force.y +
                                    dr.z * total_force.z;
                }, atom_start + 1);

                for (int32_t dx = -1; dx <= 1; ++dx)
                    for (int32_t dy = -1; dy <= 1; ++dy)
                    {
                        for (int32_t dz = -1; dz <= 1; ++dz)
                        {
                            if (dx == 0 && dy == 0 && dz == 0) continue;

                            int32_t n_ix = ix + dx;
                            int32_t n_iy = iy + dy;
                            int32_t n_iz = iz + dz;

                            n_ix = (n_ix % universe_grid.gridDimensions.x + universe_grid.gridDimensions.x) % universe_grid.gridDimensions.x;
                            n_iy = (n_iy % universe_grid.gridDimensions.y + universe_grid.gridDimensions.y) % universe_grid.gridDimensions.y;
                            n_iz = (n_iz % universe_grid.gridDimensions.z + universe_grid.gridDimensions.z) % universe_grid.gridDimensions.z;

                            int32_t neighbor_id = universe_grid.cellToIndex(n_ix, n_iy, n_iz);

                            if (!m_universe.wallcollision() && !m_universe.rooffloorcollision() && neighbor_id < -1) continue;

                            universe_grid.foreach(neighbor_id, [&](const uint32_t& j)
                            {
                                if (j <= i)
                                    return;

                                glm::vec3 dr = m_universe.minImageVec(data.positions[j] - data.positions[i]);

                                if (m_universe.areBonded(i, j))
                                    return;

                                glm::vec3 cForce = computeCoulombForce(i, j, dr);
                                glm::vec3 lForce = computeLJforce(i, j, dr);

                                glm::vec3 total_force = cForce + lForce;
                                local_forces[i] += total_force;
                                local_forces[j] -= total_force;

                                local_virial += dr.x * total_force.x +
                                                dr.y * total_force.y +
                                                dr.z * total_force.z;
                            }, atom_start, atom_end);
                        }
                    } }, atom_start, atom_end);

        total_virial.fetch_add(local_virial, std::memory_order_relaxed);
        return local_forces;
    }

    void sim_dynamics::computeUnbondedCPU()
    {
        auto &data = m_universe.getData();

        const size_t cells = universe_grid.cellOffsets.size();
        const int32_t n_threads = std::thread::hardware_concurrency();
        // const int32_t n_threads = 1;
        const int32_t threads_per_top = cells < 9 ? n_threads : std::max(1u, static_cast<uint32_t>(std::floor(static_cast<double>(n_threads / 2)))); // how many threads will run on cells with lots of work
        constexpr int32_t subdivide_top = 4;                                                                                                         // how many of the top cells to subdivide

        struct task
        {
            uint64_t work;
            int32_t cell_idx;
            int32_t atom_start;
            int32_t atom_end = -1;
        };

        struct cell_work
        {
            uint32_t cell_idx;
            uint32_t work;

            constexpr bool operator<(const cell_work &other) const
            {
                return work > other.work;
            }
        };

        std::set<cell_work> sorted_cells;

        for (int32_t c = 0; c < cells; ++c)
        {
            uint32_t offset = universe_grid.cellOffsets[c];
            uint32_t n = universe_grid.particleIndices[offset];

            if (n == 0)
                continue;

            sorted_cells.emplace(c, n * n * 9);
        }

        std::vector<task> tasks;

        int32_t rank = 0;
        for (const auto &cw : sorted_cells)
        {
            int32_t c = cw.cell_idx;
            uint32_t offset = universe_grid.cellOffsets[c];
            uint32_t n = universe_grid.particleIndices[offset];

            if (rank < subdivide_top)
            {
                int32_t slice_size = static_cast<int32_t>((n + threads_per_top) / threads_per_top);
                for (int32_t p = 0; p < threads_per_top; ++p)
                {
                    int32_t start = p * slice_size;
                    int32_t end = std::min(start + slice_size, static_cast<int32_t>(n));
                    if (start >= end)
                        break;

                    uint64_t sub_work = static_cast<uint64_t>(end - start) * n * 9ull;
                    tasks.emplace_back(sub_work, c, start, end);
                }
            }
            else
            {
                tasks.push_back({cw.work, c, 0, -1});
            }

            ++rank;
        }

        auto worker = [&](const std::vector<task> &my_tasks) -> std::vector<glm::vec3>
        {
            std::vector<glm::vec3> thread_forces(m_universe.numAtoms(), {0, 0, 0});

            for (auto &t : my_tasks)
            {
                glm::ivec3 cell = universe_grid.indexToCell(t.cell_idx);

                uint32_t offset = universe_grid.cellOffsets[t.cell_idx];
                uint32_t n = universe_grid.particleIndices[offset];

                std::vector<glm::vec3> cell_forces;

                if (t.atom_end < 0)
                {
                    cell_forces = processCellUnbonded(cell.x, cell.y, cell.z, 0, universe_grid.particleIndices[n]);
                }
                else
                {
                    cell_forces = processCellUnbonded(cell.x, cell.y, cell.z, t.atom_start, t.atom_end);
                }

                for (size_t i = 0; i < m_universe.numAtoms(); ++i)
                {
                    thread_forces[i] += cell_forces[i];
                }
            }
            return thread_forces;
        };

        int32_t n_used_threads = std::min(n_threads, static_cast<int32_t>(tasks.size()));
        std::vector<uint64_t> thread_load(n_used_threads, 0);
        std::vector<std::vector<task>> thread_tasks(n_used_threads);

        for (const auto &task : tasks)
        {
            auto min_it = std::min_element(thread_load.begin(), thread_load.end());
            int32_t t = std::distance(thread_load.begin(), min_it);

            thread_tasks[t].push_back(task);
            thread_load[t] += task.work;
        }

        std::vector<std::future<std::vector<glm::vec3>>> futures;
        futures.reserve(n_used_threads);

        for (int32_t t = 0; t < n_used_threads; ++t)
        {
            if (!thread_tasks[t].empty())
            {
                futures.push_back(
                    std::async(std::launch::async,
                               worker,
                               std::cref(thread_tasks[t])));
            }
        }

        for (auto &fut : futures)
        {
            auto local_f = fut.get();

            for (size_t i = 0; i < m_universe.numAtoms(); ++i)
                data.forces[i] += local_f[i];
        }
    }

    void sim_dynamics::computeBondedCPU()
    {
        auto &univ = m_universe;
        auto &atomData = univ.getAtomData();
        auto &data = univ.getData();
        const uint32_t n_threads = std::max(1u, std::thread::hardware_concurrency());

        std::vector<std::future<std::vector<glm::vec3>>> futures;

        auto make_task = [this, &atomData](auto &&func)
        {
            return [this, &atomData, func](int32_t start, int32_t end) -> std::vector<glm::vec3>
            {
                thread_local std::vector<glm::vec3> local_forces;
                if (local_forces.size() != atomData.atoms.size())
                {
                    local_forces.assign(atomData.atoms.size(), glm::vec3{0.0f});
                }
                std::fill(local_forces.begin(), local_forces.end(), glm::vec3{0.0f});
                func(start, end, local_forces);
                return local_forces;
            };
        };

        auto bond_func = [&](int32_t start, int32_t end, std::vector<glm::vec3> &lf)
        {
            for (int32_t i = start; i < end; ++i)
            {
                const sim::fun::bond &b = atomData.bonds[i];
                int32_t a = b.bondedAtom;
                int32_t c = b.centralAtom;

                glm::vec3 dr = univ.minImageVec(data.positions[c] - data.positions[a]);
                float len = glm::length(dr);
                if (len <= EPSILON)
                    continue;

                float dl = len - b.equilibriumLength;
                glm::vec3 force = (b.k * dl / len) * dr;

                lf[a] += force;
                lf[c] -= force;
            }
        };

        for (int32_t t = 0; t < n_threads; ++t)
        {
            int32_t start = atomData.bonds.size() * t / n_threads;
            int32_t end = atomData.bonds.size() * (t + 1) / n_threads;
            futures.emplace_back(std::async(std::launch::async, make_task(bond_func), start, end));
        }

        auto angle_func = [&](int32_t start, int32_t end, std::vector<glm::vec3> &lf)
        {
            for (int32_t a = start; a < end; ++a)
            {
                const fun::angle &ang = atomData.angles[a];
                int32_t i = ang.A, j = ang.B, k = ang.C;

                glm::vec3 r_ji = univ.minImageVec(data.positions[i] - data.positions[j]);
                glm::vec3 r_jk = univ.minImageVec(data.positions[k] - data.positions[j]);

                float len_ji = glm::length(r_ji);
                float len_jk = glm::length(r_jk);
                if (len_ji < EPSILON || len_jk < EPSILON)
                    continue;

                glm::vec3 u_ji = r_ji / len_ji;
                glm::vec3 u_jk = r_jk / len_jk;

                float cos_theta = std::clamp(glm::dot(u_ji, u_jk), -1.0f, 1.0f);
                float sin_theta = std::sqrt(std::max(1.0f - cos_theta * cos_theta, 0.0f));
                if (sin_theta < 1e-6f)
                    sin_theta = 1e-6f;

                float theta = std::acos(cos_theta);
                float delta_theta = theta - ang.rad;

                glm::vec3 dtheta_dri = (cos_theta * u_ji - u_jk) / (len_ji * sin_theta);
                glm::vec3 dtheta_drk = (cos_theta * u_jk - u_ji) / (len_jk * sin_theta);

                glm::vec3 F_i = -ang.K * delta_theta * dtheta_dri;
                glm::vec3 F_k = -ang.K * delta_theta * dtheta_drk;
                glm::vec3 F_j = -F_i - F_k;

                lf[i] += F_i;
                lf[j] += F_j;
                lf[k] += F_k;
            }
        };

        for (int32_t t = 0; t < n_threads; ++t)
        {
            int32_t start = atomData.angles.size() * t / n_threads;
            int32_t end = atomData.angles.size() * (t + 1) / n_threads;
            futures.emplace_back(std::async(std::launch::async, make_task(angle_func), start, end));
        }

        auto dihedral_func = [&](int32_t start, int32_t end, std::vector<glm::vec3> &lf)
        {
            for (int32_t d = start; d < end; ++d)
            {
                const fun::dihedral_angle &da = atomData.dihedral_angles[d];

                glm::vec3 pa = data.positions[da.A];
                glm::vec3 pb = data.positions[da.B];
                glm::vec3 pc = data.positions[da.C];
                glm::vec3 pd = data.positions[da.D];

                float phi = computeDihedral(pa, pb, pc, pd);

                float target = da.rad;
                if (target == 0.0f && da.periodicity == 1)
                {
                    int32_t chi = atomData.atoms[da.B].chirality ? atomData.atoms[da.B].chirality : atomData.atoms[da.C].chirality;
                    if (chi == 1)
                        target = (phi < M_PI) ? 1.047f : 5.236f;
                    if (chi == 2)
                        target = (phi < M_PI) ? 5.236f : 1.047f;
                }

                float diff = phi - target;
                while (diff > M_PI)
                    diff -= 2.0f * M_PI;
                while (diff < -M_PI)
                    diff += 2.0f * M_PI;

                glm::vec3 b1 = pb - pa;
                glm::vec3 b2 = pc - pb;
                glm::vec3 b3 = pd - pc;

                glm::vec3 n1 = glm::normalize(glm::cross(b1, b2));
                glm::vec3 n2 = glm::normalize(glm::cross(b2, b3));
                glm::vec3 u2 = glm::normalize(b2);

                float sin_term = std::sin(da.periodicity * phi - da.rad);
                float torque_mag = da.K * da.periodicity * sin_term;

                glm::vec3 axis = glm::normalize(pc - pb);
                glm::vec3 rA = pa - pb;
                glm::vec3 torqueA = glm::normalize(glm::cross(axis, rA)) * torque_mag;

                glm::vec3 rD = glm::normalize(pd - pc);
                glm::vec3 torqueD = glm::normalize(glm::cross(axis, rD)) * (-torque_mag);

                glm::vec3 f1 = (torque_mag / (glm::length(b1) + 1e-6f)) * glm::cross(n1, u2);
                glm::vec3 f4 = (torque_mag / (glm::length(b3) + 1e-6f)) * glm::cross(u2, n2);

                lf[da.A] += -f1;
                lf[da.D] += -f4;
                lf[da.B] += f1;
                lf[da.C] += f4;
                lf[da.D] += -torqueD * 0.5f;
            }
        };

        for (int32_t t = 0; t < n_threads; ++t)
        {
            int32_t start = atomData.dihedral_angles.size() * t / n_threads;
            int32_t end = atomData.dihedral_angles.size() * (t + 1) / n_threads;
            futures.emplace_back(std::async(std::launch::async, make_task(dihedral_func), start, end));
        }

        auto improper_func = [&](int32_t start, int32_t end, std::vector<glm::vec3> &lf)
        {
            for (int32_t d = start; d < end; ++d)
            {
                const fun::dihedral_angle &imp = atomData.improper_angles[d];

                glm::vec3 pa = data.positions[imp.A];
                glm::vec3 pb = data.positions[imp.B];
                glm::vec3 pc = data.positions[imp.C];
                glm::vec3 pd = data.positions[imp.D];

                float phi = computeDihedral(pa, pb, pc, pd);

                float diff = phi - imp.rad;
                while (diff > M_PI)
                    diff -= 2.0f * M_PI;
                while (diff < -M_PI)
                    diff += 2.0f * M_PI;

                float dE_dphi = imp.K * diff;

                glm::vec3 b1 = pb - pa;
                glm::vec3 b2 = pc - pb;
                glm::vec3 b3 = pd - pc;

                glm::vec3 n1 = glm::normalize(glm::cross(b1, b2));
                glm::vec3 n2 = glm::normalize(glm::cross(b2, b3));
                glm::vec3 u2 = glm::normalize(b2);

                float sin_term = std::sin(imp.periodicity * phi - imp.rad);
                float torque_mag = imp.K * imp.periodicity * sin_term;

                glm::vec3 axis = glm::normalize(pc - pb);
                glm::vec3 rA = pa - pb;
                glm::vec3 torqueA = glm::normalize(glm::cross(axis, rA)) * torque_mag;

                glm::vec3 rD = glm::normalize(pd - pc);
                glm::vec3 torqueD = glm::normalize(glm::cross(axis, rD)) * (-torque_mag);

                glm::vec3 f1 = (torque_mag / (glm::length(b1) + 1e-6f)) * glm::cross(n1, u2);
                glm::vec3 f4 = (torque_mag / (glm::length(b3) + 1e-6f)) * glm::cross(u2, n2);

                lf[imp.A] += -f1;
                lf[imp.D] += -f4;
                lf[imp.B] += f1 + f4;
                lf[imp.C] += -torqueD * 0.5f;
            }
        };

        for (int32_t t = 0; t < n_threads; ++t)
        {
            int32_t start = atomData.improper_angles.size() * t / n_threads;
            int32_t end = atomData.improper_angles.size() * (t + 1) / n_threads;
            futures.emplace_back(std::async(std::launch::async, make_task(improper_func), start, end));
        }

        for (auto &fut : futures)
        {
            auto local_f = fut.get();
            for (size_t i = 0; i < local_f.size(); ++i)
            {
                data.forces[i] += local_f[i];
            }
        }
    }

    glm::vec3 sim_dynamics::computeCoulombForce(uint32_t i, uint32_t j, glm::vec3 &dr_vec)
    {
        auto &data = m_universe.getData();
        float dr2 = glm::length2(dr_vec);

        if (dr2 < EPSILON || dr2 > COULOMB_CUTOFF * COULOMB_CUTOFF)
            return {0.f, 0.f, 0.f};

        float qq = data.q[i] * data.q[j];
        if (qq == 0.f)
            return {0.f, 0.f, 0.f};

        float inv_r = 1.0f / sqrtf(dr2);
        float forceMag = COULOMB_K * qq * inv_r * inv_r * inv_r;
        return forceMag * -dr_vec;
    }

    glm::vec3 sim_dynamics::computeLJforce(uint32_t i, uint32_t j, glm::vec3 &dr_vec)
    {
        auto &data = m_universe.getData();

        const uint32_t base_i = i << 1;
        const uint32_t base_j = j << 1;

        const float sigma_i = data.lj_params[base_i];
        const float epsilon_i = data.lj_params[base_i + 1];

        const float sigma_j = data.lj_params[base_j];
        const float epsilon_j = data.lj_params[base_j + 1];

        float dr2 = glm::length2(dr_vec);

        const float sigma = (sigma_i + sigma_j) * 0.5f;
        const float epsilon = sqrtf(epsilon_i * epsilon_j);

        if (dr2 < 100.f && dr2 > EPSILON)
        {
            float inv_r = 1.f / sqrtf(dr2);
            float inv_r2 = 1.f / dr2;
            float inv_r6 = inv_r2 * inv_r2 * inv_r2;
            float inv_r12 = inv_r6 * inv_r6;

            float sr6 = powf(sigma, 6.0f) * inv_r6;
            float sr12 = sr6 * sr6;

            float force_mag = 24.0f * epsilon * (2.0f * sr12 - sr6) * inv_r2;

            return force_mag * -dr_vec;
        }

        return glm::vec3{0.f, 0.f, 0.f};
    }

    float sim_dynamics::computeDihedral(const glm::vec3 &pa, const glm::vec3 &pb, const glm::vec3 &pc, const glm::vec3 &pd)
    {
        glm::vec3 v1 = pb - pa;
        glm::vec3 v2 = pc - pb;
        glm::vec3 v3 = pd - pc;

        glm::vec3 n1 = glm::cross(v1, v2);
        glm::vec3 n2 = glm::cross(v2, v3);

        n1 = glm::normalize(n1);
        n2 = glm::normalize(n2);

        float cosPhi = std::clamp(glm::dot(n1, n2), -1.f, 1.f);
        float phi = std::acos(cosPhi);

        // sign via right-hand rule
        if (glm::dot(v2, glm::cross(n1, n2)) < 0.f)
            phi = -phi;
        if (phi < 0.f)
            phi += 2.f * M_PI;
        return phi;
    }

    void sim_dynamics::zeroForces()
    {
        auto &data = m_universe.getData();
        if (m_GPU)
        {
            ssbo_force.bind();
            glClearBufferData(GL_SHADER_STORAGE_BUFFER, GL_RGBA32F, GL_RGBA, GL_FLOAT, nullptr);
            ssbo_force.unbind();
        }
        else
        {
            std::fill(data.forces.begin(), data.forces.end(), glm::vec3{0.0f});
        }
    }

    void sim_dynamics::computeExternalForces(uint32_t i)
    {
        auto& wall_charges = m_universe.getWallCharges();
        auto& data = m_universe.getData();

        for (int32_t w = 0; w < m_universe.getWallCharges().size() && m_universe.wallChargeEnabled(); ++w)
        {
            const float &wall_q = wall_charges[w];
            if (wall_q == 0.f)
                continue;

            glm::vec3 force = data.q[i] * wall_q * wall_directions[w];
            data.forces[i] += force;
        }

        if (m_universe.magneticFieldEnabled())
        {
            glm::vec3 force = data.q[i] * glm::cross(data.velocities[i], m_universe.getMagneticFieldStrength());
            data.forces[i] += force;
        }
    }

    void sim_dynamics::computeForces()
    {
        zeroForces();

        if (m_GPU)
        {
            computeBondedGPU();
            computeUnbondedGPU();
        }
        else
        {
            computeBondedCPU();
            computeUnbondedCPU();
        }
    }

    void sim_dynamics::integrate()
    {
        auto &univ = m_universe;
        auto &data = m_universe.getData();
        auto &atomData = m_universe.getAtomData();
        const size_t N = data.positions.size();

        if (N == 0)
            return;

        computeForces();

        const float effective_dt = m_dt * m_timescale;
        const float half_dt = 0.5f * effective_dt;

        for (size_t i = 0; i < N; ++i)
        {
            if (atomData.frozen_atoms[i])
            {
                data.velocities[i] = glm::vec3{0.0f};
                continue;
            }
            
            computeExternalForces(i);

            const float inv_m = 1.0f / atomData.atoms[i].mass;
            const glm::vec3 accel = data.forces[i] * inv_m;

            data.velocities[i] += accel * half_dt;
        }

        for (size_t i = 0; i < N; ++i)
        {
            if (atomData.frozen_atoms[i])
            {
                data.velocities[i] = glm::vec3{0.0f};
                continue;
            }

            data.positions[i] += data.velocities[i] * m_dt;
            univ.boundCheck(i);
        }

        computeForces();

        if (m_GPU)
        {
            ssbo_force.bind();
            glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0,
                               data.forces.size() * sizeof(glm::vec4),
                               data.forces.data());
            ssbo_force.unbind();
        }

        for (size_t i = 0; i < N; ++i)
        {
            if (atomData.frozen_atoms[i])
            {
                data.velocities[i] = glm::vec3{0.0f};
                continue;
            }

            computeExternalForces(i);

            const float inv_m = 1.0f / atomData.atoms[i].mass;
            const glm::vec3 accel = data.forces[i] * inv_m;

            data.velocities[i] += accel * half_dt;
        }
    }

    void sim_dynamics::step(float target_temp, float target_pressure)
    {
        if (m_paused)
            return;

        auto &data = m_universe.getData();

        integrate();

        if (m_step_count % GRID_REBUILD == 0)
            universe_grid.rebuild(data.positions, m_universe.boxSizes(), CELL_CUTOFF);

        setTemperature(target_temp);
        setPressure(target_pressure);
        COMDrift();

        m_step_count++;
        m_accumulated_time += m_dt * m_timescale;
    }

    float sim_dynamics::computePressure()
    {
        auto &data = m_universe.getData();
        auto &atomData = m_universe.getAtomData();

        if (atomData.atoms.empty())
            return 0.f;

        float kinetic = m_universe.calculateKineticEnergy();

        glm::vec3 box = m_universe.boxSizes();
        float volume = box.x * box.y * box.z;
        float temperature = (2.0f / 3.0f) * kinetic / atomData.atoms.size();

        float P_ideal = atomData.atoms.size() * temperature / volume;

        float P_virial = -total_virial / (3.0f * volume);

        float pressure = P_ideal + P_virial;
        return pressure;
    }

    void sim_dynamics::setPressure(float Target_P_Bar)
    {
        if (Target_P_Bar <= 0.0f)
            return;

        if (m_step_count % BAROSTAT_INTERVAL != 0)
            return;

        glm::vec3 box = m_universe.boxSizes();
        auto &data = m_universe.getData();
        auto &atomData = m_universe.getAtomData();

        constexpr float beta_T = 4.5e-5f;
        constexpr float tau_P = 1.0f;

        m_pressure = computePressure();

        float delta_P = Target_P_Bar - m_pressure;
        float mu = 1.0f - (2.f / tau_P) * beta_T * delta_P;

        mu = std::clamp(mu, 0.5f, 1.5f);

        float scale = std::cbrt(mu);

        box.z *= scale;

        for (int32_t i = 0; i < atomData.atoms.size(); ++i)
        {
            data.positions[i].z *= scale;
        }
    }

    // Bussi–Donadio–Parrinello (CSVR) velocity rescaling
    void sim_dynamics::setTemperature(float kelvin)
    {
        if (m_step_count % THERMOSTAT_INTERVAL != 0)
            return;

        auto &atomData = m_universe.getAtomData();
        auto &data = m_universe.getData();

        size_t N = atomData.atoms.size();
        if (N == 0)
            return;

        float ndof = 3.0f * static_cast<float>(N) - 3.0f;
        if (ndof <= 0.0f)
            ndof = 1.0f;

        float KE = m_universe.calculateKineticEnergy();
        float current_temp = (2.0f * KE) / (ndof * KB);
        m_temperature = current_temp;

        float target_KE = 0.5f * ndof * KB * kelvin;

        float c = target_KE / (KE + 1e-10f);

        float r = gauss_random();
        float chi = r * sqrtf(2.0f / ndof);

        float alpha = sqrtf(c * (1.0f + chi + 0.5f * chi * chi));

        for (auto &v : data.velocities)
            v *= alpha;
    }

    void sim_dynamics::COMDrift()
    {
        if (m_step_count % 1000 != 0)
            return;

        auto &atomData = m_universe.getAtomData();
        auto &data = m_universe.getData();

        glm::vec3 totalMomentum{0.f};
        float totalMass = 0.f;

        for (int32_t i = 0; i < atomData.atoms.size(); ++i)
        {
            totalMomentum += atomData.atoms[i].mass * data.velocities[i];
            totalMass += atomData.atoms[i].mass;
        }

        glm::vec3 correction = totalMomentum / totalMass;

        for (auto &v : data.velocities)
        {
            v -= correction;
        }
    }
} // namespace sim
