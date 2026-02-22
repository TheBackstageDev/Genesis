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
        auto& data = m_universe.getData();
        auto& atoms = m_universe.getAtoms();
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
        /* const uint32_t N = static_cast<uint32_t>(m_universe.numAtoms());
        const uint32_t local_size = 256;
        const uint32_t num_groups = (N + local_size - 1) / local_size;

        unbonded_program.use();
        unbonded_program.setUniform("numParticles", N);
        unbonded_program.setUniform("box", box);
        unbonded_program.setUniform("r_cut2", CELL_CUTOFF * CELL_CUTOFF);
        unbonded_program.setUniform("collide_floor_roof", roof_floor_collision);
        unbonded_program.setUniform("collide_walls", wall_collision);

        ssbo_positions.bindBase(0);
        ssbo_lj.bindBase(1);
        ssbo_q.bindBase(2);
        ssbo_force.bindBase(3);

        glDispatchCompute(num_groups, 1, 1);

        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        core::glBuffer::unbind(GL_SHADER_STORAGE_BUFFER);

        ssbo_force.bind();
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, N * sizeof(glm::vec4), data.forces.data());
        ssbo_force.unbind(); */
    }

    void sim_dynamics::computeBondedGPU()
    {
    }

    std::vector<glm::vec3> sim_dynamics::processCellUnbonded(int32_t ix, int32_t iy, int32_t iz, int32_t atom_start, int32_t atom_end)
    {
        return {};
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

        auto worker = [this](const std::vector<task> &my_tasks) -> std::vector<glm::vec3>
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
        /* auto& atoms = m_universe.getAtoms();

        const uint32_t n_threads = std::max(1u, std::thread::hardware_concurrency());
        std::vector<std::future<std::vector<sf::Vector3f>>> futures;

        auto make_task = [this](auto &&func)
        {
            return [this, func](int32_t start, int32_t end) -> std::vector<sf::Vector3f>
            {
                thread_local std::vector<sf::Vector3f> local_forces;
                if (local_forces.size() != atoms.size())
                {
                    local_forces.assign(atoms.size(), {0.0f, 0.0f, 0.0f});
                }
                std::fill(local_forces.begin(), local_forces.end(), sf::Vector3f{0, 0, 0});

                func(start, end, local_forces);

                return local_forces;
            };
        };

        auto bond_func = [this](int32_t start, int32_t end, std::vector<sf::Vector3f> &lf)
        {
            for (int32_t i = start; i < end; ++i)
            {
                const bond &b = bonds[i];
                int32_t a = b.bondedAtom;
                int32_t c = b.centralAtom;

                sf::Vector3f dr = minImageVec(pos(c) - pos(a));
                float len = dr.length();
                if (len <= EPSILON)
                    continue;

                float dl = len - b.equilibriumLength;
                sf::Vector3f force = (b.k * dl / len) * dr;

                lf[a] += force;
                lf[c] -= force;
            }
        };

        for (int32_t t = 0; t < n_threads; ++t)
        {
            int32_t start = bonds.size() * t / n_threads;
            int32_t end = bonds.size() * (t + 1) / n_threads;
            futures.emplace_back(std::async(std::launch::async, make_task(bond_func), start, end));
        }

        auto angle_func = [this](int32_t start, int32_t end, std::vector<sf::Vector3f> &lf)
        {
            for (int32_t a = start; a < end; ++a)
            {
                const angle &ang = angles[a];
                int32_t i = ang.A, j = ang.B, k = ang.C;

                sf::Vector3f r_ji = minImageVec(pos(i) - pos(j));
                sf::Vector3f r_jk = minImageVec(pos(k) - pos(j));
                float len_ji = r_ji.length();
                float len_jk = r_jk.length();
                if (len_ji < EPSILON || len_jk < EPSILON)
                    continue;

                sf::Vector3f u_ji = r_ji / len_ji;
                sf::Vector3f u_jk = r_jk / len_jk;
                float cos_theta = std::clamp(u_ji.dot(u_jk), -1.0f, 1.0f);
                float sin_theta = std::sqrt(std::max(1.0f - cos_theta * cos_theta, 0.0f));
                if (sin_theta < 1e-6f)
                    sin_theta = 1e-6f;

                float theta = std::acos(cos_theta);
                float delta_theta = theta - ang.rad;

                sf::Vector3f dtheta_dri = (cos_theta * u_ji - u_jk) / (len_ji * sin_theta);
                sf::Vector3f dtheta_drk = (cos_theta * u_jk - u_ji) / (len_jk * sin_theta);

                sf::Vector3f F_i = -ang.K * delta_theta * dtheta_dri;
                sf::Vector3f F_k = -ang.K * delta_theta * dtheta_drk;
                sf::Vector3f F_j = -F_i - F_k;

                lf[i] += F_i;
                lf[j] += F_j;
                lf[k] += F_k;
            }
        };

        for (int32_t t = 0; t < n_threads; ++t)
        {
            int32_t start = angles.size() * t / n_threads;
            int32_t end = angles.size() * (t + 1) / n_threads;
            futures.emplace_back(std::async(std::launch::async, make_task(angle_func), start, end));
        }

        auto dihedral_func = [this](int32_t start, int32_t end, std::vector<sf::Vector3f> &lf)
        {
            for (int32_t d = start; d < end; ++d)
            {
                const dihedral_angle &da = dihedral_angles[d];
                const sf::Vector3f &pa = pos(da.A);
                const sf::Vector3f &pb = pos(da.B);
                const sf::Vector3f &pc = pos(da.C);
                const sf::Vector3f &pd = pos(da.D);

                glm::vec3 pb_g = glm::vec3(pb.x, pb.y, pb.z);

                float phi = calculateDihedral(glm::vec3(pa.x, pa.y, pa.z), pb_g, glm::vec3(pc.x, pc.y, pc.z), glm::vec3(pd.x, pd.y, pd.z));

                float target = da.rad;
                if (target == 0.0f && da.periodicity == 1)
                {
                    int32_t chi = atoms[da.B].chirality ? atoms[da.B].chirality : atoms[da.C].chirality;
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

                glm::vec3 b1 = pb_g - glm::vec3(pa.x, pa.y, pa.z);
                glm::vec3 b2 = glm::vec3(pc.x, pc.y, pc.z) - pb_g;
                glm::vec3 b3 = glm::vec3(pd.x, pd.y, pd.z) - glm::vec3(pc.x, pc.y, pc.z);

                glm::vec3 n1 = glm::cross(b1, b2);
                glm::vec3 n2 = glm::cross(b2, b3);

                n1 = glm::normalize(n1);
                n2 = glm::normalize(n2);
                glm::vec3 u2 = glm::normalize(b2);

                float sin_term = std::sin(da.periodicity * phi - da.rad);
                float torque_mag = da.K * da.periodicity * sin_term;

                glm::vec3 axis = glm::normalize(glm::vec3(pc.x, pc.y, pc.z) - pb_g);

                glm::vec3 rA = glm::vec3(pa.x, pa.y, pa.z) - pb_g;
                glm::vec3 torqueA = glm::normalize(glm::cross(axis, rA)) * torque_mag;

                glm::vec3 rD = glm::normalize(glm::vec3(pd.x, pd.y, pd.z) - glm::vec3(pc.x, pc.y, pc.z));
                glm::vec3 torqueD = glm::normalize(glm::cross(axis, rD)) * (-torque_mag);

                glm::vec3 f1 = (torque_mag / glm::length(b1) + 0.0001f) * (glm::cross(n1, u2));
                glm::vec3 f4 = (torque_mag / glm::length(b3) + 0.0001f) * (glm::cross(u2, n2));

                sf::Vector3f f1_s = sf::Vector3f(f1.x, f1.y, f1.z);
                sf::Vector3f torqueD_s = sf::Vector3f(torqueD.x, torqueD.y, torqueD.z);
                sf::Vector3f f4_s = sf::Vector3f(f4.x, f4.y, f4.z);

                lf[da.A] += -f1_s;
                lf[da.D] += -f4_s;
                lf[da.B] += f1_s;
                lf[da.C] += f4_s;
                lf[da.D] += -torqueD_s * 0.5f;
            }
        };

        for (int32_t t = 0; t < n_threads; ++t)
        {
            int32_t start = dihedral_angles.size() * t / n_threads;
            int32_t end = dihedral_angles.size() * (t + 1) / n_threads;
            futures.emplace_back(std::async(std::launch::async, make_task(dihedral_func), start, end));
        }

        auto improper_func = [this](int32_t start, int32_t end, std::vector<sf::Vector3f> &lf)
        {
            for (int32_t d = start; d < end; ++d)
            {
                const dihedral_angle &imp = improper_angles[d];

                const sf::Vector3f &pa = pos(imp.A);
                const sf::Vector3f &pb = pos(imp.B);
                const sf::Vector3f &pc = pos(imp.C);
                const sf::Vector3f &pd = pos(imp.D);

                glm::vec3 pb_g(pb.x, pb.y, pb.z);

                float phi = calculateDihedral(
                    glm::vec3(pa.x, pa.y, pa.z),
                    pb_g,
                    glm::vec3(pc.x, pc.y, pc.z),
                    glm::vec3(pd.x, pd.y, pd.z));

                float diff = phi - imp.rad;
                while (diff > M_PI)
                    diff -= 2.0f * M_PI;
                while (diff < -M_PI)
                    diff += 2.0f * M_PI;

                float dE_dphi = imp.K * diff;

                glm::vec3 b1 = pb_g - glm::vec3(pa.x, pa.y, pa.z);
                glm::vec3 b2 = glm::vec3(pc.x, pc.y, pc.z) - pb_g;
                glm::vec3 b3 = glm::vec3(pd.x, pd.y, pd.z) - glm::vec3(pc.x, pc.y, pc.z);

                glm::vec3 n1 = glm::normalize(glm::cross(b1, b2));
                glm::vec3 n2 = glm::normalize(glm::cross(b2, b3));
                glm::vec3 u2 = glm::normalize(b2);

                float sin_term = std::sin(imp.periodicity * phi - imp.rad);
                float torque_mag = imp.K * imp.periodicity * sin_term;

                glm::vec3 axis = glm::normalize(glm::vec3(pc.x, pc.y, pc.z) - pb_g);

                glm::vec3 rA = glm::vec3(pa.x, pa.y, pa.z) - pb_g;
                glm::vec3 torqueA = glm::normalize(glm::cross(axis, rA)) * torque_mag;

                glm::vec3 rD = glm::normalize(glm::vec3(pd.x, pd.y, pd.z) - glm::vec3(pc.x, pc.y, pc.z));
                glm::vec3 torqueD = glm::normalize(glm::cross(axis, rD)) * (-torque_mag);

                glm::vec3 f1 = (torque_mag / (glm::length(b1) + 1e-6f)) * glm::cross(n1, u2);
                glm::vec3 f4 = (torque_mag / (glm::length(b3) + 1e-6f)) * glm::cross(u2, n2);

                sf::Vector3f f1_s(f1.x, f1.y, f1.z);
                sf::Vector3f f4_s(f4.x, f4.y, f4.z);
                sf::Vector3f torqueD_s(torqueD.x, torqueD.y, torqueD.z);

                lf[imp.A] += -f1_s;
                lf[imp.D] += -f4_s;
                lf[imp.B] += f1_s + f4_s;
                lf[imp.C] += -torqueD_s * 0.5f;
            }
        };

        for (int t = 0; t < n_threads; ++t)
        {
            int start = improper_angles.size() * t / n_threads;
            int end = improper_angles.size() * (t + 1) / n_threads;
            futures.emplace_back(std::async(std::launch::async, make_task(improper_func), start, end));
        }

        for (auto &fut : futures)
        {
            auto local_f = fut.get();
            for (int32_t i = 0; i < atoms.size(); ++i)
                add_force(i, local_f[i]);
        } */
    }

    glm::vec3 sim_dynamics::computeCoulombForce(uint32_t i, uint32_t j, glm::vec3 &dr_vec)
    {
        auto& data = m_universe.getData();
        float dr = glm::length2(dr_vec);

        if (dr < EPSILON || dr > COULOMB_CUTOFF * COULOMB_CUTOFF)
            return {0.f, 0.f, 0.f};

        float qq = data.q[i] * data.q[j];
        if (qq == 0.f)
            return {0.f, 0.f, 0.f};

        float forceMag = COULOMB_K * qq / (dr * dr * dr);
        return forceMag * -dr_vec;
    }

    glm::vec3 sim_dynamics::computeLJforce(uint32_t i, uint32_t j, glm::vec3 &dr_vec)
    {
        auto& data = m_universe.getData();

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
            float inv_r2 = 1.f / dr2;
            float inv_r6 = inv_r2 * inv_r2 * inv_r2;
            float inv_r12 = inv_r6 * inv_r6;

            float sr6 = powf(sigma, 6.0f) * inv_r6;
            float sr12 = sr6 * sr6;

            float force_mag = 24.0f * epsilon * (2.0f * sr12 - sr6) * inv_r2;

            return force_mag * dr_vec;
        }

        return glm::vec3{0.f, 0.f, 0.f};
    }

    void sim_dynamics::integrate()
    {
    
    }

    void sim_dynamics::step(float target_temp, float target_pressure)
    {
        auto &data = m_universe.getData();

        integrate();

        if (m_step_count % GRID_REBUILD == 0)
            universe_grid.rebuild(data.positions, m_universe.boxSizes(), CELL_CUTOFF);

        ++m_step_count;
    }
} // namespace sim
