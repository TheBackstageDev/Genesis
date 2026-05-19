#include "dynamics.hpp"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>

#include <iostream>

namespace sim
{
    sim_dynamics::sim_dynamics(fun::universe &u)
        : m_universe(u)
    {
        universe_verlet.cutoff = 12.f;
        universe_verlet.skin = VERLET_SKIN;

        createComputeShaders();
    }

    sim_dynamics::~sim_dynamics()
    {
        destroySSBOs();
    }

    void sim_dynamics::createComputeShaders()
    {
        const std::filesystem::path shaders_path = "resource/shaders";

        unbonded_program = core::glProgram{
            core::glShader{GL_COMPUTE_SHADER, shaders_path / "unbonded.comp"},
        };

        bonded_program = core::glProgram{
            core::glShader{GL_COMPUTE_SHADER, shaders_path / "bonded.comp"},
        };

        angles_program = core::glProgram{
            core::glShader{GL_COMPUTE_SHADER, shaders_path / "angle.comp"},
        };

        integrate_program = core::glProgram{
            core::glShader{GL_COMPUTE_SHADER, shaders_path / "integrate_program.comp"},
        };
    }

    void sim_dynamics::destroySSBOs()
    {
        if (!m_builtSSBOs) return;

        ssbo_pos.destroy();
        ssbo_lj.destroy(); 
        ssbo_q.destroy();
        ssbo_force.destroy();

        m_builtSSBOs = false;
    }

    void sim_dynamics::updateSSBOs()
    {
        auto &data = m_universe.getData();
        auto &atomData = m_universe.getAtomData();
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

        const size_t bondCount = atomData.bonds.size();
        const size_t angleCount = atomData.angles.size();
        
        if (m_lastBondCount != bondCount) 
        {
            m_lastBondCount = bondCount;

            std::vector<bufferBond> bonds;
            bonds.reserve(bondCount);

            std::vector<glm::uvec2> bondedPairs;
            bondedPairs.reserve(bondCount);

            for (const auto& bond : atomData.bonds)
            {
                uint32_t a = bond.centralAtom;
                uint32_t b = bond.bondedAtom;

                uint8_t zA = atomData.atoms[a].ZIndex;
                uint8_t zB = atomData.atoms[b].ZIndex;

                float k  = constants::getBondHarmonicConstantFromEnergy(zA, zB, bond.type);
                float r0 = constants::getBondLength(zA, zB, bond.type);

                if (m_universe.reactive())
                {
                    k  = constants::getBondHarmonicConstantFromEnergy(zA, zB, bond.order);
                    r0 = constants::getContinuousBondLength(zA, zB, bond.order);
                }

                bonds.emplace_back(b, a, k, r0);

                uint32_t i = std::min(a, b);
                uint32_t j = std::max(a, b);
                bondedPairs.emplace_back(i, j);
            }

            if (!ssbo_bonds.id() || 
                ssbo_bonds.getSize() != static_cast<GLsizeiptr>(bonds.size() * sizeof(bufferBond)))
            {
                ssbo_bonds = core::glBuffer(GL_SHADER_STORAGE_BUFFER,
                                            bonds.data(),
                                            bonds.size() * sizeof(bufferBond),
                                            GL_DYNAMIC_DRAW);
            }
            else
            {
                ssbo_bonds.update(bonds.data(), bonds.size() * sizeof(bufferBond));
            }

            if (!ssbo_bondedPairs.id() || 
                ssbo_bondedPairs.getSize() != static_cast<GLsizeiptr>(bondedPairs.size() * sizeof(glm::uvec2)))
            {
                ssbo_bondedPairs = core::glBuffer(GL_SHADER_STORAGE_BUFFER,
                                                bondedPairs.data(),
                                                bondedPairs.size() * sizeof(glm::uvec2),
                                                GL_DYNAMIC_DRAW);
            }
            else
            {
                ssbo_bondedPairs.update(bondedPairs.data(), bondedPairs.size() * sizeof(glm::uvec2));
            }
        }

        if (m_lastAngleCount != angleCount)
        {
            m_lastAngleCount = angleCount;

            std::vector<bufferAngle> angles;
            angles.reserve(angleCount);

            for (const auto& angle : atomData.angles)
            {
                uint32_t a = angle.A;
                uint32_t b = angle.B;
                uint32_t c = angle.C;

                uint8_t zA = atomData.atoms[a].ZIndex;
                uint8_t zB = atomData.atoms[b].ZIndex;
                uint8_t zC = atomData.atoms[c].ZIndex;

                angles.emplace_back(a, b, c, constants::getAngleHarmonicConstant(zA, zB, zC), angle.rad);
            }

            if (!ssbo_angles.id() || 
                ssbo_angles.getSize() != static_cast<GLsizeiptr>(angles.size() * sizeof(bufferAngle)))
            {
                ssbo_angles = core::glBuffer(GL_SHADER_STORAGE_BUFFER,
                                                angles.data(),
                                                angles.size() * sizeof(bufferAngle),
                                                GL_DYNAMIC_DRAW);
            }
            else
            {
                ssbo_angles.update(angles.data(), angles.size() * sizeof(bufferAngle));
            }
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

        m_builtSSBOs = true;
    }

    void checkGLError(const char* location) 
    {
        GLenum err = glGetError();
        if (err != GL_NO_ERROR) {
            std::cerr << "GL ERROR " << err << " at " << location << std::endl;
            if (err == GL_INVALID_VALUE) std::cerr << "  (likely bad buffer size, binding, or uniform)\n";
            if (err == GL_INVALID_OPERATION) std::cerr << "  (likely program not linked or bound wrong)\n";
        }
    }

    void sim_dynamics::computeUnbondedGPU()
    {
        auto &data = m_universe.getData();
        const uint32_t N = static_cast<uint32_t>(m_universe.numAtoms());
        if (N == 0) return;

        const uint32_t local_size = 256;
        const uint32_t num_groups = (N + local_size - 1) / local_size;

        unbonded_program.use(); checkGLError("use program");

        unbonded_program.setUniform("numParticles", N); checkGLError("set numParticles");
        unbonded_program.setUniform("numBondedPairs", (uint32_t)m_universe.numBonds()); checkGLError("set numBondedPairs");
        unbonded_program.setUniform("box", m_universe.boxSizes()); checkGLError("set box");
        unbonded_program.setUniform("r_cut2", CELL_CUTOFF * CELL_CUTOFF); checkGLError("set r_cut2");
        unbonded_program.setUniform("collide_floor_roof", m_universe.rooffloorcollision()); checkGLError("set collide_floor_roof");
        unbonded_program.setUniform("collide_walls", m_universe.wallcollision()); checkGLError("set collide_walls");

        ssbo_pos.bindBase(0);   checkGLError("bind pos");
        ssbo_lj.bindBase(1);    checkGLError("bind lj");
        ssbo_q.bindBase(2);     checkGLError("bind q");
        ssbo_force.bindBase(3); checkGLError("bind force");
        ssbo_bondedPairs.bindBase(4); checkGLError("bind bonded pairs");

        glDispatchCompute(num_groups, 1, 1); checkGLError("dispatch compute");
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT); checkGLError("memory barrier");

        core::glBuffer::unbind(GL_SHADER_STORAGE_BUFFER); checkGLError("unbind SSBO");
    }

    void sim_dynamics::computeBondedGPU()
    {
        auto &atomData = m_universe.getAtomData();
        const uint32_t numBonds = static_cast<uint32_t>(atomData.bonds.size());
        if (numBonds == 0) return;

        bonded_program.use(); checkGLError("bonded_program.use()");

        bonded_program.setUniform("numBonds", numBonds);
        bonded_program.setUniform("box", m_universe.boxSizes()); checkGLError("set box");
        bonded_program.setUniform("collide_floor_roof", m_universe.rooffloorcollision()); checkGLError("set collide_floor_roof");
        bonded_program.setUniform("collide_walls", m_universe.wallcollision()); checkGLError("set collide_walls");

        ssbo_pos.bindBase(0);
        ssbo_bonds.bindBase(2);
        ssbo_force.bindBase(3);

        uint32_t local_size = 256;
        uint32_t num_groups = (numBonds + local_size - 1) / local_size;

        glDispatchCompute(num_groups, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        core::glBuffer::unbind(GL_SHADER_STORAGE_BUFFER);
    }

    void sim_dynamics::computeAnglesGPU()
    {
        auto &atomData = m_universe.getAtomData();
        const uint32_t numAngles = static_cast<uint32_t>(atomData.angles.size());
        if (numAngles == 0) return;

        angles_program.use(); 
        checkGLError("angles_program.use()");

        angles_program.setUniform("numAngles", numAngles);
        angles_program.setUniform("box", m_universe.boxSizes());
        angles_program.setUniform("collide_floor_roof", m_universe.rooffloorcollision());
        angles_program.setUniform("collide_walls", m_universe.wallcollision());

        ssbo_pos.bindBase(0);
        ssbo_angles.bindBase(1);
        ssbo_force.bindBase(2);

        uint32_t local_size = 256;
        uint32_t num_groups = (numAngles + local_size - 1) / local_size;

        glDispatchCompute(num_groups, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        core::glBuffer::unbind(GL_SHADER_STORAGE_BUFFER);
    }

    std::vector<glm::vec3> sim_dynamics::processCellUnbonded(int32_t ix, int32_t iy, int32_t iz, int32_t atom_start, int32_t atom_end)
    {
        if (universe_verlet.verlet.empty()) return {};

        float local_virial = 0.f;
        auto &data = m_universe.getData();
        auto &atomData = m_universe.getAtomData();
        int32_t cellID = universe_grid.cellToIndex(ix, iy, iz);

        thread_local std::vector<glm::vec3> local_forces(atomData.atoms.size(), {0, 0, 0});
        if (local_forces.size() != atomData.atoms.size())
            local_forces.assign(atomData.atoms.size(), {0.0f, 0.0f, 0.0f});
        std::fill(local_forces.begin(), local_forces.end(), glm::vec3{0, 0, 0});

        universe_grid.foreach(cellID, [&](const uint32_t& i)
        {
            const glm::vec3& pi = data.positions[i];

            for (uint32_t j : universe_verlet.verlet[i])
            {
                if (j <= i) continue;

                glm::vec3 dr = m_universe.minImageVec(data.positions[j] - pi);
                float dr2 = glm::length2(dr);

                if (dr2 > CELL_CUTOFF * CELL_CUTOFF || dr2 < EPSILON)
                    continue;

                if (m_universe.areBonded(i, j))
                    continue;

                glm::vec3 cForce = computeCoulombForce(i, j, dr);
                glm::vec3 lForce = computeLJforce(i, j, dr);
                glm::vec3 total_force = cForce + lForce;

                local_forces[i] += total_force;
                local_forces[j] -= total_force;

                local_virial += glm::dot(dr, total_force);
            }
        }, atom_start, atom_end);

        total_virial.fetch_add(local_virial, std::memory_order_relaxed);
        return local_forces;
    }

    void sim_dynamics::computeUnbondedCPU()
    {
        auto &data = m_universe.getData();

        const size_t cells = universe_grid.cellOffsets.size();
        const int32_t n_threads = std::thread::hardware_concurrency();
        //const int32_t n_threads = 1;
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

                for (size_t i = 0; i < cell_forces.size(); ++i)
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
                data.forces[i] += glm::vec4(local_f[i], 0.0f);
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
                if (len <= EPSILON || len > 3.f)
                    continue;

                uint8_t aZindex = atomData.atoms[a].ZIndex, cZindex = atomData.atoms[c].ZIndex;

                float bondLength = !m_universe.reactive() ? constants::getBondLength(aZindex, cZindex, b.type) : constants::getContinuousBondLength(aZindex, cZindex, b.order);
                float dl = len - bondLength;
                
                float bond_k = !m_universe.reactive() ? constants::getBondHarmonicConstantFromEnergy(aZindex, cZindex, b.type) : constants::getBondHarmonicConstantFromEnergy(aZindex, cZindex, b.order);
                glm::vec3 force = (bond_k * dl / len) * dr;

                atomData.atoms[a].total_BO += b.order;
                atomData.atoms[c].total_BO += b.order;

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

                float K = constants::getAngleHarmonicConstant(ang.A, ang.B, ang.C);

                /* if (m_universe.reactive())
                {
                    const sim::fun::bond& AB = atomData.bonds[m_universe.getBond(i, j)];
                    const sim::fun::bond& BC = atomData.bonds[m_universe.getBond(k, j)];

                    K = constants::getContinuousAngleConstant(atomData.atoms[i].ZIndex, atomData.atoms[j].ZIndex, atomData.atoms[k].ZIndex, AB.order, BC.order);
                } */
                
                glm::vec3 F_i = -K * delta_theta * dtheta_dri;
                glm::vec3 F_k = -K * delta_theta * dtheta_drk;
                
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
                float target = 0.0f;
                float K = 6.0f;
                if (da.periodicity == 3) 
                {
                    float phi = computeDihedral(pa, pb, pc, pd);
                    K          = 1.5f;
                    target = round(phi / (2.0f * M_PI / 3.0f)) * (2.0f * M_PI / 3.0f);
                }
                else if (da.periodicity == 2) 
                {
                    K          = 6.0f;
                    target = (phi > M_PI/2.0f) ? M_PI : 0.0f;
                }

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
                data.forces[i] += glm::vec4(local_f[i], 0.0f);
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
        auto &atomData = m_universe.getAtomData();

        const uint32_t base_i = i << 1;
        const uint32_t base_j = j << 1;

        const float sigma_i = data.lj_params[base_i];
        const float epsilon_i = data.lj_params[base_i + 1];

        const float sigma_j = data.lj_params[base_j];
        const float epsilon_j = data.lj_params[base_j + 1];

        float dr2 = glm::length2(dr_vec);

        float sigma = (sigma_i + sigma_j) * 0.5f;
        float epsilon = sqrtf(epsilon_i * epsilon_j);

        if (dr2 < 100.f && dr2 > EPSILON)
        {
            float inv_r = 1.f / sqrtf(dr2);
            float inv_r2 = 1.f / dr2;
            float inv_r6 = inv_r2 * inv_r2 * inv_r2;
            float inv_r12 = inv_r6 * inv_r6;

            float force_mag = 0.f;

            float sum_cov = constants::covalent_radius[atomData.atoms[i].ZIndex] +
                            constants::covalent_radius[atomData.atoms[j].ZIndex];

            float sr6 = powf(sigma, 6.0f) * inv_r6;
            float sr12 = sr6 * sr6;

            force_mag += 24.0f * epsilon * (2.0f * sr12 - sr6) * inv_r2;
            
            if (m_universe.reactive() && dr2 < sum_cov * sum_cov)
                force_mag *= 1.f / (1e-6f + std::exp(5.f * dr2));

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
            std::fill(data.forces.begin(), data.forces.end(), glm::vec4{0.0f});
        }
    }

    void sim_dynamics::computeExternalForces(uint32_t i)
    {
        auto& wall_charges = m_universe.getWallCharges();
        auto& data = m_universe.getData();
        auto& atomData = m_universe.getAtomData();

        if (m_universe.wallChargeEnabled())
        for (int32_t w = 0; w < m_universe.getWallCharges().size(); ++w)
        {
            const float &wall_q = wall_charges[w];
            if (wall_q == 0.f)
                continue;

            glm::vec3 force = data.q[i] * wall_q * wall_directions[w];
            data.forces[i] += glm::vec4(force, 0.0f);
        }

        if (m_universe.magneticFieldEnabled())
        {
            glm::vec3 force = data.q[i] * glm::cross(data.velocities[i], m_universe.getMagneticFieldStrength());
            data.forces[i] += glm::vec4(force, 0.0f);
        }
    }

    void sim_dynamics::computeForces()
    {
        zeroForces();

        if (m_GPU)
        {
            updateSSBOs();
            computeBondedGPU();
            computeAnglesGPU();
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

        if (N == 0) return;

        const float effective_dt = m_dt * m_timescale;
        const float half_dt = 0.5f * effective_dt;

        computeForces();

        for (int32_t i = 0; i < N; ++i)
        {
            // Resets for reaction system
            atomData.atoms[i].total_BO = 0.f;

            if (atomData.frozen_atoms[i])
            {
                data.velocities[i] = glm::vec3{0.0f};
                continue;
            }

            computeExternalForces(i);

            const float inv_m = 1.0f / atomData.atoms[i].mass;
            const glm::vec3 accel = glm::vec3(data.forces[i]) * inv_m;

            data.velocities[i] += accel * half_dt;
        }

        for (int32_t i = 0; i < N; ++i)
        {
            if (atomData.frozen_atoms[i]) continue;

            data.positions[i] += data.velocities[i] * effective_dt;
            univ.boundCheck(static_cast<uint32_t>(i));
        }

        computeForces();

        if (m_GPU)
        {
            ssbo_force.bind(); checkGLError("bind force for readback");
            glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, N * sizeof(glm::vec4), data.forces.data());
            checkGLError("get buffer sub data");
            ssbo_force.unbind(); checkGLError("unbind after readback");
        }

        glm::vec3 grav_accel(0.0f, 0.0f, -m_universe.gravityMagnitude());

        for (int32_t i = 0; i < N; ++i)
        {
            if (atomData.frozen_atoms[i])
            {
                data.velocities[i] = glm::vec3{0.0f};
                continue;
            }

            computeExternalForces(i);

            const float inv_m = 1.0f / atomData.atoms[i].mass;
            const glm::vec3 accel = glm::vec3(data.forces[i]) * inv_m;

            data.velocities[i] += accel * half_dt;

            if (m_universe.gravityEnabled())
            {
                if (atomData.frozen_atoms[i]) continue;
                data.velocities[i] += grav_accel * 1e-3f * effective_dt;
            }
        }
    }

    void sim_dynamics::step(float target_temp, float target_pressure)
    {
        if (m_paused)
            return;

        auto &data = m_universe.getData();
        auto &atomData = m_universe.getAtomData();

        if (universe_verlet.verlet.size() != m_universe.numAtoms())
        {
            universe_grid.rebuild(data.positions, m_universe.boxSizes(), CELL_CUTOFF);
            universe_verlet.construct(universe_grid, m_universe);        
        }

        integrate();

        if (m_step_count % GRID_REBUILD == 0)
            universe_grid.rebuild(data.positions, m_universe.boxSizes(), CELL_CUTOFF);
        
        float verletRebuildRate = 20.f / (m_timescale + std::numeric_limits<float>::epsilon());

        if (m_step_count % std::max(1u, static_cast<uint32_t>(verletRebuildRate)) == 0)
            universe_verlet.construct(universe_grid, m_universe);

        setTemperature(target_temp);
        setPressure(target_pressure);
        COMDrift();

        m_step_count++;
        m_accumulated_time += m_dt * m_timescale;

        /* float avg_speed = 0.0f;
            for (const auto& v : data.velocities) avg_speed += glm::length(v);
            avg_speed /= atomData.atoms.size();
            std::cout << "Average molecular speed: " << avg_speed << " A/Ps\n"; */
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

    void sim_dynamics::setPressure(float target_p_kpa)
    {
        if (target_p_kpa <= 0.0f)
            return;

        if (m_step_count % BAROSTAT_INTERVAL != 0)
            return;

        glm::vec3 box = m_universe.boxSizes();
        auto &data = m_universe.getData();
        auto &atomData = m_universe.getAtomData();

        constexpr float beta_T = 4.5e-5f;
        constexpr float tau_P = 1e2;

        m_pressure = computePressure();

        float delta_P = target_p_kpa - m_pressure;
        float mu = 1.0f - (2.f / tau_P) * beta_T * delta_P;

        mu = std::clamp(mu, 0.5f, 1.5f);

        float scale = std::cbrt(mu);

        box.z *= scale;

        for (int32_t i = 0; i < atomData.atoms.size(); ++i)
        {
            data.positions[i].z *= scale;
            m_universe.setBoxSize(glm::vec3(m_universe.boxSizes().x, m_universe.boxSizes().y, m_universe.boxSizes().z * scale));
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
        if (m_step_count % 100000 != 0)
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
