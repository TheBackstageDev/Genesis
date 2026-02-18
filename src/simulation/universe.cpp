#include "universe.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <set>
#include <queue>
#include <fstream>

#include "constants.hpp"

namespace sim
{
    namespace fun
    {

        universe::universe(const universe_create_info &create_info, rendering_engine &rendering_engine)
            : box(create_info.box), gravity(create_info.has_gravity), mag_gravity(create_info.mag_gravity),
              wall_collision(create_info.wall_collision), isothermal(create_info.isothermal),
              HMassRepartitioning(create_info.HMassRepartitioning), roof_floor_collision(create_info.roof_floor_collision),
              log_flags(create_info.log_flags), rendering_eng(rendering_engine)
        {
            createComputeShaders();
        }

        universe::universe(const std::filesystem::path path, rendering_engine &rendering_engine)
            : rendering_eng(rendering_engine)
        {
            loadScene(path);
            createComputeShaders();
        }

        // Compute Shaders

        void universe::createComputeShader(const std::filesystem::path shader)
        {
            
        }

        void universe::createComputeShaders()
        {
            const std::filesystem::path shaders_path = "resource/shaders";
            createComputeShader(shaders_path / "bonded.comp");
            createComputeShader(shaders_path / "unbonded.comp");
            createComputeShader(shaders_path / "simulation.comp");
        }

        void universe::draw(sf::RenderTarget &target, rendering_info info)
        {
            rendering_simulation_info sim_info{.positions = m_displayPositions.empty() ? data.positions : m_displayPositions,
                                               .velocities = data.velocities,
                                               .q = data.q,
                                               .atoms = atoms,
                                               .bonds = bonds,
                                               .molecules = molecules,
                                               .box = box};

            if (info.flag_highlights)
            {
                info.highlight_indices = m_highlightedAtoms;
                info.highlight_bonds = m_highlightedBonds;
            }

            if (info.flag_arrows)
            {
                info.arrows = m_Arrows;
            }

            rendering_eng.draw(target, info, sim_info);
        }

        int32_t universe::createAtom(glm::vec3 p, glm::vec3 v, uint8_t ZIndex, uint8_t numNeutrons, uint8_t numElectron, int32_t chirality)
        {
            atom newAtom{};
            newAtom.ZIndex = ZIndex;

            emplace_vel(v);
            emplace_pos(p);

            std::pair<float, float> constants = constants::getAtomConstants(ZIndex);

            data.lj_params.emplace_back(constants.first);
            data.lj_params.emplace_back(constants.second);
            newAtom.radius = constants::VDW_RADII[ZIndex] * 1.5f;
            newAtom.electrons = numElectron;
            newAtom.NCount = numNeutrons == 0 ? constants::NEUTRON_COUNTS[ZIndex] : numNeutrons;
            newAtom.mass = ZIndex * MASS_PROTON + newAtom.NCount * MASS_NEUTRON + numElectron * MASS_ELECTRON;
            newAtom.chirality = chirality;
            newAtom.bondCount = 0;

            if (HMassRepartitioning)
            {
                if (ZIndex == 1)
                    newAtom.mass *= 3.f;
                if (ZIndex > 5)
                    newAtom.mass -= 3.f;
            }

            atoms.emplace_back(std::move(newAtom));
            data.q.emplace_back(ZIndex - numElectron);
            data.forces.resize(atoms.size());

            frozen_atoms.emplace_back(false);

            return atoms.size() - 1;
        }

        void universe::createBond(int32_t idx1, int32_t idx2, BondType type)
        {
            if (idx1 >= atoms.size() || idx2 >= atoms.size() || idx1 == idx2)
                return;

            if (std::find_if(bonds.begin(), bonds.end(), [&](bond &a)
                             { return a.bondedAtom == idx1 && a.centralAtom == idx2 || a.bondedAtom == idx2 && a.centralAtom == idx1; }) != bonds.end())
                return;

            bond nBond{};
            nBond.bondedAtom = idx1;
            nBond.centralAtom = idx2;
            nBond.type = type;
            nBond.equilibriumLength = constants::getBondLength(atoms[idx1].ZIndex, atoms[idx2].ZIndex, type);
            nBond.k = constants::getBondHarmonicConstantFromEnergy(atoms[idx1].ZIndex, atoms[idx2].ZIndex, type);

            int8_t bondCount = static_cast<int8_t>(type);

            atoms[idx1].bondCount += bondCount;
            atoms[idx2].bondCount += bondCount;

            float EN1 = constants::getElectronegativity(atoms[idx1].ZIndex);
            float EN2 = constants::getElectronegativity(atoms[idx2].ZIndex);
            float deltaEN = std::abs(EN1 - EN2);

            if (deltaEN > 0.1f) // Significant electronegativity difference
            {
                float charge = deltaEN * 0.31f;

                if (EN1 > EN2)
                {
                    data.q[idx1] -= charge;
                    data.q[idx2] += charge;
                }
                else
                {
                    data.q[idx2] -= charge;
                    data.q[idx1] += charge;
                }
            }

            bonds.emplace_back(std::move(nBond));
        }

        int32_t universe::createSubset(const def_subset &nSub, const int32_t baseAtom, const int32_t baseSubset)
        {
            subset nSubset{};
            nSubset.mainAtomIdx = nSub.mainAtomIdx + baseAtom;
            nSubset.bondedSubsetIdx = nSub.bondedSubset == SIZE_MAX ? SIZE_MAX : nSub.bondedSubset + baseSubset;
            nSubset.bondingSubsetIdx = nSub.bondingSubset == SIZE_MAX ? SIZE_MAX : nSub.bondingSubset + baseSubset;

            std::vector<uint8_t> neighbourZs;
            neighbourZs.reserve(nSub.connectedIdx.size() + nSub.hydrogensIdx.size());

            for (int32_t i = 0; i < nSub.connectedIdx.size(); ++i)
            {
                const int32_t bondedAtom = nSub.connectedIdx[i] + baseAtom;
                neighbourZs.emplace_back(atoms[bondedAtom].ZIndex);
            }

            for (int32_t h = 0; h < nSub.hydrogensIdx.size(); ++h)
            {
                neighbourZs.emplace_back(1);
            }

            if (nSub.hydrogensIdx.size() > 0)
            {
                nSubset.hydrogenBegin = nSub.hydrogensIdx[0] + baseAtom;
                nSubset.hydrogenCount = nSub.hydrogensIdx.size();
            }
            if (nSub.connectedIdx.size() > 0)
            {
                nSubset.connectedBegin = nSub.connectedIdx[0] + baseAtom;
                nSubset.connectedCount = nSub.connectedIdx.size();
            }

            subsets.emplace_back(std::move(nSubset));
            return subsets.size() - 1; // Index
        }

        void universe::createMolecule(molecule_structure structure, sf::Vector3f pos, sf::Vector3f vel)
        {
            int32_t baseAtomIndex = atoms.size();

            molecule nMolecule{};

            sf::Vector3f centroid{0.f, 0.f, 0.f};
            for (const auto &p : structure.positions)
            {
                centroid += p;
            }
            centroid /= static_cast<float>(structure.positions.size());

            for (auto &p : structure.positions)
            {
                p -= centroid;
            }

            for (int32_t i = 0; i < structure.atoms.size(); ++i)
            {
                const def_atom &a = structure.atoms[i];

                sf::Vector3f end_pos = structure.positions[i] + pos;
                createAtom(glm::vec3(end_pos.x, end_pos.y, end_pos.z), glm::vec3(vel.x, vel.y, vel.z), a.ZIndex, a.NIndex, a.ZIndex - a.charge, a.chirality);
                data.q[i] += structure.atoms[i].charge;
            }

            nMolecule.atomBegin = baseAtomIndex;
            nMolecule.atomCount = structure.atoms.size();

            int32_t baseBondIndex = bonds.size();
            for (int32_t b = 0; b < structure.bonds.size(); ++b)
            {
                const def_bond &db = structure.bonds[b];
                int32_t central = baseAtomIndex + db.centralAtomIdx;
                int32_t bonded = baseAtomIndex + db.bondingAtomIdx;

                createBond(bonded, central, db.type);
            }

            nMolecule.bondBegin = baseBondIndex;
            nMolecule.bondCount = structure.bonds.size();

            int32_t baseSubset = subsets.size();
            for (int32_t s = 0; s < structure.subsets.size(); ++s)
                createSubset(structure.subsets[s], baseAtomIndex, baseSubset);

            nMolecule.subsetBegin = baseSubset;
            nMolecule.subsetCount = structure.subsets.size();

            if (nMolecule.subsetCount > 0 && structure.atoms[0].charge == 0)
                balanceMolecularCharges(subsets[baseSubset]);

            rebuildBondTopology();

            int32_t baseAngle = angles.size();
            for (int32_t a = 0; a < structure.angles.size(); ++a)
            {
                angle angle = structure.angles[a];
                angle.A += baseAtomIndex;
                angle.B += baseAtomIndex;
                angle.C += baseAtomIndex;

                angles.emplace_back(std::move(angle));
            }

            nMolecule.angleBegin = baseAngle;
            nMolecule.angleCount = structure.angles.size();

            int32_t baseDihedral = dihedral_angles.size();
            for (int32_t a = 0; a < structure.dihedral_angles.size(); ++a)
            {
                dihedral_angle angle = structure.dihedral_angles[a];
                angle.A += baseAtomIndex;
                angle.B += baseAtomIndex;
                angle.C += baseAtomIndex;
                angle.D += baseAtomIndex;

                dihedral_angles.emplace_back(std::move(angle));
            }

            int32_t baseImproper = improper_angles.size();
            for (int32_t a = 0; a < structure.improper_angles.size(); ++a)
            {
                dihedral_angle angle = structure.improper_angles[a];
                angle.A += baseAtomIndex;
                angle.B += baseAtomIndex;
                angle.C += baseAtomIndex;
                angle.D += baseAtomIndex;

                improper_angles.emplace_back(std::move(angle));
            }

            nMolecule.dihedralBegin = baseDihedral;
            nMolecule.dihedralCount = structure.dihedral_angles.size();

            molecules.emplace_back(std::move(nMolecule));
        }

        void universe::balanceMolecularCharges(subset &mol)
        {
            auto it = std::find_if(subsets.begin(), subsets.end(), [&](const subset &s)
                                   { return &s == &mol; });
            if (it == subsets.end())
                return;

            int32_t startIdx = std::distance(subsets.begin(), it);

            std::vector<int32_t> atomsToBalance;
            std::unordered_set<int32_t> visited;
            std::queue<int32_t> q;
            q.push(startIdx);
            visited.insert(startIdx);

            while (!q.empty())
            {
                int32_t currentIdx = q.front();
                q.pop();
                const auto &currentSubset = subsets[currentIdx];

                atomsToBalance.push_back(currentSubset.mainAtomIdx);

                if (currentSubset.connectedCount != UINT32_MAX)
                    for (int32_t i = 0; i < currentSubset.connectedCount; ++i)
                    {
                        atomsToBalance.emplace_back(currentSubset.connectedBegin + i);
                    }

                if (currentSubset.hydrogenCount != UINT32_MAX)
                    for (int32_t i = 0; i < currentSubset.hydrogenCount; ++i)
                    {
                        atomsToBalance.emplace_back(currentSubset.hydrogenCount + i);
                    }

                int32_t nextIdx = currentSubset.bondedSubsetIdx;
                if (nextIdx < subsets.size() && nextIdx != SIZE_MAX && visited.insert(nextIdx).second)
                {
                    q.push(nextIdx);
                }
            }

            if (atomsToBalance.empty())
                return;

            std::sort(atomsToBalance.begin(), atomsToBalance.end());
            auto last = std::unique(atomsToBalance.begin(), atomsToBalance.end());
            atomsToBalance.erase(last, atomsToBalance.end());

            std::vector<std::vector<int32_t>> neighborLists(atoms.size());
            for (const auto &bond : bonds)
            {
                int32_t a = bond.bondedAtom;
                int32_t b = bond.centralAtom;
                bool aIn = std::binary_search(atomsToBalance.begin(), atomsToBalance.end(), a);
                bool bIn = std::binary_search(atomsToBalance.begin(), atomsToBalance.end(), b);
                if (aIn && bIn)
                {
                    neighborLists[a].push_back(b);
                    neighborLists[b].push_back(a);
                }
            }

            float totalCharge = std::accumulate(data.q.begin(), data.q.end(), 0.f);

            if (std::abs(totalCharge) < 1e-6f)
                return;

            std::vector<float> valenceDeficit(atomsToBalance.size(), 0.0f);
            for (int32_t i = 0; i < atomsToBalance.size(); ++i)
            {
                int32_t idx = atomsToBalance[i];
                if (idx >= atoms.size())
                    continue;

                uint8_t usualValence = constants::getUsualBonds(atoms[idx].ZIndex);
                int deficit = usualValence - static_cast<int>(atoms[idx].bondCount);
                valenceDeficit[i] = static_cast<float>(deficit);
            }

            float chargeToDistribute = -totalCharge;
            float totalDeficit = 0.0f;
            for (float d : valenceDeficit)
            {
                if (d > 0)
                    totalDeficit += d;
            }

            if (totalDeficit > 0.5f)
            {
                for (int32_t i = 0; i < atomsToBalance.size(); ++i)
                {
                    if (valenceDeficit[i] > 0)
                    {
                        float share = chargeToDistribute * (valenceDeficit[i] / totalDeficit);
                        data.q[atomsToBalance[i]] += share;
                    }
                }
            }
            else
            {
                float adjustment = chargeToDistribute / atomsToBalance.size();
                for (int32_t idx : atomsToBalance)
                {
                    data.q[idx] += adjustment;
                }
            }
        }

        void universe::boundCheck(uint32_t i)
        {
            float &x = data.positions[i].x;
            float &y = data.positions[i].y;
            float &z = data.positions[i].z;

            float &vx = data.velocities[i].x;
            float &vy = data.velocities[i].y;
            float &vz = data.velocities[i].z;

            if (roof_floor_collision)
            {
                if (z < 0.0f)
                {
                    z = 0.0f;
                    vz = -vz;
                }
                else if (z > box.z)
                {
                    z = box.z;
                    vz = -vz;
                }
                vz *= 0.99f;
            }
            else
            {
                z = std::fmod(std::fmod(z, box.z) + box.z, box.z);
            }

            if (wall_collision)
            {
                if (x < 0.0f)
                {
                    x = 0.0f;
                    vx = -vx;
                }
                else if (x > box.x)
                {
                    x = box.x;
                    vx = -vx;
                }

                if (y < 0.0f)
                {
                    y = 0.0f;
                    vy = -vy;
                }
                else if (y > box.y)
                {
                    y = box.y;
                    vy = -vy;
                }

                vx *= 0.99f;
                vy *= 0.99f;
            }
            else
            {
                x = std::fmod(std::fmod(x, box.x) + box.x, box.x);
                y = std::fmod(std::fmod(y, box.y) + box.y, box.y);
            }
        }

        float universe::ljPot(uint32_t i, uint32_t j)
        {
            float potential = 0.f;

            float dr = minImageVec(pos(i) - pos(j)).length();

            auto [sigma_i, epsilon_i] = constants::getAtomConstants(atoms[i].ZIndex);
            auto [sigma_j, epsilon_j] = constants::getAtomConstants(atoms[j].ZIndex);
            float sigma = (sigma_i + sigma_j) / 2.0f;

            if (dr < sigma * CUTOFF && dr > EPSILON)
            {
                float epsilon = sqrtf(epsilon_i * epsilon_j);
                float r6 = powf((sigma / dr), 6);
                float r12 = r6 * r6;
                potential += 4.0f * epsilon * (r12 - r6);
            }

            return potential;
        }

        glm::vec3 universe::ljForce(uint32_t i, uint32_t j)
        {
            const uint32_t base_i = i << 1;
            const uint32_t base_j = j << 1;

            const float sigma_i = data.lj_params[base_i];
            const float sigma_j = data.lj_params[base_j];
            const float epsilon_i = data.lj_params[base_i + 1];
            const float epsilon_j = data.lj_params[base_j + 1];

            sf::Vector3f dr_vec = minImageVec(pos(i) - pos(j));
            float dr = dr_vec.length();

            const float sigma = sqrtf(sigma_i * sigma_j);
            const float epsilon = sqrtf(epsilon_i * epsilon_j);

            if (dr < sigma * CUTOFF && dr > EPSILON)
            {
                float inv_r2 = 1.0f / (dr * dr);
                float inv_r6 = inv_r2 * inv_r2 * inv_r2;
                float inv_r12 = inv_r6 * inv_r6;

                float sigma6 = powf(sigma, 6);
                float sigma12 = sigma6 * sigma6;

                float force_mag = 24.0f * epsilon * inv_r2 * (2.0f * sigma12 * inv_r12 - sigma6 * inv_r6);

                return force_mag * glm::vec3(dr_vec.x, dr_vec.y, dr_vec.z);
            }

            return glm::vec3{0.f, 0.f, 0.f};
        }

        // Wolf method
        float universe::wolfForce(float r, float qi_qj)
        {
            constexpr float alpha = 0.25f; // damping
            constexpr float rc = COULOMB_CUTOFF;
            if (r >= rc)
                return 0.f;

            float ir = 1.f / r;
            float ir_c = 1.f / rc;

            float erfc_ar = std::erfc(alpha * r);
            float exp_term = std::exp(-(alpha * alpha * r * r));

            float force = COULOMB_K * qi_qj * ir * ir * (erfc_ar + 2.f * alpha / std::sqrt(M_PI) * r * exp_term);

            return force;
        }

        glm::vec3 universe::coulombForce(uint32_t i, uint32_t j, glm::vec3 &dr_vec)
        {
            float dr = glm::length(dr_vec);

            if (dr < EPSILON || dr > COULOMB_CUTOFF)
                return {0.f, 0.f, 0.f};

            float qq = data.q[i] * data.q[j];
            if (qq == 0.f)
                return {0.f, 0.f, 0.f};

            float forceMag = COULOMB_K * qq / (dr * dr * dr);
            return forceMag * -dr_vec;
        }

        float universe::calculateDihedral(const glm::vec3 &pa, const glm::vec3 &pb, const glm::vec3 &pc, const glm::vec3 &pd)
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

        std::vector<glm::vec3> universe::processCellUnbonded(int32_t ix, int32_t iy, int32_t iz, int32_t atom_start, int32_t atom_end)
        {
            float local_virial = 0.f;

            int32_t cellID = universe_grid.cellToIndex(ix, iy, iz);

            thread_local std::vector<glm::vec3> local_forces(atoms.size(), {0, 0, 0});

            if (local_forces.size() != atoms.size())
                local_forces.assign(atoms.size(), {0.0f, 0.0f, 0.0f});

            std::fill(local_forces.begin(), local_forces.end(), glm::vec3{0, 0, 0});

            universe_grid.foreach(cellID, [&](const uint32_t& i)
            {
                universe_grid.foreach(cellID, [&](const uint32_t& j)
                {
                    if (j <= i)
                        return;

                    sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                    glm::vec3 dr_g(dr.x, dr.y, dr.z);

                    if (areBonded(i, j))
                        return;

                    glm::vec3 cForce = coulombForce(i, j, dr_g);
                    glm::vec3 lForce = ljForce(i, j);

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

                            if (!wall_collision && !roof_floor_collision && neighbor_id < -1) continue;

                            universe_grid.foreach(neighbor_id, [&](const uint32_t& j)
                            {
                                if (j <= i)
                                    return;

                                sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                                glm::vec3 dr_g(dr.x, dr.y, dr.z);

                                if (areBonded(i, j))
                                    return;

                                glm::vec3 cForce = coulombForce(i, j, dr_g);
                                glm::vec3 lForce = ljForce(i, j);

                                glm::vec3 total_force = cForce + lForce;
                                local_forces[i] += total_force;
                                local_forces[j] -= total_force;

                                local_virial += dr.x * total_force.x +
                                                dr.y * total_force.y +
                                                dr.z * total_force.z;
                            }, atom_start, atom_end);
                        }
                    }
            }, atom_start, atom_end);

            total_virial.fetch_add(local_virial, std::memory_order_relaxed);
            return local_forces;
        }

        void universe::calcBondedForcesParallel()
        {
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
            }
        }

        void universe::calcUnbondedForcesParallel()
        {
            const size_t cells = universe_grid.cellOffsets.size();
            const int32_t n_threads = std::thread::hardware_concurrency();
            //const int32_t n_threads = 1;
            const int32_t threads_per_top = cells < 9 ? n_threads : std::max(1u, static_cast<uint32_t>(std::floor(static_cast<double>(n_threads / 2)))); // how many threads will run on cells with lots of work
            constexpr int32_t subdivide_top = 4;                                                                           // how many of the top cells to subdivide

            struct task
            {
                uint32_t work;
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
                std::vector<glm::vec3> thread_forces(atoms.size(), {0, 0, 0});

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

                    for (size_t i = 0; i < atoms.size(); ++i)
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
                
                for (size_t i = 0; i < atoms.size(); ++i)
                    data.forces[i] += local_f[i];
            }
        }

        void universe::update(float targetTemperature, float targetPressure)
        {
            int32_t N = atoms.size();

            if (N == 0 || m_paused)
                return;

            data.forces.assign(N, glm::vec3(0.0f));
            total_virial = 0.0f;

            calcBondedForcesParallel();
            calcUnbondedForcesParallel();

            setPressure(targetPressure);
            setTemperature(targetTemperature);

            float dt = FEMTOSECOND * m_Timescale;

            for (int32_t i = 0; i < N; ++i)
            {
                for (int32_t w = 0; w < wall_charges.size() && f_wallCharges; ++w)
                {
                    const float& wall_q = wall_charges[w];
                    if (wall_q == 0.f) continue;
         
                    glm::vec3 force = data.q[i] * wall_q * wall_directions[w];
                    data.forces[i] += force;
                }

                if (f_magneticField)
                {
                    glm::vec3 force = data.q[i] * glm::cross(data.velocities[i], magnetic_strength);
                    data.forces[i] += force;
                }
                            
                sf::Vector3f accel = force(i) / atoms[i].mass;
                add_vel(i, accel * (0.5f * dt));

                if (gravity)
                    add_vel(i, sf::Vector3f(0.f, 0.f, -mag_gravity * 0.5f * dt));
            }

            for (int32_t i = 0; i < N; ++i)
            {
                sf::Vector3f dPos = vel(i) * dt;
                add_pos(i, dPos);
                boundCheck(i);
            }

            data.forces.assign(N, glm::vec3(0.0f));

            calcBondedForcesParallel();
            calcUnbondedForcesParallel();

            for (int32_t i = 0; i < N; ++i)
            {
                for (int32_t w = 0; w < wall_charges.size() && f_wallCharges; ++w)
                {
                    const float& wall_q = wall_charges[w];
                    if (wall_q == 0.f) continue;
         
                    glm::vec3 force = data.q[i] * wall_q * wall_directions[w];
                    data.forces[i] += force;
                }

                if (f_magneticField)
                {
                    glm::vec3 force = data.q[i] * glm::cross(data.velocities[i], magnetic_strength);
                    data.forces[i] += force;
                }

                sf::Vector3f accel = force(i) / atoms[i].mass;
                add_vel(i, accel * (0.5f * dt));

                if (gravity)
                    add_vel(i, sf::Vector3f(0.f, 0.f, -mag_gravity * 0.5f * dt));
            }

            if (timeStep % GRID_REBUILD == 0)
                universe_grid.rebuild(data.positions, box, CELL_CUTOFF);

            ++timeStep;
            m_accumulatedTime += dt;

            if (timeStep % 1000 == 0)
                COMDrift(); // fixes simulation box COM drift from numerical errors
        }

        float universe::calculatePressure()
        {
            if (atoms.empty())
                return 0.f;

            float kinetic = calculateKineticEnergy();
            float volume = box.x * box.y * box.z;
            float temperature = (2.0f / 3.0f) * kinetic / atoms.size();

            float P_ideal = atoms.size() * temperature / volume;

            float P_virial = -total_virial / (3.0f * volume);

            float pressure = P_ideal + P_virial;
            return pressure;
        }

        void universe::setPressure(float Target_P_Bar)
        {
            if (Target_P_Bar <= 0.0f)
                return;
            if (timeStep % BAROSTAT_INTERVAL != 0)
                return;

            constexpr float beta_T = 4.5e-5f;
            constexpr float tau_P = 1.0f;

            pres = calculatePressure();

            float delta_P = Target_P_Bar - pres;
            float mu = 1.0f - (2.f / tau_P) * beta_T * delta_P;

            mu = std::clamp(mu, 0.5f, 1.5f);

            float scale = std::cbrt(mu);

            box.z *= scale;

            for (int32_t i = 0; i < atoms.size(); ++i)
            {
                data.positions[i].z *= scale;
            }
        }

        // Bussi–Donadio–Parrinello (CSVR) velocity rescaling
        void universe::setTemperature(float kelvin)
        {
            if (timeStep % THERMOSTAT_INTERVAL != 0)
                return;

            float d = 3.0f * atoms.size() - 3.0f;
            if (d <= 0)
                d = 1;

            float KE = calculateKineticEnergy();
            float current_temp = (2.0f * KE) / (d * KB);
            temp = current_temp;

            if (!isothermal) return;
            
            float target_KE = 0.5f * d * KB * kelvin;
            float c = target_KE / (KE + std::numeric_limits<float>::epsilon());

            float r = gauss_random();
            float chi = r * sqrtf(2.0f / d);

            float alpha = sqrtf(c * (1.0f + chi + 0.5f * chi * chi));

            // alpha = sqrtf((c + sigma * r)^2 / 2 + c * (1 - c) * chi²(d-1) / d)

            for (auto &v : data.velocities)
                v *= alpha;
        }

        // Energy Calculation
        float universe::calculateKineticEnergy()
        {
            float kinetic_energy = 0.0f;
            for (int32_t i = 0; i < atoms.size(); ++i)
                kinetic_energy += 0.5f * atoms[i].mass * vel(i).lengthSquared();

            return kinetic_energy;
        }

        float universe::calculateAtomTemperature(int32_t i)
        {
            float ke = 0.5f * atoms[i].mass * vel(i).lengthSquared();
            return (2.0f / 3.0f) * ke * KB;
        }

        void universe::COMDrift()
        {
            glm::vec3 totalMomentum{0.f};
            float totalMass = 0.f;

            for (int32_t i = 0; i < atoms.size(); ++i)
            {
                totalMomentum += atoms[i].mass * data.velocities[i];
                totalMass += atoms[i].mass;
            }

            glm::vec3 correction = totalMomentum / totalMass;

            for (auto &v : data.velocities)
            {
                v -= correction;
            }
        }

        // TO DO
        std::string universe::moleculeName(const std::vector<uint32_t> &subsetIdx)
        {
            std::map<uint8_t, int32_t> ZIndices;

            for (int32_t i = 0; i < subsetIdx.size(); ++i)
            {
                auto &s = subsets[subsetIdx[i]];
                int32_t centralAtom = s.mainAtomIdx;
                ++ZIndices[atoms[centralAtom].ZIndex];
                ZIndices[1] += s.connectedCount;
            }

            std::unordered_map<int32_t, std::string> prefixes{
                {1, "mono"},
                {2, "di"},
                {3, "tri"},
                {4, "tetra"},
                {5, "penta"},
                {6, "hexa"},
                {7, "hepta"},
                {8, "octa"},
                {9, "nona"},
                {10, "deca"}};

            std::string name("");

            bool organic = ZIndices[6] > 0;

            if (!organic)
            {
                for (int32_t i = 1; i <= COUNT_ATOMS; ++i)
                {
                    if (ZIndices[i] <= 0)
                        continue;
                    std::string name = constants::getAtomName(i);
                    std::transform(name.begin(), name.end(), name.begin(), [](uint8_t c)
                                   { return std::tolower(c); });

                    name += prefixes[ZIndices[i]] + name + " ";
                }
            }
            else // organic nomeclature
            {
            }

            return name;
        }

        void universe::initReaxParams()
        {
        }

        // loading and saving scenes

        void universe::saveFrame()
        {
            frame nFrame{};
            nFrame.positions = data.positions;
            nFrame.global_temperature = temp;

            m_frames.emplace_back(std::move(nFrame));
        }

        // Helper
        std::string formatTime()
        {
            const auto now = std::chrono::system_clock::now();

            std::string time_str = std::format("{:%Y_%m_%d%H_%M}", now);
            return time_str;
        }

        void universe::runVideo(const video &vid)
        {
        }

        video universe::saveAsVideo(const std::filesystem::path path, const std::string name)
        {
            if (m_frames.size() == 0)
                return {};

            int32_t m_framesToSaveFactor = 5;

            nlohmann::json json_video{};

            videoMetaData nMetadata{};
            nMetadata.box = box;
            nMetadata.num_atoms = atoms.size();
            size_t savedFrameCount = 0;
            for (size_t i = 0; i < m_frames.size(); ++i)
            {
                if (i % m_framesToSaveFactor == 0)
                    ++savedFrameCount;
            }
            nMetadata.num_frames = savedFrameCount;
            nMetadata.title = name.empty() ? formatTime() : name;

            json_video["metadata"] =
                {
                    {"title", nMetadata.title},
                    {"description", "Molecular dynamics trajectory"},
                    {"atoms", nMetadata.num_atoms},
                    {"frames", nMetadata.num_frames},
                    {"box", {box.x, box.y, box.z}}};

            json_video["positions"] = nlohmann::json::array();

            for (int32_t i = 0; i < m_frames.size(); ++i)
            {
                if (i % m_framesToSaveFactor != 0)
                    continue;

                const auto &frame = m_frames[i];
                if (frame.positions.size() != nMetadata.num_atoms)
                {
                    std::cerr << "[Video Save] Frame has wrong atom count! Skipping.\n";
                    continue;
                }

                for (const auto &pos : frame.positions)
                {
                    json_video["positions"].push_back(pos.x);
                    json_video["positions"].push_back(pos.y);
                    json_video["positions"].push_back(pos.z);
                }

                json_video["temperatures"].push_back(frame.global_temperature);
            }

            try
            {
                if (!std::filesystem::is_directory(path))
                    std::filesystem::create_directory(path);

                std::filesystem::path filepath = name.empty() ? path / ("video_" + formatTime() + ".json") : path / ("video_" + name + ".json");
                std::ofstream file(filepath);
                file.flush();
                file << json_video.dump();
                file.close();
            }
            catch (std::exception &e)
            {
                std::cout << "[Recorder] Failed to save trajectory: " << e.what() << std::endl;
                return video{};
            }

            std::cout << "[Recorder] Saved trajectory with " << json_video["metadata"]["frames"] << " frames and " << json_video["metadata"]["atoms"] << " atoms to " << path << '\n';

            video nVideo{};
            nVideo.frames = m_frames;
            nVideo.metadata = std::move(nMetadata);

            return nVideo;
        }

        void universe::saveScene(const std::filesystem::path path, const std::string name)
        {
            if (!std::filesystem::is_directory(path))
                std::filesystem::create_directory(path);

            nlohmann::json scene{};

            for (int32_t x = 0; x < data.positions.size(); ++x)
            {
                scene["posx"].emplace_back(data.positions[x].x);
                scene["posy"].emplace_back(data.positions[x].y);
                scene["posz"].emplace_back(data.positions[x].z);
                scene["velx"].emplace_back(data.velocities[x].x);
                scene["vely"].emplace_back(data.velocities[x].y);
                scene["velz"].emplace_back(data.velocities[x].z);

                scene["charge"].emplace_back(data.q[x]);
                scene["ZIndex"].emplace_back(atoms[x].ZIndex);
                scene["neutrons"].emplace_back(atoms[x].NCount);
                scene["electrons"].emplace_back(atoms[x].electrons);
                scene["boundCount"].emplace_back(atoms[x].bondCount);
                scene["chirality"].emplace_back(atoms[x].chirality);
            }

            for (int32_t b = 0; b < bonds.size(); ++b)
            {
                auto &bond = bonds[b];
                nlohmann::json b_json{};
                b_json["central"] = bond.centralAtom;
                b_json["bonded"] = bond.bondedAtom;
                b_json["equilibrium"] = bond.equilibriumLength;
                b_json["type"] = static_cast<uint8_t>(bond.type);
                b_json["k"] = bond.k;

                scene["bonds"].emplace_back(b_json);
            }

            for (int32_t a = 0; a < angles.size(); ++a)
            {
                auto &angle = angles[a];
                nlohmann::json a_json{};
                a_json["A"] = angle.A;
                a_json["B"] = angle.B;
                a_json["C"] = angle.C;
                a_json["K"] = angle.K;
                a_json["rad"] = angle.rad;

                scene["angles"].emplace_back(a_json);
            }

            for (int32_t d = 0; d < dihedral_angles.size(); ++d)
            {
                auto &dihedral = dihedral_angles[d];
                nlohmann::json d_json{};

                d_json["A"] = dihedral.A;
                d_json["B"] = dihedral.B;
                d_json["C"] = dihedral.C;
                d_json["D"] = dihedral.D;
                d_json["K"] = dihedral.K;
                d_json["periodicity"] = dihedral.periodicity;
                d_json["rad"] = dihedral.rad;

                scene["dihedrals"].emplace_back(d_json);
            }

            for (int32_t s = 0; s < subsets.size(); ++s)
            {
                auto &subset = subsets[s];
                nlohmann::json s_json{};

                s_json["mainAtom"] = subset.mainAtomIdx;
                s_json["bonded"] = subset.bondedSubsetIdx;
                s_json["bonding"] = subset.bondingSubsetIdx;
                s_json["connected_begin"] = subset.connectedBegin;
                s_json["connected_count"] = subset.connectedCount;
                s_json["hydrogen_begin"] = subset.hydrogenBegin;
                s_json["hydrogen_count"] = subset.hydrogenCount;

                scene["subsets"].emplace_back(s_json);
            }

            for (int32_t m = 0; m < molecules.size(); ++m)
            {
                auto &molecule = molecules[m];
                nlohmann::json m_json{};

                m_json["atomBegin"] = molecule.atomBegin;
                m_json["atomCount"] = molecule.atomCount;
                m_json["bondBegin"] = molecule.bondBegin;
                m_json["bondCount"] = molecule.bondCount;
                m_json["angleBegin"] = molecule.angleBegin;
                m_json["angleCount"] = molecule.angleCount;
                m_json["dihedralBegin"] = molecule.dihedralBegin;
                m_json["dihedralCount"] = molecule.dihedralCount;
                m_json["subsetBegin"] = molecule.subsetBegin;
                m_json["subsetCount"] = molecule.subsetCount;

                scene["molecules"].emplace_back(m_json);
            }

            scene["gravity"] = gravity;
            scene["isothermal"] = isothermal;
            scene["wall_collision"] = wall_collision;
            scene["roof_floor_collision"] = roof_floor_collision;
            scene["mag_gravity"] = mag_gravity;
            scene["boxx"] = box.x;
            scene["boxy"] = box.y;
            scene["boxz"] = box.z;

            try
            {
                std::filesystem::path filepath = name.empty() ? path / (formatTime() + ".json") : path / (name + ".json");
                std::ofstream file(filepath);
                if (!file.is_open())
                {
                    std::cerr << "[Save] Failed to create file: " << filepath << "\n";
                    return;
                }

                file << scene.dump();
                file.flush();
                file.close();

                std::cout << "[Save] Successfully saved: " << filepath << "\n";
            }
            catch (const std::exception &e)
            {
                std::cerr << "[Save] Exception while saving: " << e.what() << "\n";
            }
        }

        void universe::loadScene(const std::filesystem::path path)
        {
            std::ifstream file(path);
            if (!file.is_open() || path.extension() != ".json")
            {
                std::cerr << "[Simulation] Cannot open: " << path << '\n';
                return;
            }

            nlohmann::json scene;
            try
            {
                file >> scene;
            }
            catch (const nlohmann::json::parse_error &e)
            {
                std::cerr << "[Simulation] JSON parse error: " << e.what() << '\n';
                return;
            }

            atoms.clear();
            bonds.clear();
            angles.clear();
            dihedral_angles.clear();
            subsets.clear();
            molecules.clear();
            data.positions.clear();
            data.velocities.clear();
            data.q.clear();

            box.x = scene.value("boxx", 50.0f);
            box.y = scene.value("boxy", 50.0f);
            box.z = scene.value("boxz", 50.0f);

            gravity = scene.value("gravity", false);
            isothermal = scene.value("isothermal", true);
            wall_collision = scene.value("wall_collision", false);
            roof_floor_collision = scene.value("roof_floor_collision", false);
            mag_gravity = scene.value("mag_gravity", 9.81f);

            const auto &posx = scene["posx"];
            const auto &posy = scene["posy"];
            const auto &posz = scene["posz"];
            const auto &velx = scene["velx"];
            const auto &vely = scene["vely"];
            const auto &velz = scene["velz"];
            const auto &charges = scene["charge"];
            const auto &ZIndices = scene["ZIndex"];
            const auto &neutrons = scene["neutrons"];
            const auto &electrons = scene["electrons"];
            const auto &bondCounts = scene["boundCount"];
            const auto &chiralities = scene["chirality"];
            const auto &temperatures = scene["temperature"];
            const auto &bondorders = scene["bondOrder"];

            int32_t N = posx.size();

            atoms.reserve(N);
            data.positions.reserve(N);
            data.velocities.reserve(N);
            data.q.reserve(N);

            for (int32_t i = 0; i < N; ++i)
            {
                uint8_t Z = ZIndices[i];
                auto [sigma, epsilon] = constants::getAtomConstants(Z);
                float radius = sigma / 1.3f;
                uint8_t nNeutrons = neutrons[i];

                atom newAtom{};
                newAtom.ZIndex = Z;
                newAtom.radius = radius;
                newAtom.electrons = electrons[i];
                newAtom.NCount = nNeutrons;
                newAtom.mass = Z * MASS_PROTON + nNeutrons * MASS_NEUTRON + ((float)electrons[i]) * MASS_ELECTRON;
                newAtom.chirality = chiralities[i];
                newAtom.bondCount = bondCounts[i];

                atoms.emplace_back(std::move(newAtom));

                emplace_pos(glm::vec3(posx[i], posy[i], posz[i]));
                emplace_vel(glm::vec3(velx[i], vely[i], velz[i]));

                data.q.emplace_back(charges[i]);
            }

            data.forces.assign(N, glm::vec3(0.0f));

            if (scene.contains("bonds"))
            {
                for (const auto &b_json : scene["bonds"])
                {
                    bond b{};
                    b.centralAtom = b_json["central"];
                    b.bondedAtom = b_json["bonded"];
                    b.equilibriumLength = b_json["equilibrium"];
                    b.type = static_cast<BondType>(b_json["type"]);
                    b.k = b_json["k"];
                    bonds.emplace_back(b);
                }
            }

            if (scene.contains("angles"))
            {
                for (const auto &a_json : scene["angles"])
                {
                    angle a{};
                    a.A = a_json["A"];
                    a.B = a_json["B"];
                    a.C = a_json["C"];
                    a.K = a_json["K"];
                    a.rad = a_json["rad"];
                    angles.emplace_back(a);
                }
            }

            if (scene.contains("dihedrals"))
            {
                for (const auto &d_json : scene["dihedrals"])
                {
                    dihedral_angle d{};
                    d.A = d_json["A"];
                    d.B = d_json["B"];
                    d.C = d_json["C"];
                    d.D = d_json["D"];
                    d.K = d_json["K"];
                    d.periodicity = d_json["periodicity"];
                    d.rad = d_json["rad"];
                    dihedral_angles.emplace_back(d);
                }
            }

            if (scene.contains("subsets"))
            {
                for (const auto &s_json : scene["subsets"])
                {
                    subset s{};
                    s.mainAtomIdx = s_json["mainAtom"];
                    s.bondedSubsetIdx = s_json.value("bonded", UINT32_MAX);
                    s.bondingSubsetIdx = s_json.value("bonding", UINT32_MAX);

                    s.connectedBegin = s_json.value("connected_begin", UINT32_MAX);
                    s.connectedCount = s_json.value("connected_count", UINT32_MAX);
                    s.hydrogenCount = s_json.value("hydrogen_count", UINT32_MAX);
                    s.hydrogenBegin = s_json.value("hydrogen_begin", UINT32_MAX);

                    subsets.emplace_back(s);
                }
            }

            if (scene.contains("molecules"))
            {
                for (const auto &m_json : scene["molecules"])
                {
                    molecule m{};

                    m.atomBegin = m_json["atomBegin"];
                    m.atomCount = m_json["atomCount"];
                    m.bondBegin = m_json["bondBegin"];
                    m.bondCount = m_json["bondCount"];
                    m.angleBegin = m_json["angleBegin"];
                    m.angleCount = m_json["angleCount"];
                    m.dihedralBegin = m_json["dihedralBegin"];
                    m.dihedralCount = m_json["dihedralCount"];
                    m.subsetBegin = m_json["subsetBegin"];
                    m.subsetCount = m_json["subsetCount"];

                    molecules.emplace_back(std::move(m));
                }
            }

            rebuildBondTopology();

            std::cout << "[Simulation] Successfully loaded scene: " << path.filename()
                      << " (" << atoms.size() << " atoms, "
                      << molecules.size() << " molecules)\n";
        }

        void universe::loadFrames(const std::filesystem::path path)
        {
            if (!std::filesystem::exists(path) || path.extension() != ".json")
            {
                std::cerr << "[Video Load] Invalid or missing file: " << path << "\n";
                return;
            }

            std::ifstream file(path, std::ios::binary);
            if (!file.is_open())
            {
                std::cerr << "[Video Load] Cannot open file: " << path << "\n";
                return;
            }

            nlohmann::json j;
            try
            {
                file >> j;
            }
            catch (const nlohmann::json::parse_error &e)
            {
                std::cerr << "[Video Load] JSON parse error in " << path << ": " << e.what() << "\n";
                return;
            }

            m_frames.clear();

            if (!j.contains("metadata") || !j["metadata"].is_object())
            {
                std::cerr << "[Video Load] Missing or invalid metadata in " << path << "\n";
                return;
            }

            auto &meta = j["metadata"];
            size_t expectedAtoms = meta.value("atoms", 0ull);
            size_t expectedFrames = meta.value("frames", 0ull);
            glm::vec3 loadedBox =
                {
                    meta["box"][0].get<float>(),
                    meta["box"][1].get<float>(),
                    meta["box"][2].get<float>()};

            if (expectedAtoms == 0 || expectedFrames == 0)
            {
                std::cerr << "[Video Load] Empty or invalid frame/atom count\n";
                return;
            }

            if (glm::length(loadedBox - box) > 1e-3f)
            {
                box = loadedBox;
            }

            if (!j.contains("positions") || !j["positions"].is_array())
            {
                std::cerr << "[Video Load] Missing or invalid 'positions' array\n";
                return;
            }

            auto &posArray = j["positions"];
            /*             if (posArray.size() != expectedFrames * expectedAtoms * 3)
                        {
                            std::cerr << "[Video Load] Position array size mismatch! Expected "
                                    << expectedFrames * expectedAtoms * 3 << ", got " << posArray.size() << "\n";
                            return;
                        } */

            std::vector<float> temperatures;
            if (j.contains("temperatures") && j["temperatures"].is_array())
            {
                temperatures = j["temperatures"].get<std::vector<float>>();
                if (temperatures.size() != expectedFrames)
                {
                    std::cerr << "[Video Load] Temperatures size mismatch! Expected "
                              << expectedFrames << ", got " << temperatures.size() << "\n";
                    temperatures.clear();
                }
            }

            std::vector<float> posArrayVec = posArray.get<std::vector<float>>();
            m_frames.resize(expectedFrames);

            size_t idx = 0;
            for (size_t f = 0; f < expectedFrames; ++f) {
                auto& newFrame = m_frames[f];
                newFrame.positions.resize(expectedAtoms);

                for (size_t a = 0; a < expectedAtoms; ++a) {
                    newFrame.positions[a] = {posArrayVec[idx], posArrayVec[idx+1], posArrayVec[idx+2]};
                    idx += 3;
                }

                newFrame.global_temperature = !temperatures.empty() ? temperatures[f] : 300.0f;
            }


            std::cout << "[Video Load] Successfully loaded " << m_frames.size()
                      << " frames with " << expectedAtoms << " atoms from " << path << "\n";
        }
    } // namespace fun
} // namespace sim
