#include "universe.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <set>
#include <queue>
#include <fstream>

#include "simulation/constants.hpp"

#include <glad/glad.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/norm.hpp>

namespace sim
{
    namespace fun
    {
        universe::universe(const universe_create_info &create_info, rendering_engine &rendering_engine, parameter_table& paramTable)
            : box(create_info.box), gravity(create_info.has_gravity), mag_gravity(create_info.mag_gravity),
              wall_collision(create_info.wall_collision), isothermal(create_info.isothermal),
              HMassRepartitioning(create_info.HMassRepartitioning), roof_floor_collision(create_info.roof_floor_collision),
              log_flags(create_info.log_flags), rendering_eng(rendering_engine), m_reaction(create_info.reaction),
              m_parameterTable(paramTable)
        {
        }

        universe::universe(const std::filesystem::path path, rendering_engine &rendering_engine, parameter_table& paramTable)
            : rendering_eng(rendering_engine), m_parameterTable(paramTable)
        {
            loadScene(path);
        }

        void universe::draw(rendering_info info)
        {
            auto positionsvec = positions();
            auto velocitiesvec = velocities();

            std::vector<float> qvec{};
            float* q = m_atomStorage.qData();
            for (int32_t i = 0; i < m_atomStorage.mobileCount(); ++i)
                qvec.emplace_back(q[i]);

            rendering_simulation_info sim_info{.positions = m_displayPositions.empty() ? positionsvec : m_displayPositions,
                                               .velocities = velocitiesvec,
                                               .q = qvec,
                                               .atoms = atomData.atoms,
                                               .bonds = atomData.bonds,
                                               .molecules = atomData.molecules,
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

            rendering_eng.draw(info, sim_info);
        }

        void universe::moleculeDrawOrder()
        {
        
        }

        int32_t universe::createAtom(glm::vec3 p, glm::vec3 v, uint8_t ZIndex, uint8_t numNeutrons, uint8_t numElectron, int32_t chirality)
        {
            atom newAtom{};
            newAtom.ZIndex = ZIndex;
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

            atomData.atoms.emplace_back(std::move(newAtom));
            atomData.frozen_atoms.emplace_back(false);

            size_t atomSize = atomData.atoms.size();

            m_atomStorage.resize(atomSize);

            const LJParams ljParams = m_parameterTable.lj(ZIndex);

            sim::atomView view{};
            view.i = atomSize - 1;
            view.vx = v.x;
            view.vy = v.y;
            view.vz = v.z;
            view.fx = 0.f;
            view.fy = 0.f;
            view.fz = 0.f;
            view.px = p.x;
            view.py = p.y;
            view.pz = p.z;
            view.invMass = 1.f / newAtom.mass;
            view.q = float(ZIndex - numElectron);

            view.sigma = ljParams.sigma;            
            view.epsilon = ljParams.epsilon;

            // Morse Params (Placeholder)
            view.De = 976.f;
            view.re = 1.10f;
            view.a  = 2.3f;

            m_atomStorage.setAtom(view);

            return atomSize - 1;
        }

        void universe::createBond(int32_t idx1, int32_t idx2, BondType type)
        {
            if (idx1 >= atomData.atoms.size() || idx2 >= atomData.atoms.size() || idx1 == idx2)
                return;

            if (areBonded(idx1, idx2)) return;

            bond nBond{};
            nBond.bondedAtom = idx1;
            nBond.centralAtom = idx2;
            nBond.type = type;

            if (!m_reaction)
            {
                int8_t bondCount = static_cast<int8_t>(type);
    
                atomData.atoms[idx1].bondCount += bondCount;
                atomData.atoms[idx2].bondCount += bondCount;
            }
            
            float EN1 = constants::getElectronegativity(atomData.atoms[idx1].ZIndex);
            float EN2 = constants::getElectronegativity(atomData.atoms[idx2].ZIndex);
            float deltaEN = std::abs(EN1 - EN2);
            
            if (deltaEN > 0.1f) // Significant electronegativity difference
            {
                float charge = deltaEN * 0.31f;

                if (EN1 > EN2)
                {
                    m_atomStorage.charge(idx1) -= charge;
                    m_atomStorage.charge(idx2) += charge;
                }
                else
                {
                    m_atomStorage.charge(idx2) -= charge;
                    m_atomStorage.charge(idx1) += charge;
                }
            }
            
            nBond.order = static_cast<float>(type) + 0.1f;
            
            atomData.bonds.emplace_back(std::move(nBond));
        }

        void universe::updateAngles(int32_t idx1)
        {
            atomData.angles.erase(
                std::remove_if(atomData.angles.begin(), atomData.angles.end(),
                    [&](const angle& a) { return a.B == idx1; }),
                atomData.angles.end());

            std::vector<uint32_t> neighbors;
            for (const auto& bond : atomData.bonds) 
            {
                if (bond.bondedAtom == idx1) neighbors.push_back(bond.centralAtom);
                else if (bond.centralAtom == idx1) neighbors.push_back(bond.bondedAtom);
            }

            /* for (int32_t i = 0; i < neighbors.size(); ++i) 
            {
                for (int32_t j = i + 1; j < neighbors.size(); ++j) 
                {
                    angle newAngle{neighbors[i], (uint32_t)idx1, neighbors[j]};
                    atomData.angles.emplace_back(std::move(newAngle));
                }
            } */

            /* atomData.dihedral_angles.erase(
                std::remove_if(atomData.dihedral_angles.begin(), atomData.dihedral_angles.end(),
                    [&](const dihedral_angle& d) { return (d.B == idx1 || d.C == idx1); }),
                atomData.dihedral_angles.end());

            for (auto n : neighbors) 
            {
                for (auto n1 : neighbors) 
                {
                    if (n1 == n) continue;

                    for (const auto& bond : atomData.bonds) 
                    {
                        int32_t n2 = -1;
                        if (bond.centralAtom == n && bond.bondedAtom != idx1) n2 = bond.bondedAtom;
                        else if (bond.bondedAtom == n && bond.centralAtom != idx1) n2 = bond.centralAtom;

                        if (n2 != -1) 
                        {
                            dihedral_angle dih{n1, (uint32_t)idx1, n, (uint32_t)n2};
                            atomData.dihedral_angles.emplace_back(std::move(dih));
                        }
                    }
                }
            } */

            /* atomData.improper_angles.erase(
                std::remove_if(atomData.improper_angles.begin(), atomData.improper_angles.end(),
                    [&](const dihedral_angle& d) { return d.A == idx1 || d.B == idx1 || d.C == idx1 || d.D == idx1; }),
                atomData.improper_angles.end());

            if (neighbors.size() >= 3) 
            {
                for (int32_t i = 0; i < neighbors.size(); ++i) 
                {
                    for (int32_t j = i + 1; j < neighbors.size(); ++j) 
                    {
                        for (int32_t k = j + 1; k < neighbors.size(); ++k) 
                        {
                            dihedral_angle improper{idx1, neighbors[i], neighbors[j], neighbors[k]};
                            atomData.improper_angles.emplace_back(std::move(improper));
                        }
                    }
                }
            } */
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
                neighbourZs.emplace_back(atomData.atoms[bondedAtom].ZIndex);
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

            atomData.subsets.emplace_back(std::move(nSubset));
            return atomData.subsets.size() - 1; // Index
        }

        void universe::createMolecule(molecule_structure structure, glm::vec3 pos, glm::vec3 vel)
        {
            int32_t baseAtomIndex = atomData.atoms.size();

            molecule nMolecule{};

            glm::vec3 centroid{0.f, 0.f, 0.f};
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

                glm::vec3 end_pos = structure.positions[i] + pos;
                createAtom(end_pos, vel, a.ZIndex, a.NIndex, a.ZIndex, a.chirality);
                m_atomStorage.charge(baseAtomIndex + i) += structure.atoms[i].charge;
            }

            nMolecule.atomBegin = baseAtomIndex;
            nMolecule.atomCount = structure.atoms.size();

            int32_t baseBondIndex = atomData.bonds.size();
            for (int32_t b = 0; b < structure.bonds.size(); ++b)
            {
                const def_bond &db = structure.bonds[b];
                
                if (db.type == BondType::NONE) continue;

                int32_t central = baseAtomIndex + db.centralAtomIdx;
                int32_t bonded = baseAtomIndex + db.bondingAtomIdx;

                createBond(bonded, central, db.type);
            }

            nMolecule.bondBegin = baseBondIndex;
            nMolecule.bondCount = structure.bonds.size();

            int32_t baseSubset = atomData.subsets.size();
            for (int32_t s = 0; s < structure.subsets.size(); ++s)
                createSubset(structure.subsets[s], baseAtomIndex, baseSubset);

            nMolecule.subsetBegin = baseSubset;
            nMolecule.subsetCount = structure.subsets.size();

            if (!m_reaction && nMolecule.subsetCount > 0 && structure.atoms[0].charge == 0)
                balanceMolecularCharges(atomData.subsets[baseSubset]);

            rebuildBondTopology();

            int32_t baseAngle = atomData.angles.size();
            for (int32_t a = 0; a < structure.angles.size(); ++a)
            {
                angle angle = structure.angles[a];
                angle.A += baseAtomIndex;
                angle.B += baseAtomIndex;
                angle.C += baseAtomIndex;

                atomData.angles.emplace_back(std::move(angle));
            }

            nMolecule.angleBegin = baseAngle;
            nMolecule.angleCount = structure.angles.size();

            int32_t baseDihedral = atomData.dihedral_angles.size();
            for (int32_t a = 0; a < structure.dihedral_angles.size(); ++a)
            {
                dihedral_angle angle = structure.dihedral_angles[a];
                angle.A += baseAtomIndex;
                angle.B += baseAtomIndex;
                angle.C += baseAtomIndex;
                angle.D += baseAtomIndex;

                atomData.dihedral_angles.emplace_back(std::move(angle));
            }

            int32_t baseImproper = atomData.improper_angles.size();
            for (int32_t a = 0; a < structure.improper_angles.size(); ++a)
            {
                dihedral_angle angle = structure.improper_angles[a];
                angle.A += baseAtomIndex;
                angle.B += baseAtomIndex;
                angle.C += baseAtomIndex;
                angle.D += baseAtomIndex;

                atomData.improper_angles.emplace_back(std::move(angle));
            }

            nMolecule.dihedralBegin = baseDihedral;
            nMolecule.dihedralCount = structure.dihedral_angles.size();

            atomData.molecules.emplace_back(std::move(nMolecule));
        }

        void universe::balanceMolecularCharges(subset &mol)
        {
            /* auto it = std::find_if(atomData.subsets.begin(), atomData.subsets.end(), [&](const subset &s)
                                   { return &s == &mol; });
            if (it == atomData.subsets.end())
                return;

            int32_t startIdx = std::distance(atomData.subsets.begin(), it);

            std::vector<int32_t> atomsToBalance;
            std::unordered_set<int32_t> visited;
            std::queue<int32_t> q;
            q.push(startIdx);
            visited.insert(startIdx);

            while (!q.empty())
            {
                int32_t currentIdx = q.front();
                q.pop();
                const auto &currentSubset = atomData.subsets[currentIdx];

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
                if (nextIdx < atomData.subsets.size() && nextIdx != SIZE_MAX && visited.insert(nextIdx).second)
                {
                    q.push(nextIdx);
                }
            }

            if (atomsToBalance.empty())
                return;

            std::sort(atomsToBalance.begin(), atomsToBalance.end());
            auto last = std::unique(atomsToBalance.begin(), atomsToBalance.end());
            atomsToBalance.erase(last, atomsToBalance.end());

            std::vector<std::vector<int32_t>> neighborLists(atomData.atoms.size());
            for (const auto &bond : atomData.bonds)
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
                if (idx >= atomData.atoms.size())
                    continue;

                uint8_t usualValence = constants::getUsualBonds(atomData.atoms[idx].ZIndex);
                int deficit = usualValence - static_cast<int>(atomData.atoms[idx].bondCount);
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
            } */
        }

        void reflectVelocity(float& vx, float& vy, float& vz, const glm::vec3& normal)
        {
            glm::vec3 v(vx, vy, vz);

            float dot = glm::dot(v, normal);
            glm::vec3 reflected = v - 2.0f * dot * normal;

            vx = reflected.x;
            vy = reflected.y;
            vz = reflected.z;
        }

        void universe::boundCheck(uint32_t i)
        {
            float* __restrict x  = m_atomStorage.xData();
            float* __restrict y  = m_atomStorage.yData();
            float* __restrict z  = m_atomStorage.zData();
            float* __restrict vx = m_atomStorage.vxData();
            float* __restrict vy = m_atomStorage.vyData();
            float* __restrict vz = m_atomStorage.vzData();

            bool collided = false;

            // ====================== X Axis ======================
            if (wall_collision)
            {
                if (x[i] < 0.0f)
                {
                    x[i] = 0.0f;
                    reflectVelocity(vx[i], vy[i], vz[i], glm::vec3(1.0f, 0.0f, 0.0f));
                    collided = true;
                }
                else if (x[i] > box.x)
                {
                    x[i] = box.x;
                    reflectVelocity(vx[i], vy[i], vz[i], glm::vec3(-1.0f, 0.0f, 0.0f));
                    collided = true;
                }
            }
            else // Periodic
            {
                x[i] = std::fmod(std::fmod(x[i], box.x) + box.x, box.x);
            }

            // ====================== Y Axis ======================
            if (wall_collision)
            {
                if (y[i] < 0.0f)
                {
                    y[i] = 0.0f;
                    reflectVelocity(vx[i], vy[i], vz[i], glm::vec3(0.0f, 1.0f, 0.0f));
                    collided = true;
                }
                else if (y[i] > box.y)
                {
                    y[i] = box.y;
                    reflectVelocity(vx[i], vy[i], vz[i], glm::vec3(0.0f, -1.0f, 0.0f));
                    collided = true;
                }
            }
            else
            {
                y[i] = std::fmod(std::fmod(y[i], box.y) + box.y, box.y);
            }

            // ====================== Z Axis ======================
            if (roof_floor_collision)
            {
                if (z[i] < 0.0f)
                {
                    z[i] = 0.0f;
                    reflectVelocity(vx[i], vy[i], vz[i], glm::vec3(0.0f, 0.0f, 1.0f));
                    collided = true;
                }
                else if (z[i] > box.z)
                {
                    z[i] = box.z;
                    reflectVelocity(vx[i], vy[i], vz[i], glm::vec3(0.0f, 0.0f, -1.0f));
                    collided = true;
                }
            }
            else
            {
                z[i] = std::fmod(std::fmod(z[i], box.z) + box.z, box.z);
            }

            if (collided)
            {
                vx[i] *= 0.999f;
                vy[i] *= 0.999f;
                vz[i] *= 0.999f;
            }
        }

        // Energy Calculation
        float universe::calculateKineticEnergy()
        {
            float kinetic_energy = 0.0f;
            for (int32_t i = 0; i < atomData.atoms.size(); ++i)
                kinetic_energy += 0.5f * atomData.atoms[i].mass * glm::length2(m_atomStorage.velocity(i));
 
            return kinetic_energy;
        }

        void universe::COMDrift()
        {
            auto &storage = m_atomStorage;

            float* __restrict vx = storage.vxData();
            float* __restrict vy = storage.vyData();
            float* __restrict vz = storage.vzData();
            float* __restrict invMass = storage.invMassData();

            size_t N = storage.mobileCount();

            glm::vec3 totalMomentum{0.f};
            float totalMass = 0.f;

            for (size_t i = 0; i < N; ++i)
            {
                float mass = 1.0f / invMass[i]; 
                totalMomentum.x += mass * vx[i];
                totalMomentum.y += mass * vy[i];
                totalMomentum.z += mass * vz[i];
                totalMass += mass;
            }

            glm::vec3 correction = totalMomentum / totalMass;

            for (size_t i = 0; i < N; ++i)
            {
                vx[i] -= correction.x;
                vy[i] -= correction.y;
                vz[i] -= correction.z;
            }
        }

        // TO DO
        std::string universe::moleculeName(const std::vector<uint32_t> &subsetIdx)
        {
            std::map<uint8_t, int32_t> ZIndices;

            for (int32_t i = 0; i < subsetIdx.size(); ++i)
            {
                auto &s = atomData.subsets[subsetIdx[i]];
                int32_t centralAtom = s.mainAtomIdx;
                ++ZIndices[atomData.atoms[centralAtom].ZIndex];
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
                for (int32_t i = 1; i <= 119; ++i)
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

        // loading and saving scenes

        void universe::saveFrame()
        {
            frame nFrame{};
            nFrame.positions = positions();

            m_frames.emplace_back(std::move(nFrame));
        }

        // Helper
        std::string formatTime()
        {
            const auto now = std::chrono::system_clock::now();

            std::string time_str = std::format("{:%Y_%m_%d%H_%M}", now);
            return time_str;
        }

        video universe::saveAsVideo(const std::filesystem::path path, const std::string name)
        {
            if (m_frames.size() == 0)
                return {};

            nlohmann::json json_video{};

            videoMetaData nMetadata{};
            nMetadata.box = box;
            nMetadata.num_atoms = atomData.atoms.size();
            nMetadata.num_frames = m_frames.size();
            nMetadata.title = name.empty() ? formatTime() : name;

            json_video["metadata"] =
                {
                    {"title", nMetadata.title},
                    {"description", "Molecular dynamics trajectory"},
                    {"atomData.atoms", nMetadata.num_atoms},
                    {"frames", nMetadata.num_frames},
                    {"box", {box.x, box.y, box.z}}};

            json_video["positions"] = nlohmann::json::array();

            for (int32_t i = 0; i < m_frames.size(); ++i)
            {
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

            std::cout << "[Recorder] Saved trajectory with " << json_video["metadata"]["frames"] << " frames and " << json_video["metadata"]["atomData.atoms"] << " atomData.atoms to " << path << '\n';

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

            float* __restrict x  = m_atomStorage.xData();
            float* __restrict y  = m_atomStorage.yData();
            float* __restrict z  = m_atomStorage.zData();
            float* __restrict vx = m_atomStorage.vxData();
            float* __restrict vy = m_atomStorage.vyData();
            float* __restrict vz = m_atomStorage.vzData();
            float* __restrict q  = m_atomStorage.qData();

            for (int32_t i = 0; i < m_atomStorage.mobileCount(); ++i)
            {

                scene["posx"].emplace_back(x[i]);
                scene["posy"].emplace_back(y[i]);
                scene["posz"].emplace_back(z[i]);

                scene["velx"].emplace_back(vx[i]);
                scene["vely"].emplace_back(vy[i]);
                scene["velz"].emplace_back(vz[i]);

                scene["charge"].emplace_back(q[i]);

                scene["ZIndex"].emplace_back(atomData.atoms[i].ZIndex);
                scene["neutrons"].emplace_back(atomData.atoms[i].NCount);
                scene["electrons"].emplace_back(atomData.atoms[i].electrons);
                scene["boundCount"].emplace_back(atomData.atoms[i].bondCount);
                scene["chirality"].emplace_back(atomData.atoms[i].chirality);
            }

            for (int32_t b = 0; b < atomData.bonds.size(); ++b)
            {
                auto &bond = atomData.bonds[b];
                nlohmann::json b_json{};
                b_json["central"] = bond.centralAtom;
                b_json["bonded"] = bond.bondedAtom;
                b_json["type"] = static_cast<uint8_t>(bond.type);

                scene["bonds"].emplace_back(b_json);
            }

            for (int32_t a = 0; a < atomData.angles.size(); ++a)
            {
                auto &angle = atomData.angles[a];
                nlohmann::json a_json{};
                a_json["A"] = angle.A;
                a_json["B"] = angle.B;
                a_json["C"] = angle.C;
                a_json["K"] = angle.K;
                a_json["rad"] = angle.rad;

                scene["angles"].emplace_back(a_json);
            }

            for (int32_t d = 0; d < atomData.dihedral_angles.size(); ++d)
            {
                auto &dihedral = atomData.dihedral_angles[d];
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

            for (int32_t s = 0; s < atomData.subsets.size(); ++s)
            {
                auto &subset = atomData.subsets[s];
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

            for (int32_t m = 0; m < atomData.molecules.size(); ++m)
            {
                auto &molecule = atomData.molecules[m];
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

            clear();

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

            atomData.atoms.reserve(N);
            m_atomStorage.resize(N);

            for (int32_t i = 0; i < N; ++i)
            {
                uint8_t Z = ZIndices[i];
                auto [sigma, epsilon] = constants::getAtomConstants(Z);
                float radius = sigma / 1.3f;
                uint8_t nNeutrons = neutrons[i];

                atom newAtom{};
                newAtom.ZIndex = Z;
                newAtom.electrons = electrons[i];
                newAtom.NCount = nNeutrons;
                newAtom.mass = Z * MASS_PROTON + nNeutrons * MASS_NEUTRON + ((float)electrons[i]) * MASS_ELECTRON;
                newAtom.chirality = chiralities[i];
                newAtom.bondCount = bondCounts[i];

                atomData.atoms.emplace_back(std::move(newAtom));
                atomData.frozen_atoms.emplace_back(false);

                size_t atomSize = atomData.atoms.size();

                m_atomStorage.resize(atomSize);

                const auto ljParams = constants::getAtomConstants(Z);

                sim::atomView view{};
                view.i = i;
                view.vx = velx[i];
                view.vy = vely[i];
                view.vz = velz[i];
                view.fx = 0.f;
                view.fy = 0.f;
                view.fz = 0.f;
                view.px = posx[i];
                view.py = posy[i];
                view.pz = posz[i];
                view.invMass = 1.f / newAtom.mass;
                view.q = charges[i];
                view.sigma = ljParams.first;            
                view.epsilon = ljParams.second;            

                m_atomStorage.setAtom(view);
            }

            m_atomStorage.clearForces();

            if (scene.contains("bonds"))
            {
                for (const auto &b_json : scene["bonds"])
                {
                    bond b{};
                    b.centralAtom = b_json["central"];
                    b.bondedAtom = b_json["bonded"];
                    b.type = static_cast<BondType>(b_json["type"]);
                    atomData.bonds.emplace_back(b);
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
                    atomData.angles.emplace_back(a);
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
                    atomData.dihedral_angles.emplace_back(d);
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

                    atomData.subsets.emplace_back(s);
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

                    atomData.molecules.emplace_back(std::move(m));
                }
            }

            rebuildBondTopology();

            std::cout << "[Simulation] Successfully loaded scene: " << path.filename()
                      << " (" << atomData.atoms.size() << " atomData.atoms, "
                      << atomData.molecules.size() << " molecules)\n";
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
            /*             if (posArray.size() != expectedFrames * expectedatomData.atoms * 3)
                        {
                            std::cerr << "[Video Load] Position array size mismatch! Expected "
                                    << expectedFrames * expectedatomData.atoms * 3 << ", got " << posArray.size() << "\n";
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
                      << " frames with " << expectedAtoms << " atomData.atoms from " << path << "\n";
        }
    } // namespace fun
} // namespace sim
