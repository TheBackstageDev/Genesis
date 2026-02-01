#pragma once

#include "fundamental_structures.hpp"
#include "rendering_engine.hpp"
#include <SFML/Graphics.hpp>

#include <thread>
#include <future>
#include <queue>
#include <atomic>
#include <filesystem>
#include <json.hpp>

#ifndef _NODISCARD
#define _NODISCARD
#endif

namespace sim
{
    namespace fun
    {
        struct simData
        {
            alignas(64) std::vector<glm::vec3> positions;
            alignas(64) std::vector<glm::vec3> velocities;
            alignas(64) std::vector<glm::vec3> forces;
            alignas(64) std::vector<float> q, bond_orders, temperature;
        };

        constexpr int32_t offsets[14][3] = 
        {
            {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0},
            {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {-1,0,1},
            {1,-1,1}, {0,-1,1}, {-1,-1,1}, {-1,0,0}
        };

        struct logging_flags
        {
            bool log_reactions = true;
        };

        struct universe_create_info
        {
            bool has_gravity = false;
            bool wall_collision = false;
            bool isothermal = true;
            bool HMassRepartitioning = false;

            logging_flags log_flags;

            float mag_gravity = 9.8f;
            glm::vec3 box{CELL_CUTOFF * 2.f, CELL_CUTOFF * 2.f, CELL_CUTOFF * 2.f};
        };

        struct videoMetaData
        {
            std::string title;
            std::string description;

            glm::vec3 box;
            size_t num_atoms;
            size_t num_frames;
        };

        struct frame
        {
            alignas(64) std::vector<glm::vec3> positions{};
            std::vector<float> temperatures{};

            std::unordered_map<uint32_t, float> energyLog{};

            float global_temperature = 0.0f;
        };

        struct video
        {
            std::vector<frame> frames; // index is defined by order
            std::vector<size_t> keyFrames; // for events

            videoMetaData metadata;
        };

        class universe
        {
        public:
            universe(const universe_create_info& create_info, rendering_engine& rendering);
            universe(std::filesystem::path scene, rendering_engine& rendering);

            int32_t createAtom(glm::vec3 p, glm::vec3 v, uint8_t ZIndex = 1, uint8_t numNeutrons = 0, uint8_t numElectrons = 1, int32_t chirality = 0);
            int32_t createSubset(const def_subset& nSub, const int32_t baseAtom, const int32_t baseSubset);
            void createMolecule(molecule_structure structure, sf::Vector3f pos, sf::Vector3f vel = {0.f, 0.f, 0.f});

            void createBond(int32_t idx1, int32_t idx2, BondType type = BondType::SINGLE);
            void balanceMolecularCharges(subset& mol);
            
            void update(float targetTemperature = 1.0f, float targetPressure = 0.f);
            void draw(sf::RenderTarget& target, rendering_info info);

            void runVideo(const video& vid);

            // Scenario stuff
            void saveScene(const std::filesystem::path path, const std::string name = "");
            void loadScene(const std::filesystem::path path);

            void loadFrames(const std::vector<frame>& nFrames) { m_frames = std::move(nFrames); }
            void loadFrames(const std::filesystem::path path);
            void saveFrame();
            
            _NODISCARD video saveAsVideo(const std::filesystem::path path, const std::string name = "");
            _NODISCARD const std::vector<frame>& getFrames() const { return m_frames; }
            const frame& getFrame(size_t i) const { return m_frames[i]; }

            int32_t numFrames() { return m_frames.size(); }
            void clearFrames() { m_frames.clear(); }

            void clear()
            {
                data.positions.clear();
                data.velocities.clear();
                data.forces.clear();
                data.q.clear();
                data.temperature.clear();
                data.bond_orders.clear();

                molecules.clear();
                atoms.clear();
                bonds.clear();
                subsets.clear();
                rings.clear();

                m_highlightedAtoms.clear();
                m_highlightedBonds.clear();
                m_Arrows.clear();
            }

            void pause() { m_paused = true; }
            void unpause() { m_paused = false; }

            void highlightAtom(uint32_t index) { m_highlightedAtoms.emplace_back(index); }
            void highlightBond(uint32_t index1, uint32_t index2) { m_highlightedBonds.emplace_back(index1, index2); }
            void createArrow(uint32_t from, uint32_t to) { m_Arrows.emplace_back(from, to); }

            // Sets

            void setTimescale(float timescale = 1.0f) { m_Timescale = timescale; }

            void setDisplayPositions(const std::vector<glm::vec3>& nPositions)
            {
                m_displayPositions = nPositions;
            }

            void setPositions(const std::vector<glm::vec3>& nPositions)
            {
                data.positions = nPositions;
            }

            void setPosition(size_t i, glm::vec3 p)
            {
                data.positions[i] = p;
            }
            
            void addVelocity(size_t i, glm::vec3 v)
            {
                data.velocities[i] += v;
            }

            // Gets
            
            core::camera_t& getRenderingCamera() { return rendering_eng.camera(); }
            bool isPaused() { return m_paused; }

            const std::vector<atom>& getAtoms() const { return atoms; }
            std::vector<subset>& getSubsets() { return subsets; }
            int32_t numBonds() { return bonds.size(); }
            int32_t numAtoms() { return atoms.size(); }
            int32_t numMolecules() { return molecules.size(); }
            float getTimescale() { return m_Timescale; }
            float getEffectiveDT() { return m_Timescale * FEMTOSECOND; }
            const subset& getSubset(int32_t index) { return subsets[index]; }

            float temperature() const { return temp; }
            float pressure() { return calculatePressure(); }
            size_t timestep() const { return timeStep; }

            glm::vec3 boxSizes() { return box; }

            _NODISCARD std::vector<glm::vec3> positions() const 
            {
                return data.positions;
            }

            void clearDisplayPositions()
            {
                m_displayPositions.clear();
            }

            glm::vec3 getPosition(size_t i)
            {
                return data.positions[i];
            }

            // Helper Funcs

            sf::Vector3f minImageVec(sf::Vector3f dr)
            {
                if (wall_collision) return dr;

                dr.x -= box.x * std::round(dr.x / box.x);
                dr.y -= box.y * std::round(dr.y / box.y);
                dr.z -= box.z * std::round(dr.z / box.z);
                return dr;
            }

            glm::vec3 minImageVec(glm::vec3 dr)
            {
                if (wall_collision) return dr;

                dr.x -= box.x * std::round(dr.x / box.x);
                dr.y -= box.y * std::round(dr.y / box.y);
                dr.z -= box.z * std::round(dr.z / box.z);
                return dr;
            }

            inline bool areBonded(uint32_t i, uint32_t j) const
            {
                if (i >= bondedBits.size() || j >= bondedBits[i].size() * 64) return false;
                size_t word = j / 64;
                size_t bit  = j % 64;
                return (bondedBits[i][word] & (1ull << bit)) != 0;
            }

            // Energies
            float calculateKineticEnergy();
        private:            
            void boundCheck(uint32_t i);

            float ljPot(uint32_t i, uint32_t j);
            float wolfForce(float r, float qi_qj); 
            sf::Vector3f ljForce(uint32_t i, uint32_t j);
            sf::Vector3f coulombForce(uint32_t i, uint32_t j, sf::Vector3f& dr_vec);

            std::vector<sf::Vector3f> processCellUnbonded(int32_t ix, int32_t iy, int32_t iz);
            std::vector<sf::Vector3f> processPartialCellUnbonded(int32_t ix, int32_t iy, int32_t iz, int32_t atom_start, int32_t atom_end);
            void calcUnbondedForcesParallel();
            void calcBondedForcesParallel();

            float calculatePressure();
            void setPressure(float bar = 100);
            void setTemperature(float kelvin = 0.f);
            float calculateDihedral(const glm::vec3& pa, const glm::vec3& pb, const glm::vec3& pc, const glm::vec3& pd);

            // Energies
            float calculateAtomTemperature(int32_t i);
            float calculateBondEnergy(int32_t i, int32_t j, float bo_sigma, float bo_pi, float bo_pp);

            rendering_engine& rendering_eng;

            simData data;
            std::vector<glm::vec3> m_displayPositions{};

            std::vector<std::vector<uint64_t>> bondedBits;
            void markBonded(uint32_t i, uint32_t j)
            {
                if (i >= bondedBits.size()) bondedBits.resize(atoms.size());
                if (j >= bondedBits.size()) bondedBits.resize(atoms.size());
                while (bondedBits[i].size() * 64 <= j) bondedBits[i].push_back(0);
                while (bondedBits[j].size() * 64 <= i) bondedBits[j].push_back(0);

                size_t word_i = j / 64;
                size_t bit_i  = j % 64;
                size_t word_j = i / 64;
                size_t bit_j  = i % 64;

                bondedBits[i][word_i] |= (1ull << bit_i);
                bondedBits[j][word_j] |= (1ull << bit_j);
            }

            std::vector<atom> atoms;
            std::vector<subset> subsets;
            std::vector<molecule> molecules;
            std::vector<bool> frozen_atoms;
            
            std::vector<angle> angles;
            std::vector<dihedral_angle> dihedral_angles;
            std::vector<dihedral_angle> improper_angles;
            
            std::vector<bond> bonds;

            std::vector<std::vector<uint32_t>> rings; // for drawing on non-reactive mode;

            // CellList
            std::vector<std::vector<uint32_t>> cells;

            glm::vec3 box{20.f, 20.f, 20.f};
            uint32_t cx = 0, cy = 0, cz = 0; // cell dimensions

            void buildCells();
            int32_t getCellID(int32_t ix, int32_t iy, int32_t iz)
            {
                ix = (ix % (int32_t)cx + (int32_t)cx) % (int32_t)cx;
                iy = (iy % (int32_t)cy + (int32_t)cy) % (int32_t)cy;
                iz = (iz % (int32_t)cz + (int32_t)cz) % (int32_t)cz;
                return static_cast<int32_t>(ix + cx * (iy + cy * iz));
            }

            std::atomic<float> total_virial{0.0f};

            float temp = 0;
            float pres = 0;
            size_t timeStep = 0;

            inline bool areNearNeighbours(uint32_t i, uint32_t j, uint32_t n) const
            {
                if (i >= bondedBits.size() || j >= bondedBits.size()) return false;
                if (n == 0) return i == j;
                if (n == 1) return areBonded(i, j);

                std::queue<std::pair<uint32_t,uint32_t>> q;
                std::vector<bool> visited(bondedBits.size(), false);

                q.push({i, 0});
                visited[i] = true;

                while (!q.empty()) {
                    auto [current, depth] = q.front();
                    q.pop();

                    if (depth == n) {
                        if (current == j) return true;
                        continue;
                    }

                    for (uint32_t k = 0; k < bondedBits.size(); ++k) 
                    {
                        if (!visited[k] && areBonded(current, k)) {
                            visited[k] = true;
                            q.push({k, depth + 1});
                        }
                    }
                }

                return false;
            }

            void rebuildBondTopology()
            {
                int32_t N = atoms.size();

                for (const auto& b : bonds)
                {
                    markBonded(b.centralAtom, b.bondedAtom);
                }
            }

            uint32_t getBond(uint32_t i, uint32_t j)
            {
                for (int32_t b = 0; b < bonds.size(); ++b)
                {
                    const bond& bond = bonds[b];
                    if ((bond.bondedAtom == i && bond.centralAtom == j) || (bond.bondedAtom == j && bond.centralAtom == i))
                    {
                        return b;
                    }
                }

                return UINT32_MAX;
            }
            
            // Flags

            bool m_paused = false;

            bool gravity = false;
            bool isothermal = true;
            bool wall_collision = false;
            bool HMassRepartitioning = true;
            logging_flags log_flags;

            float mag_gravity = 9.8f;
            float m_Timescale = 1.0f;

            // Visual

            std::vector<uint32_t> m_highlightedAtoms{};
            std::vector<std::pair<uint32_t, uint32_t>> m_highlightedBonds{};
            std::vector<std::pair<uint32_t, uint32_t>> m_Arrows{};

            // Logging

            std::vector<frame> m_frames;

            // Other
            std::string moleculeName(const std::vector<uint32_t>& subsetIdx);
            void drawBox(core::window_t& window);

            inline float& bo(uint32_t i, uint32_t j) 
            {
                if (i > j) std::swap(i, j);
                return data.bond_orders[i * atoms.size() + j];
            }
            inline float bo(uint32_t i, uint32_t j) const 
            {
                if (i > j) std::swap(i, j);
                return data.bond_orders[i * atoms.size() + j];
            }
            inline sf::Vector3f pos(uint32_t i) const   { return {data.positions[i].x, data.positions[i].y, data.positions[i].z}; }
            inline sf::Vector3f vel(uint32_t i) const   { return {data.velocities[i].x, data.velocities[i].y, data.velocities[i].z}; }
            inline sf::Vector3f force(uint32_t i) const { return {data.forces[i].x, data.forces[i].y, data.forces[i].z}; }
            inline void add_force(uint32_t i, sf::Vector3f f) { data.forces[i].x += f.x, data.forces[i].y += f.y, data.forces[i].z += f.z; }
            inline void add_pos(uint32_t i, sf::Vector3f p) { data.positions[i].x += p.x, data.positions[i].y += p.y, data.positions[i].z += p.z; }
            inline void add_vel(uint32_t i, sf::Vector3f v) { data.velocities[i].x += v.x, data.velocities[i].y += v.y, data.velocities[i].z += v.z; }

            inline void emplace_vel(glm::vec3 v) { data.velocities.emplace_back(v); }
            inline void emplace_pos(glm::vec3 p) { data.positions.emplace_back(p); }

            void initReaxParams();
        };
    } // namespace fun
} // namespace sim

