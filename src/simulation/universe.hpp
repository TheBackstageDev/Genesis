#pragma once

#include "fundamental_structures.hpp"
#include "rendering_engine.hpp"
#include <SFML/Graphics.hpp>

#include <thread>
#include <future>
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
            alignas(32) std::vector<glm::vec3> positions;
            alignas(32) std::vector<glm::vec3> velocities;
            alignas(32) std::vector<glm::vec3> forces;
            alignas(32) std::vector<float> q, bond_orders, temperature;
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
            bool reactive = false;
            bool wall_collision = false;
            bool isothermal = true;
            bool render_water = true;
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
            std::vector<glm::vec3> positions;
            std::map<size_t, float> temperatures;

            float global_temperature;
        };

        struct video
        {
            std::vector<frame> frames; // index is defined by order
            std::vector<size_t> keyFrames; // for events

            std::vector<std::string> text; // changes every keyframe (or continue);

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
            void draw(sf::RenderTarget& target, const rendering_info& info);

            void runVideo(const video& vid);

            // Scenario stuff
            void saveScene(const std::filesystem::path path, const std::string name = "");
            void loadScene(const std::filesystem::path path);

            void saveAsVideo(const std::filesystem::path path, const std::string name = "");
            void saveFrame();

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
            }

            // Gets

            int32_t numBonds() { return bonds.size(); }
            int32_t numAtoms() { return atoms.size(); }
            int32_t numMolecules() { return molecules.size(); }
            const subset& getSubset(int32_t index) { return subsets[index]; }

            float temperature() const { return temp; }
            float pressure() const { return pres; }
            float timestep() const { return timeStep; }

            glm::vec3 boxSizes() { return box; }

            _NODISCARD std::vector<glm::vec3> positions() const 
            {
                return data.positions;
            }

            void setPosition(size_t i, glm::vec3 p)
            {
                data.positions[i] = p;
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
        private:
            void boundCheck(uint32_t i);

            float ljPot(uint32_t i, uint32_t j);
            float wolfForce(float r, float qi_qj); 
            sf::Vector3f ljForce(uint32_t i, uint32_t j);
            sf::Vector3f coulombForce(uint32_t i, uint32_t j, sf::Vector3f& dr_vec);

            void calcBondForces();
            void calcAngleForces();
            void calcDihedralForces();

            void calcReactiveAngleForces();
            void calcLjForces();
            void calcElectrostaticForces();

            std::vector<sf::Vector3f> processCellUnbonded(int32_t ix, int32_t iy, int32_t iz);
            std::vector<sf::Vector3f> processPartialCellUnbonded(int32_t ix, int32_t iy, int32_t iz, int32_t atom_start, int32_t atom_end);
            void calcUnbondedForcesParallel();
            void calcBondedForcesParallel();

            void calcBondedForces(); // deprecated
            void calcUnbondedForces(); // deprecated

            float calculatePressure();
            void setPressure(float bar = 100);
            void setTemperature(float kelvin = 0.f);
            float calculateDihedral(const glm::vec3& pa, const glm::vec3& pb, const glm::vec3& pc, const glm::vec3& pd);

            // Energies
            float calculateKineticEnergy();
            float calculateAtomTemperature(int32_t i);
            float calculateBondEnergy(int32_t i, int32_t j, float bo_sigma, float bo_pi, float bo_pp);

            // Reactions
            float calculateUncorrectedBondOrder(int32_t i, int32_t j);
            void processReactivePair(int32_t i, int32_t j, float cutoff = CELL_CUTOFF, float vis_thresh = 0.3f);
            void calcUnbondedForcesReactive();
            void handleReactiveForces();

            rendering_engine& rendering_eng;

            simData data;
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
            
            std::vector<bond> bonds;
            std::vector<reactive_bond> reactive_bonds;

            std::vector<std::vector<uint32_t>> rings; // for drawing on non-reactive mode;

            // CellList
            std::vector<std::vector<uint32_t>> cells;

            glm::vec3 box{0.f, 0.f, 0.f};
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

            inline bool areBonded(uint32_t i, uint32_t j) const
            {
                if (i >= bondedBits.size() || j >= bondedBits[i].size() * 64) return false;
                size_t word = j / 64;
                size_t bit  = j % 64;
                return (bondedBits[i][word] & (1ull << bit)) != 0;
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

            bool gravity = false;
            bool react = false;
            bool isothermal = true;
            bool render_water = true;
            bool wall_collision = false;
            bool HMassRepartitioning = true;
            logging_flags log_flags;

            float mag_gravity = 9.8f;

            // Logging
            struct ReactionEvent 
            {
                enum Type { BOND_FORM, BOND_BREAK, PROTON_TRANSFER } type;
                uint32_t atom1, atom2;
                float old_bo, new_bo;
                float time;
            };

            std::vector<ReactionEvent> reactionLog{};
            std::vector<std::vector<float>> positionxLog{};
            std::vector<std::vector<float>> positionyLog{};
            std::vector<std::vector<float>> positionzLog{};
            std::map<uint32_t, float> energyLog{};
            std::vector<float> temperatureLog{};
            
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

