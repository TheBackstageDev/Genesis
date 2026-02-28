#pragma once

#include "fundamental_structures.hpp"
#include "rendering_engine.hpp"
#include "core/shader.hpp"
#include "core/buffer.hpp"

#include "core/spatialgrid.hpp"

#include <SFML/Graphics.hpp>

#include <thread>
#include <future>
#include <queue>
#include <atomic>
#include <filesystem>
#include <json.hpp>
#include <random>

#ifndef _NODISCARD
#define _NODISCARD
#endif

//#define CALCULATIONS_GPU

namespace sim
{
    namespace fun
    {
        struct simData
        {
            alignas(64) std::vector<glm::vec3> positions;
            alignas(64) std::vector<glm::vec3> velocities;
            alignas(64) std::vector<glm::vec3> forces;
            alignas(64) std::vector<float> q, lj_params;
        };

        struct atomData
        {
            alignas(64) std::vector<atom> atoms;
            alignas(64) std::vector<subset> subsets;
            alignas(64) std::vector<molecule> molecules;
            alignas(64) std::vector<bool> frozen_atoms;
            
            alignas(64) std::vector<angle> angles;
            alignas(64) std::vector<dihedral_angle> dihedral_angles;
            alignas(64) std::vector<dihedral_angle> improper_angles;
            
            alignas(64) std::vector<bond> bonds;
        };

        struct logging_flags
        {
            bool log_reactions = true;
        };

        struct universe_create_info
        {
            bool has_gravity = false;
            bool wall_collision = false;
            bool roof_floor_collision = false;
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
            void createMolecule(molecule_structure structure, sf::Vector3f pos, sf::Vector3f vel =
                 {0.f, 0.f, 0.f});

            void createBond(int32_t idx1, int32_t idx2, BondType type = BondType::SINGLE);
            void balanceMolecularCharges(subset& mol);
            
            void draw(sf::RenderTarget& target, rendering_info info);

            // Scenario stuff
            void saveScene(const std::filesystem::path path, const std::string name = "");
            void loadScene(const std::filesystem::path path);

            void loadFrames(const std::vector<frame>& nFrames) { m_frames = std::move(nFrames); }
            void loadFrames(const std::filesystem::path path);
            void saveFrame();

            simData& getData() { return data; }
            atomData& getAtomData() { return atomData; }
            
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

                atomData.molecules.clear();
                atomData.atoms.clear();
                atomData.bonds.clear();
                atomData.subsets.clear();

                atomData.angles.clear();
                atomData.dihedral_angles.clear();
                atomData.improper_angles.clear();

                m_highlightedAtoms.clear();
                m_highlightedBonds.clear();
                m_Arrows.clear();
            }

            void clearArrows()
            {
                m_Arrows.clear();
            }

            void highlightAtom(uint32_t index) { m_highlightedAtoms.emplace_back(index); }
            void highlightBond(uint32_t index1, uint32_t index2) { m_highlightedBonds.emplace_back(index1, index2); }
            void createArrow(glm::vec3 from, glm::vec3 to) { m_Arrows.emplace_back(from, to); }

            // Sets

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

            void setMagneticFieldEnabled(bool state) { f_magneticField = state; }
            void setWallChargesEnabled(bool state) { f_wallCharges = state; }
            void setMagneticFieldStrength(glm::vec3 strength) { magnetic_strength = strength; }
            
            // Gets
            
            core::camera_t& getRenderingCamera() { return rendering_eng.camera(); }
            sim::rendering_engine& getRenderingEngine() { return rendering_eng; }
            bool& wallChargeEnabled() { return f_wallCharges; }
            bool& magneticFieldEnabled() { return f_magneticField; }

            bool wallcollision() { return wall_collision; }
            bool rooffloorcollision() { return roof_floor_collision; }
            
            const std::vector<atom>& getAtoms() const { return atomData.atoms; }
            std::vector<subset>& getSubsets() { return atomData.subsets; }
            int32_t numBonds() const { return atomData.bonds.size(); }
            int32_t numAtoms() const { return atomData.atoms.size(); }
            int32_t numMolecules() const { return atomData.molecules.size(); }

            std::array<float, 6>& getWallCharges() { return wall_charges; }
            glm::vec3& getMagneticFieldStrength() { return magnetic_strength; }

            glm::vec3 getForce(uint32_t i) { return data.forces[i]; }
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

            inline bool areBonded(uint32_t i, uint32_t j) const
            {
                if (i >= bondedBits.size() || j >= bondedBits[i].size() * 32) return false;
                uint32_t word = j / 32;
                uint32_t bit  = j % 32;
                return (bondedBits[i][word] & (1u << bit)) != 0;
            }

            // Energies
            float calculateKineticEnergy();

            sf::Vector3f minImageVec(sf::Vector3f dr)
            {
                bool wall_col = wallcollision();
                bool rf_col = rooffloorcollision();

                if (wall_col && rf_col) return dr;

                if (!rf_col)
                    dr.z -= box.z * std::round(dr.z / box.z);

                if (!wall_col)
                {
                    dr.x -= box.x * std::round(dr.x / box.x);
                    dr.y -= box.y * std::round(dr.y / box.y);
                }

                return dr;
            }

            glm::vec3 minImageVec(glm::vec3 dr)
            {
                bool wall_col = wallcollision();
                bool rf_col = rooffloorcollision();

                if (wall_col && rf_col) return dr;

                glm::vec3 box_sizes = boxSizes();

                if (!rf_col)
                    dr.z -= box_sizes.z * std::round(dr.z / box_sizes.z);

                if (!wall_col)
                {
                    dr.x -= box_sizes.x * std::round(dr.x / box_sizes.x);
                    dr.y -= box_sizes.y * std::round(dr.y / box_sizes.y);
                }

                return dr;
            }

            void boundCheck(uint32_t i);
        private:

            // Energies
            float calculateAtomTemperature(int32_t i);

            float gauss_random()
            {
                static std::random_device rd;
                static std::mt19937 gen(rd());
                static std::normal_distribution<float> dist(0.0f, 1.0f);

                return dist(gen);
            }

            rendering_engine& rendering_eng;

            simData data;
            std::vector<glm::vec3> m_displayPositions{};

            std::vector<std::vector<uint32_t>> bondedBits;
            void markBonded(uint32_t i, uint32_t j)
            {
                if (i >= bondedBits.size()) bondedBits.resize(atomData.atoms.size());
                if (j >= bondedBits.size()) bondedBits.resize(atomData.atoms.size());
                while (bondedBits[i].size() * 32 <= j) bondedBits[i].push_back(0);
                while (bondedBits[j].size() * 32 <= i) bondedBits[j].push_back(0);

                uint32_t word_i = j / 32;
                uint32_t bit_i  = j % 32;
                uint32_t word_j = i / 32;
                uint32_t bit_j  = i % 32;

                bondedBits[i][word_i] |= (1ull << bit_i);
                bondedBits[j][word_j] |= (1ull << bit_j);
            }

            atomData atomData{};

            glm::vec3 box{20.f, 20.f, 20.f};

            std::array<float, 6> wall_charges{0.f, 0.f, 0.f, 0.f, 0.f, 0.f}; // right, left, top, bottom, front, back
            glm::vec3 magnetic_strength{0.f, 0.f, 2.0f};

            bool f_wallCharges = false;
            bool f_magneticField = false;

            void rebuildBondTopology()
            {
                int32_t N = atomData.atoms.size();

                for (const auto& b : atomData.bonds)
                {
                    markBonded(b.centralAtom, b.bondedAtom);
                }
            }

            uint32_t getBond(uint32_t i, uint32_t j)
            {
                for (int32_t b = 0; b < atomData.bonds.size(); ++b)
                {
                    const bond& bond = atomData.bonds[b];
                    if ((bond.bondedAtom == i && bond.centralAtom == j) || (bond.bondedAtom == j && bond.centralAtom == i))
                    {
                        return b;
                    }
                }

                return UINT32_MAX;
            }

            void COMDrift();
            
            // Flags

            bool gravity = false;
            bool isothermal = true;
            bool wall_collision = false;
            bool roof_floor_collision = false;
            bool HMassRepartitioning = false;
            logging_flags log_flags;

            float mag_gravity = 9.8f;

            // Visual

            std::vector<uint32_t> m_highlightedAtoms{};
            std::vector<std::pair<uint32_t, uint32_t>> m_highlightedBonds{};
            std::vector<std::pair<glm::vec3, glm::vec3>> m_Arrows{};

            // Logging

            std::vector<frame> m_frames;

            // Other
            std::string moleculeName(const std::vector<uint32_t>& subsetIdx);
            void drawBox(core::window_t& window);

            inline void emplace_vel(glm::vec3 v) { data.velocities.emplace_back(v); }
            inline void emplace_pos(glm::vec3 p) { data.positions.emplace_back(p); }

            void initReaxParams();
        };
    } // namespace fun
} // namespace sim

