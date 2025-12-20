#pragma once

#include "fundamental_structures.hpp"
#include <SFML/Graphics.hpp>

#include <thread>
#include <future>
#include <atomic>
#include <filesystem>
#include <json.hpp>

namespace sim
{
    namespace fun
    {
        struct simData
        {
            alignas(32) std::vector<float> x, y, z;
            alignas(32) std::vector<float> vx, vy, vz;
            alignas(32) std::vector<float> fx, fy, fz;
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
            logging_flags log_flags;

            float mag_gravity = 9.8f;
            sf::Vector3f box{CELL_CUTOFF * 2.f, CELL_CUTOFF * 2.f, CELL_CUTOFF * 2.f};
        };

        class universe
        {
        public:
            universe(const universe_create_info& create_info);
            universe(std::filesystem::path scene);

            int32_t createAtom(sf::Vector3f p, sf::Vector3f v, uint8_t ZIndex = 1, uint8_t numNeutrons = 0, uint8_t numElectrons = 1, int32_t chirality = 0);
            int32_t createSubset(const def_subset& nSub, const int32_t baseAtom, const int32_t baseSubset);
            void createMolecule(molecule_structure structure, sf::Vector3f pos, sf::Vector3f vel = {0.f, 0.f, 0.f});

            void createBond(int32_t idx1, int32_t idx2, BondType type = BondType::SINGLE);
            void balanceMolecularCharges(subset& mol);

            void linkSubset(int32_t subset, int32_t subset2) { subsets[subset].bondingSubsetIdx = subset2; }

            void update(float targetTemperature = 1.0f, float targetPressure = 0.f);
            void draw(core::window_t &window, sf::RenderTarget& target, const rendering_info info);
            void drawHydrogenBond(core::window_t &window, sf::RenderTarget& target, int32_t H);
            void drawBonds(core::window_t &window, sf::RenderTarget& target, const std::vector<int32_t>& no_draw);
            void drawRings(core::window_t &window, sf::RenderTarget& target);
            void drawChargeField(core::window_t &window, sf::RenderTarget& target);

            void drawDebug(core::window_t& window);

            void insertGhostAtom(uint32_t g) { ghost_atoms.emplace_back(g); };
            void insertGhostAtom(std::vector<uint32_t>& g) { ghost_atoms.insert(ghost_atoms.end(), g.begin(), g.end()); }

            void removeGhostAtom(uint32_t g) 
            {
                auto it = std::find(ghost_atoms.begin(), ghost_atoms.end(), g);
                
                if (it != ghost_atoms.end())
                    ghost_atoms.erase(it); 
            }

            void removeGhostAtom(const std::vector<uint32_t>& g) 
            {
                for (const uint32_t& a : g)
                    removeGhostAtom(a);
            }

            // Scenario stuff
            void highlightAtom(core::window_t& window, size_t i);
            void highlightBond(core::window_t& window, size_t i, size_t j);

            void saveScene(const std::filesystem::path path, const std::string name = "");
            void loadScene(const std::filesystem::path path);

            void saveAsVideo(const std::filesystem::path path, const std::string name = "");
            void saveFrame();

            void clear()
            {
                data.x.clear();
                data.y.clear();
                data.z.clear();
                data.vx.clear();
                data.vy.clear();
                data.vz.clear();
                data.fx.clear();
                data.fy.clear();
                data.fz.clear();
                data.q.clear();
                data.temperature.clear();
                data.bond_orders.clear();

                molecules.clear();
                atoms.clear();
                bonds.clear();
                subsets.clear();
                rings.clear();
            }

            int32_t numBonds() { return bonds.size(); }
            int32_t numAtoms() { return atoms.size(); }
            int32_t numMolecules() { return molecules.size(); }
            const subset& getSubset(int32_t index) { return subsets[index]; }

            float temperature() const { return temp; }
            float pressure() const { return pres; }
            float timestep() const { return timeStep; }

            _NODISCARD std::vector<sf::Vector3f> positions() const 
            {
                std::vector<sf::Vector3f> positions;
                positions.reserve(data.x.size());

                for (int32_t i = 0; i < data.x.size(); ++i)
                {
                    positions.emplace_back(pos(i));
                }

                return positions;
            }

            void setPosition(size_t i, sf::Vector3f p)
            {
                data.x[i] = p.x;
                data.y[i] = p.y;
                data.z[i] = p.z;
            }

            void setRenderMode();

            // cam
            void handleCamera();
            core::camera_t& camera() { return cam; }

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
            void calcUnbondedForcesParallel();
            void calcBondedForcesParallel();

            void calcBondedForces(); // deprecated
            void calcUnbondedForces(); // deprecated

            float calculatePressure();
            void setPressure(float bar = 100);
            void setTemperature(float kelvin = 0.f);
            float calculateDihedral(const sf::Vector3f& pa, const sf::Vector3f& pb, const sf::Vector3f& pc, const sf::Vector3f& pd);

            // Energies
            float calculateKineticEnergy();
            float calculateAtomTemperature(int32_t i);
            float calculateBondEnergy(int32_t i, int32_t j, float bo_sigma, float bo_pi, float bo_pp);

            // Reactions
            float calculateUncorrectedBondOrder(int32_t i, int32_t j);
            void processReactivePair(int32_t i, int32_t j, float cutoff = CELL_CUTOFF, float vis_thresh = 0.3f);
            void calcUnbondedForcesReactive();
            void handleReactiveForces();

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
            std::vector<uint32_t> ghost_atoms; // for displaying
            
            std::vector<angle> angles;
            std::vector<dihedral_angle> dihedral_angles;
            
            std::vector<bond> bonds;
            std::vector<reactive_bond> reactive_bonds;

            std::vector<std::vector<uint32_t>> rings; // for drawing on non-reactive mode;

            // CellList
            std::vector<std::vector<uint32_t>> cells;

            sf::Vector3f box{0.f, 0.f, 0.f};
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

            // Camera
            ImVec2 lastMouse;
            sf::Vector2f project(core::window_t& window, const sf::Vector3f& p) const;
            core::camera_t cam;

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
            inline sf::Vector3f pos(uint32_t i) const   { return {data.x[i], data.y[i], data.z[i]}; }
            inline sf::Vector3f vel(uint32_t i) const   { return {data.vx[i], data.vy[i], data.vz[i]}; }
            inline sf::Vector3f force(uint32_t i) const { return {data.fx[i], data.fy[i], data.fz[i]}; }
            inline void add_force(uint32_t i, sf::Vector3f f) { data.fx[i] += f.x, data.fy[i] += f.y, data.fz[i] += f.z; }
            inline void add_pos(uint32_t i, sf::Vector3f p) { data.x[i] += p.x, data.y[i] += p.y, data.z[i] += p.z; }
            inline void add_vel(uint32_t i, sf::Vector3f v) { data.vx[i] += v.x, data.vy[i] += v.y, data.vz[i] += v.z; }

            inline void emplace_vel(sf::Vector3f v) { data.vx.emplace_back(v.x); data.vy.emplace_back(v.y); data.vz.emplace_back(v.z); }
            inline void emplace_pos(sf::Vector3f p) { data.x.emplace_back(p.x); data.y.emplace_back(p.y); data.z.emplace_back(p.z); }

            void initReaxParams();
        };
    } // namespace fun
} // namespace sim

