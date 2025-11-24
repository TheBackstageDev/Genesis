#pragma once

#include "fundamental_structures.hpp"
#include <SFML/Graphics.hpp>

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
            alignas(32) std::vector<float> old_fx, old_fy, old_fz;
            alignas(32) std::vector<float> q, mass;

            alignas(32) std::vector<float> x0, y0, z0;
        };

        constexpr int32_t offsets[14][3] = 
        {
            {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {-1,1,0},
            {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}, {-1,0,1},
            {1,-1,1}, {0,-1,1}, {-1,-1,1}, {-1,0,0}
        };

        class universe
        {
        public:
            universe(float universeSize = 10.f, float cell_size = CELL_CUTOFF);
            universe(std::filesystem::path scene);
            ~universe() = default;

            size_t createAtom(sf::Vector3f p, sf::Vector3f v, uint8_t ZIndex = 1, uint8_t numNeutrons = 0, uint8_t numElectrons = 1, int32_t chirality = 0);
            size_t createSubset(const def_subset& nSub, const size_t baseAtom, const size_t baseSubset);
            void createMolecule(molecule_structure& structure, sf::Vector3f pos, sf::Vector3f vel = {0.f, 0.f, 0.f});

            void createBond(size_t idx1, size_t idx2, BondType type = BondType::SINGLE);
            void balanceMolecularCharges(subset& mol);

            void linkSubset(size_t subset, size_t subset2) { subsets[subset].bondingSubsetIdx = subset2; }

            void update(float targetTemperature = 1.0f, bool reactions = true);
            void draw(core::window_t &window, bool letter = false, bool lennardBall = true);
            void drawHydrogenBond(core::window_t& window, size_t H);
            void drawBonds(core::window_t& window);
            void drawDebug(core::window_t& window);

            void saveScene(const std::filesystem::path path);
            void loadScene(const std::filesystem::path path);

            size_t numAtoms() { return atoms.size(); }
            const subset& getSubset(size_t index) { return subsets[index]; }

            float temperature() { return temp; }
            float timestep() { return timeStep; }

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

            void log(size_t step = 100);
            void handleCamera(bool leftDown, bool rightDown, const sf::Vector2i& mousePos, float wheelDelta, const std::vector<sf::Keyboard::Key>& keys);
        private:
            void boundCheck(size_t i);

            float ljPot(size_t i, float epsilon, float sigma);
            float wolfForce(float r, float qi_qj); 
            sf::Vector3f ljGrad(size_t i);
            sf::Vector3f ljForce(size_t i, size_t j);
            sf::Vector3f coulombForce(size_t i, size_t j, sf::Vector3f& dr_vec);

            void calcBondForces();
            void calcAngleForces();
            void calcDihedralForces();
            void calcLjForces();
            void calcElectrostaticForces();
            void calcBondedForces();
            void calcUnbondedForces();

            void setTemperature(float kelvin = 0.f);
            float calculateDihedral(const sf::Vector3f& pa, const sf::Vector3f& pb, const sf::Vector3f& pc, const sf::Vector3f& pd);

            // Energies
            float calculateKineticEnergy();
            float calculateAtomTemperature(size_t i);
            float calculateBondEnergy(size_t i, size_t j);
            float calculateNonBondedEnergy(size_t i, size_t j);

            // Reactions
            void handleReactions();
            void checkBonds(); // for either forming or breaking
            void breakBond(size_t atom1, size_t atom2);

            bool isReactive(size_t i) { return atoms[i].bondCount < constants::getUsualBonds(atoms[i].ZIndex) || atoms[i].bondCount > constants::getUsualBonds(atoms[i].ZIndex); }

            BondType chooseBestBondType(size_t i, size_t j);

            float boxSize = 10.f;

            simData data;
            std::vector<std::vector<uint64_t>> bondedBits;
            void markBonded(size_t i, size_t j)
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
            std::vector<bond> bonds;
            std::vector<subset> subsets;
            std::vector<molecule> molecules;

            std::vector<angle> angles;
            std::vector<dihedral_angle> dihedral_angles;
            
            // CellList
            std::vector<std::vector<size_t>> cells;

            float cell_size = 0.f;
            size_t cx = 0, cy = 0, cz = 0; // cell dimensions

            void buildCells();
            size_t getCellID(int32_t ix, int32_t iy, int32_t iz)
            {
                ix = (ix % (int32_t)cx + (int32_t)cx) % (int32_t)cx;
                iy = (iy % (int32_t)cy + (int32_t)cy) % (int32_t)cy;
                iz = (iz % (int32_t)cz + (int32_t)cz) % (int32_t)cz;
                return static_cast<size_t>(ix + cx * (iy + cy * iz));
            }

            float temp = 0;
            float neighbourInterval = 0.f;
            size_t timeStep = 0;

            // Helper Funcs
            sf::Vector3f minImageVec(sf::Vector3f dr)
            {
                dr.x -= boxSize * std::round(dr.x / boxSize);
                dr.y -= boxSize * std::round(dr.y / boxSize);
                dr.z -= boxSize * std::round(dr.z / boxSize);
                return dr;
            }

            std::vector<size_t> getAllSubsetAtoms(size_t s)
            {
                const subset& sub = subsets[s];
                std::vector<size_t> indices;
                indices.reserve(sub.connectedIdx.size() + sub.hydrogenIdx.size() + 1);
                indices.push_back(sub.mainAtomIdx);
                indices.insert(indices.end(), sub.connectedIdx.begin(), sub.connectedIdx.end());
                indices.insert(indices.end(), sub.hydrogenIdx.begin(), sub.hydrogenIdx.end());
                return indices;
            }

            inline bool areBonded(size_t i, size_t j) const
            {
                if (i >= bondedBits.size() || j >= bondedBits[i].size() * 64) return false;
                size_t word = j / 64;
                size_t bit  = j % 64;
                return (bondedBits[i][word] & (1ull << bit)) != 0;
            }
            void rebuildBondTopology()
            {
                size_t N = atoms.size();

                for (const auto& b : bonds)
                {
                    markBonded(b.centralAtom, b.bondedAtom);
                }
            }

            size_t getBond(size_t i, size_t j)
            {
                for (size_t b = 0; b < bonds.size(); ++b)
                {
                    const bond& bond = bonds[b];
                    if ((bond.bondedAtom == i && bond.centralAtom == j) || (bond.bondedAtom == j && bond.centralAtom == i))
                    {
                        return b;
                    }
                }

                return SIZE_MAX;
            }
            
            // Camera
            sf::Vector2i lastMouse;
            sf::Vector2f project(core::window_t& window, const sf::Vector3f& p) const;
            core::camera_t cam;

            // Other
            std::string moleculeName(const std::vector<size_t>& subsetIdx);
            void drawBox(core::window_t& window);

            inline sf::Vector3f neigh_pos(size_t i) const   { return {data.x0[i], data.y0[i], data.z0[i]}; }
            inline sf::Vector3f pos(size_t i) const   { return {data.x[i], data.y[i], data.z[i]}; }
            inline sf::Vector3f vel(size_t i) const   { return {data.vx[i], data.vy[i], data.vz[i]}; }
            inline sf::Vector3f old_force(size_t i) const { return {data.old_fx[i], data.old_fy[i], data.old_fz[i]}; }
            inline sf::Vector3f force(size_t i) const { return {data.fx[i], data.fy[i], data.fz[i]}; }
            inline void add_force(size_t i, sf::Vector3f f) { data.fx[i] += f.x, data.fy[i] += f.y, data.fz[i] += f.z; }
            inline void add_pos(size_t i, sf::Vector3f p) { data.x[i] += p.x, data.y[i] += p.y, data.z[i] += p.z; }
            inline void add_vel(size_t i, sf::Vector3f v) { data.vx[i] += v.x, data.vy[i] += v.y, data.vz[i] += v.z; }

            inline void emplace_vel(sf::Vector3f v) { data.vx.emplace_back(v.x); data.vy.emplace_back(v.y); data.vz.emplace_back(v.z); }
            inline void emplace_pos(sf::Vector3f p) { data.x.emplace_back(p.x); data.y.emplace_back(p.y); data.z.emplace_back(p.z); }
        };
    } // namespace fun
} // namespace sim
