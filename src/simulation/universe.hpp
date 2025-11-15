#pragma once

#include "fundamental_structures.hpp"
#include <SFML/Graphics.hpp>

namespace sim
{
    namespace fun
    {
        class universe
        {
        public:
            universe(float universeSize = 10.f);
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

            void saveScene();
            void loadScene(const std::string& path);

            std::vector<sf::Vector3f> getPositions() const { return positions; }

            size_t numAtoms() { return atoms.size(); }
            const subset& getSubset(size_t index) { return subsets[index]; }

            float temperature() { return temp; }
            float timestep() { return timeStep; }

            void log(size_t step = 100);
            void handleCamera(bool leftDown, bool rightDown, const sf::Vector2i& mousePos, float wheelDelta, const std::vector<sf::Keyboard::Key>& keys);
        private:
            void boundCheck(size_t i);

            float ljPot(size_t i, float epsilon, float sigma);
            sf::Vector3f ljGrad(size_t i);
            sf::Vector3f ljForce(size_t i, size_t j);
            sf::Vector3f coulombForce(size_t i, size_t j, sf::Vector3f& dr_vec);

            void calcBondForces();
            void calcAngleForces();
            void calcDihedralForces();
            void calcLjForces();
            void calcElectrostaticForces();

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
            void checkStrainedBonds();
            void breakBond(size_t atom1, size_t atom2);

            bool isReactive(size_t i) { return atoms[i].bondCount < constants::getUsualBonds(atoms[i].ZIndex) || atoms[i].bondCount > constants::getUsualBonds(atoms[i].ZIndex); }

            BondType chooseBestBondType(size_t i, size_t j);

            float boxSize = 10.f;
            std::vector<atom> atoms;
            std::vector<bond> bonds;
            std::vector<subset> subsets;
            std::vector<molecule> molecules;

            std::vector<angle> angles;
            std::vector<dihedral_angle> dihedral_angles;

            std::vector<sf::Vector3f> forces;
            std::vector<sf::Vector3f> positions;
            std::vector<sf::Vector3f> neighbourPos;
            std::vector<sf::Vector3f> velocities;

            std::vector<std::vector<size_t>> neighbourList; 
            void buildNeighborList();

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

            bool neighbourListNeedRebuild();
            bool areBonded(size_t i, size_t j) 
            {
                for (const auto& bond : bonds)
                {
                    if ((bond.bondedAtom == i && bond.centralAtom == j) || (bond.bondedAtom == j && bond.centralAtom == i))
                    {
                        return true;
                    }
                }

                return false;
            }
            
            // Camera
            sf::Vector2i lastMouse;
            sf::Vector2f project(core::window_t& window, const sf::Vector3f& p) const;
            core::camera_t cam;

            // Other
            std::string moleculeName(const std::vector<size_t>& subsetIdx);
            void drawBox(core::window_t& window);
        };
    } // namespace fun
} // namespace sim
