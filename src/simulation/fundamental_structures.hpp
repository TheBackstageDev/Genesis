#pragma once

#include <SFML/Graphics.hpp>

#include "core/window.hpp"
#include <vector>
#include <map>
#include <set>

#include <chrono>
#include <algorithm>

namespace sim
{
    namespace fun { enum class BondType { SINGLE, DOUBLE, TRIPLE }; };

    namespace constants
    {
        #define MASS_PROTON 1.007 // Daltons
        #define MASS_NEUTRON 1.008 // Daltons
        #define MASS_ELECTRON 1/1337 // Daltons

        #define EPSILON 0.1f 
        #define DT 0.01f // ps
        #define MULT_FACTOR 5.f
        #define ANGSTROM 1e20f
        #define PICOSECOND 1e24f
        
        #define VERLET_SKIN 2.f
        #define CUTOFF 2.5f
        #define COULOMB_CUTOFF 10.f * MULT_FACTOR
        #define ELECTROSTATIC_CUTOFF 10.f * MULT_FACTOR

        #define AVOGADRO 6.02214076e26f // conversion from Daltons to Kg
        #define BOLTZMAN_CONSTANT 1.380649e-23f // Boltzman Constant m^2 kg s^-2 K^-1
        #define KB (BOLTZMAN_CONSTANT * AVOGADRO * ANGSTROM) / PICOSECOND // A^2 D ps^-2 K^-1
        #define THERMOSTAT_INTERVAL 5

        #define COULOMB_K 1389.3546f // kJ·mol⁻¹· Å * 1/5 ·e⁻²
        #define BOND_K 1500.f // Harmonic force constant 
        #define ANGLE_K 10000.f // kJ/mol/rad² for angular potential
        #define BOND_LENGTH_FACTOR 1.0f

        #define M_PI 3.141592653589
        #define RADIAN M_PI / 180
        #define DEGREE 180 / M_PI

        inline std::pair<float, float> getAtomConstants(uint32_t ZIndex)
        {
            std::pair<float, float> constants; // Sigma (Å), Epsilon (kJ/mol)
            switch (ZIndex)
            {
                case 1:  // H
                    constants = {2.50f, 0.1255f}; break;
                case 2:  // He
                    constants = {2.58f, 0.0870f}; break;
                case 3:  // Li
                    constants = {1.82f, 0.1100f}; break;
                case 4:  // Be
                    constants = {2.75f, 0.2510f}; break;
                case 5:  // B
                    constants = {3.64f, 0.3347f}; break;
                case 6:  // C
                    constants = {3.55f, 0.2929f}; break;
                case 7:  // N
                    constants = {3.25f, 0.7113f}; break;
                case 8:  // O
                    constants = {3.03f, 0.7113f}; break;
                case 9:  // F
                    constants = {2.95f, 0.6276f}; break;
                case 10: // Ne
                    constants = {2.78f, 0.1632f}; break;
                case 11: // Na
                    constants = {2.43f, 0.4184f}; break;
                case 12: // Mg
                    constants = {3.43f, 0.4184f}; break;
                case 13: // Al
                    constants = {4.01f, 0.5021f}; break;
                case 14: // Si
                    constants = {3.83f, 0.8368f}; break;
                case 15: // P
                    constants = {3.74f, 0.8368f}; break;
                case 16: // S
                    constants = {3.56f, 1.0460f}; break;
                case 17: // Cl
                    constants = {3.47f, 1.2552f}; break;
                default:
                    constants = {2.50f, 0.1255f}; // Default to H
            }
            constants.first *= MULT_FACTOR;
            constants.second *= MULT_FACTOR;
            return constants;
        }

        inline const std::map<uint8_t, float> electronegativity = 
        {
            {1, 2.2f},  // H
            {2, 0.0f},  // He (not typically used, set to 0)
            {3, 1.0f},  // Li
            {4, 1.5f},  // Be
            {5, 2.0f},  // B
            {6, 2.5f},  // C
            {7, 3.0f},  // N
            {8, 3.5f},  // O
            {9, 4.0f},  // F
            {10, 0.0f}, // Ne (not typically used, set to 0)
            {11, 0.9f}, // Na
            {12, 1.2f}, // Mg
            {13, 1.5f}, // Al
            {14, 1.8f}, // Si
            {15, 2.1f}, // P
            {16, 2.5f}, // S
            {17, 3.0f}  // Cl
        };

        // Temporary
        inline uint8_t getValenceElectrons(uint8_t ZIndex)
        {
            switch (ZIndex)
            {
                case 1:  return 1;  // H
                case 2:  return 2;  // He
                case 3:  return 1;  // Li
                case 4:  return 2;  // Be
                case 5:  return 3;  // B
                case 6:  return 4;  // C
                case 7:  return 5;  // N
                case 8:  return 6;  // O
                case 9:  return 7;  // F
                case 10: return 8;  // Ne
                case 11: return 1;  // Na
                case 12: return 2;  // Mg
                case 13: return 3;  // Al
                case 14: return 4;  // Si
                case 15: return 5;  // P
                case 16: return 6;  // S
                case 17: return 7;  // Cl
                default: return 0;
            }
        }

        inline float getBondLength(uint8_t ZIndex1, uint8_t ZIndex2, fun::BondType type)
        {
            float base = (getAtomConstants(ZIndex1).first + getAtomConstants(ZIndex2).first) / 2.0f * BOND_LENGTH_FACTOR;
            const float baseBase = base;

            switch (type)
            {
                case fun::BondType::SINGLE:
                    if ((ZIndex1 == 6 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 6)) base = 1.29f; // C-H 
                    if ((ZIndex1 == 8 && ZIndex2 == 6) || (ZIndex1 == 6 && ZIndex2 == 8)) base = 1.12f; // C-O 
                    if ((ZIndex1 == 6 && ZIndex2 == 6)) base = 1.74f; // C-C 
                    if ((ZIndex1 == 6 && ZIndex2 == 7) || (ZIndex1 == 7 && ZIndex2 == 6)) base = 1.47f; // C-N 
                    if ((ZIndex1 == 6 && ZIndex2 == 16) || (ZIndex1 == 16 && ZIndex2 == 6)) base = 1.81f; // C-S
                    if ((ZIndex1 == 8 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 8)) base = 0.96f; // O-H 
                    if ((ZIndex1 == 16 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 16)) base = 1.34f; // S-H 
                    if ((ZIndex1 == 17 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 17)) base = 1.27f; // Cl-H 
                    if ((ZIndex1 == 7 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 7)) base = 1.01f; // N-H 
                    if ((ZIndex1 == 7 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 7)) base = 1.40f; // N-O 
                    if ((ZIndex1 == 9 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 9)) base = 0.93f; // F-H
                    break;
                case fun::BondType::DOUBLE:
                    base *= 0.8f;
                    if ((ZIndex1 == 6 && ZIndex2 == 6)) base = 1.34f; // C=C 
                    if ((ZIndex1 == 6 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 6)) base = 1.21f; // C=O
                    if ((ZIndex1 == 6 && ZIndex2 == 7) || (ZIndex1 == 7 && ZIndex2 == 6)) base = 1.38f; // C=N
                    if ((ZIndex1 == 6 && ZIndex2 == 16) || (ZIndex1 == 16 && ZIndex2 == 6)) base = 1.71f; // C=S
                    if ((ZIndex1 == 7 && ZIndex2 == 7)) base = 1.20f; // N=N 
                    if ((ZIndex1 == 16 && ZIndex2 == 16)) base = 2.05f; // S=S 
                    break;
                case fun::BondType::TRIPLE:
                    base *= 0.7f;
                    if ((ZIndex1 == 6 && ZIndex2 == 6)) base = 1.20f; // C≡C 
                    if ((ZIndex1 == 6 && ZIndex2 == 7) || (ZIndex1 == 7 && ZIndex2 == 6)) base = 1.16f; // C≡N 
                    if ((ZIndex1 == 7 && ZIndex2 == 7)) base = 1.10f; // N≡N
                    break;
                default:
                    break;
            }

            if (baseBase != base) // Has passed through one of the IF's
            {
                base *= MULT_FACTOR;
            }

            return base;
        }

        inline float getAngles(uint8_t centralZIndex, const std::vector<uint8_t>& neighborZs, const std::vector<fun::BondType>& types)
        {
            size_t bond_count = neighborZs.size();
            if (bond_count < 2) return 0.0f; // No angle if fewer than 2 neighbors

            uint32_t valence = getValenceElectrons(centralZIndex);
            uint32_t bonding_electrons = 0;
            for (const auto& type : types)
            {
                bonding_electrons += (type == fun::BondType::SINGLE ? 2 : (type == fun::BondType::DOUBLE ? 4 : 6));
            }
            uint32_t lone_pairs = (valence - (bonding_electrons / 2)) / 2; // Approximate lone pairs
            uint32_t total_pairs = bond_count + lone_pairs;
            float ideal_angle = 0.0f;

            switch (total_pairs)
            {
                case 2: ideal_angle = M_PI; break; // Linear
                case 3: ideal_angle = 120.0f * RADIAN; break; // Trigonal planar
                case 4:
                    ideal_angle = 109.5f * RADIAN; // Tetrahedral base
                    if (lone_pairs == 1 && bond_count == 3) ideal_angle = 107.0f * RADIAN; 
                    if (lone_pairs == 2 && bond_count == 2) ideal_angle = 104.5f * RADIAN; 
                    break;
                case 5: ideal_angle = 120.0f * RADIAN; break; // Trigonal bipyramidal (equatorial)
                case 6: ideal_angle = 90.0f * RADIAN; break; // Octahedral
            }

            return ideal_angle;
        }

        inline std::string getAtomLetter(uint32_t ZIndex)
        {
            std::string name = "H";
            switch (ZIndex)
            {
                case 1:  name = "H";  break;
                case 2:  name = "He"; break;
                case 3:  name = "Li"; break;
                case 4:  name = "Be"; break;
                case 5:  name = "B";  break;
                case 6:  name = "C";  break;
                case 7:  name = "N";  break;
                case 8:  name = "O";  break;
                case 9:  name = "F";  break;
                case 10: name = "Ne"; break;
                case 11: name = "Na"; break;
                case 12: name = "Mg"; break;
                case 13: name = "Al"; break;
                case 14: name = "Si"; break;
                case 15: name = "P";  break;
                case 16: name = "S";  break;
                case 17: name = "Cl"; break;
                default: name = "H";  break;
            }
            return name;
        }
    };

    namespace fun
    {
        struct bond {
            uint32_t bondedAtom; 
            uint32_t centralAtom; 
            
            double equilibriumLength = 0.f;
            BondType type;
        };

        struct subset
        {
            size_t mainAtomIdx;              
            size_t bondedSubsetIdx = SIZE_MAX;       // Index of the last subset (optional, max(size_t) if none)
            size_t bondingSubsetIdx = SIZE_MAX;      // Index of the next subset (optional, max(size_t) if none)
            std::vector<size_t> connectedIdx; // connected to the main atom
            std::vector<BondType> bondTypes;
            float idealAngle;                // In radians
        };

        struct atom
        {
            float mass;            
            float radius;

            float epsilon; // LJ 
            float sigma; // LJ

            float charge;
            int8_t electrons;

            uint8_t ZIndex;
            uint8_t NCount; // neutrons
            int8_t bondCount;

            void draw(sf::Vector2f& pos, core::window_t &window, bool letterMode);
        };

        // Molecule-Structure related

        struct def_atom
        {
            uint8_t ZIndex;
            uint8_t NIndex;
            int8_t charge = 0; // non-Zero for ions 
        };

        struct def_bond {
            size_t bondingAtomIdx;   
            size_t centralAtomIdx;  
        };

        struct def_subset {
            size_t mainAtomIdx;      
            size_t bondedSubset;
            size_t bondingSubset;
            std::vector<size_t> connectedIdx; 
        };

        struct molecule_structure 
        {
            std::vector<def_atom> atoms;
            std::vector<def_subset> subsets;
            std::vector<def_bond> bonds;
            std::vector<BondType> bondTypes;
        };

        class universe
        {
        public:
            universe(float universeSize = 10.f);
            ~universe() = default;

            size_t createAtom(sf::Vector2f p, sf::Vector2f v, uint8_t ZIndex = 1, uint8_t numNeutrons = 0, uint8_t numElectrons = 1);
            size_t createSubset(const size_t central, const size_t subsetNext, const size_t subsetLast, const std::vector<size_t>& bonds, const std::vector<fun::BondType>& bondTypes);
            void createMolecule(const molecule_structure& structure, sf::Vector2f pos);

            void createBond(size_t idx1, size_t idx2, BondType type = BondType::SINGLE);
            void balanceMolecularCharges(subset& mol);

            void linkSubset(size_t subset, size_t subset2) { subsets[subset].bondingSubsetIdx = subset2; }

            void update(float targetTemperature = 1.0f);
            void draw(core::window_t &window, bool letter = false);
            void drawDebug(core::window_t& window);

            size_t numAtoms() { return atoms.size(); }
            const subset& getSubset(size_t index) { return subsets[index]; }

            float temperature() { return temp; }
            float timestep() { return timeStep; }

            void log(size_t step = 100);
        private:
            void boundCheck(size_t i);

            float ljPot(size_t i, float epsilon, float sigma);
            sf::Vector2f ljGrad(size_t i);
            sf::Vector2f ljForce(size_t i, size_t j);
            sf::Vector2f coulombForce(size_t i, size_t j, sf::Vector2f& dr_vec);

            void calcBondForces();
            void calcAngleForces();
            void calcLjForces();
            void calcElectrostaticForces();
            
            // Molecule Creation
            void organizeMolecule(const molecule_structure& structure, const sf::Vector2f& initPos);
            void positionMolecule(size_t firstSubsetIndex);

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

            bool areOnSameSubset(size_t i, size_t j)
            {
                for (const auto& subset : subsets)
                {
                    bool i_onSub = std::find(subset.connectedIdx.begin(), subset.connectedIdx.end(), i) != subset.connectedIdx.end();
                    bool j_onSub = std::find(subset.connectedIdx.begin(), subset.connectedIdx.end(), j) != subset.connectedIdx.end();
                    if (i_onSub && j_onSub)
                        return true;

                    if (i_onSub && subset.mainAtomIdx == j || j_onSub && subset.mainAtomIdx == i)
                        return true;
                }

                return false;
            }

            float boxSize = 10.f;
            std::vector<atom> atoms;
            std::vector<subset> subsets;

            std::vector<bond> bonds;

            std::vector<sf::Vector2f> forces;
            std::vector<sf::Vector2f> positions;
            std::vector<sf::Vector2f> velocities;
            std::vector<sf::Vector2f> prev_positions;

            float temp = 0;
            size_t timeStep = 0;

            // Helper Funcs
            sf::Vector2f minImageVec(sf::Vector2f dr)
            {
                dr.x -= boxSize * std::round(dr.x / boxSize);
                dr.y -= boxSize * std::round(dr.y / boxSize);
                return dr;
            }
        };
    } // namespace fun
} // namespace sim
