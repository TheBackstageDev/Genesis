#pragma once

#include <SFML/Graphics.hpp>

#include "core/window.hpp"
#include <vector>
#include <map>
#include <set>

#include <chrono>
#include <algorithm>

#include "smiles_parser.hpp"

namespace sim
{
    namespace constants
    {
        #define MASS_PROTON 1.0007 // Daltons
        #define MASS_NEUTRON 1.0008 // Daltons
        #define MASS_ELECTRON 1/1337 // Daltons

        #define EPSILON 0.1f 
        #define DT 0.01f // ps
        #define MULT_FACTOR 5.f
        #define ANGSTROM 1e20f
        #define PICOSECOND 1e24f
        
        #define VERLET_SKIN 2.f
        #define CUTOFF 2.5f
        #define COULOMB_CUTOFF 10.f * MULT_FACTOR

        #define AVOGADRO 6.02214076e26f // conversion from Daltons to Kg
        #define BOLTZMAN_CONSTANT 1.380649e-23f // Boltzman Constant m^2 kg s^-2 K^-1
        #define KB (BOLTZMAN_CONSTANT * AVOGADRO * ANGSTROM) / PICOSECOND // A^2 D ps^-2 K^-1
        #define THERMOSTAT_INTERVAL 5

        #define COULOMB_K 1389.3546f // kJ·mol⁻¹· Å * 1/5 ·e⁻²
        #define BOND_K 1000.f // Harmonic force constant 
        #define ANGLE_K 15000.f // kJ/mol/rad² for angular potential
        #define BOND_LENGTH_FACTOR 1.0f

        #define M_PI 3.1415926535
        #define RADIAN M_PI / 180
        #define DEGREE 180 / M_PI

        inline std::pair<float, float> getAtomConstants(uint32_t ZIndex)
        {
            std::pair<float, float> constants; // {sigma (Å), epsilon (kJ/mol)}
            switch (ZIndex)
            {
                case 1:  constants = {2.50f, 0.1255f}; break; // H
                case 2:  constants = {2.58f, 0.0870f}; break; // He
                case 3:  constants = {1.82f, 0.1100f}; break; // Li
                case 4:  constants = {2.75f, 0.2510f}; break; // Be
                case 5:  constants = {3.64f, 0.3347f}; break; // B
                case 6:  constants = {3.55f, 0.2929f}; break; // C
                case 7:  constants = {3.25f, 0.7113f}; break; // N
                case 8:  constants = {3.03f, 0.7113f}; break; // O
                case 9:  constants = {2.95f, 0.6276f}; break; // F
                case 10: constants = {2.78f, 0.1632f}; break; // Ne
                case 11: constants = {2.43f, 0.4184f}; break; // Na
                case 12: constants = {3.43f, 0.4184f}; break; // Mg
                case 13: constants = {4.01f, 0.5021f}; break; // Al
                case 14: constants = {3.83f, 0.8368f}; break; // Si
                case 15: constants = {3.74f, 0.8368f}; break; // P
                case 16: constants = {3.56f, 1.0460f}; break; // S
                case 17: constants = {3.47f, 1.2552f}; break; // Cl
                case 18: constants = {3.30f, 0.4184f}; break; // Ar
                case 19: constants = {2.75f, 0.4184f}; break; // K
                case 20: constants = {3.50f, 0.5021f}; break; // Ca
                case 21: constants = {3.80f, 0.7113f}; break; // Sc
                case 22: constants = {3.75f, 0.8368f}; break; // Ti
                case 23: constants = {3.70f, 0.9205f}; break; // V
                case 24: constants = {3.65f, 1.0460f}; break; // Cr
                case 25: constants = {3.60f, 1.1715f}; break; // Mn
                case 26: constants = {3.55f, 1.2970f}; break; // Fe
                default: constants = {2.50f, 0.1255f}; break; // H fallback
            }

            constants.first  *= MULT_FACTOR;  
            constants.second *= MULT_FACTOR; 
            return constants;
        }

        inline const std::map<uint8_t, float> electronegativity =
        {
            {1,  2.20f}, // H
            {2,  0.00f}, // He
            {3,  0.98f}, // Li
            {4,  1.57f}, // Be
            {5,  2.04f}, // B
            {6,  2.55f}, // C
            {7,  3.04f}, // N
            {8,  3.44f}, // O
            {9,  3.98f}, // F
            {10, 0.00f}, // Ne
            {11, 0.93f}, // Na
            {12, 1.31f}, // Mg
            {13, 1.61f}, // Al
            {14, 1.90f}, // Si
            {15, 2.19f}, // P
            {16, 2.58f}, // S
            {17, 3.16f}, // Cl
            {18, 0.00f}, // Ar
            {19, 0.82f}, // K
            {20, 1.00f}, // Ca
            {21, 1.36f}, // Sc
            {22, 1.54f}, // Ti
            {23, 1.63f}, // V
            {24, 1.66f}, // Cr
            {25, 1.55f}, // Mn
            {26, 1.83f}  // Fe
        };

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
                case 18: return 8;  // Ar
                case 19: return 1;  // K
                case 20: return 2;  // Ca
                case 21: return 3;  // Sc (3d¹ 4s² → 3)
                case 22: return 4;  // Ti (3d² 4s² → 4)
                case 23: return 5;  // V  (3d³ 4s² → 5)
                case 24: return 6;  // Cr (3d⁵ 4s¹ → 6)
                case 25: return 7;  // Mn (3d⁵ 4s² → 7)
                case 26: return 8;  // Fe (3d⁶ 4s² → 8)
                default: return 0;
            }
        }

        inline std::string getAtomLetter(uint32_t ZIndex)
        {
            switch (ZIndex)
            {
                case 1:  return "H";  // Hydrogen
                case 2:  return "He"; // Helium
                case 3:  return "Li"; // Lithium
                case 4:  return "Be"; // Beryllium
                case 5:  return "B";  // Boron
                case 6:  return "C";  // Carbon
                case 7:  return "N";  // Nitrogen
                case 8:  return "O";  // Oxygen
                case 9:  return "F";  // Fluorine
                case 10: return "Ne"; // Neon
                case 11: return "Na"; // Sodium
                case 12: return "Mg"; // Magnesium
                case 13: return "Al"; // Aluminum
                case 14: return "Si"; // Silicon
                case 15: return "P";  // Phosphorus
                case 16: return "S";  // Sulfur
                case 17: return "Cl"; // Chlorine
                case 18: return "Ar"; // Argon
                case 19: return "K";  // Potassium
                case 20: return "Ca"; // Calcium
                case 21: return "Sc"; // Scandium
                case 22: return "Ti"; // Titanium
                case 23: return "V";  // Vanadium
                case 24: return "Cr"; // Chromium
                case 25: return "Mn"; // Manganese
                case 26: return "Fe"; // Iron
                default: return "H";
            }
        }

        inline float getBondLength(uint8_t ZIndex1, uint8_t ZIndex2, fun::BondType type)
        {
            float base = (getAtomConstants(ZIndex1).first / MULT_FACTOR + getAtomConstants(ZIndex2).first / MULT_FACTOR) / 2.0f * BOND_LENGTH_FACTOR;
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
                    base *= 0.7f;
                    if ((ZIndex1 == 6 && ZIndex2 == 6)) base = 1.34f; // C=C 
                    if ((ZIndex1 == 6 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 6)) base = 1.21f; // C=O
                    if ((ZIndex1 == 6 && ZIndex2 == 7) || (ZIndex1 == 7 && ZIndex2 == 6)) base = 1.38f; // C=N
                    if ((ZIndex1 == 6 && ZIndex2 == 16) || (ZIndex1 == 16 && ZIndex2 == 6)) base = 1.71f; // C=S
                    if ((ZIndex1 == 7 && ZIndex2 == 7)) base = 1.20f; // N=N 
                    if ((ZIndex1 == 16 && ZIndex2 == 16)) base = 2.05f; // S=S 
                    break;
                case fun::BondType::TRIPLE:
                    base *= 0.6f;
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

            void draw(sf::Vector3f& pos, core::window_t &window, bool letterMode);
        };

        class universe
        {
        public:
            universe(float universeSize = 10.f);
            ~universe() = default;

            size_t createAtom(sf::Vector3f p, sf::Vector3f v, uint8_t ZIndex = 1, uint8_t numNeutrons = 0, uint8_t numElectrons = 1);
            size_t createSubset(const size_t central, const size_t subsetNext, const size_t subsetLast, 
                const size_t mainNext, const size_t mainLast, const std::vector<size_t>& bonds, const std::vector<fun::BondType>& bondTypes);
            void createMolecule(const molecule_structure& structure, sf::Vector3f pos);

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
            sf::Vector3f ljGrad(size_t i);
            sf::Vector3f ljForce(size_t i, size_t j);
            sf::Vector3f coulombForce(size_t i, size_t j, sf::Vector3f& dr_vec);

            void calcBondForces();
            void calcAngleForces();
            void calcLjForces();
            void calcElectrostaticForces();

            void setTemperature(float kelvin = 0.f);
            
            // Molecule Creation
            void organizeMolecule(const molecule_structure& structure, const sf::Vector3f& initPos);
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

            // Others

            float calculateKineticEnergy();

            // Reactions

            void handleReactions();

            float boxSize = 10.f;
            std::vector<atom> atoms;
            std::vector<subset> subsets;

            std::vector<bond> bonds;

            std::vector<sf::Vector3f> forces;
            std::vector<sf::Vector3f> positions;
            std::vector<sf::Vector3f> velocities;

            float temp = 0;
            size_t timeStep = 0;

            // Helper Funcs
            sf::Vector3f minImageVec(sf::Vector3f dr)
            {
                dr.x -= boxSize * std::round(dr.x / boxSize);
                dr.y -= boxSize * std::round(dr.y / boxSize);
                return dr;
            }
        };
    } // namespace fun
} // namespace sim
