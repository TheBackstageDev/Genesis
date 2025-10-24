#pragma once

#include <SFML/Graphics.hpp>

#include "core/window.hpp"
#include <vector>
#include <map>

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
        #define DT 0.005f // ps
        #define ANGSTROM 1e20f
        #define PICOSECOND 1e24f

        #define LJ_EPSILON_H 4.36f
        #define LJ_SIGMA_H 0.74f 

        #define MULT_FACTOR 5.f
        
        #define CUTOFF 2.5f
        #define COULOMB_CUTOFF 15.f * MULT_FACTOR
        #define ELECTROSTATIC_CUTOFF 10.0f * MULT_FACTOR

        #define AVOGADRO 6.02214076e26f // conversion from Daltons to Kg
        #define BOLTZMAN_CONSTANT 1.380649e-23f // Boltzman Constant m^2 kg s^-2 K^-1
        #define KB (BOLTZMAN_CONSTANT * AVOGADRO * ANGSTROM) / PICOSECOND // A^2 D ps^-2 K^-1
        #define THERMOSTAT_INTERVAL 5

        #define COULOMB_K 1389.3546f // kJ·mol⁻¹·Å·e⁻²
        
        #define BOND_K 2000.f // Harmonic force constant 
        #define ANGLE_K 12000.0f // kJ/mol/rad² for angular potential
        #define BOND_LENGTH_FACTOR 1.0f

        #define M_PI 3.1415926535897932
        #define RADIAN M_PI / 180


        inline std::pair<float, float> getAtomConstants(uint32_t ZIndex)
        {
            std::pair<float, float> constants; // Sigma (Å), Epsilon (kJ/mol)
            switch (ZIndex)
            {
                case 1:  // H
                    constants = {LJ_SIGMA_H, LJ_EPSILON_H}; break;
                case 2:  // He
                    constants = {0.80f, 1.20f}; break;
                case 3:  // Li
                    constants = {1.20f, 6.02f}; break;
                case 4:  // Be
                    constants = {1.10f, 4.82f}; break;
                case 5:  // B
                    constants = {1.50f, 3.61f}; break;
                case 6:  // C
                    constants = {1.70f, 3.01f}; break;
                case 7:  // N
                    constants = {1.60f, 2.71f}; break;
                case 8:  // O
                    constants = {1.50f, 2.89f}; break;
                case 9:  // F
                    constants = {1.40f, 3.13f}; break;
                case 10: // Ne
                    constants = {0.90f, 1.81f}; break;
                case 11: // Na
                    constants = {2.00f, 0.52f}; break;
                case 12: // Mg
                    constants = {1.80f, 1.10f}; break;
                case 13: // Al
                    constants = {2.10f, 2.45f}; break;
                case 14: // Si
                    constants = {2.10f, 2.49f}; break;
                case 15: // P
                    constants = {1.80f, 2.10f}; break;
                case 16: // S
                    constants = {1.80f, 2.08f}; break;
                case 17: // Cl
                    constants = {1.75f, 2.50f}; break;
                default:
                    constants = {LJ_SIGMA_H, LJ_EPSILON_H}; // Default to H
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

            switch (type)
            {
                case fun::BondType::SINGLE:
                    if ((ZIndex1 == 6 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 6)) base = 1.09f; // C-H
                    if ((ZIndex1 == 8 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 8)) base = 0.96f; // O-H
                    if ((ZIndex1 == 6 && ZIndex2 == 6)) base = 1.54f; // C-C
                    if ((ZIndex1 == 7 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 7)) base = 1.01f; // N-H
                    if ((ZIndex1 == 16 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 16)) base = 1.34f; // S-H
                    if ((ZIndex1 == 17 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 17)) base = 1.27f; // Cl-H
                    if ((ZIndex1 == 17 && ZIndex2 == 11) || (ZIndex1 == 11 && ZIndex2 == 17)) base = 2.82f; // Na-Cl
                    break;
                case fun::BondType::DOUBLE:
                    base *= 0.8f;
                    if ((ZIndex1 == 6 && ZIndex2 == 6)) base = 1.34f; // C=C
                    if ((ZIndex1 == 6 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 6)) base = 1.21f; // C=O
                    if ((ZIndex1 == 7 && ZIndex2 == 7)) base = 1.20f; // N=N
                    if ((ZIndex1 == 16 && ZIndex2 == 16)) base = 2.05f; // S=S
                    break;
                case fun::BondType::TRIPLE:
                    base *= 0.7f;
                    if ((ZIndex1 == 6 && ZIndex2 == 6)) base = 1.20f; // C≡C
                    if ((ZIndex1 == 7 && ZIndex2 == 7)) base = 1.10f; // N≡N
                    break;
                default:
                    break;
            }

            return base * MULT_FACTOR;
        }

        inline std::pair<float, float> getPartialCharges(uint8_t Z1, uint8_t Z2, fun::BondType type)
        {
            if ((Z1 == 8 && Z2 == 1) || (Z1 == 1 && Z2 == 8)) // O-H bond
            {
                return (type == fun::BondType::SINGLE) ? std::make_pair(-1.2f, 0.6f) : std::make_pair(0.0f, 0.0f); // O: -1.2e, H: +0.6e
            }
            return std::make_pair(0.0f, 0.0f);
        }

        inline float getAngles(uint8_t centralZIndex, const std::vector<uint8_t>& neighborZs, const std::vector<fun::BondType>& types)
        {
            size_t bond_count = neighborZs.size();
            if (bond_count < 2) return 0.0f; // No angle if fewer than 2 neighbors

            uint8_t valence = getValenceElectrons(centralZIndex);
            uint8_t bonding_electrons = 0;
            for (const auto& type : types)
            {
                bonding_electrons += (type == fun::BondType::SINGLE ? 2 : (type == fun::BondType::DOUBLE ? 4 : 6));
            }
            uint8_t lone_pairs = (valence - (bonding_electrons / 2)) / 2; // Approximate lone pairs
            uint8_t total_pairs = bond_count + lone_pairs;
            float ideal_angle = 0.0f;

            switch (total_pairs)
            {
                case 2: ideal_angle = 180.0f * RADIAN; break; // Linear
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
            uint32_t atom1; 
            uint32_t atom2; // bonded to
            
            double equilibriumLength = 0.f;
            BondType type;
        };

        struct subset
        {
            size_t mainAtomIdx;              
            std::vector<size_t> connectedIdx; // connected to the main atom
            float idealAngle;                // In radians
            size_t connectedSubsetIdx;       // Index of the next subset (optional, max(size_t) if none)
            subset(size_t mainIdx, const std::vector<size_t>& connIdx, float angle, size_t connSubset = -1)
                : mainAtomIdx(mainIdx), connectedIdx(connIdx), idealAngle(angle), connectedSubsetIdx(connSubset) {}
        };

        struct atom
        {
            sf::Vector2f position; 
            sf::Vector2f velocity; 
            float mass;            
            float radius;

            float epsilon; // LJ 
            float sigma; // LJ

            float charge;
            int8_t electrons;

            uint8_t ZIndex;
            uint8_t NCount; // neutrons
            int8_t bondCount;

            void draw(core::window_t &window, bool letterMode);
        };

        class universe
        {
        public:
            universe(float universeSize = 10.f);
            ~universe() = default;

            size_t createAtom(sf::Vector2f p, sf::Vector2f v, uint8_t ZIndex = 1, uint8_t numNeutrons = 0, uint8_t numElectrons = 1);
            void createBond(size_t idx1, size_t idx2, BondType type = BondType::SINGLE);

            void balanceMolecularCharges();

            void update(float targetTemperature = 1.0f);
            void draw(core::window_t &window, bool letter = false);
            void drawDebug(core::window_t& window);

            size_t numAtoms() { return atoms.size(); }

            float temperature() { return temp; }
            float timestep() { return timeStep; }

        private:
            void boundCheck(atom &a);

            float ljPot(size_t i, float epsilon, float sigma);
            sf::Vector2f ljGrad(size_t i);
            sf::Vector2f ljForce(size_t i, size_t j);
            sf::Vector2f coulombForce(size_t i, size_t j, sf::Vector2f& dr_vec);

            void calcBondForces();
            void calcAngleForces();
            void calcLjForces();
            void calcElectrostaticForces();

            float boxSize = 10.f;
            std::vector<atom> atoms;
            std::vector<bond> bonds;

            std::vector<sf::Vector2f> forces;
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
