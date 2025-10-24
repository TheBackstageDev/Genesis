#pragma once

#include <SFML/Graphics.hpp>

#include "core/window.hpp"
#include <vector>

#include <algorithm>

namespace sim
{
    namespace fun { enum class BondType { SINGLE, DOUBLE, TRIPLE }; };

    namespace constants
    {
        #define MASS_PROTON 1 // Daltons
        #define MASS_ELECTRON 1/1337 // Daltons
        #define EPSILON 0.1f 
        #define DT 0.01f // ps
        #define ANGSTROM 1e20f
        #define PICOSECOND 1e24f

        #define LJ_EPSILON_H 4.36f
        #define LJ_SIGMA_H 0.74f 

        #define CUTOFF 2.5f
        #define COULOMB_CUTOFF 15.f
        #define ELECTROSTATIC_CUTOFF 10.0f

        #define AVOGADRO 6.02214076e26f // conversion from Daltons to Kg
        #define BOLTZMAN_CONSTANT 1.380649e-23f // Boltzman Constant m^2 kg s^-2 K^-1
        #define KB (BOLTZMAN_CONSTANT * AVOGADRO * ANGSTROM) / PICOSECOND // A^2 D ps^-2 K^-1
        #define THERMOSTAT_INTERVAL 5

        #define COULOMB_K (138.935455f / ANGSTROM) // kJ·mol⁻¹·Å·e⁻² (1/(4πε₀) in Å units)
        
        #define BOND_K 3000.f // Harmonic force constant 
        #define ANGLE_K 13000.0f // kJ/mol/rad² for angular potential
        #define BOND_LENGTH_FACTOR 1.0f

        #define M_PI 3.1415926535897932
        #define RADIAN M_PI / 180

        #define MULT_FACTOR 5.f

        inline std::pair<float, float> getAtomConstants(uint32_t ZIndex)
        {
            std::pair<float, float> constants; // Sigma (Å), Epsilon (kJ/mol)
            switch (ZIndex)
            {
            case 1: // H
                constants = {LJ_SIGMA_H, LJ_EPSILON_H}; break;
            case 2: // He
                constants = {0.80f, 1.20f}; break; 
            case 3: // Li
                constants = {1.20f, 6.02f}; break; 
            case 4: // Be
                constants = {1.10f, 4.82f}; break; 
            case 5: // B
                constants = {1.50f, 3.61f}; break; 
            case 6: // C
                constants = {1.70f, 3.01f}; break; 
            case 7: // N
                constants = {1.60f, 2.71f}; break;
            case 8: // O
                constants = {1.50f, 2.89f}; break; 
            case 9: // F
                constants = {1.40f, 3.13f}; break; 
            case 10: // Ne
                constants = {0.90f, 1.81f}; break; 
            default:
                constants = {LJ_SIGMA_H, LJ_EPSILON_H};
            }
            constants.first *= MULT_FACTOR;
            constants.second *= MULT_FACTOR;
            return constants;
        }

        // Temporary
        inline uint8_t getValenceElectrons(uint8_t ZIndex)
        {
            switch (ZIndex)
            {
                case 1: return 1;  // H
                case 2: return 2;  // He
                case 3: return 1;  // Li
                case 4: return 2;  // Be
                case 5: return 3;  // B
                case 6: return 4;  // C
                case 7: return 5;  // N
                case 8: return 6;  // O
                case 9: return 7;  // F
                case 10: return 8; // Ne
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
                break; 
            case fun::BondType::DOUBLE:
                base *= 0.8;
                if ((ZIndex1 == 6 && ZIndex2 == 6)) base = 1.34f; // C=C
                if ((ZIndex1 == 6 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 6)) base = 1.21f; // C=O
                break;
            case fun::BondType::TRIPLE:
                base *= 0.7f; 
                if ((ZIndex1 == 6 && ZIndex2 == 6)) base = 1.20f; // C≡C
                break;
            default:
                break;
            }

            return base * MULT_FACTOR;
        }

        inline float getAngles(uint8_t centralZIndex, const std::vector<uint8_t>& neighborZs, const std::vector<fun::BondType>& types)
        {
            size_t bond_count = neighborZs.size();
            if (bond_count < 2) return 0.0f; // No angle if fewer than 2 neighbors

            uint8_t valence = getValenceElectrons(centralZIndex);
            uint8_t bonding_electrons = bond_count * 2; // Assuming single bonds
            uint8_t lone_pairs = (valence - bonding_electrons) / 2;
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
            switch(ZIndex)
            {
            case 1: name = "H"; break;
            case 2: name = "He"; break;
            case 3: name = "Li"; break;
            case 4: name = "Be"; break;
            case 5: name = "B"; break;
            case 6: name = "C"; break;
            case 7: name = "N"; break;
            case 8: name = "O"; break;
            case 9: name = "F"; break;
            case 10: name = "Ne"; break;
            default:
                name = "H";
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

            int8_t charge;
            int8_t electrons;

            uint8_t ZIndex;
            int8_t bondCount;

            void draw(core::window_t &window, bool letterMode);
        };

        class universe
        {
        public:
            universe(float universeSize = 10.f);
            ~universe() = default;

            void createAtom(sf::Vector2f p, sf::Vector2f v, uint8_t ZIndex = 1, uint8_t numElectron = 1);
            void createBond(size_t idx1, size_t idx2, BondType type = BondType::SINGLE);

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
