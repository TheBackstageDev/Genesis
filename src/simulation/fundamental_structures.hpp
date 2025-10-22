#pragma once

#include <SFML/Graphics.hpp>

#include "core/window.hpp"
#include <vector>

#include <algorithm>

namespace sim
{
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

        #define AVOGADRO 6.02214076e26f // conversion from Daltons to Kg
        #define BOLTZMAN_CONSTANT 1.380649e-23f // Boltzman Constant m^2 kg s^-2 K^-1
        #define KB (BOLTZMAN_CONSTANT * AVOGADRO * ANGSTROM) / PICOSECOND // A^2 D ps^-2 K^-1
        #define THERMOSTAT_INTERVAL 5

        #define BOND_K 16.0f // Harmonic force constant (kJ/mol/Å²)
        #define BOND_LENGTH_FACTOR 0.8f

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
            constants.first *= 5.f;
            constants.second *= 5.f;
            return constants;
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
        enum class BondType { SINGLE, DOUBLE, TRIPLE };

        struct bond {
            uint32_t atom1; 
            uint32_t atom2; // bonded to
            
            double equilibriumLength = 0.f;
            double forceConstant = 0.f;
            BondType type;
        };

        struct atom
        {
            sf::Vector2f position; 
            sf::Vector2f velocity; 
            float mass;            
            float radius;

            float epsilon; // LJ 
            float sigma; // LJ

            uint32_t ZIndex;

            void draw(core::window_t &window, bool letterMode);
        };

        class universe
        {
        public:
            universe(float universeSize = 10.f);
            ~universe() = default;

            void createAtom(sf::Vector2f p, sf::Vector2f v, uint32_t ZIndex = 1);
            void createBond(uint32_t idx1, uint32_t idx2, BondType type = BondType::SINGLE);

            void update(float targetTemperature = 1.0f);
            void draw(core::window_t &window, bool letter = false);
            void drawDebug(core::window_t& window);

            size_t numAtoms() { return atoms.size(); }

            float temperature() { return temp; };
            float timestep() { return timeStep; }

        private:
            void boundCheck(atom &a);

            float ljPot(size_t i, float epsilon, float sigma);
            sf::Vector2f ljGrad(size_t i, float epsilon, float sigma);

            void calcBondForces(bond& bond);

            float boxSize = 10.f;
            std::vector<atom> atoms;
            std::vector<bond> bonds;

            std::vector<sf::Vector2f> forces;
            std::vector<sf::Vector2f> prev_positions;

            float temp = 0;
            size_t timeStep = 0;
        };
    } // namespace fun
} // namespace sim
