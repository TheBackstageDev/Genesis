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
        #define DT 0.001f // ps
        #define ANGSTROM 1e20f
        #define PICOSECOND 1e24f

        #define LJ_EPSILON 318.f 
        #define LJ_SIGMA 2.5f 

        #define CUTOFF 15.0f // A

        #define AVOGADRO 6.02214076e26f // conversion from Daltons to Kg
        #define BOLTZMAN_CONSTANT 1.380649e-23f // Boltzman Constant m^2 kg s^-2 K^-1
        #define KB (BOLTZMAN_CONSTANT * AVOGADRO * ANGSTROM) / PICOSECOND // A^2 D ps^-2 K^-1
        #define TARGET_TEMP 15 // K
    };

    namespace fun
    {
        struct atom
        {
            sf::Vector2f position; 
            sf::Vector2f velocity; 
            float mass;            
            float radius;

            float epsilon; // LJ 
            float sigma; // LJ

            void draw(core::window_t &window);
        };

        class universe
        {
        public:
            universe(float universeSize = 10.f);
            ~universe() = default;

            void createAtom(sf::Vector2f p, sf::Vector2f v);

            void update();
            void draw(core::window_t &window);
            void drawDebug(core::window_t& window);

        private:
            void boundCheck(atom &a);

            float ljPot(float i, float epsilon, float sigma);
            sf::Vector2f ljGrad(float i, float epsilon, float sigma);

            float boxSize = 10.f;
            std::vector<atom> atoms;
            std::vector<std::pair<size_t, size_t>> bonds;

            size_t timeStep = 0;
        };
    } // namespace fun
} // namespace sim
