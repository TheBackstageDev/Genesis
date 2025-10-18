#pragma once

#include <SFML/Graphics.hpp>

#include "core/window.hpp"
#include <vector>

#include <algorithm>

namespace sim
{
    namespace constants
    {
        #define MASS_PROTON 1
        #define MASS_ELECTRON 1/1337
        #define DT 0.01f 
        #define EPSILON 0.1f 
        #define H_RADIUS 0.31f

        #define LJ_EPSILON 1.0f 
        #define LJ_SIGMA 1.2f 

        #define CUTOFF 15.0f // A

        #define KB 1.380649e-23f // Boltzman Constant
        #define TARG_TEMP 100 // K
    };

    namespace fun
    {
        struct atom
        {
            sf::Vector2f position; 
            sf::Vector2f velocity; 
            float mass;            
            float radius;

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
