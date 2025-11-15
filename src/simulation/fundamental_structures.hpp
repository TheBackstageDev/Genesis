#pragma once

#include <SFML/Graphics.hpp>

#include "core/window.hpp"
#include "core/camera.hpp"
#include <vector>
#include <map>
#include <unordered_map>
#include <set>

#include <chrono>
#include <algorithm>

#include "constants.hpp"
#include "smiles_parser.hpp"

namespace sim
{
    namespace fun
    {
        // for vizualization and debugging
        struct hydrogen_bond
        {
            size_t H, A, D;
            float strength = 0.f;
        };

        struct bond 
        {
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
            
            std::vector<size_t> hydrogenIdx;
            std::vector<size_t> connectedIdx; // non hydrogen connected to the main atom
        };

        struct molecule
        {
            std::string name = "Unknown"; // For Visualization
            std::vector<size_t> subsetIdx;
            std::vector<size_t> atomIdx;
            std::vector<size_t> bondIdx;

            std::vector<size_t> dihedrals;
            std::vector<size_t> angles;
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

            int32_t chirality;
            void draw(float temperature, sf::Vector2f& pos, float camDistance, core::window_t &window, bool letterMode, bool lennardBall);
        };
    } // namespace fun
} // namespace sim
