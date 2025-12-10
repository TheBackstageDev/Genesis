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

namespace sim
{
    namespace fun
    {
        struct angle
        {
            size_t A, B, C; // atom idx

            float K; // spring constant
            float rad; // angle in radians
        };

        struct dihedral_angle
        {
            size_t A, B, C, D; // atom idx

            float K; // spring constant
            float rad; // angle in radians
            int32_t periodicity; // 1, 2, 3
        };

        struct def_atom
        {
            uint8_t ZIndex;
            uint8_t NIndex;
            int8_t charge = 0; // non-Zero for ions 
            int32_t nBonds = 0;
            int32_t chirality = 0;

            int32_t nHydrogens = 0;
            bool aromatic = false;
            bool hydrogenize = true; // add hydrogens
        };

        struct def_bond 
        {
            size_t bondingAtomIdx;   
            size_t centralAtomIdx;  

            sim::fun::BondType type;
        };

        struct def_subset 
        {
            size_t mainAtomIdx;      
            size_t bondedSubset = SIZE_MAX;
            size_t bondingSubset = SIZE_MAX;
            float idealAngle;

            std::vector<size_t> hydrogensIdx;
            std::vector<size_t> connectedIdx; // any connection other than hydrogen 
        };

        struct molecule_structure 
        {
            std::vector<def_atom> atoms;
            std::vector<def_subset> subsets;
            std::vector<def_bond> bonds;
            std::vector<angle> angles;
            std::vector<dihedral_angle> dihedral_angles;
            std::vector<sf::Vector3f> positions;
        };

        // for vizualization and debugging
        struct hydrogen_bond
        {
            size_t H, A, D;
            float strength = 0.f;
        };

        struct reactive_bond
        {
            size_t i, j;
            float bo;
            BondType type = BondType::NONE;
        };

        struct bond 
        {
            uint32_t bondedAtom; 
            uint32_t centralAtom; 
            
            float equilibriumLength = 0.f;
            float k = 100.f;
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

            bool water = false;
        };

        struct atom
        {
            float mass;            
            float radius;

            float epsilon; // LJ 
            float sigma; // LJ

            int8_t electrons;

            uint8_t ZIndex;
            uint8_t NCount; // neutrons
            int8_t bondCount;

            int32_t chirality;
            void draw(float temperature, sf::Vector2f& pos, float camDistance, float q, core::window_t &window, bool letterMode, bool lennardBall);
        };
    } // namespace fun
} // namespace sim
