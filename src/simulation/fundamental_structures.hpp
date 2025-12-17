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
            uint32_t A, B, C; // atom idx

            float K; // spring constant
            float rad; // angle in radians

            int32_t padding[3];
        };

        struct dihedral_angle
        {
            uint32_t A, B, C, D; // atom idx

            float K; // spring constant
            float rad; // angle in radians
            int32_t periodicity; // 1, 2, 3

            int32_t padding;
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
            uint32_t bondingAtomIdx;   
            uint32_t centralAtomIdx;  

            sim::fun::BondType type;
            
            int32_t padding;
        };

        struct def_subset 
        {
            uint32_t mainAtomIdx;      
            uint32_t bondedSubset = UINT32_MAX;
            uint32_t bondingSubset = UINT32_MAX;
            float idealAngle;

            std::vector<uint32_t> hydrogensIdx;
            std::vector<uint32_t> connectedIdx; // any connection other than hydrogen 
        };

        struct molecule_structure 
        {
            std::vector<def_atom> atoms;
            std::vector<def_subset> subsets;
            std::vector<def_bond> bonds;
            std::vector<angle> angles;
            std::vector<std::vector<uint32_t>> rings_aromatic; // for drawing
            std::vector<dihedral_angle> dihedral_angles;
            std::vector<sf::Vector3f> positions;
        };

        // for vizualization and debugging
        struct hydrogen_bond
        {
            uint32_t H, A, D;
            float strength = 0.f;
        };

        struct reactive_bond
        {
            uint32_t i, j;
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

        struct alignas(32) subset
        {
            uint32_t mainAtomIdx;              
            uint32_t bondedSubsetIdx   = UINT32_MAX;       // Index of the last subset (optional, max(uint32_t) if none)
            uint32_t bondingSubsetIdx  = UINT32_MAX;      // Index of the next subset (optional, max(uint32_t) if none)
            
            uint32_t hydrogenBegin     = UINT32_MAX;
            uint32_t hydrogenCount     = UINT32_MAX;
            
            uint32_t connectedBegin    = UINT32_MAX;
            uint32_t connectedCount    = UINT32_MAX;

            int32_t padding;
        };

        struct alignas(64) molecule
        {
            uint32_t subsetBegin = UINT32_MAX;
            uint32_t subsetCount = UINT32_MAX;

            uint32_t atomBegin = UINT32_MAX;
            uint32_t atomCount = UINT32_MAX;

            uint32_t bondBegin = UINT32_MAX;
            uint32_t bondCount = UINT32_MAX;

            uint32_t angleBegin = UINT32_MAX;
            uint32_t angleCount = UINT32_MAX;

            uint32_t dihedralBegin = UINT32_MAX;
            uint32_t dihedralCount = UINT32_MAX;

            int32_t padding[6];

            bool exclude = false;
        };

        struct atom
        {
            float mass;            
            float radius;

            float epsilon; // LJ 
            float sigma; // LJ

            uint8_t electrons;
            uint8_t ZIndex;
            uint8_t NCount; // neutrons
            uint8_t bondCount;

            int32_t chirality;
            int8_t padding[12];
            void draw(float temperature, sf::Vector2f& pos, float camDistance, float q, core::window_t &window, sf::RenderTarget& target, bool letterMode, bool lennardBall, bool spaceFilling);
        };
    } // namespace fun
} // namespace sim
