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

#define GLM_FORCE_RADIANS
#define GLM_FORCE_INTRINSICS
#define GLM_FORCE_ALIGNED_GENTYPES

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

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
            std::vector<dihedral_angle> improper_angles;
            std::vector<sf::Vector3f> positions;
        };

        // for vizualization and debugging
        struct hydrogen_bond
        {
            uint32_t H, A, D;
            float strength = 0.f;
        };

        struct bond 
        {
            uint32_t bondedAtom; 
            uint32_t centralAtom; 
            
            float equilibriumLength = 0.f;
            float k = 100.f;
            BondType type;
            int32_t padding;
        };

        struct subset
        {
            uint32_t mainAtomIdx;              
            uint32_t bondedSubsetIdx   = UINT32_MAX;       // Index of the last subset (optional, max(uint32_t) if none)
            uint32_t bondingSubsetIdx  = UINT32_MAX;      // Index of the next subset (optional, max(uint32_t) if none)
            
            uint32_t hydrogenBegin     = UINT32_MAX;
            uint32_t hydrogenCount     = 0;
            
            uint32_t connectedBegin    = UINT32_MAX;
            uint32_t connectedCount    = 0;

            uint32_t pad;
        };

        struct molecule
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
        };

        enum color_rendering_mode
        {
            COLOR,
            VELOCITY,
            CHARGE,
            COUNT
        };
        
        struct rendering_info
        {
            bool lennardBall = true; 
            bool licorice = false; 
            bool hyperBalls = false;
            bool spaceFilling = false; 
            bool universeBox = true;
            bool renderWater = true;

            color_rendering_mode color_mode{color_rendering_mode::COLOR};

            float opacity = 1.0f;
            ImVec4 color_addition{0.0f, 0.0f, 0.0f, 0.0f};

            bool flag_highlights = true;
            bool flag_arrows = true;
            std::vector<uint32_t> highlight_indices;
            std::vector<std::pair<uint32_t, uint32_t>> highlight_bonds; // first and second are the indices of the atoms
            std::vector<std::pair<glm::vec3, glm::vec3>> arrows;
        };

        enum class compound_type : uint32_t
        {
            ORGANIC,
            BIOMOLECULE,
            INORGANIC,
            ION,
            NANOMATERIAL,
            ELEMENTS,
            COUNT
        };

        struct compound_preset_info
        {
            std::string name;
            std::string formula;
            float molecular_weight = 0.0f;

            compound_type type;
            sim::fun::molecule_structure structure;
            uint32_t id = 0;
        };

        struct rendering_simulation_info
        {
            std::vector<glm::vec3>& positions;
            std::vector<glm::vec3>& velocities;

            std::vector<float>& q;
            std::vector<atom>& atoms;
            std::vector<bond>& bonds;
            std::vector<molecule>& molecules;
            glm::vec3 box;
        };
    } // namespace fun
} // namespace sim
