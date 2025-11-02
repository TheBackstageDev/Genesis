#pragma once

#include <vector>
#include <map>
#include <string>

#include <SFML/Graphics.hpp>

namespace sim
{
    namespace fun { enum class BondType { SINGLE, DOUBLE, TRIPLE }; };

    namespace constants
    {
        #define M_PI 3.1415926535
        #define RADIAN M_PI / 180
        #define DEGREE 180 / M_PI

        inline uint8_t symbolToZ(const std::string& sym)
        {
            static const std::map<std::string, uint8_t> elements = {
                {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6}, {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10},
                {"Na", 11}, {"Mg", 12}, {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18},
                {"K", 19}, {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24}, {"Mn", 25}, {"Fe", 26},
            };
            auto it = elements.find(sym);
            return it != elements.end() ? it->second : 6; // default C
        }
    
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
                case 21: return 3;  // Sc 
                case 22: return 4;  // Ti 
                case 23: return 5;  // V  
                case 24: return 6;  // Cr 
                case 25: return 7;  // Mn 
                case 26: return 8;  // Fe 
                default: return 0;
            }
        }

        inline uint8_t getUsualBonds(uint8_t ZIndex)
        {
            uint8_t valenceShell = getValenceElectrons(ZIndex);

            // special cases
            if (ZIndex == 1) return 1;   
            if (valenceShell == 8) return 0; // noble gases

            // octet rule
            uint8_t bonds = 8 - valenceShell;

            // exceptions
            switch (ZIndex)
            {
                case 4:  return 2;   // Be  
                case 5:  return 3;   // B   
                case 15: return 3;   // P   
                case 16: return 2;   // S   
                case 21: return 3;   // Sc  
                case 22: return 4;   // Ti 
                case 23: return 5;   // V  
                case 24: return 3;   // Cr  
                case 25: return 2;   // Mn  
                case 26: return 2;   // Fe 
            }

            return bonds;
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
    }; // namespace constants

    struct def_atom
    {
        uint8_t ZIndex;
        uint8_t NIndex;
        int8_t charge = 0; // non-Zero for ions 
        int32_t nBonds = 0;
        bool aromatic = false;
    };

    struct def_bond {
        size_t bondingAtomIdx;   
        size_t centralAtomIdx;  
        sim::fun::BondType type;
    };

    struct def_subset {
        size_t mainAtomIdx;      
        size_t bondedSubset;
        size_t bondingSubset;
        std::vector<size_t> connectedIdx; 
    };

    struct molecule_structure 
    {
        std::vector<def_atom> atoms;
        std::vector<def_subset> subsets;
        std::vector<def_bond> bonds;
        std::vector<sf::Vector3f> positons{};
    };

    molecule_structure parseSMILES(const std::string& molecule);
    std::vector<float> getAngles(const std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds);
    void organizeSubsets(std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds);
    void addImplicitHydrogens(std::vector<def_atom>& nAtoms, std::vector<def_bond>& nBonds);
    void positionAtoms(const std::vector<def_atom>& nAtoms, const std::vector<float>& angles, std::vector<sf::Vector3f>& positions);
} // namespace sim
