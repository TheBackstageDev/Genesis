#pragma once

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

#include <SFML/Graphics.hpp>

#include "constants.hpp"

namespace sim
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
        std::vector<size_t> branches; 
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

    molecule_structure parseSMILES(const std::string& molecule, bool implicitHydrogens = true);
    void organizeAngles(std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds, 
                        std::vector<dihedral_angle>& dihedral_angles, std::vector<angle>& angles);
    void organizeSubsets(std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds);
    void addImplicitHydrogens(std::vector<def_atom>& nAtoms, std::vector<def_bond>& nBonds);
    void positionAtoms(const std::string& SMILES, std::vector<def_bond>& nBonds, const std::vector<std::vector<size_t>>& rings, 
        const std::vector<def_atom>& nAtoms, const std::vector<def_subset>& nSubsets, std::vector<sf::Vector3f>& positions);

    sf::Vector3f rotateDirection(const sf::Vector3f& v, const sf::Vector3f& axis, float angle);
} // namespace sim
