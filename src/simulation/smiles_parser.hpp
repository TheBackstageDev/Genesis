#pragma once

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

#include <SFML/Graphics.hpp>

#include "fundamental_structures.hpp"
#include "constants.hpp"

using namespace sim::fun;
namespace sim
{
    molecule_structure parseSMILES(const std::string& molecule, bool implicitHydrogens = true);
    void organizeAngles(std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds, 
                        std::vector<dihedral_angle>& dihedral_angles, std::vector<angle>& angles);
    void organizeSubsets(std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds);
    void addImplicitHydrogens(std::vector<def_atom>& nAtoms, std::vector<def_bond>& nBonds);
    void positionAtoms(const std::string& SMILES, std::vector<def_bond>& nBonds, const std::vector<std::vector<size_t>>& rings, 
        const std::vector<def_atom>& nAtoms, const std::vector<def_subset>& nSubsets, std::vector<sf::Vector3f>& positions);

    sf::Vector3f rotateDirection(const sf::Vector3f& v, const sf::Vector3f& axis, float angle);
} // namespace sim
