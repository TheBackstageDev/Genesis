#pragma once

#include <vector>
#include <map>
#include <unordered_map>
#include <string>

#include "simulation/fundamental_structures.hpp"
#include "simulation/constants.hpp"

using namespace sim::fun;
namespace sim
{
    molecule_structure parseSMILES(const std::string& molecule, bool implicitHydrogens = true);
    void organizeAngles(std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds, 
                        std::vector<dihedral_angle>& dihedral_angles, std::vector<dihedral_angle>& improper_angles, std::vector<angle>& angles, bool carborane = false);
    void organizeSubsets(std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds);
    void addImplicitHydrogens(std::vector<def_atom>& nAtoms, std::vector<def_bond>& nBonds);

    glm::vec3 rotateDirection(const glm::vec3& v, const glm::vec3& axis, float angle);
} // namespace sim
