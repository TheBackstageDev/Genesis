#include "smiles_parser.hpp"

#include <iostream>

namespace sim
{
    molecule_structure sim::parseSMILES(const std::string& molecule)
    {
        molecule_structure nStructure{};
        std::vector<def_bond> nBonds{};
        std::vector<def_subset> nSubsets{};
        std::vector<def_atom> nAtoms{};

        for (size_t i = 0; i < molecule.size(); ++i)
        {
            const uint8_t c = molecule[i];
            
        }

        return {};
    }
} // namespace sim
