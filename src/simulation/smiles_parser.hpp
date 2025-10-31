#pragma once

#include <vector>

namespace sim
{
    namespace fun { enum class BondType { SINGLE, DOUBLE, TRIPLE }; };

    struct def_atom
    {
        uint8_t ZIndex;
        uint8_t NIndex;
        int8_t charge = 0; // non-Zero for ions 
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
    };

    class smiles_parser
    {
    public:
        smiles_parser();
        ~smiles_parser();

    private:
    };
} // namespace sim
