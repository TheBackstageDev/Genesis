#pragma once

#include <vector>
#include <map>
#include <string>

#include <SFML/Graphics.hpp>

namespace sim
{
    namespace fun { enum class BondType { SINGLE, DOUBLE, TRIPLE }; };

    inline std::map<std::string, uint8_t> AtomZTable = 
    {
        {"H",1}, {"He",2}, {"Li",3}, {"B",5}, {"C",6}, {"N",7}, {"O",8},
        {"F",9}, {"Ne",10}, {"Na",11}, {"Mg",12}, {"Al",13}, {"Si",14},
        {"P",15}, {"S",16}, {"Cl",17}, {"K",19}, {"Ca",20}, {"Fe",26}
    };

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

    molecule_structure parseSMILES(const std::string& molecule);
} // namespace sim
