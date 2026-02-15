#include "format_loader.hpp"

#include <fstream>
#include <iostream>

namespace sim::io
{
    inline fun::BondType estimateBondType(uint8_t z1, uint8_t z2, float distance)
    {
        if (z1 == 1 || z2 == 1) 
        {
            return fun::BondType::SINGLE;
        }

    #define BONDED(a, b) ((z1 == (a) && z2 == (b)) || (z1 == (b) && z2 == (a)))

        if (BONDED(6, 6))
        {
            if (distance < 1.24f) return fun::BondType::TRIPLE;     // C≡C     
            if (distance < 1.38f) return fun::BondType::DOUBLE;     // C=C     
            return fun::BondType::SINGLE;                           // C–C    
        }

        if (BONDED(6, 8))
        {
            if (distance < 1.24f) return fun::BondType::DOUBLE;     // C=O
            if (distance < 1.17f) return fun::BondType::TRIPLE;     // C#O
            return fun::BondType::SINGLE;                           // C–O
        }

        if (BONDED(6, 7))
        {
            if (distance < 1.18f) return fun::BondType::TRIPLE;
            if (distance < 1.35f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE; 
        }

        if (BONDED(6, 17))
        {
            if (distance < 1.69f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(7, 8))
        {
            if (distance < 1.23f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(14, 8))
        {
            if (distance < 1.51f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(16, 6))
        {
            if (distance < 1.68f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(16, 16))
        {
            if (distance < 1.90) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(16, 8))
        {
            if (distance < 1.45f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(15, 8))
        {
            if (distance < 1.51f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(15, 7))
        {
            if (distance < 1.55f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(15, 16))
        {
            if (distance < 1.95f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(7, 7))
        {
            if (distance < 1.10f) return fun::BondType::TRIPLE;
            if (distance < 1.24f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(8, 8))
        {
            if (distance < 1.21f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(5, 6))   
        {
            if (distance < 1.4f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(5, 17))
        {
            if (distance < 1.68f) return fun::BondType::DOUBLE;
            if (distance < 1.88f) return fun::BondType::SINGLE;
            return fun::BondType::NONE;
        }

        if (BONDED(5, 5))   
        {
            if (distance < 1.5f) return fun::BondType::DOUBLE;
            return fun::BondType::SINGLE;
        }

        if (BONDED(1, 5))  
        {
            if (distance < 1.20f) return fun::BondType::SINGLE;
            return fun::BondType::NONE;
        }

        if (BONDED(1, 6))  
        {
            if (distance < 1.15f) return fun::BondType::SINGLE; 
            return fun::BondType::NONE;
        }

        if (BONDED(6, 6))
        {
            if (distance < 1.45f) return fun::BondType::DOUBLE;
            if (distance < 1.58f) return fun::BondType::SINGLE;
        }

        return fun::BondType::NONE;
    #undef BONDED
    }

    bool loadXYZ(const std::filesystem::path path, std::vector<fun::def_atom>& atoms, std::vector<fun::def_bond>& bonds, std::vector<sf::Vector3f>& positions)
    {
        std::ifstream file(path);
        if (!file.is_open() || path.extension() != ".xyz") 
        {
            std::cerr << "[XYZ] Cannot open: " << path << '\n';
            return false;
        }

        int32_t natoms = 0;
        std::string nAtoms_s;
        std::getline(file, nAtoms_s);
        natoms = std::atoi( nAtoms_s.c_str() );
        if (natoms <= 0) return false;
        
        atoms.clear();
        positions.clear();
        atoms.reserve(natoms);
        positions.reserve(natoms);

        std::string line;
        std::string title;
        std::getline(file, title);

        for (size_t i = 0; i < natoms; ++i)
        {
            std::getline(file, line);
            std::istringstream iss(line);
            std::string symbol;
            float x, y, z, q = 0.0f;

            if (!(iss >> symbol >> x >> y >> z)) 
            {
                std::cerr << "[XYZ] Parse error at line " << i + 3 << '\n';
                continue;
            }

            iss >> q; // custom charge, if it is unable to read, it stays at 0.0f

            fun::def_atom nAtom{};
            nAtom.aromatic = false;
            nAtom.charge = q;
            nAtom.chirality = 0;
            nAtom.hydrogenize = false;
            nAtom.nBonds = 0;
            nAtom.nHydrogens = 0;
            nAtom.ZIndex = constants::symbolToZ(symbol);
            nAtom.NIndex = constants::NEUTRON_COUNTS[nAtom.ZIndex];

            atoms.emplace_back(std::move(nAtom));
            positions.emplace_back(x, y, z);
        }

        constexpr float tolerance = 0.25f;

        bonds.clear();
        for (size_t i = 0; i < atoms.size(); ++i)
        {
            uint8_t z1 = atoms[i].ZIndex;
            for (size_t j = i + 1; j < atoms.size(); ++j)
            {
                uint8_t z2 = atoms[j].ZIndex;
                float dist = (positions[i] - positions[j]).length();
                float r_cov = (constants::covalent_radius[z1] + constants::covalent_radius[z2]);

                if (dist < r_cov + tolerance) 
                {
                    fun::def_bond nBond{};
                    
                    if (z1 == 1) 
                    {
                        nBond.bondingAtomIdx = i;
                        nBond.centralAtomIdx = j;
                    } 
                    else if (z2 == 1) 
                    {
                        nBond.bondingAtomIdx = j;
                        nBond.centralAtomIdx = i;
                    } 
                    else 
                    {
                        nBond.centralAtomIdx = std::min(i, j);
                        nBond.bondingAtomIdx = std::max(i, j);
                    }
                    nBond.type = estimateBondType(z1, z2, dist);

                    ++atoms[i].nBonds;
                    ++atoms[j].nBonds;
                    bonds.emplace_back(std::move(nBond));
                }
            }
        }
        return true;
    }
} // namespace sim::io
