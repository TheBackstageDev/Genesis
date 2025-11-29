#include "format_loader.hpp"

#include <fstream>
#include <iostream>
#include <strstream>

namespace sim::io
{
    inline fun::BondType estimateBondType(uint8_t z1, uint8_t z2, float distance)
    {
        if (z1 > z2) std::swap(z1, z2);

        if (z1 == 1 || z2 == 1) return fun::BondType::SINGLE;

        if (z1 == 6 && z2 == 6) 
        {
            if (distance < 1.26f) return fun::BondType::TRIPLE;     // C≡C
            if (distance < 1.4f) return fun::BondType::DOUBLE;     // C=C

            return fun::BondType::SINGLE;
        }

        if ((z1 == 6 && z2 == 8) || (z1 == 8 && z2 == 6)) 
        {
            if (distance < 1.30f) return fun::BondType::DOUBLE;     // C=O
            return fun::BondType::SINGLE;
        }

        if ((z1 == 6 && z2 == 7) || (z1 == 7 && z2 == 6)) 
        {
            if (distance < 1.26f) return fun::BondType::TRIPLE;     // C≡N
            if (distance < 1.32f) return fun::BondType::DOUBLE;     // C=N
            return fun::BondType::SINGLE;
        }

        if ((z1 == 7 && z2 == 8) || (z1 == 8 && z2 == 7)) 
        {
            if (distance < 1.30f) return fun::BondType::DOUBLE;     // N=O
            return fun::BondType::SINGLE;
        }

        if ((z1 == 6 && z2 == 16) || (z1 == 16 && z2 == 6)) 
        {
            if (distance < 1.70f) return fun::BondType::DOUBLE;     // C=S
            return fun::BondType::SINGLE;
        }

        float r_single = constants::getBondLength(z1, z2, fun::BondType::SINGLE);
        float r_double = constants::getBondLength(z1, z2, fun::BondType::DOUBLE);
        float r_triple = constants::getBondLength(z1, z2, fun::BondType::TRIPLE);

        // Fallback
        if (r_single <= 0.1f) r_single = (constants::covalent_radius[z1] + constants::covalent_radius[z2]) * 1.05f;
        if (r_double <= 0.1f) r_double = r_single * 0.87f;
        if (r_triple <= 0.1f) r_triple = r_single * 0.78f;

        constexpr float pauling_const = 0.30f;
        if (distance >= r_single + 0.3f) return fun::BondType::SINGLE;
        if (distance <= r_triple - 0.1f) return fun::BondType::TRIPLE;

        float order = std::exp((r_single - distance) / pauling_const);

        order = std::clamp(order, 1.0f, 3.0f);
        if (order > 2.7f)  return fun::BondType::TRIPLE;
        if (order > 1.7f)  return fun::BondType::DOUBLE;

        return fun::BondType::SINGLE;
    }

    bool io::loadXYZ(const std::filesystem::path path, std::vector<fun::def_atom>& atoms, std::vector<fun::def_bond>& bonds, std::vector<sf::Vector3f>& positions)
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
            double x, y, z;

            if (!(iss >> symbol >> x >> y >> z)) 
            {
                std::cerr << "[XYZ] Parse error at line " << i + 3 << '\n';
                continue;
            }

            fun::def_atom nAtom{};
            nAtom.aromatic = false;
            nAtom.charge = 0.f;
            nAtom.chirality = 0;
            nAtom.hydrogenize = false;
            nAtom.nBonds = 0;
            nAtom.nHydrogens = 0;
            nAtom.ZIndex = constants::symbolToZ(symbol);
            nAtom.NIndex = nAtom.ZIndex;

            atoms.emplace_back(std::move(nAtom));
            positions.emplace_back(x, y, z);
        }

        constexpr float tolerance = 0.45f;

        bonds.clear();
        for (size_t i = 0; i < atoms.size(); ++i)
        {
            uint8_t z1 = atoms[i].ZIndex;
            for (size_t j = i + 1; j < atoms.size(); ++j)
            {
                uint8_t z2 = atoms[j].ZIndex;
                double dist = (positions[i] - positions[j]).length();
                float r_cov = (constants::covalent_radius[z1] + constants::covalent_radius[z2]);

                if (dist < r_cov + tolerance) 
                {
                    fun::def_bond nBond{};
                    
                    if (z1 == 1) 
                    {
                        nBond.bondingAtomIdx = i;
                        nBond.centralAtomIdx = j;
                    } 
                    else if (z2 == 1) {
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
