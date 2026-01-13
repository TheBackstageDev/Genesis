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
            if (distance < 1.24f) return fun::BondType::TRIPLE;     // C≡C     ~1.20–1.23 Å
            if (distance < 1.39f) return fun::BondType::DOUBLE;     // C=C     ~1.33–1.38 Å (alkene/aromatic)
            return fun::BondType::SINGLE;                           // C–C     ~1.47–1.54 Å
        }

        if (BONDED(6, 8))
        {
            if (distance < 1.24f) return fun::BondType::DOUBLE;     // C=O     ~1.20–1.23 Å (carbonyl)
            return fun::BondType::SINGLE;                           // C–O     ~1.36–1.43 Å
        }

        if (BONDED(6, 7))
        {
            if (distance < 1.18f) return fun::BondType::TRIPLE;     // C≡N     ~1.15–1.17 Å (nitrile)
            if (distance < 1.32f) return fun::BondType::DOUBLE;     // C=N     ~1.27–1.31 Å (imine)
            return fun::BondType::SINGLE;                           // C–N     ~1.45–1.47 Å
        }

        if (BONDED(7, 8))
        {
            if (distance < 1.30f) return fun::BondType::DOUBLE;     // N=O     ~1.22–1.28 Å (nitro, oxime)
            return fun::BondType::SINGLE;                           // N–O     ~1.40–1.46 Å
        }

        if (BONDED(6, 16))
        {
            if (distance < 1.68f) return fun::BondType::DOUBLE;     // C=S     ~1.60–1.67 Å (thio-carbonyl)
            return fun::BondType::SINGLE;                           // C–S     ~1.78–1.82 Å
        }

        if (BONDED(15, 8))
        {
            if (distance < 1.51f) return fun::BondType::DOUBLE;     // P=O     ~1.44–1.51 Å
            return fun::BondType::SINGLE;                           // P–O     ~1.58–1.65 Å
        }

        if (BONDED(15, 7))
        {
            if (distance < 1.55f) return fun::BondType::DOUBLE;     // P=N     ~1.50–1.57 Å
            return fun::BondType::SINGLE;                           // P–N     ~1.65–1.80 Å
        }

        if (BONDED(15, 16))
        {
            if (distance < 1.95f) return fun::BondType::DOUBLE;     // P=S     ~1.90–1.94 Å
            return fun::BondType::SINGLE;                           // P–S     ~2.05–2.12 Å
        }

        if (BONDED(7, 7))
        {
            if (distance < 1.12f) return fun::BondType::TRIPLE;     // N≡N     ~1.10 Å
            if (distance < 1.25f) return fun::BondType::DOUBLE;     // N=N     ~1.24 Å (azo)
            return fun::BondType::SINGLE;                           // N–N     ~1.45 Å
        }

        if (BONDED(8, 8))
        {
            if (distance < 1.23f) return fun::BondType::DOUBLE;     // O=O     ~1.21 Å
            return fun::BondType::SINGLE;                           // O–O     ~1.48 Å (peroxide)
        }

        float r_single = constants::getBondLength(z1, z2, fun::BondType::SINGLE);
        float r_double = constants::getBondLength(z1, z2, fun::BondType::DOUBLE);
        float r_triple = constants::getBondLength(z1, z2, fun::BondType::TRIPLE);

        if (r_single < 0.3f || r_single > 5.0f) 
        {
            r_single = (constants::covalent_radius[z1] + constants::covalent_radius[z2]) * 1.10f;
        }
        if (r_double < 0.3f || r_double > 5.0f) 
        {
            r_double = r_single * 0.88f;
        }
        if (r_triple < 0.3f || r_triple > 5.0f) 
        {
            r_triple = r_single * 0.78f;
        }

        if (distance <= r_triple + 0.05f) 
        {
            return fun::BondType::TRIPLE;
        }

        if (distance >= r_single + 0.35f) 
        {
            return fun::BondType::SINGLE;
        }

        constexpr float pauling_const = 0.30f;
        float delta = r_single - distance;
        float bond_order = std::exp(delta / pauling_const);

        bond_order = std::clamp(bond_order, 1.0f, 3.2f);

        if (bond_order >= 2.6f) 
        {
            return fun::BondType::TRIPLE;
        }
        else if (bond_order >= 1.65f) 
        {
            return fun::BondType::DOUBLE;
        }

        return fun::BondType::SINGLE;
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

        constexpr float tolerance = 0.35f;

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
