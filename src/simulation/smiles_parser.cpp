#include "smiles_parser.hpp"

#include <stack>
#include <iostream>

namespace sim
{
    using namespace sim::fun;
    molecule_structure parseSMILES(const std::string &molecule, bool implicitHydrogens)
    {
        molecule_structure nStructure{};
        std::vector<def_bond> nBonds{};
        std::vector<def_subset> nSubsets{};
        std::vector<def_atom> nAtoms{};
        std::vector<angle> nAngles{};
        std::vector<dihedral_angle> nDihedralAngles{};
        std::vector<dihedral_angle> nImproperAngles{};
        std::vector<glm::vec3> npositions{};

        if (molecule.empty())
            return nStructure;

        uint32_t prevAtom = UINT32_MAX;
        BondType bondType = BondType::SINGLE; // Current Bond

        std::stack<uint32_t> branchStack;
        std::vector<uint32_t> branchAtoms;

        std::unordered_map<int32_t, uint32_t> ringOpen;                // ring number > atom index
        std::unordered_map<uint32_t, std::vector<uint32_t>> ringIndices; // atom indices
        std::unordered_map<uint32_t, std::vector<uint32_t>> ringBondIndices;

        std::vector<std::vector<uint32_t>> rings;

        std::string currentMolecule = molecule;
        bool newMolecule = false;

        for (uint32_t i = 0; i < currentMolecule.size(); ++i)
        {
            const uint8_t c = currentMolecule[i];

            if (isspace(c))
            {
                ++i;
                continue;
            }

            if (c == '=' || c == '#')
            {
                bondType = c == '=' ? BondType::DOUBLE : BondType::TRIPLE;
                continue;
            }

            // is Ring
            if (std::isdigit(c) || c == '%')
            {
                int ringID = 0;
                if (c == '%')
                {
                    i++;

                    while (i < currentMolecule.length() && std::isdigit(currentMolecule[i]))
                    {
                        ringID = ringID * 10 + (currentMolecule[i] - '0');
                        i++;
                    }
                }
                else
                {
                    ringID = c - '0';
                }

                if (ringOpen.count(ringID))
                {
                    uint32_t openAtom = ringOpen[ringID];
                    uint32_t closeAtom = prevAtom;
                    def_bond nBond{};
                    nBond.centralAtomIdx = closeAtom;
                    nBond.bondingAtomIdx = openAtom;
                    nBond.type = (nAtoms[openAtom].aromatic && nAtoms[closeAtom].aromatic)
                                ? BondType::SINGLE : bondType; // AROMATIC for later double and singles
                    nBonds.emplace_back(std::move(nBond));
                    nAtoms[openAtom].nBonds += static_cast<uint32_t>(nBond.type);
                    nAtoms[closeAtom].nBonds += static_cast<uint32_t>(nBond.type);

                    ringIndices[ringID].emplace_back(closeAtom);
                    rings.emplace_back(ringIndices[ringID]);

                    ringOpen.erase(ringID);
                }
                else
                {
                    ringOpen[ringID] = prevAtom;

                    auto& indices = ringIndices[ringID];
                    if (indices.empty())
                        indices.push_back(prevAtom);
                }
                bondType = BondType::SINGLE;
                continue;
            }

            if (c == '(')
            {
                if (prevAtom != UINT32_MAX)
                {
                    branchStack.emplace(prevAtom);
                    branchAtoms.emplace_back(prevAtom);
                }
                continue;
            }
            else if (c == ')')
            {
                if (!branchStack.empty())
                {
                    prevAtom = branchStack.top();
                    branchStack.pop();

                    branchAtoms.erase(
                        std::remove(branchAtoms.begin(), branchAtoms.end(), prevAtom),
                        branchAtoms.end()
                    );
                }
                continue;
            }

            if (c == '[')
            {
                uint32_t j = i + 1;
                std::string atomStr;
                while (j < currentMolecule.size() && currentMolecule[j] != ']') atomStr += currentMolecule[j++];
                if (j >= currentMolecule.size() || currentMolecule[j] != ']') { i = j; continue; }

                std::string symbol;
                int isotope = 0, charge = 0;
                uint32_t chirality = 0;
                uint32_t k = 0;

                if (std::isdigit(atomStr[k]))
                    isotope = std::stoi(atomStr.substr(k)) - 1;

                if (k >= atomStr.size() || !std::isalpha(atomStr[k])) { i = j; continue; }
                symbol += atomStr[k++];

                if (k < atomStr.size() && std::islower(atomStr[k]) &&
                    !std::isdigit(atomStr[k]) && atomStr[k] != '@' &&
                    atomStr[k] != '+' && atomStr[k] != '-' && atomStr[k] != ']')
                {
                    symbol += atomStr[k++];
                }

                if (symbol.empty()) { i = j; continue; }

                bool aromatic = (symbol.size() == 1 && std::islower(symbol[0]));
                if (aromatic) symbol[0] = std::toupper(symbol[0]);

                if (k < atomStr.size() && atomStr[k] == '@')
                {
                    chirality = 1;
                    if (k < atomStr.size() && atomStr[k] == '@') { chirality = 2; ++k; }
                }

                uint32_t Hcount = 0;
                if (k < atomStr.size() && atomStr[k] == 'H')
                {
                    ++k;
                    if (k < atomStr.size() && std::isdigit(atomStr[k]))
                        Hcount = atomStr[k] - '0';
                    else
                        Hcount = 1;
                }

                if (k < atomStr.size() && (atomStr[k] == '+' || atomStr[k] == '-'))
                {
                    char sign = atomStr[k];
                    int32_t mult = 1;
                    if (k < atomStr.size() && std::isdigit(atomStr[k]))
                        mult = atomStr[k] - '0';
                    charge = (sign == '+') ? mult : -mult;
                }

                uint8_t Z = constants::symbolToZ(symbol);
                if (Z == 0) { i = j; continue; }

                def_atom nAtom{};
                nAtom.charge    = charge;
                nAtom.NIndex    = isotope ? isotope : Z;
                nAtom.ZIndex    = Z;
                nAtom.aromatic  = aromatic;
                nAtom.chirality = chirality;
                nAtom.hydrogenize = Hcount > 0;
                nAtom.nHydrogens = Hcount;
                nAtoms.emplace_back(std::move(nAtom));

                uint32_t currentAtom = nAtoms.size() - 1;

                if (prevAtom != UINT32_MAX)
                {
                    def_bond nb{};
                    nb.centralAtomIdx = prevAtom;
                    nb.bondingAtomIdx = currentAtom;
                    nb.type           = bondType;
                    nBonds.emplace_back(std::move(nb));
                    nAtoms[prevAtom].nBonds += static_cast<uint32_t>(bondType);
                    nAtoms[currentAtom].nBonds += static_cast<uint32_t>(bondType);
                }

                prevAtom   = currentAtom;
                bondType   = BondType::SINGLE;

                i = j;
                continue;
            }

            // Different molecule
            if (c == '.')
            {
                prevAtom = UINT32_MAX;
                bondType = BondType::SINGLE;
                newMolecule = true;
                continue;
            }

            if (isalpha(c))
            {
                bool aromatic = islower(c);
                std::string sym(1, toupper(c));
                uint8_t ZIndex = constants::symbolToZ(sym);

                def_atom nAtom{};
                nAtom.charge = 0;
                nAtom.NIndex = 0;
                nAtom.ZIndex = ZIndex;
                nAtom.aromatic = aromatic;

                if (i + 1 < currentMolecule.size())
                {
                    char next = currentMolecule[i + 1];
                    if (isupper(c) && isalpha(next) && !isupper(next))
                    {
                        std::string two = sym + next;

                        uint8_t newZ = constants::symbolToZ(two);
                        if (newZ != 0)      
                        {
                            sym = two;
                            nAtom.ZIndex = newZ;
                            nAtom.NIndex = 0;
                            ++i;                            
                            aromatic = false;                 
                        }
                    }
                }

                nAtoms.emplace_back(std::move(nAtom));

                uint32_t currentAtom = nAtoms.size() - 1;

                if (prevAtom != UINT32_MAX)
                {   
                    def_bond nBond{};
                    nBond.bondingAtomIdx = currentAtom;
                    nBond.centralAtomIdx = prevAtom;
                    nBond.type = bondType;
                    nBonds.emplace_back(std::move(nBond));

                    nAtoms[prevAtom].nBonds += static_cast<uint32_t>(bondType);
                    nAtoms[currentAtom].nBonds += static_cast<uint32_t>(bondType);
                }
                else if (newMolecule)
                {
                    nAtoms[currentAtom].hydrogenize = false;
                }

                int32_t currentRingID = -1;
                for (const auto& [id, openAtom] : ringOpen)
                {
                    if (openAtom == prevAtom) 
                    {
                        currentRingID = id;
                        break;
                    }
                }

                for (const auto& [ringID, openAtom] : ringOpen)
                {
                    auto& idxList = ringIndices[ringID];
                    if (std::find(idxList.begin(), idxList.end(), prevAtom) == idxList.end() &&
                        std::find(branchAtoms.begin(), branchAtoms.end(), prevAtom) == branchAtoms.end())
                    {
                        idxList.push_back(prevAtom);
                        ringBondIndices[ringID].push_back(nBonds.size() - 1);
                    }
                }

                prevAtom = currentAtom;
                bondType = BondType::SINGLE;
            }
        }

        for (uint32_t ringIdx = 0; ringIdx < rings.size(); ++ringIdx) {
            const auto& ring = rings[ringIdx];
            if (ring.size() < 3) continue;

            bool isAromatic = true;
            for (uint32_t a : ring)
                if (!nAtoms[a].aromatic) { isAromatic = false; break; }
            if (!isAromatic) continue;
            
            bool doubleBond = true;
            for (uint32_t i = 0; i < ring.size(); ++i) {
                uint32_t a = ring[i];
                uint32_t b = ring[(i + 1) % ring.size()];
                
                auto bondIt = std::find_if(nBonds.begin(), nBonds.end(),
                [&](const def_bond& bb) {
                    return (bb.centralAtomIdx == a && bb.bondingAtomIdx == b) ||
                    (bb.centralAtomIdx == b && bb.bondingAtomIdx == a);
                });
                
                if (bondIt != nBonds.end()) 
                  bondIt->type = doubleBond ? BondType::DOUBLE : BondType::SINGLE;
                
                if (doubleBond)
                {
                    ++nAtoms[a].nBonds;
                    ++nAtoms[b].nBonds;
                }
                
                doubleBond = !doubleBond;

                nStructure.rings_aromatic.emplace_back(std::move(ring));
            }
        }

        if (implicitHydrogens)
            addImplicitHydrogens(nAtoms, nBonds);

        npositions.resize(nAtoms.size());
        organizeSubsets(nSubsets, nAtoms, nBonds);
        organizeAngles(nSubsets, nAtoms, nBonds, nDihedralAngles, nImproperAngles, nAngles);

        nStructure.atoms = std::move(nAtoms);
        nStructure.bonds = std::move(nBonds);
        nStructure.subsets = std::move(nSubsets);
        nStructure.angles = std::move(nAngles);
        nStructure.dihedral_angles = std::move(nDihedralAngles);
        nStructure.improper_angles = std::move(nImproperAngles);
        nStructure.positions = std::move(npositions);

        return nStructure;
    }

    struct pair_hash {
        uint32_t operator() (const std::pair<int64_t, int64_t>& p) const {
            return (uint32_t)(p.first << 32) | p.second;
        }
    };

    void organizeAngles(std::vector<def_subset> &nSubsets, const std::vector<def_atom> &nAtoms, const std::vector<def_bond> &nBonds,
                            std::vector<dihedral_angle>& dihedral_angles, std::vector<dihedral_angle>& improper_angles, std::vector<angle>& angles, bool carborane)
    {
        std::unordered_map<std::pair<uint32_t, uint32_t>, def_bond, pair_hash> bond_map;
        for (const auto& b : nBonds) 
        {
            uint32_t i = b.centralAtomIdx;
            uint32_t j = b.bondingAtomIdx;
            if (i > j) std::swap(i,j);
            bond_map[{i, j}] = b;
        }

        for (def_subset& sub : nSubsets) 
        {
            const uint32_t B = sub.mainAtomIdx;

            if (sub.connectedIdx.size() + sub.hydrogensIdx.size() < 2) continue;

            std::vector<uint8_t> Z(sub.connectedIdx.size() + sub.hydrogensIdx.size());
            std::vector<BondType> type(sub.connectedIdx.size() + sub.hydrogensIdx.size(), BondType::SINGLE);

            std::vector<uint32_t> neigh = sub.connectedIdx;
            neigh.insert(neigh.end(), sub.hydrogensIdx.begin(), sub.hydrogensIdx.end());
            for (uint32_t i = 0; i < neigh.size(); ++i) 
            {
                const uint32_t A = neigh[i];
                Z[i] = nAtoms[A].ZIndex;

                uint32_t a = B, b = A;
                if (a > b) std::swap(a,b);
                auto it = bond_map.find({a,b});
                if (it != bond_map.end()) type[i] = it->second.type;
            }

            for (uint32_t i = 0; i < neigh.size(); ++i) 
            {
                for (uint32_t j = i + 1; j < neigh.size(); ++j) 
                {
                    angle ang;
                    ang.A = neigh[i];
                    ang.B = B;
                    ang.C = neigh[j];
                    ang.rad = constants::getAngles(nAtoms[B].ZIndex, Z, type, carborane);
                    ang.K   = constants::getAngleHarmonicConstant(ang.A, ang.B, ang.C);
                    angles.emplace_back(std::move(ang));
                }
            }

            for (uint32_t c_idx = 0; c_idx < neigh.size(); ++c_idx) 
            {
                const uint32_t C = neigh[c_idx];

                uint32_t a = B, b = C;
                if (a > b) std::swap(a, b);
                const def_bond* bc_bond = nullptr;
                auto it = bond_map.find({a, b});
                if (it != bond_map.end()) bc_bond = &it->second;
                if (!bc_bond) continue;

                const auto& neigh_B = neigh; 

                std::vector<uint32_t> neigh_C;
                for (const auto& bond : nBonds) 
                {
                    if (bond.centralAtomIdx == C) neigh_C.push_back(bond.bondingAtomIdx);
                    else if (bond.bondingAtomIdx == C) neigh_C.push_back(bond.centralAtomIdx);
                }

                if (sub.connectedIdx.size() + sub.hydrogensIdx.size() == 3 || nAtoms[sub.mainAtomIdx].aromatic)
                {
                    uint32_t B = sub.mainAtomIdx;
                    uint8_t centralZ = nAtoms[B].ZIndex;

                    std::vector<uint32_t> bonded3;
                    bonded3.reserve(3);

                    for (uint32_t idx : sub.connectedIdx)   bonded3.push_back(idx);
                    for (uint32_t idx : sub.hydrogensIdx)   bonded3.push_back(idx);

                    if (bonded3.size() != 3) continue;

                    dihedral_angle imp{};
                    imp.A = bonded3[0];
                    imp.B = B;
                    imp.C = bonded3[1];
                    imp.D = bonded3[2];

                    imp.K           = 10.0f;
                    imp.periodicity = 2.0f;
                    imp.rad         = M_PI;

                    if (nAtoms[B].aromatic)
                        imp.K = 20.0f;

                    improper_angles.emplace_back(std::move(imp));
                }
                
                if (nAtoms[a].aromatic || nAtoms[b].aromatic) continue;

                for (uint32_t A : neigh_B) 
                {
                    if (A == C) continue;

                    for (uint32_t D : neigh_C) 
                    {
                        if (D == B) continue;
                        if (nAtoms[A].ZIndex == 1 && nAtoms[D].ZIndex == 1) continue;
                        if (nAtoms[A].aromatic && nAtoms[D].aromatic) continue;

                        dihedral_angle dh{};
                        dh.A = A; dh.B = B; dh.C = C; dh.D = D;

                        dh.rad      = M_PI / 3.0f;      // 60°
                        dh.K          = 1.5f;
                        dh.periodicity = 3;

                        if (bc_bond->type == BondType::DOUBLE) 
                        {
                            dh.rad      = M_PI;        // trans
                            dh.K          = 6.0f;
                            dh.periodicity = 2;
                        }
                        else if (nAtoms[B].ZIndex == 7 || nAtoms[C].ZIndex == 7) 
                        {
                            dh.rad      = M_PI;
                            dh.K          = 4.0f;
                            dh.periodicity = 2;
                        }

                        char chi = nAtoms[B].chirality ? nAtoms[B].chirality : nAtoms[C].chirality;
                        if (chi != 0) {
                            dh.rad      = 0.0f;
                            dh.K          = 12.0f;
                            dh.periodicity = 1;
                        }

                        dihedral_angles.emplace_back(std::move(dh));
                    }
                }
            }
        }
    }

    void organizeSubsets(std::vector<def_subset> &nSubsets, const std::vector<def_atom> &nAtoms, const std::vector<def_bond> &nBonds)
    {
        std::unordered_map<uint32_t, std::vector<uint32_t>> atomBonds; // atomIndex > atomsIndices

        for (uint32_t b = 0; b < nBonds.size(); ++b)
        {
            const def_bond &bond = nBonds[b];

            atomBonds[bond.centralAtomIdx].emplace_back(bond.bondingAtomIdx);
            atomBonds[bond.bondingAtomIdx].emplace_back(bond.centralAtomIdx);
        }

        if (nAtoms.size() == 2 && nAtoms[0].ZIndex == 1 && nAtoms[1].ZIndex == 1) // H2
        {
            def_subset sub{};
            sub.mainAtomIdx = 0;
            sub.hydrogensIdx = {1};

            nSubsets.emplace_back(std::move(sub));
            return;
        }

        if (nAtoms.size() == 1 && nAtoms[0].ZIndex == 1) // H
        {
            def_subset sub{};
            sub.mainAtomIdx = 0;
            sub.hydrogensIdx = {};

            nSubsets.emplace_back(std::move(sub));
            return;
        }

        nSubsets.reserve(nAtoms.size()); // upper bound
        for (uint32_t a = 0; a < nAtoms.size(); ++a)
        {
            const def_atom &atom = nAtoms[a];

            if (atom.ZIndex == 1)
                continue; // Hydrogen is not a central atom other than in H2

            def_subset nSubset{};
            nSubset.mainAtomIdx = a;

            std::vector<uint32_t> connected = atomBonds[a];

            for (uint32_t i = 0; i < connected.size(); ++i)
            {
                if (nAtoms[connected[i]].ZIndex == 1)
                    nSubset.hydrogensIdx.emplace_back(connected[i]);
                else
                    nSubset.connectedIdx.emplace_back(connected[i]);
            }

            nSubsets.emplace_back(std::move(nSubset));
        }

        std::unordered_map<uint32_t, uint32_t> atomToSubset;
        for (uint32_t s = 0; s < nSubsets.size(); ++s)
            atomToSubset[nSubsets[s].mainAtomIdx] = s;

        for (uint32_t s = 0; s < nSubsets.size(); ++s)
        {
            uint32_t main = nSubsets[s].mainAtomIdx;

            for (const auto &bond : nBonds)
            {
                if (bond.centralAtomIdx != main)
                    continue;
                uint32_t neigh = bond.bondingAtomIdx;

                if (nAtoms[neigh].ZIndex == 1)
                    continue;

                auto it = atomToSubset.find(neigh);
                if (it == atomToSubset.end())
                    continue;

                uint32_t neighSubset = it->second;

                nSubsets[s].bondingSubset = neighSubset;
                nSubsets[neighSubset].bondedSubset = s;
            }
        }
    }

    void addImplicitHydrogens(std::vector<def_atom> &nAtoms, std::vector<def_bond> &nBonds)
    {
        uint32_t originalSize = nAtoms.size();

        for (uint32_t i = 0; i < originalSize; ++i)
        {
            const def_atom &a = nAtoms[i];
            if (a.ZIndex == 1)
                continue;

            if (!a.hydrogenize) continue;

            int8_t bonds = constants::getUsualBonds(a.ZIndex);
            if (bonds == 0)
                continue; // Noble Gas
            int32_t numHydrogen = bonds - a.nBonds;
            if (numHydrogen <= 0)
                continue;

            numHydrogen = a.nHydrogens == 0 ? numHydrogen : a.nHydrogens;

            for (int32_t h = 0; h < numHydrogen; ++h)
            {
                nAtoms.emplace_back(def_atom{1, 0, 0, 1});
                uint32_t hIdx = nAtoms.size() - 1;

                def_bond hBond{};
                hBond.centralAtomIdx = i;
                hBond.bondingAtomIdx = hIdx;
                hBond.type = BondType::SINGLE;
                nBonds.emplace_back(std::move(hBond));

                nAtoms[i].nBonds += 1;
            }
        }
    }

    sf::Vector3f rotateDirection(const sf::Vector3f &v, const sf::Vector3f &axis, float angle)
    {
        const float c = std::cos(angle);
        const float s = std::sin(angle);
        const float t = 1.0f - c;

        const float x = axis.x, y = axis.y, z = axis.z;
        const float xx = x * x, yy = y * y, zz = z * z;
        const float xy = x * y, xz = x * z, yz = y * z;

        const sf::Vector3f row1(t * xx + c, t * xy - s * z, t * xz + s * y);
        const sf::Vector3f row2(t * xy + s * z, t * yy + c, t * yz - s * x);
        const sf::Vector3f row3(t * xz - s * y, t * yz + s * x, t * zz + c);

        return sf::Vector3f(
            v.x * row1.x + v.y * row1.y + v.z * row1.z,
            v.x * row2.x + v.y * row2.y + v.z * row2.z,
            v.x * row3.x + v.y * row3.y + v.z * row3.z);
    }
} // namespace sim
