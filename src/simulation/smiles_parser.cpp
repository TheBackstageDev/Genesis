#include "smiles_parser.hpp"

#include <stack>
#include <iostream>

namespace sim
{
    using namespace sim::fun;
    molecule_structure sim::parseSMILES(const std::string &molecule, bool implicitHydrogens)
    {
        molecule_structure nStructure{};
        std::vector<def_bond> nBonds{};
        std::vector<def_subset> nSubsets{};
        std::vector<def_atom> nAtoms{};
        std::vector<angle> nAngles{};
        std::vector<dihedral_angle> nDihedralAngles{};
        std::vector<sf::Vector3f> npositions{};

        if (molecule.empty())
            return nStructure;

        uint32_t prevAtom = SIZE_MAX;
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
                if (prevAtom != SIZE_MAX)
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
                char chirality = 0;
                uint32_t k = 0;

                if (std::isdigit(atomStr[k]))
                    isotope = std::stoi(atomStr.substr(k), &k) - 1;

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
                    chirality = atomStr[k++];
                    if (k < atomStr.size() && atomStr[k] == '@') { chirality = '@@'; ++k; }
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

                if (prevAtom != SIZE_MAX)
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
                prevAtom = SIZE_MAX;
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
                nAtom.NIndex = ZIndex;
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
                            nAtom.NIndex = newZ; // Usually yes
                            ++i;                            
                            aromatic = false;                 
                        }
                    }
                }

                nAtoms.emplace_back(std::move(nAtom));

                uint32_t currentAtom = nAtoms.size() - 1;

                if (prevAtom != SIZE_MAX)
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
        organizeAngles(nSubsets, nAtoms, nBonds, nDihedralAngles, nAngles);
        positionAtoms(molecule, nBonds, rings, nAtoms, nSubsets, npositions, nAngles);

        nStructure.atoms = std::move(nAtoms);
        nStructure.bonds = std::move(nBonds);
        nStructure.subsets = std::move(nSubsets);
        nStructure.angles = std::move(nAngles);
        nStructure.dihedral_angles = std::move(nDihedralAngles);
        nStructure.positions = std::move(npositions);

        return nStructure;
    }

    struct pair_hash {
        uint32_t operator() (const std::pair<int64_t, int64_t>& p) const {
            return (uint32_t)(p.first << 32) | p.second;
        }
    };

    void sim::organizeAngles(std::vector<def_subset> &nSubsets, const std::vector<def_atom> &nAtoms, const std::vector<def_bond> &nBonds,
                            std::vector<dihedral_angle>& dihedral_angles, std::vector<angle>& angles)
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
                    ang.rad = constants::getAngles(nAtoms[B].ZIndex, Z, type);
                    ang.K   = constants::getAngleHarmonicConstant(ang.A, ang.B, ang.C);
                    angles.push_back(ang);
                }
            }

            for (uint32_t c_idx = 0; c_idx < neigh.size(); ++c_idx) {
                const uint32_t C = neigh[c_idx];

                uint32_t a = B, b = C;
                if (a > b) std::swap(a, b);
                const def_bond* bc_bond = nullptr;
                auto it = bond_map.find({a, b});
                if (it != bond_map.end()) bc_bond = &it->second;
                if (!bc_bond) continue;

                if (nAtoms[a].aromatic || nAtoms[b].aromatic) continue;

                const auto& neigh_B = neigh; 

                std::vector<uint32_t> neigh_C;
                for (const auto& bond : nBonds) 
                {
                    if (bond.centralAtomIdx == C) neigh_C.push_back(bond.bondingAtomIdx);
                    else if (bond.bondingAtomIdx == C) neigh_C.push_back(bond.centralAtomIdx);
                }

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

                        dihedral_angles.push_back(dh);
                    }
                }
            }
        }
    }

    void sim::organizeSubsets(std::vector<def_subset> &nSubsets, const std::vector<def_atom> &nAtoms, const std::vector<def_bond> &nBonds)
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

    void sim::addImplicitHydrogens(std::vector<def_atom> &nAtoms, std::vector<def_bond> &nBonds)
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

    void sim::positionAtoms(const std::string &SMILES, std::vector<def_bond>& nBonds, const std::vector<std::vector<uint32_t>> &rings,
                            const std::vector<def_atom> &nAtoms, const std::vector<def_subset> &nSubsets,
                            std::vector<sf::Vector3f> &positions, std::vector<angle>& angles)
    {
        positions.assign(nAtoms.size(), sf::Vector3f(0, 0, 0));
        if (nAtoms.empty())
            return;

        std::unordered_map<uint32_t, uint32_t> atomToSubset;
        for (uint32_t s = 0; s < nSubsets.size(); ++s)
            atomToSubset[nSubsets[s].mainAtomIdx] = s;

        std::unordered_map<uint32_t, float> ringDir;
        std::unordered_map<int32_t, std::vector<uint32_t>> ringBonds; 
        std::unordered_map<int32_t, uint32_t> ringOpen; 
        std::vector<uint32_t> branchStack;

        BondType curBond = BondType::SINGLE;
        sf::Vector3f pos(0,0,0), dir(1,0,0);
        std::stack<sf::Vector3f> posStk, dirStk;
        std::stack<uint32_t> atomStk;     
        
        std::unordered_map<uint32_t,bool> flipMap;
        uint32_t atomIdx = 0;             
        uint32_t prevAtom = SIZE_MAX;

        auto place = [&](uint32_t idx, const sf::Vector3f& p){
            positions[idx] = p; pos = p;
            pos.z += 0.001f * idx;
        };
        auto rotate2D = [&](float rad){
            float c = std::cos(rad), s = std::sin(rad);
            dir = sf::Vector3f(dir.x*c - dir.y*s,
                            dir.x*s + dir.y*c, 0.f);
            if (dir.length() > 0.f)
                dir = dir.normalized();
        };

        for (const auto &ring : rings)
        {
            if (ring.size() < 3) continue;
            float turn = 360.f / ring.size() * RADIAN;
            for (uint32_t a : ring) ringDir[a] = turn;
        }

        for (uint32_t i = 0; i < SMILES.size(); ++i)
        {
            char c = SMILES[i];
            if (isspace(c)) continue;
            if (c=='='||c=='#') { curBond = (c=='=')?BondType::DOUBLE:BondType::TRIPLE; continue; }

            float ideal = 0.f;
            if (prevAtom != SIZE_MAX)
            {
                auto it = atomToSubset.find(prevAtom);
                if (it != atomToSubset.end()) ideal = nSubsets[it->second].idealAngle;
            }

            // ----- cis / trans -----
            if (c=='/') { rotate2D(-M_PI/3); continue; }
            if (c=='\\'){ rotate2D( M_PI/3); continue; }

            // ----- branches -----
            if (c=='(')
            {
                if (prevAtom != UINT32_MAX)
                {
                    atomStk.push(prevAtom);
                    posStk.push(pos);
                    dirStk.push(dir);

                    uint32_t A = prevAtom;
                    uint32_t B = UINT32_MAX;
                    uint32_t C = UINT32_MAX;

                    for (const auto& bond : nBonds)
                    {
                        if (bond.centralAtomIdx == A && bond.bondingAtomIdx != prevAtom)
                        {
                            if (B == UINT32_MAX) B = bond.bondingAtomIdx;
                            else if (C == UINT32_MAX) { C = bond.bondingAtomIdx; break; }
                        }
                        else if (bond.bondingAtomIdx == A && bond.centralAtomIdx != prevAtom)
                        {
                            if (B == UINT32_MAX) B = bond.centralAtomIdx;
                            else if (C == UINT32_MAX) { C = bond.centralAtomIdx; break; }
                        }
                    }

                    if (B != UINT32_MAX && C != UINT32_MAX)
                    {
                        sf::Vector3f v1 = positions[A] - positions[B];
                        sf::Vector3f v2 = positions[A] - positions[C];
                        sf::Vector3f normal = v1.cross(v2);
                        if (normal.length() > 1e-6f)
                        {
                            normal = normal.normalized();

                            sf::Vector3f newDir = rotateDirection(dir, normal, -30.0f * RADIAN);
                            dir = newDir;
                        }
                    }
                    else if (B != UINT32_MAX)
                    {
                        if (ringDir.count(A))
                        {
                            sf::Vector3f ringNormal(0, 0, 1);
                            sf::Vector3f axis = positions[A] - positions[B];
                            if (axis.length() > 1e-6f) axis = axis.normalized();
                            dir = rotateDirection(axis, ringNormal, -30.0f * RADIAN);
                        }
                    }
                }
                continue;
            }
            if (c==')')
            {
                if (!atomStk.empty())
                {
                    prevAtom = atomStk.top(); atomStk.pop();
                    pos = posStk.top(); posStk.pop();
                    dir = dirStk.top(); dirStk.pop();

                    float ideal = 0.f;
                    auto it = atomToSubset.find(prevAtom);
                    if (it != atomToSubset.end()) ideal = nSubsets[it->second].idealAngle;

                    if (ideal > 0.f)
                    {
                        bool &fl = flipMap[prevAtom];
                        rotate2D((fl ? 1.f : -1.f) * ideal * 0.5f);
                        fl = !fl;
                    }
                    else
                    {
                        rotate2D(M_PI / 2.f); 
                    }
                }
                continue;
            }

            if (isalpha(c) || c == '[')
            {
                if (atomIdx >= nAtoms.size()) break;
                const def_atom& a = nAtoms[atomIdx];

                for (const auto& [ringID, openAtom] : ringOpen)
                {
                    if (rings.size() > 0 && std::find(rings[ringID - 1].begin(),
                                rings[ringID - 1].end(),
                                prevAtom) == rings[ringID - 1].end() &&
                        std::find(branchStack.begin(), branchStack.end(), prevAtom) == branchStack.end())
                    {
                        ringBonds[ringID - 1].push_back(nBonds.size() - 1);
                    }
                }

                float len = 2.0f * MULT_FACTOR;
                if (prevAtom != SIZE_MAX && a.ZIndex != 1)
                {
                    len = constants::getBondLength(nAtoms[prevAtom].ZIndex, a.ZIndex, curBond);

                    if (ringDir.count(prevAtom))
                    {
                        auto ring_it = std::find_if(rings.begin(), rings.end(), [&](const std::vector<uint32_t>& atoms) {
                            return std::find(atoms.begin(), atoms.end(), prevAtom) != atoms.end();
                        });

                        int32_t ringID = std::distance(rings.begin(), ring_it);

                        float sign =  ringID % 2 == 0 ? 1.f : -1.f;
                        rotate2D(ringDir[prevAtom] * sign);
                    }
                    else if (ideal > 0.f)
                    {
                        bool& fl = flipMap[prevAtom];
                        rotate2D((fl ? 1.f : -1.f) * ideal * 0.3f);
                        fl = !fl;
                    }
                }

                place(atomIdx, pos + dir * len);

                prevAtom = atomIdx;
                ++atomIdx;
                curBond = BondType::SINGLE;
                continue;
            }

            // ----- ring numbers -----
            if (std::isdigit(c) || c == '%')
            {
                int32_t ringID = (c == '%') ? (SMILES[i+1]-'0')*10+(SMILES[i+2]-'0') : (c-'0');
                if (c=='%') i+=2;

                if (ringOpen.count(ringID))
                {
                    uint32_t openAtom  = ringOpen[ringID];
                    
                    bool isFused = false;
                    for (uint32_t i = 1; i <= rings.size(); ++i)
                    {
                        if (ringID != i)
                        {
                            for (const uint32_t& atomID : rings[i - 1])
                            {
                                if (atomID == prevAtom)
                                {
                                    isFused = true;
                                    break; 
                                }
                            }
                        }
                    }
                        
                    if (!isFused)
                    {
                        sf::Vector3f bondDir = (positions[prevAtom] - positions[ringOpen[ringID] + 1]).normalized();
                        pos = positions[prevAtom] + bondDir * 3.f;
                        rotate2D(M_PI);
                    }
                    else
                    {
                        dir = (positions[openAtom] - positions[prevAtom - 1]).normalized();
                        pos = positions[prevAtom] - dir;
                        rotate2D(M_PI);
                    }

                    ringOpen.erase(ringID);
                }
                else
                {
                    ringOpen[ringID] = prevAtom;
                }
                curBond = BondType::SINGLE;
            }
        }

        // Hydrogen Placement
        if (nSubsets.size() > 1)
        for (uint32_t s = 0; s < nSubsets.size(); ++s)
        {
            const def_subset &sub = nSubsets[s];
            const std::vector<uint32_t> hydrogens = sub.hydrogensIdx;
            const uint32_t centralAtom = sub.mainAtomIdx;

            const def_atom& atom = nAtoms[centralAtom];
            const uint8_t Z = atom.ZIndex;
            const char chirality = atom.chirality;

            sf::Vector3f axis(1,0,0);
            if (sub.bondedSubset != SIZE_MAX && sub.bondingSubset != SIZE_MAX)
            {
                const uint32_t c1 = nSubsets[sub.bondedSubset].mainAtomIdx;
                const uint32_t c2 = nSubsets[sub.bondingSubset].mainAtomIdx;
                const sf::Vector3f v1 = (positions[centralAtom] - positions[c1]).normalized();
                const sf::Vector3f v2 = (positions[centralAtom] - positions[c2]).normalized();
                axis = (v1 + v2).length() == 0.f ? axis : (v1 + v2).normalized();
            }
            else if (sub.bondingSubset != SIZE_MAX)
            {
                const uint32_t prev = nSubsets[sub.bondingSubset].mainAtomIdx;
                axis = (positions[centralAtom] - positions[prev]).normalized();
            }
            else if (sub.bondedSubset != SIZE_MAX)
            {
                const uint32_t next = nSubsets[sub.bondedSubset].mainAtomIdx;
                axis = (positions[centralAtom] - positions[next]).normalized();
            }

            sf::Vector3f baseDir = axis;
            auto turnIt = ringDir.find(centralAtom);
            if (turnIt != ringDir.end())
            {
                float turnRad = turnIt->second;
                const sf::Vector3f ringNormal(0, 0, 1);

                baseDir = rotateDirection(axis, ringNormal, -30.0f * RADIAN);
            }

            std::vector<sf::Vector3f> neighDirs;
            for (uint32_t neighIdx : sub.connectedIdx)
            {
                sf::Vector3f dir = (positions[neighIdx] - positions[centralAtom]).normalized();
                neighDirs.push_back(dir);
            }

            const uint32_t nH = sub.hydrogensIdx.size();
            for (uint32_t h = 0; h < nH; ++h)
            {
                const uint32_t hIdx = sub.hydrogensIdx[h];
                positions[hIdx] = positions[centralAtom];

                sf::Vector3f dir = baseDir; 
                const sf::Vector3f ringNormal(0, 0, 1);

                if (nH == 3)
                    dir = rotateDirection(dir, axis, -90.f * RADIAN);

                if (sub.bondedSubset != SIZE_MAX && sub.bondingSubset != SIZE_MAX && ringDir.count(sub.mainAtomIdx) > 0)
                    dir = rotateDirection(dir, axis, 120.f * h * RADIAN);
                else if (nH == 1 && chirality != 0)
                {
                    // Chirality
                    const float cosT = -1.0f/3.0f;
                    const float sinT = std::sqrt(1.0f - cosT*cosT);
                    sf::Vector3f X = (std::abs(axis.x) < 0.9f) ? sf::Vector3f(0,1,0) : sf::Vector3f(1,0,0);
                    X = X == axis ? axis : axis.cross(X).normalized();
                    const sf::Vector3f Y = axis.cross(X);

                    sf::Vector3f tet[4] = {
                        cosT*axis + sinT*X,
                        cosT*axis - sinT*X,
                        cosT*axis + sinT*Y,
                        cosT*axis - sinT*Y
                    };

                    std::vector<sf::Vector3f> assigned = {neighDirs[0], neighDirs[1], neighDirs[2]};
                    sf::Vector3f Hdir = tet[3];

                    std::vector<std::pair<int32_t, uint32_t>> priority;
                    for (uint32_t i = 0; i < 3; ++i)
                    {
                        int pri = nAtoms[sub.connectedIdx[i]].ZIndex;
                        priority.emplace_back(pri, i);
                    }
                    std::sort(priority.rbegin(), priority.rend());

                    sf::Vector3f ordered[3];
                    for (int i = 0; i < 3; ++i)
                        ordered[i] = assigned[priority[i].second];

                    sf::Vector3f v1 = ordered[0] - Hdir;
                    sf::Vector3f v2 = ordered[1] - Hdir;
                    sf::Vector3f v3 = ordered[2] - Hdir;
                    sf::Vector3f cross = v1.cross(v2);
                    float dot = cross.dot(v3);

                    bool isClockwise = (dot > 0);
                    bool wantClockwise = (chirality == '@');

                    if (isClockwise != wantClockwise)
                        std::swap(ordered[1], ordered[2]);

                    tet[0] = ordered[0];
                    tet[1] = ordered[1];
                    tet[2] = ordered[2];
                    tet[3] = -(tet[0] + tet[1] + tet[2]); // H opposite

                    dir = tet[3].normalized();
                }
                else if (nH == 2) 
                    dir = rotateDirection(dir, ringNormal, M_PI * h);
                else if (nH == 3) 
                    dir = rotateDirection(dir, ringNormal, 90.f * h * RADIAN);
                else if (nH == 4)
                {
                    const float cosT = -1.0f / 3.0f;
                    const float sinT = std::sqrt(1.0f - cosT * cosT);

                    sf::Vector3f X = (std::abs(ringNormal.x) < 0.9f) ? sf::Vector3f(0, 1, 0) : sf::Vector3f(1, 0, 0);
                    X = ringNormal.cross(X).normalized();
                    const sf::Vector3f Y = ringNormal.cross(X);

                    sf::Vector3f tet[4] = {
                        cosT * ringNormal + sinT * X,
                        cosT * ringNormal - sinT * X,
                        cosT * ringNormal + sinT * Y,
                        cosT * ringNormal - sinT * Y};

                    if (baseDir.dot(tet[0]) < 0.999f)
                    {
                        sf::Vector3f rotate2DAxis = tet[0].cross(baseDir);
                        if (rotate2DAxis.length() > 1e-6f)
                        {
                            rotate2DAxis = rotate2DAxis.normalized();
                            const float ang = std::acos(std::clamp(tet[0].dot(baseDir), -1.0f, 1.0f));
                            for (int k = 0; k < 4; ++k)
                                tet[k] = rotateDirection(tet[k], rotate2DAxis, ang);
                        }
                    }
                    dir = tet[h];
                }
                else 
                {
                    dir = rotateDirection(dir, ringNormal, 90 * h * RADIAN);
                }

                positions[hIdx] += dir;
            }
        }

        const uint32_t nH = nSubsets[0].hydrogensIdx.size();
        const uint32_t nS = nSubsets.size();

        for (uint32_t i = 0; i < nH; ++i)
        {
            sf::Vector3f dir;

            if (nH == 4)
            {
                const sf::Vector3f tetra[4] = {
                    sf::Vector3f( 1,  1,  1).normalized(),
                    sf::Vector3f( 1, -1, -1).normalized(),
                    sf::Vector3f(-1,  1, -1).normalized(),
                    sf::Vector3f(-1, -1,  1).normalized()
                };
                dir = tetra[i];
            }
            else
            {
                sf::Vector3f ref = (nH == 2) ? sf::Vector3f(1,0,0) : sf::Vector3f(0,0,1);

                float angle = 2.0f * M_PI * i / nH;
                dir = rotateDirection(ref, sf::Vector3f(0,0,1), angle);

                float cosTheta = (nH == 2) ? std::cos(109.471f * RADIAN / 2.0f)
                                        : std::cos(90.0f * RADIAN);
                dir.z = (nH == 2) ? -std::sqrt(1.0f - cosTheta*cosTheta) : dir.z;
                dir = dir.normalized();
            }

            positions[nSubsets[0].hydrogensIdx[i]] = positions[nSubsets[0].mainAtomIdx] + dir;
        }

        if (nAtoms.size() <= 3) return;

        const int32_t max_iters = 30 * static_cast<int32_t>(nAtoms.size());
        constexpr float   dt          = 0.01f;
        constexpr float   repulse     = 2.2f;
        constexpr float   k_spring    = 5.0f;
        constexpr float   convergence = 0.01f;

        for (int32_t it = 0; it < max_iters; ++it)
        {
            std::vector<sf::Vector3f> forces(positions.size(), sf::Vector3f(0.f, 0.f, 0.f));
            float maxF = 0.f;

            for (uint32_t i = 0; i < positions.size(); ++i)
            for (uint32_t j = i + 1; j < positions.size(); ++j)
            {
                sf::Vector3f dr = positions[j] - positions[i];
                float dist = dr.length();

                if (dist < EPSILON) 
                {
                    positions[i] += sf::Vector3f( 0.2f,  0.1f, 0.0f);
                    positions[j] += sf::Vector3f(-0.2f, -0.1f, 0.0f);
                    continue;
                }

                float r0 = constants::getBondLength(nAtoms[i].ZIndex, nAtoms[j].ZIndex, BondType::SINGLE) * 1.4f;
                if (dist < r0)
                {
                    float f = repulse * (r0 - dist) / dist;
                    sf::Vector3f F = dr * (f / dist);
                    forces[i] -= F;
                    forces[j] += F;
                }

                if (dist < COULOMB_CUTOFF)
                {
                    float qq = nAtoms[i].charge * nAtoms[j].charge;
                    if (qq != 0.f)
                    {
                        float factor = COULOMB_K * qq / (dist * dist * dist);
                        sf::Vector3f F = -dr * factor;
                        forces[i] += F;
                        forces[j] -= F;
                    }
                }
            }

            for (const auto& b : nBonds)
            {
                uint32_t c = b.centralAtomIdx;
                uint32_t a = b.bondingAtomIdx;

                sf::Vector3f dr = positions[a] - positions[c];
                float dist = dr.length();

                if (dist < EPSILON) continue;

                float target = constants::getBondLength(
                                nAtoms[c].ZIndex,
                                nAtoms[a].ZIndex,
                                b.type);

                float f = k_spring * (dist - target);
                sf::Vector3f F = dr * (f / dist);

                forces[c] += F;
                forces[a] -= F;
            }

            for (uint32_t i = 0; i < positions.size(); ++i)
            {
                sf::Vector3f delta = forces[i] * dt;
                positions[i] += delta;

                maxF = std::max(maxF, delta.length());
            }

            if (maxF <= convergence) 
                break;
        }
    }

    sf::Vector3f sim::rotateDirection(const sf::Vector3f &v, const sf::Vector3f &axis, float angle)
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
