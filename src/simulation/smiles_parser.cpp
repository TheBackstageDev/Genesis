#include "smiles_parser.hpp"

#include <stack>
#include <iostream>

namespace sim
{
    using namespace sim::fun;
    molecule_structure sim::parseSMILES(const std::string& molecule, bool implicitHydrogens)
    {
        molecule_structure nStructure{};
        std::vector<def_bond> nBonds{};
        std::vector<def_subset> nSubsets{};
        std::vector<def_atom> nAtoms{};
        std::vector<sf::Vector3f> nPositons{};

        if (molecule.empty()) return nStructure;

        size_t prevAtom = SIZE_MAX;
        BondType bondType = BondType::SINGLE; // Current Bond
        
        std::stack<size_t> branchStack;

        std::unordered_map<int32_t, size_t> ringOpen;  // ring number > atom index
        std::unordered_map<size_t, std::vector<size_t>> ringIndices;  // atom indices
        std::unordered_map<size_t, std::vector<size_t>> ringBondIndices; 

        std::vector<std::vector<size_t>> rings;

        for (size_t i = 0; i < molecule.size(); ++i)
        {
            const uint8_t c = molecule[i];
            
            if (isspace(c)) { ++i; continue; }

            if (c == '=' || c == '#')
            {
                bondType = c == '=' ? BondType::DOUBLE : BondType::TRIPLE;
                continue;
            }

            // is Ring
            if (isdigit(c) || c == '%')
            {
                int32_t ringID = c - '0';

                // Ring number > 9
                if (c == '%')
                {
                    ringID = (molecule[i+1]-'0') * 10 + (molecule[i+2]-'0'); 
                    i += 2;
                }

                if (ringOpen.count(ringID))
                {
                    size_t openAtom = ringOpen[ringID];
                    size_t closeAtom = nAtoms.size() - 1;

                    bool isAromatic = nAtoms[openAtom].aromatic && nAtoms[closeAtom].aromatic;

                    nAtoms[openAtom].nBonds += static_cast<uint32_t>(bondType);
                    nAtoms[closeAtom].nBonds += static_cast<uint32_t>(bondType);

                    if (isAromatic && ringBondIndices[ringID].size() >= 5 && ringIndices[ringID].size() <= 7)
                    {
                        std::vector<size_t> ringBondIdx;

                        for (size_t k = 0; k < ringIndices[ringID].size(); ++k)
                        {
                            size_t a1 = ringIndices[ringID][k];
                            size_t a2 = ringIndices[ringID][(k+1) % ringIndices[ringID].size()];

                            for (size_t b = 0; b < nBonds.size(); ++b)
                            {
                                const def_bond& bb = nBonds[b];
                                if ((bb.centralAtomIdx == a1 && bb.bondingAtomIdx == a2) ||
                                    (bb.centralAtomIdx == a2 && bb.bondingAtomIdx == a1))
                                {
                                    ringBondIdx.push_back(b);
                                    break;
                                }
                            }
                        }

                        bool makeDouble = false;
                        for (size_t bIdx : ringBondIdx)
                        {
                            if (makeDouble)
                            {
                                nBonds[bIdx].type = BondType::DOUBLE;
                                nAtoms[nBonds[bIdx].centralAtomIdx].nBonds += 1;
                                nAtoms[nBonds[bIdx].bondingAtomIdx].nBonds += 1;
                            }
                            makeDouble = !makeDouble;
                        }

                        bondType = !makeDouble ? BondType::DOUBLE : BondType::SINGLE;
                        nAtoms[openAtom].nBonds += static_cast<uint32_t>(bondType);
                        nAtoms[closeAtom].nBonds += static_cast<uint32_t>(bondType);
                    }

                    def_bond nBond{};
                    nBond.centralAtomIdx = closeAtom;
                    nBond.bondingAtomIdx = openAtom;
                    nBond.type = bondType;

                    ringIndices[ringID].insert(ringIndices[ringID].end(), closeAtom);
                    rings.emplace_back(ringIndices[ringID]);

                    nBonds.emplace_back(std::move(nBond));
                }
                else
                {
                    ringOpen[ringID] = nAtoms.size() - 1;
                }

                bondType = BondType::SINGLE;
                continue;
            }

            if (c == '(')
            {
                branchStack.emplace(prevAtom);
                continue;
            }
            else if (c == ')')
            {
                if (!branchStack.empty())
                {
                    prevAtom = branchStack.top();
                    branchStack.pop();
                }
                continue;
            }

            if (c == '[')
            {
                
            }
            else if (c == ']')
            {
            
            }

            // Different molecule
            if (c == '.')
            {

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

                nAtoms.emplace_back(std::move(nAtom));

                size_t currentAtom = nAtoms.size() - 1;

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
                
                if (ringOpen.size() > 0 && prevAtom != SIZE_MAX && 
                    std::find(ringIndices[ringOpen.size()].begin(), ringIndices[ringOpen.size()].end(), prevAtom) 
                    == ringIndices[ringOpen.size()].end())
                {
                    ringBondIndices[ringOpen.size()].emplace_back(nBonds.size() - 1);
                    ringIndices[ringOpen.size()].emplace_back(prevAtom);
                }

                prevAtom = currentAtom;
                bondType = BondType::SINGLE;
            }
        }
        
        if (implicitHydrogens)
            addImplicitHydrogens(nAtoms, nBonds);
            
        nPositons.resize(nAtoms.size());
        organizeSubsets(nSubsets, nAtoms, nBonds);
        getAngles(nSubsets, nAtoms, nBonds);    
        positionAtoms(molecule, rings, nAtoms, nSubsets, nPositons);

        nStructure.atoms = std::move(nAtoms);
        nStructure.bonds = std::move(nBonds);
        nStructure.subsets = std::move(nSubsets);
        nStructure.positons = std::move(nPositons);

        return nStructure;
    }

    std::vector<float> sim::getAngles(std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds)
    {
        std::vector<float> angles(nSubsets.size());

        for (size_t s = 0; s < nSubsets.size(); ++s)
        {
            def_subset& sub = nSubsets[s];

            std::vector<size_t> connected = sub.connectedIdx;
            connected.insert(connected.end(), sub.hydrogensIdx.begin(), sub.hydrogensIdx.end());

            if (sub.connectedIdx.size() < 2 && sub.hydrogensIdx.size() < 2) continue; 

            std::vector<uint8_t> ZIndices;
            std::vector<BondType> bondTypes;
            ZIndices.reserve(connected.size());
            bondTypes.reserve(connected.size());

            // no angle for one bonds;

            for (size_t b = 0; b < connected.size(); ++b)
            {
                auto bondIt = std::find_if(nBonds.begin(), nBonds.end(), 
                    [&](const def_bond& a){ return a.bondingAtomIdx == connected[b] && a.centralAtomIdx == sub.mainAtomIdx; });

                if (bondIt != nBonds.end()) 
                {
                    const def_bond& bond = *bondIt;
                    bondTypes.emplace_back(bond.type);
                }

                ZIndices.emplace_back(nAtoms[connected[b]].ZIndex);
            }

            float idealAngle = constants::getAngles(nAtoms[sub.mainAtomIdx].ZIndex, ZIndices, bondTypes);
            sub.idealAngle = idealAngle;
            angles[s] = idealAngle;
        }

        return angles;
    }
    
    void sim::organizeSubsets(std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds)
    {
        std::unordered_map<size_t, std::vector<size_t>> atomBonds; // atomIndex > atomsIndices

        for (size_t b = 0; b < nBonds.size(); ++b)
        {
            const def_bond& bond = nBonds[b];

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

        nSubsets.reserve(nAtoms.size());  // upper bound
        for (size_t a = 0; a < nAtoms.size(); ++a)
        {
            const def_atom& atom = nAtoms[a];

            if (atom.ZIndex == 1) continue; // Hydrogen is not a central atom other than in H2

            def_subset nSubset{};
            nSubset.mainAtomIdx = a;
            
            std::vector<size_t> connected = atomBonds[a];

            for (size_t i = 0; i < connected.size(); ++i)
            {
                if (nAtoms[connected[i]].ZIndex == 1)
                    nSubset.hydrogensIdx.emplace_back(connected[i]);
                else
                    nSubset.connectedIdx.emplace_back(connected[i]);
            }

            nSubsets.emplace_back(std::move(nSubset));
        }

        std::unordered_map<size_t, size_t> atomToSubset;
        for (size_t s = 0; s < nSubsets.size(); ++s)
            atomToSubset[nSubsets[s].mainAtomIdx] = s;

        for (size_t s = 0; s < nSubsets.size(); ++s)
        {
            size_t main = nSubsets[s].mainAtomIdx;

            for (const auto& bond : nBonds)
            {
                if (bond.centralAtomIdx != main) continue;
                size_t neigh = bond.bondingAtomIdx;

                if (nAtoms[neigh].ZIndex == 1) continue;

                auto it = atomToSubset.find(neigh);
                if (it == atomToSubset.end()) continue;

                size_t neighSubset = it->second;

                nSubsets[s].bondingSubset = neighSubset;   
                nSubsets[neighSubset].bondedSubset = s; 
            }
        }
    }

    void sim::addImplicitHydrogens(std::vector<def_atom>& nAtoms, std::vector<def_bond>& nBonds)
    {
        size_t originalSize = nAtoms.size(); 

        for (size_t i = 0; i < originalSize; ++i)
        {
            const def_atom& a = nAtoms[i];
            if (a.ZIndex == 1) continue;

            int8_t bonds = constants::getUsualBonds(a.ZIndex);
            int32_t numHydrogen = bonds - a.nBonds;
            if (numHydrogen <= 0) continue;

            for (int32_t h = 0; h < numHydrogen; ++h)
            {
                nAtoms.emplace_back(def_atom{1, 1, 0, 1});
                size_t hIdx = nAtoms.size() - 1;

                def_bond hBond{};
                hBond.centralAtomIdx = i;
                hBond.bondingAtomIdx = hIdx;
                hBond.type = BondType::SINGLE;
                nBonds.emplace_back(std::move(hBond));

                nAtoms[i].nBonds += 1;
            }
        }
    }

    void sim::positionAtoms(const std::string& SMILES, const std::vector<std::vector<size_t>>& rings, 
                            const std::vector<def_atom>& nAtoms, const std::vector<def_subset>& nSubsets, 
                            std::vector<sf::Vector3f>& positions)
    {
        positions.assign(nAtoms.size(), sf::Vector3f(0, 0, 0));
        if (nAtoms.empty()) return;

        std::unordered_map<size_t, size_t> atomToSubset;
        for (size_t s = 0; s < nSubsets.size(); ++s)
            atomToSubset[nSubsets[s].mainAtomIdx] = s;

        std::unordered_map<size_t, float> ringDir;

        for (const auto& ring : rings) 
        {
            if (ring.size() < 3) continue;

            float turnDeg = 360.f / static_cast<float>(ring.size());
            float turnRad = turnDeg * RADIAN;

            for (size_t i = 0; i < ring.size(); ++i) {
                size_t a = ring[i];                       

                ringDir[a] = turnRad;                         
            }
        }

        BondType bondType = BondType::SINGLE;
        sf::Vector3f currentPos(0, 0, 0);
        sf::Vector3f currentDir(1, 0, 0);

        std::stack<sf::Vector3f> branchPosStack;
        std::stack<sf::Vector3f> branchDirStack;
        std::stack<size_t>       branchAtomStack;

        std::map<int32_t, size_t> ringOpen;   
        std::unordered_map<size_t, bool> flipMap;
        
        size_t atomIdx = 0;
        size_t prevAtom = SIZE_MAX; 

        auto placeAtom = [&](size_t idx, const sf::Vector3f& pos) {
            positions[idx] = pos;
            currentPos = positions[idx];
        };

        auto rotateDirection = [&](float radians) {
            float c = std::cos(radians);
            float s = std::sin(radians);
            currentDir = sf::Vector3f(
                currentDir.x * c - currentDir.y * s,
                currentDir.x * s + currentDir.y * c,
                0.0f
            ).normalized();
        };

        for (size_t i = 0; i < SMILES.size(); ++i)
        {
            const uint8_t c = SMILES[i];
            
            if (isspace(c)) { ++i; continue; }

            if (c == '=' || c == '#')
            {
                bondType = c == '=' ? BondType::DOUBLE : BondType::TRIPLE;
                continue;
            }

            if (c == '[')
            {
                
            }
            else if (c == ']')
            {
            
            }

            float ideal = 0.f;
            auto subIt = atomToSubset.find(prevAtom);
            if (subIt != atomToSubset.end()) 
            {
                ideal = nSubsets[subIt->second].idealAngle;
            }

            // Direction up or down
            if (c == '/') {

                if (prevAtom != SIZE_MAX) 
                {
                    rotateDirection(-M_PI / 3);  // 60° up
                }
                continue;
            }
            if (c == '\\') {
                if (prevAtom != SIZE_MAX) 
                {
                    rotateDirection(M_PI / 3);  // 60° down
                }
                continue;
            }
            
            if (c == '@')
            {
                if (SMILES[i + 1] == '@')
                {
                    continue;
                }

                continue;
            }

            if (c == '(') {
                if (prevAtom != SIZE_MAX) {
                    branchAtomStack.push(prevAtom);
                    branchPosStack.push(currentPos);
                    branchDirStack.push(currentDir);
                }
                continue;
            }
            if (c == ')') {
                if (!branchAtomStack.empty()) {
                    prevAtom   = branchAtomStack.top();   branchAtomStack.pop();
                    currentPos = branchPosStack.top();    branchPosStack.pop();
                    currentDir = branchDirStack.top();    branchDirStack.pop();

                    if (ideal > 0.f) {
                        bool& flip = flipMap[prevAtom];
                        float sign = flip ? 1.f : -1.f;
                        rotateDirection(sign * ideal);
                        flip = !flip;
                    } else {
                        rotateDirection(M_PI);  
                    }
                }
                continue;
            }

            if (std::isdigit(c) || c == '%') 
            {
                int32_t ringID = c - '0';
                if (c == '%') {
                    ringID = (SMILES[i + 1] - '0') * 10 + (SMILES[i + 2] - '0');
                    i += 2;
                }

                if (ringOpen.count(ringID)) 
                {                    
                    prevAtom = atomIdx;                 
                    ringOpen.erase(ringOpen[ringID]);
                } 
                else 
                {                
                    rotateDirection(M_PI * 1.5f);  
                    ringOpen[ringID] = prevAtom;                          
                }
                bondType = BondType::SINGLE;
                continue;
            }

            if (isalpha(c))
            {
                if (atomIdx >= nAtoms.size()) break;

                float bondLength = 2.f * MULT_FACTOR;
                
                if (prevAtom != SIZE_MAX && nAtoms[atomIdx].ZIndex != 1 && c != 'H')
                {
                    const def_atom& a1 = nAtoms[prevAtom];
                    const def_atom& a2 = nAtoms[atomIdx];
                    bondLength = constants::getBondLength(a1.ZIndex, a2.ZIndex, bondType);

                    static bool flip = false;

                    auto ringIt = ringDir.find(prevAtom);
                    if (ringIt != ringDir.end()) 
                    {
                        rotateDirection(ringIt->second);         
                    }
                    else if (ideal > 0.f)
                    {
                        float sign = flip ? 1.f : -1.f;
                        rotateDirection(sign * ideal * 0.5f);
                        flip = !flip;
                    }
                    
                    bondType = BondType::SINGLE;
                }
                
                placeAtom(atomIdx, currentPos + currentDir * bondLength);
                ++atomIdx;
                prevAtom = atomIdx - 1;
            }
        }

        // Hydrogen Placement
        for (size_t s = 0; s < nSubsets.size(); ++s)
        {
            const def_subset& sub = nSubsets[s];
            const std::vector<size_t> hydrogens = sub.hydrogensIdx;
            const size_t centralAtom = sub.mainAtomIdx;

            for (size_t h = 0; h < hydrogens.size(); ++h)
            {
                const size_t hydrogenIdx = hydrogens[h];
                float bondLength = 1.3f * MULT_FACTOR;

                positions[hydrogenIdx] = positions[centralAtom];

                if (sub.bondedSubset != SIZE_MAX && sub.bondingSubset != SIZE_MAX)
                {
                    const size_t central1 = nSubsets[sub.bondedSubset].mainAtomIdx;
                    const size_t central2 = nSubsets[sub.bondingSubset].mainAtomIdx;
                    
                    currentDir = (positions[central2] - positions[central1]).cross(sf::Vector3f(0, 0, 1)).normalized();

                    if (hydrogens.size() > 1 && ringDir.find(centralAtom) != ringDir.end())
                    {
                        rotateDirection(-20 * RADIAN);
                        rotateDirection(60 * h * RADIAN);
                    }
                    else if (hydrogens.size() > 1)
                    {
                        rotateDirection(180 * h * RADIAN);
                    }
                    else if (ringDir.find(centralAtom) != ringDir.end())
                    {
                        rotateDirection(-90 * RADIAN);
                        rotateDirection(sub.idealAngle);
                    }
                }
                else if (sub.bondedSubset == SIZE_MAX && sub.bondingSubset != SIZE_MAX) // Terminal Atom
                {
                    const size_t centrallast = nSubsets[sub.bondingSubset].mainAtomIdx;

                    currentDir = (positions[hydrogenIdx] - positions[centrallast]).normalized();

                    if (hydrogens.size() > 1)
                        rotateDirection(60 * h * RADIAN);
                }
                else if (sub.bondingSubset == SIZE_MAX && sub.bondedSubset != SIZE_MAX) // Initial Atom
                {
                    const size_t centralNext = nSubsets[sub.bondedSubset].mainAtomIdx;

                    currentDir = (positions[hydrogenIdx] - positions[centralNext]).normalized();

                    if (hydrogens.size() > 1)
                        rotateDirection(60 * h * RADIAN);
                }
                else if (hydrogens.size() > 1)
                {
                    rotateDirection(sub.idealAngle * h);   
                }

                positions[hydrogenIdx] += currentDir * bondLength;
            }
        }
    }
} // namespace sim
