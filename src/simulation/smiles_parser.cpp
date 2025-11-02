#include "smiles_parser.hpp"

#include <stack>
#include <iostream>

namespace sim
{
    using namespace sim::fun;
    molecule_structure sim::parseSMILES(const std::string& molecule)
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

        std::map<int32_t, size_t> ringOpen;  // ring number > atom index
        std::vector<size_t> ringIndices;  // atom indices
        std::vector<size_t> ringBondIndices; 

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

                    nAtoms[openAtom].nBonds += static_cast<uint32_t>(bondType) + 1;
                    nAtoms[closeAtom].nBonds += static_cast<uint32_t>(bondType) + 1;

                    if (isAromatic && ringBondIndices.size() >= 5 && ringIndices.size() <= 7)
                    {
                        std::vector<size_t> ringBondIdx;

                        for (size_t k = 0; k < ringIndices.size(); ++k)
                        {
                            size_t a1 = ringIndices[k];
                            size_t a2 = ringIndices[(k+1) % ringIndices.size()];

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
                        nAtoms[closeAtom].nBonds += static_cast<uint32_t>(bondType) + 1;
                        nAtoms[openAtom].nBonds += static_cast<uint32_t>(bondType) + 1;
                    }

                    def_bond nBond{};
                    nBond.centralAtomIdx = openAtom;
                    nBond.bondingAtomIdx = closeAtom;
                    nBond.type = bondType;

                    ringOpen.erase(ringID);
                    ringIndices.clear();
                    ringBondIndices.clear();

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

            if (c == '@')
            {
                if (molecule[i + 1] == '@')
                {
                    continue;
                }

                continue;
            }

            // Different molecule
            if (c == '.')
            {

            }

            // Direction up or down
            if (c == '\\')
            {
            
            }
            else if (c == '/')
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

                    if (branchStack.empty() && aromatic)
                    {
                        ringBondIndices.emplace_back(nBonds.size() - 1);
                        ringIndices.emplace_back(prevAtom);
                    }

                    nAtoms[prevAtom].nBonds += static_cast<uint32_t>(bondType) + 1;
                    nAtoms[currentAtom].nBonds += static_cast<uint32_t>(bondType) + 1;
                }
                
                prevAtom = currentAtom;
                bondType = BondType::SINGLE;
            }
        }
        
        if (ringOpen.size() > 0)
            throw std::runtime_error("Unclosed rings!");

        addImplicitHydrogens(nAtoms, nBonds);
        organizeSubsets(nSubsets, nAtoms, nBonds);
        positionAtoms(nAtoms, getAngles(nSubsets, nAtoms, nBonds), nPositons);

        nStructure.atoms = std::move(nAtoms);
        nStructure.bonds = std::move(nBonds);
        nStructure.subsets = std::move(nSubsets);
        nStructure.positons = std::move(nPositons);

        return nStructure;
    }

    std::vector<float> sim::getAngles(const std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds)
    {
        std::vector<float> angles(nSubsets.size());

        for (const def_subset& sub : nSubsets)
        {
            std::vector<size_t> connected = sub.connectedIdx;

            if (sub.bondingSubset != SIZE_MAX) connected.emplace_back(sub.bondingSubset);
            if (sub.bondedSubset != SIZE_MAX) connected.emplace_back(sub.bondedSubset);

            std::vector<uint8_t> ZIndices;
            std::vector<BondType> bondTypes;
            ZIndices.reserve(sub.connectedIdx.size());
            bondTypes.reserve(sub.connectedIdx.size());

            // no angle for one bonds;
            if (sub.connectedIdx.size() < 2) continue; 

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

            angles.emplace_back(constants::getAngles(sub.mainAtomIdx, ZIndices, bondTypes));
        }

        return angles;
    }
    
    // Returns atom and it's bonds for later stages
    void sim::organizeSubsets(std::vector<def_subset>& nSubsets, const std::vector<def_atom>& nAtoms, const std::vector<def_bond>& nBonds)
    {
        std::map<size_t, std::vector<size_t>> atomBonds; // atomIndex > atomsIndices

        for (size_t b = 0; b < nBonds.size(); ++b)
        {
            const def_bond& bond = nBonds[b];

            atomBonds[bond.centralAtomIdx].emplace_back(bond.bondingAtomIdx);
            atomBonds[bond.bondingAtomIdx].emplace_back(bond.centralAtomIdx);
        }

        for (size_t a = 0; a < nAtoms.size(); ++a)
        {
            
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

    void sim::positionAtoms(const std::vector<def_atom>& nAtoms, const std::vector<float>& angles, std::vector<sf::Vector3f>& positions)
    {
        for (size_t i = 0; i < nAtoms.size(); ++i)
        {
            const def_atom& a = nAtoms[i];

        }
    }
} // namespace sim
