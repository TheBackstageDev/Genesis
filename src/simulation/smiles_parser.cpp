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
        std::vector<sf::Vector3f> npositions{};

        if (molecule.empty())
            return nStructure;

        size_t prevAtom = SIZE_MAX;
        BondType bondType = BondType::SINGLE; // Current Bond

        std::stack<size_t> branchStack;

        std::unordered_map<int32_t, size_t> ringOpen;                // ring number > atom index
        std::unordered_map<size_t, std::vector<size_t>> ringIndices; // atom indices
        std::unordered_map<size_t, std::vector<size_t>> ringBondIndices;

        std::vector<std::vector<size_t>> rings;

        for (size_t i = 0; i < molecule.size(); ++i)
        {
            const uint8_t c = molecule[i];

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
            if (isdigit(c) || c == '%')
            {
                int32_t ringID = c - '0';

                // Ring number > 9
                if (c == '%')
                {
                    ringID = (molecule[i + 1] - '0') * 10 + (molecule[i + 2] - '0');
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
                            size_t a2 = ringIndices[ringID][(k + 1) % ringIndices[ringID].size()];

                            for (size_t b = 0; b < nBonds.size(); ++b)
                            {
                                const def_bond &bb = nBonds[b];
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

                if (i + 1 < molecule.size())
                {
                    char next = molecule[i + 1];
                    if (isalpha(next))
                    {
                        std::string two = sym + next;

                        two[0] = std::toupper(two[0]);
                        if (two.size() > 1) two[1] = std::tolower(two[1]);

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
                    std::find(ringIndices[ringOpen.size()].begin(), ringIndices[ringOpen.size()].end(), prevAtom) == ringIndices[ringOpen.size()].end())
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

        npositions.resize(nAtoms.size());
        organizeSubsets(nSubsets, nAtoms, nBonds);
        getAngles(nSubsets, nAtoms, nBonds);
        positionAtoms(molecule, rings, nAtoms, nSubsets, npositions);

        nStructure.atoms = std::move(nAtoms);
        nStructure.bonds = std::move(nBonds);
        nStructure.subsets = std::move(nSubsets);
        nStructure.positions = std::move(npositions);

        return nStructure;
    }

    std::vector<float> sim::getAngles(std::vector<def_subset> &nSubsets, const std::vector<def_atom> &nAtoms, const std::vector<def_bond> &nBonds)
    {
        std::vector<float> angles(nSubsets.size());

        for (size_t s = 0; s < nSubsets.size(); ++s)
        {
            def_subset &sub = nSubsets[s];

            std::vector<size_t> connected = sub.connectedIdx;
            connected.insert(connected.end(), sub.hydrogensIdx.begin(), sub.hydrogensIdx.end());

            if (connected.size() < 2)
                continue;

            std::vector<uint8_t> ZIndices;
            std::vector<BondType> bondTypes;
            ZIndices.reserve(connected.size());
            bondTypes.reserve(connected.size());

            // no angle for one bonds;

            for (size_t b = 0; b < connected.size(); ++b)
            {
                auto bondIt = std::find_if(nBonds.begin(), nBonds.end(),
                                           [&](const def_bond &a)
                                           { return a.bondingAtomIdx == connected[b] && a.centralAtomIdx == sub.mainAtomIdx; });

                if (bondIt != nBonds.end())
                {
                    const def_bond &bond = *bondIt;
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

    void sim::organizeSubsets(std::vector<def_subset> &nSubsets, const std::vector<def_atom> &nAtoms, const std::vector<def_bond> &nBonds)
    {
        std::unordered_map<size_t, std::vector<size_t>> atomBonds; // atomIndex > atomsIndices

        for (size_t b = 0; b < nBonds.size(); ++b)
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

        nSubsets.reserve(nAtoms.size()); // upper bound
        for (size_t a = 0; a < nAtoms.size(); ++a)
        {
            const def_atom &atom = nAtoms[a];

            if (atom.ZIndex == 1)
                continue; // Hydrogen is not a central atom other than in H2

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

            for (const auto &bond : nBonds)
            {
                if (bond.centralAtomIdx != main)
                    continue;
                size_t neigh = bond.bondingAtomIdx;

                if (nAtoms[neigh].ZIndex == 1)
                    continue;

                auto it = atomToSubset.find(neigh);
                if (it == atomToSubset.end())
                    continue;

                size_t neighSubset = it->second;

                nSubsets[s].bondingSubset = neighSubset;
                nSubsets[neighSubset].bondedSubset = s;
            }
        }
    }

    void sim::addImplicitHydrogens(std::vector<def_atom> &nAtoms, std::vector<def_bond> &nBonds)
    {
        size_t originalSize = nAtoms.size();

        for (size_t i = 0; i < originalSize; ++i)
        {
            const def_atom &a = nAtoms[i];
            if (a.ZIndex == 1)
                continue;

            int8_t bonds = constants::getUsualBonds(a.ZIndex);
            if (bonds == 0)
                continue; // Noble Gas
            int32_t numHydrogen = bonds - a.nBonds;
            if (numHydrogen <= 0)
                continue;

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

    void sim::positionAtoms(const std::string &SMILES, const std::vector<std::vector<size_t>> &rings,
                            const std::vector<def_atom> &nAtoms, const std::vector<def_subset> &nSubsets,
                            std::vector<sf::Vector3f> &positions)
    {
        positions.assign(nAtoms.size(), sf::Vector3f(0, 0, 0));
        if (nAtoms.empty())
            return;

        std::unordered_map<size_t, size_t> atomToSubset;
        for (size_t s = 0; s < nSubsets.size(); ++s)
            atomToSubset[nSubsets[s].mainAtomIdx] = s;

        std::unordered_map<size_t, float> ringDir;

        for (const auto &ring : rings)
        {
            if (ring.size() < 3)
                continue;

            float turnDeg = 360.f / static_cast<float>(ring.size());
            float turnRad = turnDeg * RADIAN;

            for (size_t i = 0; i < ring.size(); ++i)
            {
                size_t a = ring[i];

                ringDir[a] = turnRad;
            }
        }

        BondType bondType = BondType::SINGLE;
        sf::Vector3f currentPos(0, 0, 0);
        sf::Vector3f currentDir(1, 0, 0);

        std::stack<sf::Vector3f> branchPosStack;
        std::stack<sf::Vector3f> branchDirStack;
        std::stack<size_t> branchAtomStack;

        std::map<int32_t, size_t> ringOpen;
        std::unordered_map<size_t, bool> flipMap;

        size_t atomIdx = 0;
        size_t prevAtom = SIZE_MAX;

        auto placeAtom = [&](size_t idx, const sf::Vector3f &pos)
        {
            positions[idx] = pos;
            currentPos = positions[idx];
        };

        auto rot = [&](float radians)
        {
            float c = std::cos(radians);
            float s = std::sin(radians);
            currentDir = sf::Vector3f(
                             currentDir.x * c - currentDir.y * s,
                             currentDir.x * s + currentDir.y * c,
                             0.0f)
                             .normalized();
        };

        for (size_t i = 0; i < SMILES.size(); ++i)
        {
            const uint8_t c = SMILES[i];

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
            if (c == '/')
            {

                if (prevAtom != SIZE_MAX)
                {
                    rot(-M_PI / 3); // 60° up
                }
                continue;
            }
            if (c == '\\')
            {
                if (prevAtom != SIZE_MAX)
                {
                    rot(M_PI / 3); // 60° down
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

            if (c == '(')
            {
                if (prevAtom != SIZE_MAX)
                {
                    branchAtomStack.push(prevAtom);
                    branchPosStack.push(currentPos);
                    branchDirStack.push(currentDir);
                }
                continue;
            }
            if (c == ')')
            {
                if (!branchAtomStack.empty())
                {
                    prevAtom = branchAtomStack.top();
                    branchAtomStack.pop();
                    currentPos = branchPosStack.top();
                    branchPosStack.pop();
                    currentDir = branchDirStack.top();
                    branchDirStack.pop();

                    if (ideal > 0.f)
                    {
                        bool &flip = flipMap[prevAtom];
                        float sign = flip ? 1.f : -1.f;
                        rot(sign * ideal * 0.5f);
                        flip = !flip;
                    }
                    else
                    {
                        rot(M_PI);
                    }
                }
                continue;
            }

            if (std::isdigit(c) || c == '%')
            {
                int32_t ringID = c - '0';
                if (c == '%')
                {
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
                    rot(M_PI * 1.5f);
                    ringOpen[ringID] = prevAtom;
                }
                bondType = BondType::SINGLE;
                continue;
            }

            if (isalpha(c))
            {
                if (atomIdx >= nAtoms.size())
                    break;

                float bondLength = 2.f * MULT_FACTOR;

                if (prevAtom != SIZE_MAX && nAtoms[atomIdx].ZIndex != 1 && c != 'H')
                {
                    const def_atom &a1 = nAtoms[prevAtom];
                    const def_atom &a2 = nAtoms[atomIdx];
                    bondLength = constants::getBondLength(a1.ZIndex, a2.ZIndex, bondType);

                    static bool flip = false;

                    auto ringIt = ringDir.find(prevAtom);
                    if (ringIt != ringDir.end())
                    {
                        rot(ringIt->second);
                    }
                    else if (ideal > 0.f)
                    {
                        float sign = flip ? 1.f : -1.f;
                        rot(sign * ideal * 0.3f);
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
            const def_subset &sub = nSubsets[s];
            const std::vector<size_t> hydrogens = sub.hydrogensIdx;
            const size_t centralAtom = sub.mainAtomIdx;

            sf::Vector3f axis(0,0,1);
            if (sub.bondedSubset != SIZE_MAX && sub.bondingSubset != SIZE_MAX)
            {
                const size_t c1 = nSubsets[sub.bondedSubset].mainAtomIdx;
                const size_t c2 = nSubsets[sub.bondingSubset].mainAtomIdx;
                const sf::Vector3f v1 = (positions[centralAtom] - positions[c1]).normalized();
                const sf::Vector3f v2 = (positions[centralAtom] - positions[c2]).normalized();
                axis = (v1 + v2).length() == 0.f ? axis : (v1 + v2).normalized();
            }
            else if (sub.bondingSubset != SIZE_MAX)
            {
                const size_t prev = nSubsets[sub.bondingSubset].mainAtomIdx;
                axis = (positions[centralAtom] - positions[prev]).normalized();
            }
            else if (sub.bondedSubset != SIZE_MAX)
            {
                const size_t next = nSubsets[sub.bondedSubset].mainAtomIdx;
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

            const size_t nH = sub.hydrogensIdx.size();
            for (size_t h = 0; h < nH; ++h)
            {
                const size_t hIdx = sub.hydrogensIdx[h];
                positions[hIdx] = positions[centralAtom];

                sf::Vector3f dir = baseDir; 
                const sf::Vector3f ringNormal(1, 0, 1);

                if (nH == 3)
                    dir = rotateDirection(dir, axis, -90.f * RADIAN);

                if (sub.bondedSubset != SIZE_MAX && sub.bondingSubset != SIZE_MAX && ringDir.count(sub.mainAtomIdx) > 0)
                    dir = rotateDirection(dir, ringNormal, 120.f * h * RADIAN);
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
                        sf::Vector3f rotAxis = tet[0].cross(baseDir);
                        if (rotAxis.length() > 1e-6f)
                        {
                            rotAxis = rotAxis.normalized();
                            const float ang = std::acos(std::clamp(tet[0].dot(baseDir), -1.0f, 1.0f));
                            for (int k = 0; k < 4; ++k)
                                tet[k] = rotateDirection(tet[k], rotAxis, ang);
                        }
                    }
                    dir = tet[h];
                }
                else 
                {
                    dir = rotateDirection(dir, ringNormal, 90 * h * RADIAN);
                }

                positions[hIdx] += dir * 0.9f;
            }
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
