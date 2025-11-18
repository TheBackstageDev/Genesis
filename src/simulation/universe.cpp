#include "universe.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <queue>

#include "constants.hpp"

namespace sim
{
    namespace fun
    {
        
        universe::universe(float universeSize)
            : boxSize(universeSize)
        {
        }

        void universe::drawBox(core::window_t &window)
        {
            const std::array<sf::Vector3f, 8> corners =
                {{{0.f, 0.f, 0.f},
                  {boxSize, 0.f, 0.f},
                  {boxSize, boxSize, 0.f},
                  {0.f, boxSize, 0.f},
                  {0.f, 0.f, boxSize},
                  {boxSize, 0.f, boxSize},
                  {boxSize, boxSize, boxSize},
                  {0.f, boxSize, boxSize}}};

            const std::array<std::pair<int, int>, 12> edges =
                {{
                    {0, 1}, {1, 2}, {2, 3}, {3, 0}, // bottom face
                    {4, 5},
                    {5, 6},
                    {6, 7},
                    {7, 4}, // top face
                    {0, 4},
                    {1, 5},
                    {2, 6},
                    {3, 7} // vertical pillars
                }};

            sf::VertexArray lines(sf::PrimitiveType::Lines, edges.size() * 2);
            sf::Color edgeColor(200, 200, 200, 180);

            size_t idx = 0;
            for (const auto &[i, j] : edges)
            {
                sf::Vector2f p1 = project(window, corners[i]);
                sf::Vector2f p2 = project(window, corners[j]);

                const sf::Vector3f &a = corners[i];
                const sf::Vector3f &b = corners[j];

                sf::Vector3f camToA = a - cam.eye();
                sf::Vector3f camToB = b - cam.eye();

                float dot = camToA.normalized().dot(camToB.normalized());
                const float THRESHOLD = 0.3f;
                if (dot < THRESHOLD)
                    continue;

                lines[idx++] = sf::Vertex(p1, edgeColor);
                lines[idx++] = sf::Vertex(p2, edgeColor);
            }

            window.getWindow().draw(lines);
        }

        void universe::draw(core::window_t &window, bool letter, bool lennardBall)
        {
            std::vector<size_t> drawOrder(atoms.size());
            std::iota(drawOrder.begin(), drawOrder.end(), 0);
            sf::Vector3f eye = cam.eye();

            std::sort(drawOrder.begin(), drawOrder.end(),
                      [&](size_t a, size_t b)
                      { return (pos(a) - eye).lengthSquared() > (pos(b) - eye).lengthSquared(); });

            drawBox(window);
            drawBonds(window);

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector2f p2 = project(window, pos(drawOrder[i]));
                if (p2.x < -1000)
                    continue;

                float temp = calculateAtomTemperature(drawOrder[i]);
                if (atoms[drawOrder[i]].ZIndex == 1)
                    drawHydrogenBond(window, drawOrder[i]);

                float camDistance = (pos(drawOrder[i]) - cam.eye()).length();

                float T = calculateAtomTemperature(i);
                atoms[drawOrder[i]].draw(T, p2, camDistance, window, letter, lennardBall);
            }
        }

        void universe::drawBonds(core::window_t &window)
        {
            sf::Vector2f dimensions = window.getWindow().getView().getSize();

            for (size_t b = 0; b < bonds.size(); ++b)
            {
                const bond &bond = bonds[b];
                const sf::Vector3f &pCentral = pos(bond.centralAtom); // pB
                const sf::Vector3f &pBonded = pos(bond.bondedAtom);   // pA

                sf::Vector2f s1 = cam.project(pCentral, dimensions.x, dimensions.y); // center
                sf::Vector2f s2 = cam.project(pBonded, dimensions.x, dimensions.y);  // end

                if (s1.x <= -9999 || s2.x <= -9999)
                    continue;
                if ((pCentral - pBonded).length() > 20.f)
                    continue;
                if ((pCentral - cam.eye()).length() > 200.f)
                    continue;

                sf::Vector2f dir = s2 - s1;
                float len = dir.length();
                if (len < 2.0f)
                    continue;
                dir /= len;

                sf::Vector2f perp{-dir.y, dir.x};

                uint8_t lines = static_cast<uint8_t>(bond.type);
                float shrink = 2.f;   // pixels
                float spacing = 10.f; // pixels between lines

                std::vector<sf::Vertex> verts;
                verts.reserve(lines * 2);

                for (int i = 0; i < lines; ++i)
                {
                    float offset = spacing * (i - (lines - 1) * 0.5f);
                    sf::Vector2f off = perp * offset / (pCentral - cam.eye()).length();

                    sf::Vector2f start = s1 + dir * shrink + off;
                    sf::Vector2f end = s2 - dir * shrink + off;

                    verts.emplace_back(start, sf::Color::White);
                    verts.emplace_back(end, sf::Color::White);
                }

                window.getWindow().draw(verts.data(), verts.size(), sf::PrimitiveType::Lines);
            }
        }

        void universe::drawHydrogenBond(core::window_t &window, size_t H)
        {
            const auto& atomH = atoms[H];
            if (neighbourList.size() != atoms.size())
                return;
            if (atomH.charge > -0.2f && atomH.charge < 0.2f)
                return;

            size_t D = SIZE_MAX;
            for (const auto &b : bonds)
            {
                if (b.bondedAtom == H)
                {
                    D = b.centralAtom;
                    break;
                }
            }
            if (D == SIZE_MAX)
                return;

            sf::Vector2f pH = project(window, pos(H));

            constexpr float MAX_HA_DISTANCE = 2.8f;
            constexpr float MIN_COS_ANGLE   = 0.5f;
            constexpr float STRONG_ENERGY   = -40.0f;
            constexpr float WEAK_ENERGY     = -5.0f;
            
            for (size_t A : neighbourList[H]) {
                if (A == D || A == H) continue;

                const auto& atomA = atoms[A];
                if (std::abs(atomA.charge) < 0.15f) continue;
                if (atomH.charge * atomA.charge > 0.0f) continue;

                sf::Vector3f rHD = minImageVec(pos(D) - pos(H));
                sf::Vector3f rHA = minImageVec(pos(A) - pos(H));
                float dHD = rHD.length();
                float dHA = rHA.length();

                if (dHA > MAX_HA_DISTANCE || dHA < 0.7f) continue;
                if (dHD < 0.5f) continue;

                float cos_angle = rHD.dot(rHA) / (dHD * dHA);
                if (cos_angle < MIN_COS_ANGLE) continue;

                float energy_kJ = COULOMB_K * atomH.charge * atomA.charge / dHA;
                if (energy_kJ > WEAK_ENERGY) continue;

                float strength = std::clamp((energy_kJ - WEAK_ENERGY) / (STRONG_ENERGY - WEAK_ENERGY), 0.0f, 1.0f);
                if (strength < 0.15f) continue;

                uint8_t alpha = static_cast<uint8_t>(255 * strength);

                sf::Vector2f pA = project(window, pos(A));
                sf::Vector2f dir = (pA - pH).normalized();
                float dist2D = (pA - pH).length();
                if (dist2D < 1e-3f) continue;

                constexpr int   NUM_DASHES = 8;
                constexpr float DASH_FRACTION = 0.55f;
                constexpr float GAP_FRACTION  = 0.45f;

                float dash_len = dist2D * DASH_FRACTION / NUM_DASHES;
                float gap_len  = dist2D * GAP_FRACTION  / NUM_DASHES;

                sf::Vertex line[2];
                line[0].color = line[1].color = sf::Color(255, 255, 255, alpha);

                for (int i = 0; i < NUM_DASHES; ++i) {
                    float start = (dash_len + gap_len) * i;
                    float end   = start + dash_len;

                    if (end > dist2D) end = dist2D;
                    if (start >= dist2D) break;

                    line[0].position = pH + dir * start;
                    line[1].position = pH + dir * end;
                    window.getWindow().draw(line, 2, sf::PrimitiveType::Lines);
                }
            }
        }

        void universe::drawDebug(core::window_t &window)
        {
        }

        void universe::log(size_t step)
        {
        }

        size_t universe::createAtom(sf::Vector3f p, sf::Vector3f v, uint8_t ZIndex, uint8_t numNeutrons, uint8_t numElectron, int32_t chirality)
        {
            atom newAtom{};
            newAtom.ZIndex = ZIndex;

            emplace_vel(v);
            emplace_pos(p);

            std::pair<float, float> constants = constants::getAtomConstants(ZIndex);

            newAtom.sigma = constants.first;
            newAtom.epsilon = constants.second;
            newAtom.radius = constants.first;
            newAtom.electrons = numElectron;
            newAtom.NCount = numNeutrons;
            newAtom.charge = ZIndex - numElectron;
            newAtom.mass = ZIndex * MASS_PROTON + numNeutrons * MASS_NEUTRON + numElectron * MASS_ELECTRON;
            newAtom.chirality = chirality;
            newAtom.bondCount = 0;

            atoms.emplace_back(std::move(newAtom));
            data.fx.resize(atoms.size());
            data.fy.resize(atoms.size());
            data.fz.resize(atoms.size());

            return atoms.size() - 1;
        }

        void universe::createBond(size_t idx1, size_t idx2, BondType type)
        {
            if (idx1 >= atoms.size() || idx2 >= atoms.size() || idx1 == idx2)
                return;

            if (std::find_if(bonds.begin(), bonds.end(), [&](bond &a)
                             { return a.bondedAtom == idx1 && a.centralAtom == idx2 || a.bondedAtom == idx2 && a.centralAtom == idx1; }) != bonds.end())
                return;

            bond nBond{};
            nBond.bondedAtom = idx1;
            nBond.centralAtom = idx2;
            nBond.type = type;
            nBond.equilibriumLength = constants::getBondLength(atoms[idx1].ZIndex, atoms[idx2].ZIndex, type);

            int8_t bondCount = static_cast<int8_t>(type);

            atoms[idx1].bondCount += bondCount;
            atoms[idx2].bondCount += bondCount;

            float EN1 = constants::electronegativity.at(atoms[idx1].ZIndex);
            float EN2 = constants::electronegativity.at(atoms[idx2].ZIndex);
            float deltaEN = std::abs(EN1 - EN2);

            if (deltaEN > 0.4f) // Significant electronegativity difference
            {
                float charge = deltaEN;

                if (EN1 > EN2)
                {
                    atoms[idx1].charge -= charge;
                    atoms[idx2].charge += charge;
                }
                else
                {
                    atoms[idx2].charge -= charge;
                    atoms[idx1].charge += charge;
                }
            }

            bonds.emplace_back(std::move(nBond));
        }

        size_t universe::createSubset(const def_subset &nSub, const size_t baseAtom, const size_t baseSubset)
        {
            subset nSubset{};
            nSubset.mainAtomIdx = nSub.mainAtomIdx + baseAtom;
            nSubset.bondedSubsetIdx = nSub.bondedSubset == SIZE_MAX ? SIZE_MAX : nSub.bondedSubset + baseSubset;
            nSubset.bondingSubsetIdx = nSub.bondingSubset == SIZE_MAX ? SIZE_MAX : nSub.bondingSubset + baseSubset;

            std::vector<size_t> connected{nSub.connectedIdx};
            std::vector<size_t> hydrogens{nSub.hydrogensIdx};
            std::vector<uint8_t> neighbourZs;
            neighbourZs.reserve(nSub.connectedIdx.size() + nSub.hydrogensIdx.size());

            for (size_t i = 0; i < connected.size(); ++i)
            {
                const size_t bondedAtom = connected[i];

                connected[i] += baseAtom;
                neighbourZs.emplace_back(atoms[bondedAtom].ZIndex);
            }

            for (size_t h = 0; h < hydrogens.size(); ++h)
            {
                const size_t bondedAtom = hydrogens[h];

                hydrogens[h] += baseAtom;
                neighbourZs.emplace_back(1);
            }

            nSubset.connectedIdx = std::move(connected);
            nSubset.hydrogenIdx = std::move(hydrogens);
            subsets.emplace_back(std::move(nSubset));

            return subsets.size() - 1; // Index
        }

        void universe::createMolecule(molecule_structure &structure, sf::Vector3f pos, sf::Vector3f vel)
        {
            size_t baseAtomIndex = atoms.size();

            molecule nMolecule{};

            for (size_t i = 0; i < structure.atoms.size(); ++i)
            {
                const def_atom &a = structure.atoms[i];
                if (pos.z != 0)
                    structure.positions[i].z += 0.01f * i;
                createAtom(structure.positions[i] + pos, vel, a.ZIndex, a.NIndex, a.ZIndex - a.charge, a.chirality);

                nMolecule.atomIdx.emplace_back(baseAtomIndex + i);
            }

            size_t baseBondIndex = bonds.size();
            for (size_t b = 0; b < structure.bonds.size(); ++b)
            {
                const def_bond &db = structure.bonds[b];
                size_t central = baseAtomIndex + db.centralAtomIdx;
                size_t bonded = baseAtomIndex + db.bondingAtomIdx;

                createBond(bonded, central, db.type);
                nMolecule.bondIdx.emplace_back(baseBondIndex + b);
            }

            size_t baseSubset = subsets.size();
            for (size_t s = 0; s < structure.subsets.size(); ++s)
            {
                createSubset(structure.subsets[s], baseAtomIndex, baseSubset);
                nMolecule.subsetIdx.emplace_back(baseSubset + s);
            }

            size_t baseAngle = angles.size();
            for (size_t a = 0; a < structure.angles.size(); ++a)
            {
                angles.emplace_back(structure.angles[a]);
                nMolecule.angles.emplace_back(baseAngle + a);
            }

            size_t baseDihedral = dihedral_angles.size();
            for (size_t a = 0; a < structure.dihedral_angles.size(); ++a)
            {
                dihedral_angles.emplace_back(structure.dihedral_angles[a]);
                nMolecule.dihedrals.emplace_back(baseDihedral + a);
            }

            nMolecule.name = moleculeName(nMolecule.subsetIdx);
            balanceMolecularCharges(subsets[baseSubset]);

            rebuildBondTopology();
        }

        void universe::balanceMolecularCharges(subset &mol)
        {
            auto it = std::find_if(subsets.begin(), subsets.end(), [&](const subset &s) 
            {
                return &s == &mol;
            });
            if (it == subsets.end()) return;

            size_t startIdx = std::distance(subsets.begin(), it);

            std::vector<size_t> atomsToBalance;
            std::unordered_set<size_t> visited;
            std::queue<size_t> q;
            q.push(startIdx);
            visited.insert(startIdx);

            while (!q.empty()) {
                size_t currentIdx = q.front(); q.pop();
                const auto &currentSubset = subsets[currentIdx];

                atomsToBalance.push_back(currentSubset.mainAtomIdx);
                atomsToBalance.insert(atomsToBalance.end(),
                    currentSubset.connectedIdx.begin(), currentSubset.connectedIdx.end());
                atomsToBalance.insert(atomsToBalance.end(),
                    currentSubset.hydrogenIdx.begin(), currentSubset.hydrogenIdx.end());

                size_t nextIdx = currentSubset.bondedSubsetIdx;
                if (nextIdx < subsets.size() && nextIdx != SIZE_MAX && visited.insert(nextIdx).second) 
                {
                    q.push(nextIdx);
                }
            }

            if (atomsToBalance.empty()) return;

            std::sort(atomsToBalance.begin(), atomsToBalance.end());
            auto last = std::unique(atomsToBalance.begin(), atomsToBalance.end());
            atomsToBalance.erase(last, atomsToBalance.end());

            std::vector<std::vector<size_t>> neighborLists(atoms.size());
            for (const auto &bond : bonds) {
                size_t a = bond.bondedAtom;
                size_t b = bond.centralAtom;
                bool aIn = std::binary_search(atomsToBalance.begin(), atomsToBalance.end(), a);
                bool bIn = std::binary_search(atomsToBalance.begin(), atomsToBalance.end(), b);
                if (aIn && bIn) {
                    neighborLists[a].push_back(b);
                    neighborLists[b].push_back(a);
                }
            }

            float totalCharge = 0.0f;
            for (size_t idx : atomsToBalance) 
            {
                totalCharge += atoms[idx].charge;
            }

            if (std::abs(totalCharge) < 1e-6f) return;

            std::vector<float> valenceDeficit(atomsToBalance.size(), 0.0f);
            for (size_t i = 0; i < atomsToBalance.size(); ++i) 
            {
                size_t idx = atomsToBalance[i];
                uint8_t usualValence = constants::getUsualBonds(atoms[idx].ZIndex);
                int deficit = usualValence - static_cast<int>(atoms[idx].bondCount);
                valenceDeficit[i] = static_cast<float>(deficit);
            }

            float chargeToDistribute = -totalCharge;
            float totalDeficit = 0.0f;
            for (float d : valenceDeficit) 
            {
                if (d > 0) totalDeficit += d;
            }

            if (totalDeficit > 0.5f) 
            {
                for (size_t i = 0; i < atomsToBalance.size(); ++i) 
                {
                    if (valenceDeficit[i] > 0) {
                        float share = chargeToDistribute * (valenceDeficit[i] / totalDeficit);
                        atoms[atomsToBalance[i]].charge += share;
                    }
                }
            } else 
            {
                float adjustment = chargeToDistribute / atomsToBalance.size();
                for (size_t idx : atomsToBalance) {
                    atoms[idx].charge += adjustment;
                }
            }
        }

        void universe::boundCheck(size_t i)
        {
            data.x[i] = std::fmod(data.x[i] + boxSize, boxSize);
            data.y[i] = std::fmod(data.y[i] + boxSize, boxSize);
            data.z[i] = std::fmod(data.z[i] + boxSize, boxSize);
        }

        float universe::ljPot(size_t i, float epsilon_i, float sigma_i)
        {
            float potential = 0.f;

            for (size_t j = 0; j < neighbourList.size(); ++j)
            {
                float dr = minImageVec(pos(i) - pos(neighbourList[i][j])).length();

                auto [sigma_j, epsilon_j] = constants::getAtomConstants(atoms[neighbourList[i][j]].ZIndex);
                float sigma = (sigma_i + sigma_j) / 2.0f;

                if (dr < sigma * CUTOFF && dr > EPSILON)
                {
                    float epsilon = sqrtf(epsilon_i * epsilon_j);
                    float r6 = powf((sigma / dr), 6);
                    float r12 = r6 * r6;
                    potential += 4.0f * epsilon * (r12 - r6);
                }
            }

            return potential;
        }

        sf::Vector3f universe::ljForce(size_t i, size_t j)
        {
            const atom &a1 = atoms[i];
            const atom &a2 = atoms[j];

            const float sigma_i = a1.sigma;
            const float sigma_j = a2.sigma;
            const float epsilon_i = a1.epsilon;
            const float epsilon_j = a2.epsilon;

            sf::Vector3f dr_vec = minImageVec(pos(i) - pos(j));
            float dr = dr_vec.length();

            const float sigma = (sigma_i + sigma_j) / 2.0f;

            if (dr < sigma * CUTOFF && dr > EPSILON)
            {
                float sigma6 = powf(sigma, 6);
                float sigma12 = sigma6 * sigma6;

                float epsilon = sqrtf(epsilon_i * epsilon_j);
                float r8 = sigma6 / powf(dr, 8);
                float r14 = 2.f * sigma12 / powf(dr, 14);
                float du_dr = 24.0f * epsilon * (r14 - r8);
                return (du_dr / dr) * dr_vec;
            }

            return sf::Vector3f{0.f, 0.f, 0.f};
        }

        sf::Vector3f universe::ljGrad(size_t i)
        {
            sf::Vector3f gradient({0.f, 0.f, 0.f});

            for (size_t j = 0; j < atoms.size(); ++j)
            {
                if (j == i)
                    continue;

                gradient += ljForce(i, j);
            }

            return gradient;
        }

        // Wolf method
        float universe::wolfForce(float r, float qi_qj)
        {
            constexpr float alpha = 0.25f;           // damping
            constexpr float rc = COULOMB_CUTOFF;
            if (r >= rc) return 0.f;

            float ir = 1.f/r;
            float ir_c = 1.f/rc;

            float erfc_ar = std::erfc(alpha*r);
            float exp_term = std::exp(-(alpha*alpha*r*r));

            float force = COULOMB_K * qi_qj * ir*ir * (
                erfc_ar + 2.f*alpha/std::sqrt(M_PI)*r*exp_term
            );

            return force;
        }

        sf::Vector3f universe::coulombForce(size_t i, size_t j, sf::Vector3f &dr_vec)
        {
            const atom &a1 = atoms[i];
            const atom &a2 = atoms[j];

            float dr = dr_vec.length();

            if (dr < EPSILON || dr > COULOMB_CUTOFF)
                return {0.f, 0.f, 0.f};

            float qq = a1.charge * a2.charge;
            if (qq == 0.f)
                return {0.f, 0.f, 0.f};

            float forceMag = COULOMB_K * a1.charge * a2.charge / dr;
            return -forceMag * dr_vec / dr;
        }

        void universe::calcBondForces()
        {
            for (size_t i = 0; i < bonds.size(); ++i)
            {
                bond &bond = bonds[i];

                size_t idx1 = bond.bondedAtom;
                size_t idx2 = bond.centralAtom;
                sf::Vector3f r_vec = minImageVec(pos(idx2) - pos(idx1));
                float dr = r_vec.length();
                if (dr <= EPSILON)
                    continue;

                float delta_r = dr - bond.equilibriumLength;

                /* if (delta_r > bond.equilibriumLength * 2.f)
                {
                    breakBond(idx1, idx2);
                    continue;
                } */

                float force_magnitude = BOND_K * delta_r;

                sf::Vector3f force_dir = r_vec / dr;
                sf::Vector3f force = force_magnitude * force_dir;

                add_force(idx1, force);
                add_force(idx2, -force);
            }
        }

        void universe::calcAngleForces()
        {
            //size_t count = 0;

            for (const angle& ang : angles) 
            {
                size_t i = ang.A;   // left atom
                size_t j = ang.B;   // central atom
                size_t k = ang.C;   // right atom

                sf::Vector3f r_ji = minImageVec(pos(i) - pos(j));
                sf::Vector3f r_jk = minImageVec(pos(k) - pos(j)); 

                float len_ji = r_ji.length();
                float len_jk = r_jk.length();

                if (len_ji < EPSILON || len_jk < EPSILON) continue;

                sf::Vector3f u_ji = r_ji / len_ji;
                sf::Vector3f u_jk = r_jk / len_jk;

                float cos_theta = std::clamp(u_ji.dot(u_jk), -1.0f, 1.0f);
                float theta = std::acos(cos_theta);

                float sin_theta = std::sqrt(1.0f - cos_theta * cos_theta);
                if (sin_theta < 1e-6f) sin_theta = 1e-6f;

                float delta_theta = theta - ang.rad;
                float k_angle = 2.0f * ang.K;

                sf::Vector3f dtheta_dri = (cos_theta * u_ji - u_jk) / (len_ji * sin_theta);
                sf::Vector3f dtheta_drk = (cos_theta * u_jk - u_ji) / (len_jk * sin_theta);

                sf::Vector3f F_i = -k_angle * delta_theta * dtheta_dri;
                sf::Vector3f F_k = -k_angle * delta_theta * dtheta_drk;
                sf::Vector3f F_j = -F_i - F_k;  // Newton’s 3rd law

                add_force(i, F_i);
                add_force(j, F_j);
                add_force(k, F_k);

                //++count;
            }

            //std::cout << "angle forces applied: " << count << std::endl;
        }

        void universe::calcLjForces()
        {
            // size_t count = 0;

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                std::vector<size_t> &neighbours = neighbourList[i];

                for (size_t j = 0; j < neighbours.size(); ++j)
                {
                    size_t index_j = neighbours[j];

                    const float sigma = (atoms[i].sigma + atoms[index_j].sigma) / 2.0f;

                    //++count;

                    sf::Vector3f force = ljForce(i, index_j);

                    add_force(i, force);
                    add_force(index_j, -force); // Every action creates an Equal and opposite reaction - Sir Isaac Newton
                }
            }

            // std::cout << "count: " << count << std::endl;
        }

        float universe::calculateDihedral(const sf::Vector3f &pa, const sf::Vector3f &pb, const sf::Vector3f &pc, const sf::Vector3f &pd)
        {
            sf::Vector3f v1 = pb - pa;
            sf::Vector3f v2 = pc - pb;
            sf::Vector3f v3 = pd - pc;

            sf::Vector3f n1 = v1.cross(v2).normalized();
            sf::Vector3f n2 = v2.cross(v3).normalized();

            float cosPhi = std::clamp(n1.dot(n2), -1.f, 1.f);
            float phi    = std::acos(cosPhi);

            // sign via right-hand rule
            if (v2.dot(n1.cross(n2)) < 0.f) phi = -phi;
            if (phi < 0.f) phi += 2.f * M_PI;
            return phi;
        }

        void universe::calcDihedralForces()
        {
            for (size_t d = 0; d < dihedral_angles.size(); ++d)
            {
                const dihedral_angle& d_angle = dihedral_angles[d];

                const sf::Vector3f& pa = pos(d_angle.A);
                const sf::Vector3f& pb = pos(d_angle.B); 
                const sf::Vector3f& pc = pos(d_angle.C); 
                const sf::Vector3f& pd = pos(d_angle.D); 

                float phi = calculateDihedral(pa, pb, pc, pd);

                float target = d_angle.rad;
                if (target == 0.f && d_angle.periodicity == 1) 
                {
                    int32_t chi = atoms[d_angle.B].chirality ? atoms[d_angle.B].chirality : atoms[d_angle.C].chirality;
                    if (chi == 1) target = (phi < static_cast<float>(M_PI)) ? 1.047f : 5.236f;   // ~60° / ~300°
                    else if (chi == 2) target = (phi < static_cast<float>(M_PI)) ? 5.236f : 1.047f;
                    else continue;
                }

                float diff = phi - target;
                while (diff > M_PI) diff -= 2.f * M_PI;
                while (diff < -M_PI) diff += 2.f * M_PI;

                float forceMag = -d_angle.K * float(d_angle.periodicity) *
                                std::sin(float(d_angle.periodicity) * phi);

                sf::Vector3f axis = (pc - pb).normalized();

                sf::Vector3f rA = pa - pb;
                sf::Vector3f torqueA = axis.cross(rA).normalized() * forceMag;

                sf::Vector3f rD = pd - pc;
                sf::Vector3f torqueD = axis.cross(rD).normalized() * (-forceMag);

                add_force(d_angle.A, torqueA * 0.5f);
                add_force(d_angle.B, torqueA * 0.5f - torqueD * 0.5f);
                add_force(d_angle.C, torqueD * 0.5f);
                add_force(d_angle.D, -torqueD * 0.5f);
            }
        }

        void universe::calcElectrostaticForces()
        {
            float cutoff2 = COULOMB_CUTOFF * COULOMB_CUTOFF;
            for (size_t i = 0; i < atoms.size(); ++i)
            {
                if (atoms[i].charge == 0.f) continue;

                for (size_t j : neighbourList[i])
                {
                    if (atoms[j].charge == 0.f) continue;

                    sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                    float dr_m = dr.length();
                    if (dr_m > cutoff2 || dr_m < EPSILON)
                        continue;

                    sf::Vector3f force = coulombForce(i, j, dr);
                    add_force(i, force);
                    add_force(j, -force);
                }
            }
        }

        void universe::update(float targetTemperature, bool reactions)
        {
            data.old_fx.resize(atoms.size());
            data.old_fy.resize(atoms.size());
            data.old_fz.resize(atoms.size());
            std::swap(data.fx, data.old_fx);
            std::swap(data.fy, data.old_fy);
            std::swap(data.fz, data.old_fz);

            data.fx.resize(atoms.size());
            data.fy.resize(atoms.size());
            data.fz.resize(atoms.size());

            if (neighbourList.empty() || neighbourListNeedRebuild())
                buildNeighborList();

            calcLjForces();
            calcBondForces();
            calcAngleForces();
            calcDihedralForces();
            calcElectrostaticForces();

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector3f a_old = old_force(i) / atoms[i].mass;
                sf::Vector3f a_new = force(i) / atoms[i].mass;
                add_pos(i, vel(i) * DT + 0.5f * (a_old + a_new) * DT * DT);
                boundCheck(i);
            }

            std::fill(data.fx.begin(), data.fx.end(), 0.f);
            std::fill(data.fy.begin(), data.fy.end(), 0.f);
            std::fill(data.fz.begin(), data.fz.end(), 0.f);

            calcLjForces();
            calcBondForces();
            calcAngleForces();
            calcDihedralForces();
            calcElectrostaticForces();

            for (size_t i = 0; i < atoms.size(); ++i)
                add_vel(i, (old_force(i) + force(i)) * (0.5f * DT / atoms[i].mass));

            if (reactions)
                handleReactions();

            setTemperature(targetTemperature);

            ++timeStep;
        }

        void universe::setTemperature(float kelvin)
        {
            if (timeStep % THERMOSTAT_INTERVAL == 0)
            {
                float avg_KE = calculateKineticEnergy() / atoms.size();
                temp = (2.0f / 3.0f) * (avg_KE * KB);
                float lambda = sqrtf(kelvin / temp);

                if (timeStep > 100)
                {
                    constexpr float tau = 0.6f;
                    lambda = 1.0f + (lambda - 1.0f) / tau;
                }

                for (size_t i = 0; i < data.vx.size(); ++i)
                {
                    data.vx[i] *= lambda;
                    data.vy[i] *= lambda;
                    data.vz[i] *= lambda;
                }
            }
        }

        // Energy Calculation
        float universe::calculateKineticEnergy()
        {
            float kinetic_energy = 0.0f;
            for (size_t i = 0; i < atoms.size(); ++i)
                kinetic_energy += 0.5f * atoms[i].mass * vel(i).lengthSquared();

            return kinetic_energy;
        }

        void universe::buildNeighborList()
        {
            const float rskin = COULOMB_CUTOFF + VERLET_SKIN;
            
            if (neighbourList.size() != atoms.size())
                neighbourList.assign(atoms.size(), {});

            for (auto& list : neighbourList)
                list.clear();

            data.x0 = data.x;
            data.y0 = data.y;
            data.z0 = data.z;

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                for (size_t j = i + 1; j < atoms.size(); ++j)
                {
                    if (areBonded(i, j))
                        continue;

                    sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                    if (dr.length() < rskin)
                    {
                        neighbourList[i].push_back(j);
                        neighbourList[j].push_back(i);
                    }
                }
            }
        }

        bool universe::neighbourListNeedRebuild()
        {
            if (neighbourList.size() != atoms.size())
                return true;

            const float threshold2 = REBUILD_THRESHOLD * REBUILD_THRESHOLD;

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector3f d = pos(i) - neigh_pos(i);
                if (d.lengthSquared() > threshold2)
                    return true;
            }
            return false;
        }

        // Reactions

        float universe::calculateAtomTemperature(size_t i)
        {
            float ke = 0.5f * atoms[i].mass * vel(i).lengthSquared();
            return (2.0f / 3.0f) * ke * KB;
        }

        float universe::calculateBondEnergy(size_t i, size_t j)
        {
            float r_ij = (pos(i) - pos(j)).length();
            const bond& b = bonds[getBond(i, j)];

            float bondStretchingEnergy;
            float angleEnergy;
            float dihedralEnergy;

            float coulombPotential = COULOMB_K * (atoms[i].charge * atoms[j].charge) / r_ij; // Nitroglycerin is a great example

            return coulombPotential + bondStretchingEnergy + angleEnergy + dihedralEnergy;
        }

        float universe::calculateNonBondedEnergy(size_t i, size_t j)
        {
            float r_ij = (pos(i) - pos(j)).length();
            float coulombPotential = COULOMB_K * (atoms[i].charge * atoms[j].charge) / r_ij;
            float lennardjonesPotential = ljPot(i, atoms[i].epsilon, atoms[i].sigma);

            return lennardjonesPotential + coulombPotential;
        }

        void universe::breakBond(size_t atom1, size_t atom2)
        {
            auto it = std::find_if(bonds.begin(), bonds.end(), [&](const bond &bond)
                                   { return bond.centralAtom == atom1 && bond.bondedAtom == atom2 || bond.centralAtom == atom2 && bond.bondedAtom == atom1; });

            BondType type = it->type;

            if (it != bonds.end())
            {
                bonds.erase(it);
                atoms[atom1].bondCount -= static_cast<uint8_t>(type);
                atoms[atom2].bondCount -= static_cast<uint8_t>(type);
            }

            if (atoms[atom1].ZIndex == 1 || atoms[atom2].ZIndex == 1)
            {
                size_t H = (atoms[atom1].ZIndex == 1 ? atom1 : atom2);
                size_t C = (atoms[atom1].ZIndex == 1 ? atom2 : atom1);

                auto sub_it = std::find_if(subsets.begin(), subsets.end(),
                                           [&](const subset &s)
                                           {
                                               return s.mainAtomIdx == C &&
                                                      std::find(s.hydrogenIdx.begin(),
                                                                s.hydrogenIdx.end(), H) != s.hydrogenIdx.end();
                                           });

                if (sub_it != subsets.end())
                {
                    // balanceMolecularCharges(*sub_it);
                    sub_it->hydrogenIdx.erase(
                        std::remove(sub_it->hydrogenIdx.begin(),
                                    sub_it->hydrogenIdx.end(), H),
                        sub_it->hydrogenIdx.end());
                }
                return;
            }

            auto subset_it1 = std::find_if(subsets.begin(), subsets.end(),
                                           [&](const subset &s)
                                           { return s.mainAtomIdx == atom1; });

            auto subset_it2 = std::find_if(subsets.begin(), subsets.end(),
                                           [&](const subset &s)
                                           { return s.mainAtomIdx == atom2; });

            if (subset_it1 != subsets.end())
            {
                auto &s1 = *subset_it1;
                s1.connectedIdx.erase(
                    std::remove(s1.connectedIdx.begin(), s1.connectedIdx.end(), atom2),
                    s1.connectedIdx.end());

                if (s1.bondedSubsetIdx < subsets.size() && subsets[s1.bondedSubsetIdx].mainAtomIdx == atom2)
                    s1.bondedSubsetIdx = SIZE_MAX;
                if (s1.bondingSubsetIdx < subsets.size() && subsets[s1.bondingSubsetIdx].mainAtomIdx == atom2)
                    s1.bondingSubsetIdx = SIZE_MAX;

                balanceMolecularCharges(s1);
            }

            if (subset_it2 != subsets.end())
            {
                auto &s2 = *subset_it2;
                s2.connectedIdx.erase(
                    std::remove(s2.connectedIdx.begin(), s2.connectedIdx.end(), atom1),
                    s2.connectedIdx.end());

                if (s2.bondedSubsetIdx < subsets.size() && subsets[s2.bondedSubsetIdx].mainAtomIdx == atom1)
                    s2.bondedSubsetIdx = SIZE_MAX;
                if (s2.bondingSubsetIdx < subsets.size() && subsets[s2.bondingSubsetIdx].mainAtomIdx == atom1)
                    s2.bondingSubsetIdx = SIZE_MAX;

                balanceMolecularCharges(s2);
            }
        }

        BondType universe::chooseBestBondType(size_t i, size_t j)
        {
            if (areBonded(i, j))
                return BondType::NONE;

            uint8_t Z1 = atoms[i].ZIndex, Z2 = atoms[j].ZIndex;
            float r = minImageVec(pos(i) - pos(j)).length();

            uint8_t val1 = atoms[i].bondCount;
            uint8_t val2 = atoms[j].bondCount;
            uint8_t max1 = constants::getUsualBonds(Z1);
            uint8_t max2 = constants::getUsualBonds(Z2);

            for (int8_t order = 3; order >= 1; --order)
            {
                BondType type = static_cast<BondType>(order - 1);
                float r0 = constants::getBondLength(Z1, Z2, type);

                if (r0 <= 0.0f)
                    continue;
                if (r > r0 * 1.3f || r < r0 * 0.5f)
                    continue;

                if (val1 + order > max1 || val2 + order > max2)
                    continue;

                return type;
            }

            return BondType::NONE;
        }

        void universe::handleReactions()
        {
            std::vector<std::pair<size_t, size_t>> candidates;

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                if (!isReactive(i))
                    continue;

                for (size_t j : neighbourList[i])
                {
                    if (!isReactive(j))
                        continue;
                    if (areBonded(i, j))
                        continue;

                    BondType bestType = chooseBestBondType(i, j);
                    float r = minImageVec(pos(i) - pos(j)).length();
                    float r0 = constants::getBondLength(atoms[i].ZIndex, atoms[j].ZIndex, bestType);
                    if (r0 > 0 && r < r0 * 1.3f)
                    {
                        createBond(i, j, bestType);
                        auto subset_it1 = std::find_if(subsets.begin(), subsets.end(),
                                                       [&](const subset &s)
                                                       { return s.mainAtomIdx == i; });

                        auto subset_it2 = std::find_if(subsets.begin(), subsets.end(),
                                                       [&](const subset &s)
                                                       { return s.mainAtomIdx == j; });

                        if (subset_it1 != subsets.end() && subset_it2 != subsets.end())
                        {
                            subset_it1->bondingSubsetIdx = std::distance(subsets.begin(), subset_it2);
                            subset_it2->bondingSubsetIdx = std::distance(subsets.begin(), subset_it1);
                        }
                    }
                }
            }
        }

        // Camera
        sf::Vector2f universe::project(core::window_t &window, const sf::Vector3f &p) const
        {
            auto size = window.getWindow().getView().getSize();
            return cam.project(p, size.x, size.y);
        }

        void universe::handleCamera(bool leftDown, bool rightDown, const sf::Vector2i &mousePos, float wheelDelta, const std::vector<sf::Keyboard::Key> &keys)
        {
            if (leftDown)
            {
                sf::Vector2i delta = mousePos - lastMouse;
                cam.azimuth -= delta.x * 0.3f;
                cam.elevation = std::clamp(cam.elevation - delta.y * 0.3f,
                                           -89.f, 89.f);
            }

            sf::Vector3f forward = cam.distance == 0.f ? sf::Vector3f(1.f, 0.f, 0.f) : (cam.target - cam.eye()).normalized();
            sf::Vector3f right = forward.cross({0, 0, 1}).normalized();
            sf::Vector3f up = right.cross(forward);

            if (rightDown)
            {
                sf::Vector2i delta = mousePos - lastMouse;
                cam.target += sf::Vector3f{right.x * delta.x, right.y, right.z} * cam.distance * 0.002f +
                              sf::Vector3f{up.x, up.y * delta.y, up.z} * cam.distance * 0.002f;
            }

            if (wheelDelta != 0.0f)
            {
                float base_speed = 1.5f;
                float speed = base_speed * (cam.distance / 50.0f);

                if (wheelDelta > 0)
                {
                    cam.distance *= std::pow(0.9f, speed * wheelDelta);
                }
                else
                {
                    cam.distance *= std::pow(1.1f, speed * -wheelDelta);
                }
                cam.distance = std::clamp(cam.distance, 10.0f, 500.0f);
            }

            for (auto k : keys)
            {
                if (k == sf::Keyboard::Key::W)
                    cam.target += forward * 0.5f;
                if (k == sf::Keyboard::Key::S)
                    cam.target -= forward * 0.5f;
                if (k == sf::Keyboard::Key::A)
                    cam.target -= right * 0.5f;
                if (k == sf::Keyboard::Key::D)
                    cam.target += right * 0.5f;
                if (k == sf::Keyboard::Key::Q)
                    cam.target -= up * 0.5f; // Down
                if (k == sf::Keyboard::Key::E)
                    cam.target += up * 0.5f; // Up
            }

            lastMouse = mousePos;
        }

        // TO DO
        std::string universe::moleculeName(const std::vector<size_t> &subsetIdx)
        {
            std::map<uint8_t, size_t> ZIndices;

            for (int32_t i = 0; i < subsetIdx.size(); ++i)
            {
                size_t centralAtom = subsets[subsetIdx[i]].mainAtomIdx;
                ++ZIndices[atoms[centralAtom].ZIndex];
            }

            std::unordered_map<size_t, std::string> prefixes{
                {1, "mono"},
                {2, "di"},
                {3, "tri"},
                {4, "tetra"},
                {5, "penta"},
                {6, "hexa"},
                {7, "hepta"},
                {8, "octa"},
                {9, "nona"},
                {10, "deca"}};

            std::string name("");

            bool organic = ZIndices[6] > 0;

            if (!organic)
            {
                for (int32_t i = 1; i <= COUNT_ATOMS; ++i)
                {
                    if (ZIndices[i] <= 0)
                        continue;
                    name += prefixes[ZIndices[i]] + constants::getAtomName(i) + " ";
                }
            }
            else // organic nomeclature
            {
            }

            return name;
        }
    } // namespace fun
} // namespace sim
