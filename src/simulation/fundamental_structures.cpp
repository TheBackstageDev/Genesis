#include "fundamental_structures.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <unordered_set>

#include "constants.hpp"

namespace sim
{
    constexpr float pixel_per_A = 1.0f; // 50 pixels/Å

    sf::Color lerpColor(sf::Color a, sf::Color b, float t)
    {
        return sf::Color(
            static_cast<uint8_t>(a.r + (b.r - a.r) * sin(t)),
            static_cast<uint8_t>(a.g + (b.g - a.g) * sin(t)),
            static_cast<uint8_t>(a.b + (b.b - a.b) * sin(t)),
            255
        );
    }

    namespace fun
    {
        void atom::draw(float temperature, sf::Vector3f& pos, core::window_t &window, float UniverseSize, bool letter)
        {
            sf::Vector2f posVec{pos.x, pos.y};

            float z = pos.z; 
            float depthFactor = 1.0f - std::clamp((z + UniverseSize/2) / UniverseSize, 0.f, 1.f); 
            float scale = 0.7f + 0.6f * depthFactor; 

            sf::Color atomColor = lerpColor(sf::Color::Blue, sf::Color::Red, temperature / 100.f * RADIAN);

            const float cloudRadius = radius * 0.5f * scale;
            sf::CircleShape cloud(cloudRadius);
            cloud.setOrigin({cloudRadius, cloudRadius});
            cloud.setPosition(posVec);
            cloud.setFillColor(sf::Color(atomColor.r, atomColor.g, atomColor.b, 120));
            window.draw(cloud);

            if (letter)
            {
                auto& font = window.getFont();
                std::string name = constants::getAtomLetter(ZIndex);
                sf::Text name_text;
                name_text.setFont(font);
                name_text.setString(name);
                name_text.setCharacterSize(static_cast<unsigned>(20.f * scale));
                name_text.setScale({0.02f * radius * scale, 0.02f * radius * scale});
                name_text.setFillColor({
                    255,
                    255,
                    255,
                    120
                });

                sf::FloatRect bounds = name_text.getLocalBounds();
                name_text.setOrigin(bounds.getCenter());
                name_text.setPosition(posVec);
                window.draw(name_text);

                // Charge
                if (charge != 0.f && std::abs(charge) > 0.5f)
                {
                    sf::Text ion_text;
                    ion_text.setFont(font);
                    ion_text.setString(charge > 0.f ? "+" : "-");
                    ion_text.setCharacterSize(static_cast<unsigned>(16.f * scale));
                    ion_text.setScale({0.02f * radius * scale, 0.02f * radius * scale});
                    ion_text.setFillColor(sf::Color::White);
                    sf::FloatRect bounds = ion_text.getLocalBounds();
                    ion_text.setOrigin(bounds.getCenter());
                    ion_text.setPosition(posVec + sf::Vector2f(0.f, -radius * 0.3f * scale));
                    window.draw(ion_text);
                }
            }
        }

        void atom::drawProjected(sf::Vector2f screenPos,
                                core::window_t& window,
                                bool letter,
                                sf::Color baseColor)
        {       
            sf::CircleShape cloud(1.5f);
            cloud.setOrigin({1.5f, 1.5f});
            cloud.setPosition(screenPos);
            cloud.setFillColor(sf::Color(baseColor.r, baseColor.g, baseColor.b, 150));
            window.draw(cloud);

            if (letter)
            {
                auto& font = window.getFont();
                std::string name = constants::getAtomLetter(ZIndex);
                sf::Text name_text;
                name_text.setFont(font);
                name_text.setString(name);
                name_text.setCharacterSize(static_cast<unsigned>(20.f));
                name_text.setScale({0.02f * radius, 0.02f * radius});
                name_text.setFillColor({
                    255,
                    255,
                    255,
                    120
                });

                sf::FloatRect b = name_text.getLocalBounds();
                name_text.setOrigin(b.getCenter());
                name_text.setPosition(screenPos);
                window.draw(name_text);

                if (std::abs(charge) > 0.5f)
                {
                    sf::Text ion_text;
                    ion_text.setFont(font);
                    ion_text.setString(charge > 0.f ? "+" : "-");
                    ion_text.setCharacterSize(16);
                    ion_text.setScale({0.02f * radius, 0.02f * radius});
                    ion_text.setFillColor(sf::Color::White);

                    sf::FloatRect bounds = ion_text.getLocalBounds();
                    ion_text.setOrigin(bounds.getCenter());
                    ion_text.setPosition(screenPos + sf::Vector2f(8.f, -8.f));
                    window.draw(ion_text);
                }
            }
        }
        
        universe::universe(float universeSize)
            : boxSize(universeSize)
        {
        }

        void universe::draw(core::window_t& window, bool letter)
        {
            std::vector<size_t> drawOrder(atoms.size());
            std::iota(drawOrder.begin(), drawOrder.end(), 0);
            sf::Vector3f eye = cam.eye();
            
            std::sort(drawOrder.begin(), drawOrder.end(),
            [&](size_t a,size_t b){ return (positions[a]-eye).lengthSquared()
                                     < (positions[b]-eye).lengthSquared(); });

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector2f p2 = project(window, positions[i]);
                if (p2.x < -1000) continue;

                float temp = calculateAtomTemperature(drawOrder[i]);
                if (atoms[drawOrder[i]].ZIndex == 1)
                    drawHydrogenBond(window, drawOrder[i]);

                float T  = calculateAtomTemperature(i);
                atoms[drawOrder[i]].draw(temp, positions[drawOrder[i]], window, boxSize, letter);
            }

            for (size_t b = 0; b < bonds.size(); ++b)
            {
                const bond& bond = bonds[b];
                const sf::Vector3f& pA = positions[bond.bondedAtom];   
                const sf::Vector3f& pB = positions[bond.centralAtom]; 

                if ((pA - pB).length() > 5.f * MULT_FACTOR) continue;

                sf::Vector3f dir  = pA - pB == sf::Vector3f(0.f, 0.f, 0.f) ? sf::Vector3f(0.f, 0.f, 0.f) : (pA - pB).normalized();
                sf::Vector2f perp{-dir.y, dir.x};       

                int8_t lines = static_cast<int8_t>(bond.type);
                const float shrink = 0.5f * MULT_FACTOR;   

                std::vector<sf::Vertex> vertices;
                vertices.reserve(lines * 2);

                const float offsetStep = -0.2f * MULT_FACTOR;               
                for (int8_t i = 0; i < lines; ++i)
                {
                    float signedOffset = offsetStep * (i - (lines - 1) * 0.5f);
                    sf::Vector2f offset = perp * signedOffset;

                    sf::Vector2f start = sf::Vector2f{pB.x, pB.y} + sf::Vector2f{dir.x, dir.y} * shrink + offset;
                    sf::Vector2f end   = sf::Vector2f{pA.x, pA.y} - sf::Vector2f{dir.x, dir.y} * shrink + offset;

                    vertices.emplace_back(sf::Vector2f{start.x, start.y});
                    vertices.emplace_back(sf::Vector2f{end.x, end.y});
                }

                window.getWindow().draw(vertices.data(), vertices.size(), sf::PrimitiveType::Lines);
            }

            // Universe Edges

            std::array<sf::Vertex, 2> down = {
                sf::Vertex({boxSize, boxSize}, sf::Color::White),
                sf::Vertex({0.f, boxSize}, sf::Color::White)
            };

            std::array<sf::Vertex, 2> up = {
                sf::Vertex({boxSize, 0.f}, sf::Color::White),
                sf::Vertex({0.f, 0.f}, sf::Color::White)
            };

            std::array<sf::Vertex, 2> left = {
                sf::Vertex({0.f, boxSize}, sf::Color::White),
                sf::Vertex({0.f, 0.f}, sf::Color::White)
            };

            std::array<sf::Vertex, 2> right = {
                sf::Vertex({boxSize, 0.f}, sf::Color::White),
                sf::Vertex({boxSize, boxSize}, sf::Color::White)
            };

            window.getWindow().draw(up.data(), up.size(), sf::PrimitiveType::Lines);
            window.getWindow().draw(down.data(), down.size(), sf::PrimitiveType::Lines);
            window.getWindow().draw(left.data(), left.size(), sf::PrimitiveType::Lines);
            window.getWindow().draw(right.data(), right.size(), sf::PrimitiveType::Lines);
        }

        void universe::drawHydrogenBond(core::window_t& window, size_t H)
        {
            if (neighbourList.size() != atoms.size()) return;
            if (atoms[H].charge > -0.2f && atoms[H].charge < 0.2f) return;

            size_t D = SIZE_MAX;
            for (const auto& b : bonds)
            {
                if (b.bondedAtom == H) { D = b.centralAtom; break; }
                if (b.centralAtom == H) { D = b.bondedAtom; break; }
            }
            if (D == SIZE_MAX) return;

            for (size_t j : neighbourList[H])
            {
                if (j == D) continue;
                uint8_t ZA = atoms[j].ZIndex;
                
                if (ZA != 7 && ZA != 8 && ZA != 9) continue;  // N, O, F
                if (atoms[j].charge < 0.2f && atoms[j].charge > -0.2f) continue;

                sf::Vector3f r_HD = minImageVec(positions[D] - positions[H]);
                sf::Vector3f r_HA = minImageVec(positions[j] - positions[H]);
                float d_HD = r_HD.length();
                float d_HA = r_HA.length();

                if (d_HD < 0.3f * MULT_FACTOR || d_HD > 2.f * MULT_FACTOR) continue;

                float cos_angle = r_HD.dot(r_HA) / (d_HD * d_HA);
                if (cos_angle < 0.5f) continue;  // >120°
                
                sf::Vector2f pH = {positions[H].x, positions[H].y};
                sf::Vector2f pA = {positions[j].x, positions[j].y};
                
                sf::Vector2f d = pA - pH;

                if (d.length() <= 0.f) continue;

                sf::Vector2f dir =  d.normalized();
                float len = (pA - pH).length();

                if (len > 5.f * MULT_FACTOR) continue;
                float energy_kJ = COULOMB_K * atoms[H].charge * atoms[j].charge / d_HA;

                float maxEnergy = -20.0f;
                float minEnergy =  -2.0f;
                float norm = std::clamp((energy_kJ - minEnergy) / (maxEnergy - minEnergy), 0.f, 1.f);
                uint8_t alpha = static_cast<uint8_t>(norm * 255.f);

                if (norm < 0.2f) continue;

                for (float t = 0.3f; t < 0.8f; t += 0.1f)
                {
                    sf::Vector2f start = pH + dir * len * t;
                    sf::Vector2f end   = pH + dir * len * (t + 0.06f);
                    if (end.x < 0 || end.y < 0 || start.x > boxSize || start.y > boxSize) continue;

                    sf::Vertex line[] = {
                        sf::Vertex(start, sf::Color(255, 255, 255, alpha)),
                        sf::Vertex(end,   sf::Color(255, 255, 255, alpha))
                    };
                    window.getWindow().draw(line, 2, sf::PrimitiveType::Lines);
                }
            }
        }

        void universe::drawDebug(core::window_t& window)
        {

        }

        void universe::log(size_t step)
        {

        }

        size_t universe::createAtom(sf::Vector3f p, sf::Vector3f v, uint8_t ZIndex, uint8_t numNeutrons, uint8_t numElectron)
        {
            atom newAtom{};
            newAtom.ZIndex = ZIndex;

            velocities.emplace_back(v);
            positions.emplace_back(p);

            std::pair<float, float> constants = constants::getAtomConstants(ZIndex);
            
            newAtom.sigma = constants.first;
            newAtom.epsilon = constants.second;
            newAtom.radius = constants.first;
            newAtom.electrons = numElectron;
            newAtom.NCount = numNeutrons;
            newAtom.charge = ZIndex - numElectron;
            newAtom.mass = ZIndex * MASS_PROTON + numNeutrons * MASS_NEUTRON + numElectron * MASS_ELECTRON;
            newAtom.bondCount = 0;

            atoms.emplace_back(std::move(newAtom));
            forces.resize(atoms.size());
 
            return atoms.size() - 1;
        }

        void universe::createBond(size_t idx1, size_t idx2, BondType type)
        {
            if (idx1 >= atoms.size() || idx2 >= atoms.size() || idx1 == idx2)
                return;

            if (std::find_if(bonds.begin(), bonds.end(), [&](bond& a){
                return a.bondedAtom == idx1 && a.centralAtom == idx2 || a.bondedAtom == idx2 && a.centralAtom == idx1;
            }) != bonds.end()) return; 

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

        size_t universe::createSubset(const def_subset& nSub, const size_t baseAtom, const size_t baseSubset)
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
            nSubset.idealAngle = nSub.idealAngle;
            subsets.emplace_back(std::move(nSubset));

            return subsets.size() - 1; // Index
        }

        // TO DO: Support ionic-molecules 
        void universe::createMolecule(molecule_structure& structure, sf::Vector3f pos)
        {
            size_t baseAtomIndex = atoms.size();
            
            for (size_t i = 0; i < structure.atoms.size(); ++i)
            {
                const def_atom& a = structure.atoms[i];
                if (pos.z != 0)
                    structure.positions[i] += sf::Vector3f(0.f, 0.f, 0.01f * i);
                
                createAtom(structure.positions[i] + pos, {0.f, 0.f, 0.f}, a.ZIndex, a.NIndex, a.ZIndex - a.charge);
            }

            for (size_t b = 0; b < structure.bonds.size(); ++b)
            {
                const def_bond& db = structure.bonds[b];
                size_t central = baseAtomIndex + db.centralAtomIdx; 
                size_t bonded  = baseAtomIndex + db.bondingAtomIdx;

                createBond(bonded, central, db.type);
            }

            size_t baseSubset = subsets.size();
            for (const def_subset& s : structure.subsets)
            {
                createSubset(s, baseAtomIndex, baseSubset);
            }

            balanceMolecularCharges(subsets[baseSubset]);
        }

        void universe::balanceMolecularCharges(subset& mol)
        {
            auto it = std::find_if(subsets.begin(), subsets.end(), [&](const subset& s) { return &s == &mol; });
            if (it == subsets.end()) return; 

            size_t startIdx = std::distance(subsets.begin(), it);
            std::vector<size_t> atomsToBalance;
            std::unordered_set<size_t> visited;
            size_t currentIdx = startIdx;

            while (currentIdx < subsets.size()) {
                if (!visited.insert(currentIdx).second) 
                    break;

                const auto& currentSubset = subsets[currentIdx];
                atomsToBalance.push_back(currentSubset.mainAtomIdx);
                atomsToBalance.insert(atomsToBalance.end(), currentSubset.connectedIdx.begin(), currentSubset.connectedIdx.end());
                atomsToBalance.insert(atomsToBalance.end(), currentSubset.hydrogenIdx.begin(), currentSubset.hydrogenIdx.end());

                size_t nextIdx = currentSubset.bondedSubsetIdx;
                if (nextIdx >= subsets.size() || nextIdx == SIZE_MAX) 
                    break; 
                currentIdx = nextIdx;
            }

            if (atomsToBalance.empty()) return;

            std::vector<std::vector<size_t>> neighborLists(atoms.size());
            for (const auto& bond : bonds)
            {
                if (std::find(atomsToBalance.begin(), atomsToBalance.end(), bond.bondedAtom) != atomsToBalance.end() &&
                    std::find(atomsToBalance.begin(), atomsToBalance.end(), bond.centralAtom) != atomsToBalance.end())
                {
                    neighborLists[bond.bondedAtom].push_back(bond.centralAtom);
                    neighborLists[bond.centralAtom].push_back(bond.bondedAtom);
                }
            }

            for (size_t idx : atomsToBalance)
            {
                if (atoms[idx].bondCount == 0) continue;
                const auto& neighbors = neighborLists[idx];
                if (neighbors.empty()) continue;

                float total_charge = atoms[idx].charge;
                for (size_t j : neighbors)
                    total_charge += atoms[j].charge;

                if (std::abs(total_charge) > 1e-6f)
                {
                    uint8_t valence = constants::getUsualBonds(atoms[idx].ZIndex);
                    float adjustment = -total_charge / (neighbors.size() + 1);

                    atoms[idx].charge -= adjustment;
                    atoms[idx].charge += valence - atoms[idx].bondCount;  

                    for (size_t j : neighbors)
                        atoms[j].charge += adjustment;
                }
            }
        }

        void universe::boundCheck(size_t i)
        {
            positions[i].x = std::fmod(positions[i].x + boxSize, boxSize);
            positions[i].y = std::fmod(positions[i].y + boxSize, boxSize);
        }

        float universe::ljPot(size_t i, float epsilon_i, float sigma_i)
        {
            float potential = 0.f;

            for (size_t j = 0; j < atoms.size(); ++j)
            {
                if (j == i) continue;
                
                float dr = minImageVec(positions[i] - positions[j]).length();

                auto [sigma_j, epsilon_j] = constants::getAtomConstants(atoms[j].ZIndex);
                float sigma = (sigma_i + sigma_j) / 2.0f;

                if (dr < sigma * CUTOFF && dr > EPSILON)
                {
                    float epsilon = sqrtf(epsilon_i * epsilon_j);
                    float r6 = powf((sigma / dr), 6); 
                    float r12 = r6 * r6; 
                    potential += 4.0f * epsilon * (r12-r6);
                }
            }

            return potential;
        }

        sf::Vector3f universe::ljForce(size_t i, size_t j)
        {
            const atom& a1 = atoms[i];
            const atom& a2 = atoms[j];

            const float sigma_i = a1.sigma;
            const float sigma_j = a2.sigma;
            const float epsilon_i = a1.epsilon;
            const float epsilon_j = a2.epsilon;

            sf::Vector3f dr_vec = minImageVec(positions[i] - positions[j]);
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
                if (j == i) continue;

                gradient += ljForce(i, j);
            }

            return gradient;
        }

        sf::Vector3f universe::coulombForce(size_t i, size_t j, sf::Vector3f& dr_vec)
        {
            const atom& a1 = atoms[i];
            const atom& a2 = atoms[j];

            float dr = dr_vec.length();

            if (dr < EPSILON || dr > COULOMB_CUTOFF)
                return {0.f, 0.f, 0.f};

            float qq = a1.charge * a2.charge;
            if (qq == 0.f) return {0.f, 0.f, 0.f};

            float forceMag = COULOMB_K * a1.charge * a2.charge / dr;
            return -forceMag * dr_vec / dr; 
        }

        void universe::calcBondForces()
        {
            for (size_t i = 0; i < bonds.size(); ++i)
            {
                bond& bond = bonds[i];

                size_t idx1 = bond.bondedAtom;
                size_t idx2 = bond.centralAtom;
                sf::Vector3f r_vec = minImageVec(positions[idx2] - positions[idx1]);
                float dr = r_vec.length();
                if (dr <= EPSILON) return;

                float delta_r = dr - bond.equilibriumLength;

                if (delta_r > bond.equilibriumLength * 1.5f)
                {
                    breakBond(idx1, idx2);
                    continue;
                }

                float force_magnitude = BOND_K * delta_r;

                sf::Vector3f force_dir = r_vec / dr;
                sf::Vector3f force = force_magnitude * force_dir;

                forces[idx1] += force;
                forces[idx2] -= force;
            }
        }

        void universe::calcAngleForces()
        {
            for (size_t s = 0; s < subsets.size(); ++s)
            {
                const subset& current_subset = subsets[s];
                auto neighbours = current_subset.connectedIdx;
                neighbours.insert(neighbours.begin(), current_subset.hydrogenIdx.begin(), current_subset.hydrogenIdx.end());
                size_t mainAtom = current_subset.mainAtomIdx;
                float ideal_theta = current_subset.idealAngle;

                if (neighbours.size() < 2) continue; // need > 2 for angle

                for (size_t i = 0; i < neighbours.size(); ++i)
                {
                    for (size_t k = i + 1; k < neighbours.size(); ++k)
                    {
                        size_t idx_i = neighbours[i];
                        size_t idx_k = neighbours[k];

                        sf::Vector3f r_ji = minImageVec(positions[idx_i] - positions[mainAtom]);
                        sf::Vector3f r_jk = minImageVec(positions[idx_k] - positions[mainAtom]);
                        float r_ji_len = r_ji.length();
                        float r_jk_len = r_jk.length();
                        if (r_ji_len <= EPSILON || r_jk_len <= EPSILON) continue;

                        float inv_denominator = (r_ji_len * r_jk_len);
                        float cos_theta = r_ji.dot(r_jk) / inv_denominator;
                        cos_theta = std::clamp(cos_theta, -1.0f, 1.0f);
                        float theta = std::acos(cos_theta);

                        float delta_theta = theta - ideal_theta;
                        float force_magnitude = ANGLE_K * delta_theta;

                        sf::Vector3f f_i = (force_magnitude / inv_denominator) * (r_jk - cos_theta * r_ji);
                        sf::Vector3f f_k = (force_magnitude / inv_denominator) * (r_ji - cos_theta * r_jk);
                        sf::Vector3f f_j = f_i - f_k;

                        forces[idx_i] += f_i;
                        forces[mainAtom] += f_j;
                        forces[idx_k] += f_k;
                    }
                }
            }
        }

        void universe::calcLjForces()
        {
            //size_t count = 0;

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                std::vector<size_t>& neighbours = neighbourList[i];

                for (size_t j = 0; j < neighbours.size(); ++j)
                {
                    size_t index_j = neighbours[j];

                    const float sigma = (atoms[i].sigma + atoms[index_j].sigma) / 2.0f;
                    if (areBonded(i, j)) continue;

                    //++count;
                    
                    sf::Vector3f force = ljForce(i, index_j);

                    forces[i] += force;
                    forces[index_j] -= force; // Every action creates an Equal and opposite reaction - Sir Isaac Newton
                }
            }

            //std::cout << "count: " << count << std::endl;
        }

        void universe::calcElectrostaticForces()
        {
            float cutoff2 = COULOMB_CUTOFF * COULOMB_CUTOFF;
            for (size_t i = 0; i < atoms.size(); ++i)
            {
                for (size_t j : neighbourList[i])  
                {
                    sf::Vector3f dr = minImageVec(positions[j] - positions[i]);
                    float dr2 = dr.lengthSquared();
                    if (dr2 > cutoff2 || dr2 < EPSILON * EPSILON) continue;

                    sf::Vector3f force = coulombForce(i, j, dr);
                    forces[i] += force;
                    forces[j] -= force;
                }
            }
        }

        void universe::update(float targetTemperature, bool reactions)
        {
            std::vector<sf::Vector3f> old_forces = forces;
            std::fill(forces.begin(), forces.end(), sf::Vector3f{0.f, 0.f, 0.f});
            
            if (neighbourList.empty() || neighbourListNeedRebuild()) 
                buildNeighborList();  

            calcLjForces();    
            calcBondForces();
            calcAngleForces();
            calcElectrostaticForces();

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector3f a_old = old_forces[i] / atoms[i].mass;
                sf::Vector3f a_new = forces[i] / atoms[i].mass;
                positions[i] += velocities[i] * DT + 0.5f * (a_old + a_new) * DT * DT;
                boundCheck(i);
            }

            std::fill(forces.begin(), forces.end(), sf::Vector3f{0.f, 0.f, 0.f});
            calcLjForces();    
            calcBondForces();
            calcAngleForces();
            calcElectrostaticForces();

            for (size_t i = 0; i < atoms.size(); ++i)
                velocities[i] += (old_forces[i] + forces[i]) * (0.5f * DT / atoms[i].mass);

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

                for (size_t i = 0; i < velocities.size(); ++i) 
                    velocities[i] *= lambda;
            }
        }

        // Energy Calculation
        float universe::calculateKineticEnergy()
        {
            float kinetic_energy = 0.0f;
            for (size_t i = 0; i < atoms.size(); ++i) 
                kinetic_energy += 0.5f * atoms[i].mass * velocities[i].lengthSquared();
            
            return kinetic_energy;
        }

        void universe::buildNeighborList()
        {
            if (neighbourList.size() != atoms.size())
                neighbourList.resize(atoms.size());

            const float rskin = COULOMB_CUTOFF + VERLET_SKIN;
            const float rskin2 = rskin * rskin;

            for (auto& list : neighbourList) list.clear();

            for (size_t i = 0; i < atoms.size(); ++i) {
                for (size_t j = i + 1; j < atoms.size(); ++j) {
                    if (areBonded(i, j)) continue;

                    sf::Vector3f dr = minImageVec(positions[i] - positions[j]);
                    if (dr.lengthSquared() < rskin2) {
                        neighbourList[i].push_back(j);
                    }
                }
            }

            neighbourInterval = 0.0f;
            neighbourPos = positions;
        }

        bool universe::neighbourListNeedRebuild()
        {
            if (neighbourList.size() != atoms.size())
                return true;

            const float threshold2 = REBUILD_THRESHOLD * REBUILD_THRESHOLD; 

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector3f d = positions[i] - neighbourPos[i];
                if (d.lengthSquared() > threshold2)
                    return true;
            }
            return false;
        }

        // Reactions

        float universe::calculateAtomTemperature(size_t i)
        {
            float ke = 0.5f * atoms[i].mass * velocities[i].lengthSquared();
            return (2.0f / 3.0f) * ke * KB;
        }

        float universe::calculateBondEnergy(size_t i, size_t j)
        {
            return 0.f;
        }

        float universe::calculateNonBondedEnergy(size_t i, size_t j)
        {
            return 0.f;
        }

        void universe::breakBond(size_t atom1, size_t atom2)
        {
            auto it = std::find_if(bonds.begin(), bonds.end(), [&](const bond& bond)
            {
                return bond.centralAtom == atom1 && bond.bondedAtom == atom2 
                    || bond.centralAtom == atom2 && bond.bondedAtom == atom1;
            });

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
                    [&](const subset& s) {
                        return s.mainAtomIdx == C &&
                            std::find(s.hydrogenIdx.begin(),
                                        s.hydrogenIdx.end(), H) != s.hydrogenIdx.end();
                    });

                if (sub_it != subsets.end())
                {
                    //balanceMolecularCharges(*sub_it);
                    sub_it->hydrogenIdx.erase(
                        std::remove(sub_it->hydrogenIdx.begin(),
                                    sub_it->hydrogenIdx.end(), H),
                        sub_it->hydrogenIdx.end());
                }
                return;
            }

            auto subset_it1 = std::find_if(subsets.begin(), subsets.end(),
                [&](const subset& s) { return s.mainAtomIdx == atom1; });

            auto subset_it2 = std::find_if(subsets.begin(), subsets.end(),
                [&](const subset& s) { return s.mainAtomIdx == atom2; });

            if (subset_it1 != subsets.end())
            {
                auto& s1 = *subset_it1;
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
                auto& s2 = *subset_it2;
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
            if (areBonded(i, j)) return BondType::NONE;

            uint8_t Z1 = atoms

            [i].ZIndex, Z2 = atoms[j].ZIndex;
                float r = minImageVec(positions[i] - positions[j]).length();

            uint8_t val1 = atoms[i].bondCount;
            uint8_t val2 = atoms[j].bondCount;
            uint8_t max1 = constants::getUsualBonds(Z1);
            uint8_t max2 = constants::getUsualBonds(Z2);

            for (int8_t order = 3; order >= 1; --order)
            {
                BondType type = static_cast<BondType>(order - 1);
                float r0 = constants::getBondLength(Z1, Z2, type);

                if (r0 <= 0.0f) continue;
                if (r > r0 * 1.3f || r < r0 * 0.5f) continue;

                if (val1 + order > max1 || val2 + order > max2) continue;

                return type; 
            }

            return BondType::NONE; 
        }

        void universe::checkStrainedBonds()
        {
        
        }

        void universe::handleReactions()
        {
            std::vector<std::pair<size_t, size_t>> candidates;

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                if (!isReactive(i)) continue;

                for (size_t j : neighbourList[i])
                {
                    if (!isReactive(j)) continue;
                    if (areBonded(i, j)) continue;

                    BondType bestType = chooseBestBondType(i, j);
                    float r = minImageVec(positions[i] - positions[j]).length();
                    float r0 = constants::getBondLength(atoms[i].ZIndex, atoms[j].ZIndex, bestType);
                    if (r0 > 0 && r < r0 * 1.3f) 
                    {
                        createBond(i, j, bestType);
                        auto subset_it1 = std::find_if(subsets.begin(), subsets.end(),
                            [&](const subset& s) { return s.mainAtomIdx == i; });

                        auto subset_it2 = std::find_if(subsets.begin(), subsets.end(),
                            [&](const subset& s) { return s.mainAtomIdx == j; });

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
        sf::Vector2f universe::project(core::window_t& window, const sf::Vector3f& p) const
        {
            auto size = window.getWindow().getView().getSize();
            return cam.project(p, size.x, size.y);
        }

        void universe::handleCamera(bool leftDown, bool rightDown, const sf::Vector2i& mousePos, float wheelDelta, const std::vector<sf::Keyboard::Key>& keys)
        {
            if (leftDown)
            {
                sf::Vector2i delta = mousePos - lastMouse;
                cam.azimuth   -= delta.x * 0.3f;
                cam.elevation  = std::clamp(cam.elevation - delta.y*0.3f,
                                            -89.f, 89.f);
            }

            if (rightDown)
            {
                sf::Vector2i delta = mousePos - lastMouse;
                sf::Vector3f forward = (cam.target-cam.eye()).normalized();
                sf::Vector3f right   = forward.cross({0,0,1}).normalized();
                sf::Vector3f up      = right.cross(forward);
                cam.target += sf::Vector3f{right.x * delta.x, right.y, right.z} * cam.distance * 0.002f +
                              sf::Vector3f{up.x, up.y * delta.y, up.z} * cam.distance * 0.002f;
            }

            if (wheelDelta != 0.f)
            {
                cam.distance *= (wheelDelta>0) ? 0.9f : 1.1f;
                cam.distance = std::max(50.f, std::min(600.f, cam.distance));
            }

            for (auto k : keys)
            {
                if (k == sf::Keyboard::Key::X) { cam.azimuth=0.f;   cam.elevation=0.f; }
                if (k == sf::Keyboard::Key::Y) { cam.azimuth=90.f;  cam.elevation=0.f; }
                if (k == sf::Keyboard::Key::Z) { cam.azimuth=45.f;  cam.elevation=90.f;}
                if (k == sf::Keyboard::Key::R) { cam = core::camera_t(); cam.distance=180.f; }
            }

            lastMouse = mousePos;
        }
    } // namespace fun 
} // namespace sim
