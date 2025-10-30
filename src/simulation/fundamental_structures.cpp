#include "fundamental_structures.hpp"

#include <iostream>
#include <algorithm>
//#include <thread>

namespace sim
{
    constexpr float pixel_per_A = 1.0f; // 50 pixels/Å

    namespace fun
    {
        void atom::draw(sf::Vector2f& pos, core::window_t& window, bool letter)
        {
            sf::Vector2f posVec{pos.x, pos.y};

            if (letter)
            {
                auto& font = window.getFont();
                std::string name = constants::getAtomLetter(ZIndex);
                sf::Text name_text;
                name_text.setFont(font);
                name_text.setString(name);
                name_text.setCharacterSize(20.f);
                name_text.setScale({0.02f * radius, 0.02f * radius});

                sf::FloatRect bounds = name_text.getLocalBounds();
                name_text.setOrigin(bounds.getCenter());  
                name_text.setPosition(pos);
                window.draw(name_text);
            }

            if (charge != 0.f && std::abs(charge) > 0.5f)
            {
                bool cation = charge > 0.f;
                sf::Text ion_text;
                auto& font = window.getFont();

                ion_text.setFont(font);
                ion_text.setString(cation ? "+" : "-");
                ion_text.setCharacterSize(16.f);
                ion_text.setScale({0.02f * radius, 0.02f * radius});

                sf::FloatRect bounds = ion_text.getLocalBounds();
                ion_text.setOrigin(bounds.getCenter());
                ion_text.setPosition(pos + sf::Vector2f(0.f, -radius * 0.3f));
                window.draw(ion_text);
            }

            const float cloudRadius = radius * 0.5f;

            sf::CircleShape cloud(cloudRadius);
            cloud.setOrigin({cloudRadius, cloudRadius});             
            cloud.setPosition(pos);
            cloud.setFillColor(sf::Color(50, 50, 255, 80));
            window.draw(cloud);
        }

        universe::universe(float universeSize)
            : boxSize(universeSize)
        {
        }

        void universe::draw(core::window_t& window, bool letter)
        {
            for (size_t i = 0; i < atoms.size(); ++i)
            {
                atoms[i].draw(positions[i], window, letter);
            }

            for (size_t b = 0; b < bonds.size(); ++b)
            {
                const bond& bond = bonds[b];
                const sf::Vector2f& pA = positions[bond.bondedAtom];   
                const sf::Vector2f& pB = positions[bond.centralAtom]; 

                if ((pA - pB).length() > 5.f * MULT_FACTOR) continue;

                sf::Vector2f dir  = (pA - pB).normalized();
                sf::Vector2f perp = -dir.perpendicular();       

                int8_t lines = 1;
                switch (bond.type)
                {
                    case BondType::DOUBLE: lines = 2; break;
                    case BondType::TRIPLE: lines = 3; break;
                    default:               lines = 1; break;
                }
                const float shrink = 0.5f * MULT_FACTOR;   

                std::vector<sf::Vertex> vertices;
                vertices.reserve(lines * 2);

                const float offsetStep = 0.2f * MULT_FACTOR;               
                for (int8_t i = 0; i < lines; ++i)
                {
                    float signedOffset = offsetStep * (i - (lines - 1) * 0.5f);
                    sf::Vector2f offset = perp * signedOffset;

                    sf::Vector2f start = pB + dir * shrink + offset;
                    sf::Vector2f end   = pA - dir * shrink + offset;

                    vertices.emplace_back(start);
                    vertices.emplace_back(end);
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

        void universe::drawDebug(core::window_t& window)
        {

        }

        void universe::log(size_t step)
        {

        }

        size_t universe::createAtom(sf::Vector2f p, sf::Vector2f v, uint8_t ZIndex, uint8_t numNeutrons, uint8_t numElectron)
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

            int8_t bondCount = static_cast<int8_t>(type) + 1;

            atoms[idx1].charge = 0.f;
            atoms[idx2].charge = 0.f;

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

            uint8_t valence1 = constants::getValenceElectrons(atoms[idx1].ZIndex);
            uint8_t valence2 = constants::getValenceElectrons(atoms[idx2].ZIndex);

            atoms[idx1].charge += valence1 - atoms[idx1].bondCount;
            atoms[idx2].charge += valence2 - atoms[idx2].bondCount;

            bonds.emplace_back(std::move(nBond));
        }

        size_t universe::createSubset(const size_t central, const size_t subsetNext, const size_t subsetLast, const size_t mainNext, const size_t mainLast, 
            const std::vector<size_t>& bonds, const std::vector<fun::BondType>& bondTypes)
        {
            subset nSubset{};
            nSubset.mainAtomIdx = central;
            nSubset.bondedSubsetIdx = subsetNext;
            nSubset.bondingSubsetIdx = subsetLast;
            
            std::vector<size_t> connected;
            std::vector<uint8_t> neighbourZs;
            connected.reserve(bonds.size());
            neighbourZs.reserve(bonds.size() + 2);

            for (size_t i = 0; i < bonds.size(); ++i)
            {
                const size_t bondedAtom = bonds[i];

                if (bondedAtom == mainLast) createBond(bondedAtom, central, bondTypes[i - 2]);
                else createBond(bondedAtom, central, bondTypes[i]);

                if (bondedAtom != mainNext && bondedAtom != mainLast)
                    connected.emplace_back(bondedAtom);

                neighbourZs.emplace_back(atoms[bondedAtom].ZIndex);
            }

            nSubset.connectedIdx = std::move(connected);
            nSubset.idealAngle = constants::getAngles(atoms[central].ZIndex, neighbourZs, bondTypes);
            nSubset.bondTypes = std::move(bondTypes);
            
            subsets.emplace_back(std::move(nSubset));

            return subsets.size() - 1; // Index
        }

        // TO DO: Support ionic-molecules 
        void universe::createMolecule(const molecule_structure& structure, sf::Vector2f pos)
        {
            size_t baseAtomIndex = atoms.size();
            
            sf::Vector2f current_pos = pos;
            for (const def_atom& a : structure.atoms)
            {
                createAtom(current_pos, {1.f, -1.f}, a.ZIndex, a.NIndex, a.ZIndex - a.charge);
            }

            std::map<size_t, std::vector<BondType>> bondOrder{};

            for (size_t b = 0; b < structure.bonds.size(); ++b)
            {
                const def_bond& db = structure.bonds[b];

                size_t globalIdx1 = db.bondingAtomIdx;
                size_t globalIdx2 = db.centralAtomIdx;
                bondOrder[globalIdx2].emplace_back(db.type);
            }

            size_t baseSubset = subsets.size();
            for (const def_subset& s : structure.subsets)
            {
                std::vector<size_t> mainAtomBonds{};
                mainAtomBonds.reserve(s.connectedIdx.size() + 2);
                
                for (size_t b = 0; b < s.connectedIdx.size(); ++b)
                {
                    mainAtomBonds.emplace_back(s.connectedIdx[b] + baseAtomIndex);
                }

                if (s.bondedSubset != SIZE_MAX)
                    mainAtomBonds.emplace_back(structure.subsets[s.bondedSubset].mainAtomIdx + baseAtomIndex);
                if (s.bondingSubset != SIZE_MAX)
                    mainAtomBonds.emplace_back(structure.subsets[s.bondingSubset].mainAtomIdx + baseAtomIndex);

                createSubset(s.mainAtomIdx + baseAtomIndex, s.bondedSubset + baseSubset, s.bondingSubset + baseSubset, structure.subsets[s.bondedSubset].mainAtomIdx, 
                    structure.subsets[s.bondingSubset].mainAtomIdx, mainAtomBonds, bondOrder[s.mainAtomIdx]);
            }
            
            if (structure.subsets.empty())
                organizeMolecule(structure, pos); // Creates subsets, automatically

            positionMolecule(baseSubset);
            balanceMolecularCharges(subsets[baseSubset]);
        }

        void universe::organizeMolecule(const molecule_structure& structure, const sf::Vector2f& initPos)
        {
            
        }

        void universe::positionMolecule(size_t firstSubsetIndex)
        {
            const float mult_factor = MULT_FACTOR;
            sf::Vector2f lastSubsetPos{positions[subsets[firstSubsetIndex].mainAtomIdx]};
            
            for (size_t s = firstSubsetIndex; s < subsets.size(); ++s)
            {
                const subset& sub = subsets[s];
                size_t mainAtomIdx = sub.mainAtomIdx;
                
                sf::Vector2f mainAtomPos{positions[mainAtomIdx].x + 0.001f, positions[mainAtomIdx].y}; 
                sf::Vector2f direction = (mainAtomPos - lastSubsetPos).normalized();
                sf::Vector2f perpendicular = -direction.perpendicular();
                
                positions[mainAtomIdx] = mainAtomPos;

                std::vector<size_t> neighbours = sub.connectedIdx;
                neighbours.emplace_back(sub.bondedSubsetIdx);
                neighbours.emplace_back(sub.bondingSubsetIdx);

                float angleOffset = 0.0f;
                for (size_t n = 0; n < neighbours.size(); ++n)
                {
                    size_t neighbour = neighbours[n];
                    if (neighbour == SIZE_MAX) continue;
                    
                    angleOffset += sub.idealAngle;

                    sf::Vector2f neighborDir = (direction * cos(angleOffset) + perpendicular * sin(angleOffset)).normalized() * mult_factor;
                    positions[neighbour] = mainAtomPos + neighborDir;
                }
                
                if (sub.bondedSubsetIdx != SIZE_MAX)
                {
                    size_t nextIdx = sub.bondedSubsetIdx;
                    const size_t start = nextIdx;
                    
                    sf::Vector2f prevPos = lastSubsetPos;
                    sf::Vector2f currentToPrev = (positions[subsets[nextIdx].mainAtomIdx] - mainAtomPos);
                    
                    lastSubsetPos = mainAtomPos;
                    float angle = sub.idealAngle;

                    while (nextIdx != SIZE_MAX)
                    {
                        const subset& nextSub = subsets[nextIdx];
                        size_t nextMainAtomIdx = nextSub.mainAtomIdx;

                        sf::Vector2f directionToNext = sf::Vector2f(
                            cos(angle) * currentToPrev.x - sin(angle) * currentToPrev.y,
                            sin(angle) * currentToPrev.x + cos(angle) * currentToPrev.y
                        ).normalized() * mult_factor;
                        positions[nextMainAtomIdx] = mainAtomPos + directionToNext * 1.5f;

                        sf::Vector2f perpendicular = -directionToNext.perpendicular();

                        angleOffset = 0.f;
                        for (size_t n = 0; n < nextSub.connectedIdx.size(); ++n)
                        {
                            size_t neighbour = nextSub.connectedIdx[n];
                            angleOffset += nextSub.idealAngle;

                            sf::Vector2f neighborDir = (directionToNext * cos(angleOffset) + perpendicular * sin(angleOffset)).normalized() * mult_factor;
                            positions[neighbour] = positions[nextMainAtomIdx] + neighborDir + perpendicular * 1.2f;;
                        }

                        nextIdx = nextSub.bondedSubsetIdx;

                        mainAtomPos = positions[nextMainAtomIdx];
                        currentToPrev = (mainAtomPos - prevPos).normalized();

                        if (nextIdx == start) break;
                    }
                }

                lastSubsetPos = mainAtomPos;
            }
        }

        void universe::balanceMolecularCharges(subset& mol)
        {
            std::vector<size_t> atomsToBalance;

            auto it = std::find_if(subsets.begin(), subsets.end(), [&](const subset& s) { return &s == &mol; });
            if (it == subsets.end()) return; 
            size_t molIndex = std::distance(subsets.begin(), it);
            const size_t initialMol = molIndex;

            while (molIndex != static_cast<size_t>(-1)) 
            {
                if (molIndex >= subsets.size()) break; 
                atomsToBalance.push_back(subsets[molIndex].mainAtomIdx);
                atomsToBalance.insert(atomsToBalance.end(), subsets[molIndex].connectedIdx.begin(), subsets[molIndex].connectedIdx.end());
                molIndex = subsets[molIndex].bondedSubsetIdx;

                if (molIndex == initialMol)
                    break;
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
                    uint8_t valence1 = constants::getValenceElectrons(atoms[idx].ZIndex);

                    float adjustment = -total_charge / (neighbors.size() + 1);
                    atoms[idx].charge += adjustment - atoms[idx].bondCount + valence1;
                    
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

        sf::Vector2f universe::ljForce(size_t i, size_t j)
        {
            const atom& a1 = atoms[i];
            const atom& a2 = atoms[j];

            const float sigma_i = a1.sigma;
            const float sigma_j = a2.sigma;
            const float epsilon_i = a1.epsilon;
            const float epsilon_j = a2.epsilon;

            sf::Vector2f dr_vec = minImageVec(positions[i] - positions[j]);
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

            return sf::Vector2f{0.f, 0.f};
        }

        sf::Vector2f universe::ljGrad(size_t i)
        {
            sf::Vector2f gradient({0.f, 0.f});

            for (size_t j = 0; j < atoms.size(); ++j)
            {
                if (j == i) continue;

                gradient += ljForce(i, j);
            }

            return gradient;
        }

        sf::Vector2f universe::coulombForce(size_t i, size_t j, sf::Vector2f& dr_vec)
        {
            const atom& a1 = atoms[i];
            const atom& a2 = atoms[j];

            float dr = dr_vec.length();

            if (dr < EPSILON * EPSILON || dr > COULOMB_CUTOFF * COULOMB_CUTOFF)
                return {0.f, 0.f};

            float qq = a1.charge * a2.charge;
            if (qq == 0.f) return {0.f, 0.f};

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
                sf::Vector2f r_vec = minImageVec(positions[idx2] - positions[idx1]);
                float dr = r_vec.length();
                if (dr <= EPSILON) return;

                float delta_r = dr - bond.equilibriumLength;
                float force_magnitude = BOND_K * delta_r;

                sf::Vector2f force_dir = r_vec / dr;
                sf::Vector2f force = force_magnitude * force_dir;

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
                size_t mainAtom = current_subset.mainAtomIdx;
                float ideal_theta = current_subset.idealAngle;

                if (neighbours.size() < 2) continue; // need > 2 for angle

                for (size_t i = 0; i < neighbours.size(); ++i)
                {
                    for (size_t k = i + 1; k < neighbours.size(); ++k)
                    {
                        size_t idx_i = neighbours[i];
                        size_t idx_k = neighbours[k];

                        sf::Vector2f r_ji = minImageVec(positions[idx_i] - positions[mainAtom]);
                        sf::Vector2f r_jk = minImageVec(positions[idx_k] - positions[mainAtom]);
                        float r_ji_len = r_ji.length();
                        float r_jk_len = r_jk.length();
                        if (r_ji_len <= EPSILON || r_jk_len <= EPSILON) continue;

                        float inv_denominator = (r_ji_len * r_jk_len);
                        float cos_theta = r_ji.dot(r_jk) / inv_denominator;
                        cos_theta = std::clamp(cos_theta, -1.0f, 1.0f);
                        float theta = std::acos(cos_theta);

                        float delta_theta = theta - ideal_theta;
                        float force_magnitude = ANGLE_K * delta_theta;

                        sf::Vector2f f_i = (force_magnitude / inv_denominator) * (r_jk - cos_theta * r_ji);
                        sf::Vector2f f_k = (force_magnitude / inv_denominator) * (r_ji - cos_theta * r_jk);
                        sf::Vector2f f_j = f_i - f_k;

                        forces[idx_i] += f_i;
                        forces[mainAtom] += f_j;
                        forces[idx_k] += f_k;
                    }
                }
            }
        }

        void universe::calcLjForces()
        {
            for (size_t i = 0; i < atoms.size(); ++i)
            {
                for (size_t j = i + 1; j < atoms.size(); ++j)
                {
                    const float sigma = (atoms[i].sigma + atoms[j].sigma) / 2.0f;
                    sf::Vector2f r = minImageVec(positions[i] - positions[j]);
                    if (areBonded(i, j)) continue;

                    sf::Vector2f force = ljForce(i, j);

                    forces[i] += force;
                    forces[j] -= force; // Every action creates an Equal and opposite reaction - Sir Isaac Newton
                }
            }
        }

        void universe::calcElectrostaticForces()
        {
            float cutoff2 = COULOMB_CUTOFF * COULOMB_CUTOFF;
            for (size_t i = 0; i < atoms.size(); ++i)
            {
                for (size_t j = i + 1; j < atoms.size(); ++j)
                {
                    if (areBonded(i, j)) continue; // Skip electrostatics for bonded atoms

                    sf::Vector2f dr = minImageVec(positions[j] - positions[i]);
                    float dr2 = dr.lengthSquared();

                    if (dr2 > cutoff2 || dr2 < EPSILON * EPSILON) continue;

                    sf::Vector2f force = coulombForce(i, j, dr);
                    forces[i] += force;      
                    forces[j] -= force;      
                }
            }
        }

        void universe::update(float targetTemperature)
        {
            std::fill(forces.begin(), forces.end(), sf::Vector2f{0.f, 0.f});
            calcLjForces();    
            calcBondForces();
            calcAngleForces();
            calcElectrostaticForces();

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                atom& a = atoms[i];
                if (a.mass <= 0) continue;
                sf::Vector2f acc = forces[i] / a.mass;

                positions[i] += velocities[i] * DT + 0.5f * acc * DT * DT;
                boundCheck(i);
            }

            std::fill(forces.begin(), forces.end(), sf::Vector2f{0.f, 0.f});
            calcLjForces();    
            calcBondForces();
            calcAngleForces();
            calcElectrostaticForces();

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                atom& a = atoms[i];
                if (a.mass <= 0) continue;
                sf::Vector2f acc = forces[i] / a.mass;

                velocities[i] += 0.5f * (acc + forces[i] / a.mass) * DT; 
            }

            if (timeStep % THERMOSTAT_INTERVAL == 0) {
                float kinetic_energy = 0.0f;
                for (size_t i = 0; i < atoms.size(); ++i) {
                    float v_squared = velocities[i].lengthSquared();
                    kinetic_energy += 0.5f * atoms[i].mass * v_squared;
                }

                float avg_KE = kinetic_energy / atoms.size();
                temp = (2.f / 3.f) * avg_KE * KB;
                float lambda = sqrtf(targetTemperature / temp);
                // lambda = (lambda - 1.0f) * 0.5f + 1.0f; // update slower

                for (size_t i = 0; i < velocities.size(); ++i) {
                    velocities[i] *= lambda;
                }
            }
 
            ++timeStep;
        }
    } // namespace fun 
} // namespace sim
