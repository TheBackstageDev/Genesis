#include "fundamental_structures.hpp"

#include <iostream>
#include <algorithm>

namespace sim
{
    constexpr float pixel_per_A = 1.0f; // 50 pixels/Å

    namespace fun
    {
        void atom::draw(core::window_t& window, bool letter)
        {
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
                name_text.setOrigin(bounds.size / 2.0f);

                name_text.setPosition(position * pixel_per_A);

                window.draw(name_text);

                return;
            }

            sf::CircleShape cloud(radius / 1.2);
            cloud.setPosition({position.x - radius, position.y - radius});
            cloud.setFillColor(sf::Color(50, 50, 255, 80));
            window.draw(cloud);

            float nucleusRadius = 0.01f;
            sf::CircleShape nucleus(nucleusRadius);
            nucleus.setPosition(position);
            nucleus.setFillColor(sf::Color::Red);

            window.draw(nucleus);
        }

        universe::universe(float universeSize)
            : boxSize(universeSize)
        {
        }

        void universe::draw(core::window_t& window, bool letter)
        {
            for (auto& atom : atoms)
            {
                atom.draw(window, letter);
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

        void universe::createAtom(sf::Vector2f p, sf::Vector2f v, uint8_t ZIndex, uint8_t numElectron)
        {
            atom newAtom{};
            newAtom.ZIndex = ZIndex;
            newAtom.position = p;
            newAtom.velocity = v;

            std::pair<float, float> constants = constants::getAtomConstants(ZIndex);
            
            newAtom.sigma = constants.first;
            newAtom.epsilon = constants.second;
            newAtom.radius = constants.first;
            newAtom.electrons = numElectron;
            newAtom.charge = ZIndex - numElectron;
            newAtom.mass = ZIndex * MASS_PROTON + numElectron * MASS_ELECTRON;
            newAtom.bondCount = 0;

            atoms.emplace_back(std::move(newAtom));

            prev_positions.resize(atoms.size());
            forces.resize(atoms.size());
        }

        void universe::createBond(size_t idx1, size_t idx2, BondType type)
        {
            if (idx1 >= atoms.size() || idx2 >= atoms.size() || idx1 == idx2)
                return;

            // NOTE: Handle bond type transition
            if (std::find_if(bonds.begin(), bonds.end(), [&](bond& a){
                return a.atom1 == idx1 && a.atom2 == idx2;
            }) != bonds.end()) return; // bond alreadly exists within the two atoms

            bond nBond{};
            nBond.atom1 = idx1;
            nBond.atom2 = idx2;
            nBond.type = type;
            nBond.equilibriumLength = constants::getBondLength(atoms[idx1].ZIndex, atoms[idx2].ZIndex, type);

            ++atoms[idx1].bondCount;
            ++atoms[idx2].bondCount;

            bonds.emplace_back(std::move(nBond));
        }

        void universe::boundCheck(atom& a)
        {
            if (a.position.x < 0.f)
            {
                a.position.x = 0.f;
                a.velocity.x = -a.velocity.x;
            }
            else if (a.position.x > boxSize)
            {
                a.position.x = boxSize;
                a.velocity.x = -a.velocity.x;
            }

            if (a.position.y < 0.f)
            {
                a.position.y = 0.f;
                a.velocity.y = -a.velocity.y;
            }
            else if (a.position.y > boxSize)
            {
                a.position.y = boxSize;
                a.velocity.y = -a.velocity.y;
            }
            /* a.position.x = std::fmod(a.position.x + boxSize, boxSize);
            a.position.y = std::fmod(a.position.y + boxSize, boxSize); */
        }

        float universe::ljPot(size_t i, float epsilon_i, float sigma_i)
        {
            float potential = 0.f;

            for (size_t j = 0; j < atoms.size(); ++j)
            {
                if (j == i) continue;
                
                float dr = (atoms[i].position - atoms[j].position).length();

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

            sf::Vector2f dr_vec = a1.position - a2.position;
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

        void universe::calcBondForces()
        {
            for (size_t i = 0; i < bonds.size(); ++i)
            {
                bond& bond = bonds[i];

                size_t idx1 = bond.atom1;
                size_t idx2 = bond.atom2;
                sf::Vector2f r_vec = atoms[idx2].position - atoms[idx1].position;
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
            for (size_t j = 0; j < atoms.size(); ++j)
            {
                std::vector<size_t> connectedIdxs;
                std::vector<BondType> bondTypes;
                for (const auto& bond : bonds)
                {
                    if (bond.atom1 == j && std::find(connectedIdxs.begin(), connectedIdxs.end(), bond.atom2) == connectedIdxs.end())
                    {
                        connectedIdxs.push_back(bond.atom2);
                        bondTypes.push_back(bond.type);
                    }
                    else if (bond.atom2 == j && std::find(connectedIdxs.begin(), connectedIdxs.end(), bond.atom1) == connectedIdxs.end())
                    {
                        connectedIdxs.push_back(bond.atom1);
                        bondTypes.push_back(bond.type);
                    }
                }

                if (connectedIdxs.size() < 2) continue; // Need at least 2 neighbors for an angle

                for (size_t i = 0; i < connectedIdxs.size(); ++i)
                {
                    for (size_t k = i + 1; k < connectedIdxs.size(); ++k)
                    {
                        size_t idx_i = connectedIdxs[i];
                        size_t idx_k = connectedIdxs[k];

                        sf::Vector2f r_ji = atoms[idx_i].position - atoms[j].position;
                        sf::Vector2f r_jk = atoms[idx_k].position - atoms[j].position;
                        float r_ji_len = r_ji.length();
                        float r_jk_len = r_jk.length();
                        if (r_ji_len <= EPSILON || r_jk_len <= EPSILON) continue;

                        float cos_theta = r_ji.dot(r_jk) / (r_ji_len * r_jk_len);
                        cos_theta = std::clamp(cos_theta, -1.0f, 1.0f);
                        float theta = std::acos(cos_theta);

                        float ideal_theta = constants::getAngles(atoms[j].ZIndex, {atoms[idx_i].ZIndex, atoms[idx_k].ZIndex}, bondTypes);
                        if (ideal_theta == 0.0f) continue;

                        float delta_theta = theta - ideal_theta;
                        float force_magnitude = ANGLE_K * delta_theta;

                        sf::Vector2f f_i = (force_magnitude / (r_ji_len * r_jk_len)) * (r_jk - cos_theta * r_ji);
                        sf::Vector2f f_k = (force_magnitude / (r_ji_len * r_jk_len)) * (r_ji - cos_theta * r_jk);
                        sf::Vector2f f_j = -f_i - f_k;

                        forces[idx_i] += f_i;
                        forces[j] += f_j;
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
                    bool stopCalc{false};

                    for (const bond& bond : bonds)
                    {
                        if ((bond.atom1 == j && bond.atom2 == i) || (bond.atom1 == i && bond.atom2 == j))
                        {
                            stopCalc = true;
                        }
                    }

                    if (stopCalc)
                        continue;

                    sf::Vector2f force = ljForce(i, j);

                    forces[i] += force;
                    forces[j] -= force; // Every action creates an Equal and opposite reaction - Sir Isaac Newton
                }
            }
        }

        void universe::update(float targetTemperature)
        {
            std::fill(forces.begin(), forces.end(), sf::Vector2f{0.f, 0.f});
            calcLjForces();    
            calcBondForces();
            calcAngleForces();

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                atom& a = atoms[i];
                if (a.mass <= 0) continue;
                sf::Vector2f acc = forces[i] / a.mass;
                a.position += a.velocity * DT + acc * DT;
                boundCheck(a);
                sf::Vector2f new_acc = forces[i] / a.mass;
                a.velocity += 0.5f * DT * (acc + new_acc);

                prev_positions[i] = a.position;
            }

            if (timeStep % THERMOSTAT_INTERVAL == 0) {
                float kinetic_energy = 0.0f;
                for (const auto& a : atoms) {
                    float v_squared = a.velocity.lengthSquared();
                    kinetic_energy += 0.5f * a.mass * v_squared;
                }

                float avg_KE = kinetic_energy / atoms.size();
                temp = (2.f / 3.f) * avg_KE * KB;
                float lambda = sqrtf(targetTemperature / temp);
                // lambda = (lambda - 1.0f) * 0.5f + 1.0f; // update slower

                for (auto& atom : atoms) {
                    atom.velocity *= lambda;
                }
            }
 
            ++timeStep;
        }
    } // namespace fun 
} // namespace sim
