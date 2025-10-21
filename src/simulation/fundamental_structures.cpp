#include "fundamental_structures.hpp"

#include <iostream>

namespace sim
{
    namespace fun
    {
        void atom::draw(core::window_t& window, bool letter)
        {
            if (letter)
            {
                auto& font = window.getFont();

                std::string name = constants::getAtomLetter(ZIndex);

                sf::Text name_text;
                name_text.setFont(window.getFont());
                name_text.setString(name);
                name_text.setCharacterSize(20.f);
                name_text.setScale({0.1f * radius, 0.1f * radius});
                name_text.setPosition(position);

                window.draw(name_text);

                return;
            }

            sf::CircleShape cloud(radius);
            cloud.setPosition({position.x - radius, position.y - radius});
            cloud.setFillColor(sf::Color(50, 50, 255, 80));
            window.draw(cloud);

            float nucleusRadius = 0.01f; // Å
            sf::CircleShape nucleus(nucleusRadius);
            nucleus.setPosition({position.x - nucleusRadius, position.y - nucleusRadius});
            nucleus.setFillColor(sf::Color::Red);
            window.draw(nucleus);
        }

        universe::universe(float universeSize)
            : boxSize(universeSize)
        {
        }

        void universe::draw(core::window_t& window, bool letter)
        {
            for (auto& bond : bonds)
            {
                std::array<sf::Vertex, 2> line = {
                    sf::Vertex(atoms[bond.first].position, sf::Color::White),
                    sf::Vertex(atoms[bond.second].position, sf::Color::White)
                };

                window.getWindow().draw(line.data(), line.size(), sf::PrimitiveType::Lines);
            }

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

        void universe::createAtom(sf::Vector2f p, sf::Vector2f v, uint32_t ZIndex)
        {
            atom newAtom{};
            newAtom.ZIndex = ZIndex;
            newAtom.position = p;
            newAtom.velocity = v;

            std::pair<float, float> constants = constants::getAtomConstants(ZIndex);
            
            newAtom.sigma = constants.first;
            newAtom.epsilon = constants.second;
            newAtom.radius = constants.first / 2.f;
            newAtom.mass = ZIndex * MASS_PROTON;

            atoms.emplace_back(std::move(newAtom));

            forces.resize(atoms.size());
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
        }

        float universe::ljPot(size_t i, float epsilon, float sigma_i)
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
                    float r6 = powf((sigma / dr), 6); 
                    float r12 = r6 * r6; 
                    potential += 4.0f * epsilon * (r12-r6);
                }
            }

            return potential;
        }

        sf::Vector2f universe::ljGrad(size_t i, float epsilon, float sigma_i)
        {
            sf::Vector2f gradient({0.f, 0.f});

            for (size_t j = 0; j < atoms.size(); ++j)
            {
                if (j == i) continue;

                sf::Vector2f dr_vec = atoms[i].position - atoms[j].position;
                float dr = dr_vec.length();

                auto [sigma_j, epsilon_j] = constants::getAtomConstants(atoms[j].ZIndex);
                float sigma = (sigma_i + sigma_j) / 2.0f;

                if (dr < sigma * CUTOFF && dr > EPSILON)
                {

                    float r8 = powf(sigma, 6) / powf((dr), 8);
                    float r14 = 2.f * powf(sigma, 12) / powf((dr), 14);
                    float du_dr = 24.0f * epsilon * (r14 - r8) / dr;
                    sf::Vector2f force = (du_dr / dr) * dr_vec;
                    gradient += force;

                    if (j < forces.size()) // Every action creates an Equal and opposite reaction - Sir Isaac Newton
                        forces[j] -= force;
                }
            }

            return gradient;
        }

        void universe::update(float targetTemperature)
        {
            bonds.clear();

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                forces[i] = ljGrad(i, atoms[i].sigma, atoms[i].epsilon);
            }

            float kinetic_energy = 0.0f;
            for (size_t i = 0; i < atoms.size(); ++i)
            {
                atom& a = atoms[i];
                if (a.mass <= 0) continue;
                sf::Vector2f acc = forces[i] / a.mass;
                a.position += (a.velocity + acc) * DT;
                boundCheck(a);
                sf::Vector2f new_acc = forces[i] / a.mass;
                a.velocity += 0.5f * (acc + new_acc) * DT;

                float v_squared = a.velocity.lengthSquared();
                kinetic_energy += 0.5f * a.mass * v_squared;
            }

            float avg_KE = kinetic_energy / atoms.size();
            temp = (2.f/3.f) * avg_KE * KB; 
            float lambda = sqrtf(targetTemperature/temp);
            lambda = (lambda - 1.0f) * 0.5f + 1.0f; // update slower

            for (auto& atom : atoms)
            {
                atom.velocity *= lambda;
            }
 
            ++timeStep;
        }
    } // namespace fun 
} // namespace sim
