#include "fundamental_structures.hpp"

#include <iostream>

namespace sim
{
    namespace fun
    {
        void atom::draw(core::window_t& window)
        {
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

        void universe::draw(core::window_t& window)
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
                atom.draw(window);
            }
        }

        void universe::drawDebug(core::window_t& window)
        {

        }

        void universe::createAtom(sf::Vector2f p, sf::Vector2f v)
        {
            atom newAtom{};
            newAtom.position = p;
            newAtom.velocity = v;
            newAtom.radius = LJ_SIGMA;
            newAtom.mass = MASS_PROTON + MASS_ELECTRON;

            atoms.emplace_back(std::move(newAtom));
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

        float universe::ljPot(float i, float epsilon, float sigma)
        {
            float potential = 0.f;

            for (size_t j = 0; j < atoms.size(); ++j)
            {
                if (j == i) continue;

                float dr = (atoms[i].position - atoms[j].position).length();
                float r6 = powf((sigma / dr), 6); 
                float r12 = r6 * r6; 

                potential += 4.0f * epsilon * (r12-r6);
            }

            return potential;
        }

        sf::Vector2f universe::ljGrad(float i, float epsilon, float sigma)
        {
            sf::Vector2f gradient({0.f, 0.f});

            for (size_t j = 0; j < atoms.size(); ++j)
            {
                if (j == i) continue;

                sf::Vector2f dr_vec = atoms[i].position - atoms[j].position;
                float dr = dr_vec.length();
                float r8 = powf(sigma, 6) / powf((dr), 8);
                float r14 = 2.f * powf(sigma, 12) / powf((dr), 14);
                float du_dr = 24.0f * epsilon * (r14 - r8) / dr;

                gradient += (du_dr / dr) * dr_vec;
            }

            return gradient;
        }

        void universe::update()
        {
            bonds.clear();
            std::vector<sf::Vector2f> forces(atoms.size(), {0.0f, 0.0f});

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector2f grad = ljGrad(i, LJ_EPSILON, LJ_SIGMA);
                forces[i] += grad;
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
            float current_temp = (2.f/3.f) * avg_KE * KB; 
            float lambda = sqrtf(TARGET_TEMP/current_temp);
            lambda = (lambda - 1.0f) * 0.5f + 1.0f; // update slower

            for (auto& atom : atoms)
            {
                atom.velocity *= lambda;
            }
 
            std::cout << "Temperature: " << current_temp << " Kelvin \n";

            ++timeStep;
        }
    } // namespace fun 
} // namespace sim
