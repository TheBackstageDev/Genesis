#include "core/window.hpp"
#include "simulation/fundamental_structures.hpp"
#include "simulation/smiles_parser.hpp"

#include <iostream>
#include <random>
#include <chrono>

#include <algorithm>

sf::Text timeStepText;
sf::Text tempText;

void displayUI(core::window_t& window, sim::fun::universe &universe)
{
    auto& uiview = window.getuiview();
    window.getWindow().setView(uiview);

    sf::Vector2u windowSize = window.getWindow().getSize();

    uiview.setSize({static_cast<float>(windowSize.x), static_cast<float>(windowSize.y)});
    uiview.setCenter({windowSize.x / 2.0f, windowSize.y / 2.0f});
    uiview.setViewport(sf::FloatRect({0.f, 0.f}, {1.f, 1.f}));

    std::ostringstream tempStream;
    tempStream << std::fixed << std::setprecision(2) << "Temp: " << universe.temperature() << " K";
    tempText.setString(tempStream.str());

    std::ostringstream timeStream;
    timeStream << std::fixed << std::setprecision(1) << "Time: " << static_cast<float>(universe.timestep() * DT) << " ps";
    timeStepText.setString(timeStream.str());

    window.draw(timeStepText);
    window.draw(tempText);
}

void setupUI(core::window_t& window)
{
    timeStepText.setFont(window.getFont());
    timeStepText.setCharacterSize(25);

    tempText.setFont(window.getFont());
    tempText.setCharacterSize(25);
    tempText.setPosition({0.f, 26.f});
}

int main()
{    
    core::window_t window(500, 500, "Genesis Engine");
    size_t universeSize = 50.f;

    sim::fun::universe universe(universeSize);
    
    setupUI(window);

    window.setCameraCallback([&](bool left, bool right, const sf::Vector2i& mouse, float wheel, const std::vector<sf::Keyboard::Key>& keys)
    {
        universe.handleCamera(left, right, mouse, wheel, keys);
    });

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis(2.f, universeSize - 2.f); 
    float targetTemp = 100.f;

    auto water = sim::parseSMILES("O");  
    auto hy = sim::parseSMILES("C1CNCCN1c(c2)c(F)cc3c2N(C4CC4)C=C(C3=O)C(=O)O");  

    universe.createMolecule(hy, {30, 30, 1});
    //universe.createMolecule(water, {100, 40, 0});

    size_t count = 0;
    float minDistance = 3.f;

    std::vector<sf::Vector3f> centers{{150.f, 150.f, 0.f}};
    centers.reserve(count + 1);

    const float minDistSq = minDistance * minDistance;

    for (size_t i = 0; i < count; ++i)
    {
        sf::Vector3f pos;
        bool valid = false;

        for (int32_t attempt = 0; attempt < 100 && !valid; ++attempt)
        {
            pos.x = dis(gen);
            pos.y = dis(gen);
            pos.z = dis(gen);

            valid = true;
            for (const auto& c : centers)
            {
                float dx = pos.x - c.x;
                float dy = pos.y - c.y;
                float dz = pos.z - c.z;
                if (dx*dx + dy*dy + dz*dz < minDistSq)
                {
                    valid = false;
                    break;
                }
            }
        }

        if (valid)
        {
            centers.push_back(pos);
            universe.createMolecule(water, pos);
        }
    }

    while (window.isOpen())
    {
        window.pollEvents(); 
        window.refresh();
        
        if (!window.isPaused())
        {
            auto start_time = std::chrono::steady_clock::now();
            universe.update(targetTemp);
            auto end_time = std::chrono::steady_clock::now();
            auto duration = end_time - start_time;
            double delta_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            std::cout << delta_time_ms << "ms" << std::endl;
        }

        if (window.stepFrame())
        {
            universe.update(targetTemp, false);
            targetTemp += 1.f;
        }

        window.clear();
        universe.drawDebug(window);
        universe.draw(window, true);
        displayUI(window, universe);
        window.display();
    }
}