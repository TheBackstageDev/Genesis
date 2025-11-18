#include "core/window.hpp"
#include "simulation/universe.hpp"
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
    size_t universeSize = 100.f;

    sim::fun::universe universe(universeSize);
    
    setupUI(window);

    window.setCameraCallback([&](bool left, bool right, const sf::Vector2i& mouse, float wheel, const std::vector<sf::Keyboard::Key>& keys)
    {
        universe.handleCamera(left, right, mouse, wheel, keys);
    });

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis(2.f, universeSize - 2.f); 
    std::uniform_real_distribution<> ve(-5.f, 5.f); 
    float targetTemp = 200.f;

    auto sac = sim::parseSMILES("HH");
    auto sac2 = sim::parseSMILES("O=O");

    auto water = sim::parseSMILES("O");  

    auto stuff = sim::parseSMILES("CC(C)CC(C(=O)NC(C)C(=O)NC(CC1=CN=CN1)C(=O)NC(C(C)C)C(=O)NC(CCCCN)C(=O)NC(C)C(=O)NC(CC2=CC=CC=C2)C(=O)NC(C)C(=O)NC(CC3=CC=CC=C3)C(=O)NC(CCCCN)C(=O)NC(C(C)O)C(=O)NC(CCCCN)C(=O)O)N");
    auto stuff2 = sim::parseSMILES("CC(C)CC(C(=O)NC(C)C(=O)NC(CO)C(=O)NC(C(C)O)C(=O)NC(CC1=CN=CN1)C(=O)NC(CCCCN)C(=O)NC(CC2=CC=CC=C2)C(=O)NC(C)C(=O)NC(CC3=CC=CC=C3)C(=O)NC(C(C)C)C(=O)NC(CCCCN)C(=O)NC(C)C(=O)NC(CC4=CC=CC=C4)C(=O)NC(C)C(=O)NC(CCCCN)C(=O)NC(C(C)O)C(=O)NC(CCCCN)C(=O)O)N");
    //universe.createMolecule(sac, {10, 10, 15});
    universe.createMolecule(stuff, {30, 30, 30});
    universe.createMolecule(stuff2, {30, 50, 30});

    //universe.createMolecule(methane, {30, 30, 30});

    size_t count = 150;
    size_t count2 = 75;
    float minDistance = 2.f;

    std::vector<sf::Vector3f> centers{};

    const float minDistSq = minDistance * minDistance;

    for (size_t t = 0; t < 0; ++t)
    for (size_t i = 0; i < (t == 0 ? count : count2); ++i)
    {
        sf::Vector3f pos;
        sf::Vector3f vel;
        bool valid = false;

        for (int32_t attempt = 0; attempt < 100 && !valid; ++attempt)
        {
            pos.x = dis(gen);
            pos.y = dis(gen);
            pos.z = dis(gen);

            vel.x = ve(gen);
            vel.y = ve(gen);
            vel.z = ve(gen);

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
            universe.createMolecule(t == 0 ? sac : sac2, pos);
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