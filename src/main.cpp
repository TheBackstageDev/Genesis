#include "core/window.hpp"
#include "simulation/universe.hpp"
#include "simulation/smiles_parser.hpp"

#include <imgui/imgui.h>
#include <imgui/imgui-SFML.h>

#include <iostream>
#include <random>
#include <chrono>

#include <algorithm>

float targetTemp = 100.f;

void displayUI(core::window_t& window, sim::fun::universe &universe)
{
    ImGui::Begin("Info");

    ImGui::DragFloat("Target Temp", &targetTemp, 0.01f, 0.01f, 0.f);
    ImGui::Text("Current Temp: %.2fK", universe.temperature());
    ImGui::Text("Current Time: %.3fps", universe.timestep() * DT);

    ImGui::End();
}

int main()
{    
    core::window_t window(500, 500, "Genesis Engine");
    if (!ImGui::SFML::Init(window.getWindow()))
    {
        throw std::runtime_error("Failed to init imgui!");
    }

    size_t universeSize = 30.f;
    sim::fun::universe universe(universeSize);

    window.setCameraCallback([&](bool left, bool right, const sf::Vector2i& mouse, float wheel, const std::vector<sf::Keyboard::Key>& keys)
    {
        universe.handleCamera(left, right, mouse, wheel, keys);
    });

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis(2.f, universeSize - 2.f); 
    std::uniform_real_distribution<> ve(-5.f, 5.f); 

    auto water = sim::parseSMILES("O");  
    auto stuff = sim::parseSMILES("CC[C@@H](C)[C@H]1C(=O)N[C@H](C(=O)NC(C(=O)N[C@H](CSSC(C(=O)N[C@H](C(=O)N1)CC2=CC=C(C=C2)O)N)C(=O)N3CCC[C@@H]3C(=O)N[C@H](CC(C)C)C(=O)O)CC(=O)N)CCC(=O)N");
 
    universe.createMolecule(stuff, {15, 15, 15});

    size_t count = 70;
    float minDistance = 3.f;

    std::vector<sf::Vector3f> centers{universe.positions()};

    const float minDistSq = minDistance * minDistance;

    for (size_t i = 0; i < count; ++i)
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
            universe.createMolecule(water, pos);
        }
    }

    sf::Clock deltaClock;
    while (window.isOpen())
    {
        window.pollEvents(); 
        ImGui::SFML::Update(window.getWindow(), deltaClock.restart());

        window.refresh();
        
        if (!window.isPaused())
        {
            auto start_time = std::chrono::steady_clock::now();
            universe.update(targetTemp, false);
            auto end_time = std::chrono::steady_clock::now();
            auto duration = end_time - start_time;
            double delta_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            std::cout << delta_time_ms << "ms" << std::endl;
        }

        if (window.stepFrame())
            universe.update(targetTemp, false);

        window.clear();
        //universe.drawDebug(window);
        universe.draw(window, true);

        displayUI(window, universe);
        ImGui::SFML::Render(window.getWindow());
        window.display();
    }
}