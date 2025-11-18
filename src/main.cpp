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

    size_t universeSize = 50.f;
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
    auto sac = sim::parseSMILES("HH");
    auto sac2 = sim::parseSMILES("O=O");

    auto stuff = sim::parseSMILES("C1C=C(C(=CC1(NC2=NC(=NC(=N2)N(CCO)CCO)N(CCO)CCO)NC3=NC(=NC(=N3)N(CCO)CCO)N(CCO)CCO)S(=O)(=O)[O-])/C=C/C4=CC=CC=C4S(=O)(=O)[O-].C(CO)NCCO.[Na+].[K+]");
    auto stuff2 = sim::parseSMILES("CC(C)CC(C(=O)NC(C)C(=O)NC(CO)C(=O)NC(C(C)O)C(=O)NC(CC1=CN=CN1)C(=O)NC(CCCCN)C(=O)NC(CC2=CC=CC=C2)C(=O)NC(C)C(=O)NC(CC3=CC=CC=C3)C(=O)NC(C(C)C)C(=O)NC(CCCCN)C(=O)NC(C)C(=O)NC(CC4=CC=CC=C4)C(=O)NC(C)C(=O)NC(CCCCN)C(=O)NC(C(C)O)C(=O)NC(CCCCN)C(=O)O)N");
    //universe.createMolecule(sac, {10, 10, 15});
    universe.createMolecule(stuff, {30, 30, 30});
    //universe.createMolecule(stuff2, {30, 50, 30});

    //universe.createMolecule(methane, {30, 30, 30});

    size_t count = 100;
    size_t count2 = 50;
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

    sf::Clock deltaClock;
    while (window.isOpen())
    {
        window.pollEvents(); 
        ImGui::SFML::Update(window.getWindow(), deltaClock.restart());

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
            universe.update(targetTemp, false);

        window.clear();
        universe.drawDebug(window);
        universe.draw(window, true);

        displayUI(window, universe);
        ImGui::SFML::Render(window.getWindow());
        window.display();
    }
}