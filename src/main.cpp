#include "core/window.hpp"
#include "simulation/universe.hpp"
#include "simulation/smiles_parser.hpp"
#include "simulation/format_loader.hpp"

#include <imgui/imgui.h>
#include <imgui/imgui-SFML.h>

#include <iostream>
#include <random>
#include <chrono>

#include <algorithm>

float targetTemp = 300.f;
float targetPressure = 30.f;

void displayUI(core::window_t& window, sim::fun::universe& universe)
{
    ImGui::Begin("Genesis Engine");
    ImGui::DragFloat("Target Temp (K)", &targetTemp, 0.1f, 0.001f, 12000.f);
    ImGui::DragFloat("Target Pressure (Bar)", &targetPressure, 0.1f, 0.f, 12000.f);
    ImGui::Text("Current Pressure: %.1f Bar", universe.pressure());
    ImGui::Text("Current Temp: %.1f K", universe.temperature());
    ImGui::Text("Atoms: %zu", universe.numAtoms());
    ImGui::Text("Bonds: %zu", universe.numBonds());
    ImGui::Text("Molecules: %zu", universe.numMolecules());
    ImGui::Text("Time: %.3f ps", universe.timestep() * DT);
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
    sf::Vector3f box(universeSize, universeSize, universeSize);

    sim::fun::logging_flags log_flags{};
    log_flags.log_reactions = true;

    sim::fun::universe_create_info create_info{};
    create_info.has_gravity = false;
    create_info.reactive = false;
    create_info.wall_collision = false;
    create_info.render_water = true;
    create_info.log_flags = log_flags;
    create_info.mag_gravity = 9.8f;
    create_info.box = box;

    //sim::fun::universe universe(create_info);
    sim::fun::universe universe("src/scenes/2025_12_1219_46.json");

    bool drawCharge = false;
    sf::Clock deltaClock;
    while (window.isOpen())
    {
        window.pollEvents(); 
        ImGui::SFML::Update(window.getWindow(), deltaClock.restart());

        window.refresh();
        
        if (!window.isPaused())
        {
            auto start_time = std::chrono::steady_clock::now();
            universe.update(targetTemp, targetPressure);
            auto end_time = std::chrono::steady_clock::now();
            auto duration = end_time - start_time;
            double delta_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            std::cout << delta_time_ms << "ms" << std::endl;
        }

        if (window.stepFrame())
        {
            universe.update(targetTemp, targetPressure);
            universe.saveScene("src/scenes");
            universe.saveAsVideo("src/scenes/videos");
        }

        window.clear();
        //universe.drawDebug(window);
        
        universe.draw(window, true);
        if (drawCharge) universe.drawChargeField(window);

        displayUI(window, universe);
        universe.handleCamera();

        ImGui::SFML::Render(window.getWindow());
        window.display();
    }
}