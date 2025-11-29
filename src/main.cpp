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

float targetTemp = 100.0f;

void displayUI(core::window_t& window, sim::fun::universe& universe)
{
    ImGui::Begin("Genesis Engine");
    ImGui::DragFloat("Target Temp (K)", &targetTemp, 0.1f, 1.f, 12000.f);
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

    size_t universeSize = 20.f;
    sim::fun::universe universe(universeSize, 13.f, false);

    window.setCameraCallback([&](bool left, bool right, const sf::Vector2i& mouse, float wheel, const std::vector<sf::Keyboard::Key>& keys)
    {
        universe.handleCamera(left, right, mouse, wheel, keys);
    });

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis(2.f, universeSize - 2.f); 
    std::uniform_real_distribution<> ve(-5.f, 5.f); 

    auto water = sim::parseSMILES("O");  
    auto HydrochloricAcid = sim::parseSMILES("[Cl-][H+]");  
    auto stuff = sim::parseSMILES("O=P(O)(O)OP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n2cnc1c(ncnc12)N)[C@H](O)[C@@H]3O");
    auto O_Radical = sim::parseSMILES("O", false);  

    //universe.createMolecule(HydrochloricAcid, {10, 10, 10}, {0.f, 0.f, 0.0f});
    universe.createMolecule(water, {10, 15, 10}, {0.f, -0.4f, 0});
    //universe.createMolecule(stuff, {10, 5, 5}, {0.1f, 0.f, 0.f});

    /* molecule_structure mol{};
    auto dna = sim::io::loadXYZ("src/resource/protein.xyz", mol.atoms, mol.bonds, mol.positions);
    sim::organizeSubsets(mol.subsets, mol.atoms, mol.bonds);
    sim::organizeAngles(mol.subsets, mol.atoms, mol.bonds, mol.dihedral_angles, mol.angles);

    universe.createMolecule(mol, {40, 40, 40}, {0.f, 0.f, 0.0f}); */

    size_t count = 0;
    float minDistance = 6.f;

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
            universe.update(targetTemp);
            auto end_time = std::chrono::steady_clock::now();
            auto duration = end_time - start_time;
            double delta_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
            std::cout << delta_time_ms << "ms" << std::endl;
        }

        if (window.stepFrame())
            universe.update(targetTemp);

        window.clear();
        //universe.drawDebug(window);
        universe.draw(window, true);
        universe.drawChargeField(window);

        displayUI(window, universe);
        ImGui::SFML::Render(window.getWindow());
        window.display();
    }
}