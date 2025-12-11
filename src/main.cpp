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

    size_t universeSize = 30.f;
    sf::Vector3f box(universeSize, universeSize, universeSize);

    sim::fun::universe_create_info create_info{};
    create_info.has_gravity = true;
    create_info.reactive = false;
    create_info.wall_collision = true;
    create_info.render_water = true;
    create_info.mag_gravity = 9.8f;
    create_info.box = box;

    sim::fun::universe universe(create_info);

    window.setCameraCallback([&](bool left, bool right, const sf::Vector2i& mouse, float wheel, const std::vector<sf::Keyboard::Key>& keys)
    {
        universe.handleCamera(left, right, mouse, wheel, keys);
    });

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> disx(3.f, box.x - 3.f); 
    std::uniform_real_distribution<> disy(3.f, box.y - 3.f); 
    std::uniform_real_distribution<> disz(3.f, box.z / 2.f);
    std::uniform_real_distribution<> disz_gas(box.z / 2.f + 2.f, box.z - 2.f);

    std::uniform_real_distribution<> atmosphere(0.f, 1.f);

    std::uniform_real_distribution<> ve(-0.f, 0.f); 

    auto water = sim::parseSMILES("O"); 
    auto benzene = sim::parseSMILES("c1ccccc1");  
    auto HydrochloricAcid = sim::parseSMILES("[Cl-].[H+]");  
    auto Amonia = sim::parseSMILES("N");  
    auto Methane = sim::parseSMILES("C");  
    auto stuff = sim::parseSMILES("N[C@@H](CCCNC(N)=N)C(=O)N1CCC[C@H]1C(=O)N2CCC[C@H]2C(=O)NC(CC3=CC=CC=C3)C(=O)NC(CC4=CC=CC=C4)C(=O)N[C@@H](CO)C(=O)N5CCC[C@H]5C(=O)N[C@@H](CC6=CC=CC=C6)C(=O)N[C@@H](CCCNC(N)=N)C(=O)O");
    auto stuff2 = sim::parseSMILES("CCCCCCCCCCCCOS(=O)(=O)[O-]");  
    auto stuff3 = sim::parseSMILES("C1=CC2=CC3=CC4=C5C6=C3C7=C2C(=C1)C8=C9C7=C1C6=C2C3=C6C1=C1C9=C(C=CC1=CC6=CC1=C3C3=C6C(=C1)C=C1C=CC=C7C1=C6C1=C6C3=C2C5=C2C6=C3C(=CC2=C4)C=CC2=C3C1=C7C=C2)C=C8");

    auto H2 = sim::parseSMILES("HH"); 
    auto N2 = sim::parseSMILES("N#N"); 
    auto O2 = sim::parseSMILES("O=O"); 
    auto Ar = sim::parseSMILES("Ar"); 
    auto CO2 = sim::parseSMILES("O=C=O"); 
    auto CO = sim::parseSMILES("[C-]#[O+]"); 

    //universe.createMolecule(water, {30, 50, 7}, {0.f, 0.0f, 0.f});
    //universe.createMolecule(benzene, {5, 8, 20}, {0.f, 0.f, 0.f});
    //universe.createMolecule(HydrochloricAcid, {5, 8, 5}, {0.f, 0.f, 0.f});
    //universe.createMolecule(Amonia, {5, 12, 5}, {0.f, -0.1f, 0.f});
    //universe.createMolecule(Methane, {5, 12, 5}, {0.f, -0.1f, 0.f});
    //universe.createMolecule(O2, {5, 5, 5}, {0.f, 0.f, 0.f});
    //universe.createMolecule(H2, {3, 7, 5}, {0.f, -0.1f, 0.f});
    //universe.createMolecule(H2, {7, 7, 5}, {0.f, -0.1f, 0.f});
    //universe.createMolecule(stuff, {5, 15, 15}, {0.f, 0.f, 0.f});
    //universe.createMolecule(stuff2, {5, 5, 10}, {0.f, 0.f, 0.f});
    //universe.createMolecule(stuff2, {5, 9, 10}, {0.f, 0.f, 0.f});
    //universe.createMolecule(stuff2, {5, 13, 10}, {0.f, 0.f, 0.f});
    //universe.createMolecule(stuff2, {5, 17, 10}, {0.f, 0.f, 0.f});
    //universe.createMolecule(stuff3, {5, 17, 3}, {0.f, 0.f, 0.f});
    //universe.createMolecule(CO, {5, 5, 5}, {0.f, 0.f, 0.f});

    /* molecule_structure mol{};
    auto dna = sim::io::loadXYZ("src/resource/protein.xyz", mol.atoms, mol.bonds, mol.positions);
    sim::organizeSubsets(mol.subsets, mol.atoms, mol.bonds);
    sim::organizeAngles(mol.subsets, mol.atoms, mol.bonds, mol.dihedral_angles, mol.angles);

    universe.createMolecule(mol, {20, 20, 20}, {0.f, 0.f, 0.0f}); */

    size_t count = 300;
    float minDistance = 3.5f;
    float minDistanceAtm = 3.5f;

    std::vector<sf::Vector3f> centers{universe.positions()};

    for (size_t i = 0; i < count; ++i)
    {
        sf::Vector3f pos;
        sf::Vector3f vel;
        bool valid = false;

        for (int32_t attempt = 0; attempt < 100 && !valid; ++attempt)
        {
            pos.x = disx(gen);
            pos.y = disy(gen);
            pos.z = disz(gen);

            vel.x = ve(gen);
            vel.y = ve(gen);
            vel.z = ve(gen);

            valid = true;
            for (const auto& c : centers)
            {
                float dx = pos.x - c.x;
                float dy = pos.y - c.y;
                float dz = pos.z - c.z;
                if (dx*dx + dy*dy + dz*dz < minDistance*minDistance)
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

    for (size_t i = 0; i < count; ++i)
    {
        sf::Vector3f pos;
        sf::Vector3f vel;
        bool valid = false;

        for (int32_t attempt = 0; attempt < 100 && !valid; ++attempt)
        {
            pos.x = disx(gen);
            pos.y = disy(gen);
            pos.z = disz_gas(gen);

            vel.x = ve(gen);
            vel.y = ve(gen);
            vel.z = ve(gen);

            valid = true;
            for (const auto& c : centers)
            {
                float dx = pos.x - c.x;
                float dy = pos.y - c.y;
                float dz = pos.z - c.z;
                if (dx*dx + dy*dy + dz*dz < minDistanceAtm*minDistanceAtm)
                {
                    valid = false;
                    break;
                }
            }
        }

        if (valid)
        {
            molecule_structure air_particle; 

            float atm_chance = atmosphere(gen);  // 0.0 – 1.0

            if (atm_chance < 0.0004)
                air_particle = Methane;
            else if (atm_chance < 0.0093)
                air_particle = CO2;
            else if (atm_chance < 0.0188)
                air_particle = Ar;
            else if (atm_chance < 0.2093)     
                air_particle = O2;
            else
                air_particle = N2;

            centers.push_back(pos);
            universe.createMolecule(air_particle, pos);
        }
    }

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
            drawCharge = !drawCharge;
        }

        window.clear();
        //universe.drawDebug(window);
        
        universe.draw(window, true);
        if (drawCharge) universe.drawChargeField(window);

        displayUI(window, universe);
        ImGui::SFML::Render(window.getWindow());
        window.display();
    }
}