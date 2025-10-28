#include "core/window.hpp"
#include "simulation/fundamental_structures.hpp"

#include <iostream>
#include <random>
#include <chrono>

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
    size_t universeSize = 400.f;

    sim::fun::universe universe(universeSize);
    
    setupUI(window);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis(10.0f, universeSize - 10.f); 
    std::uniform_real_distribution<> vel(-5.f, 5.f); 
    
    std::vector<sf::Vector2f> oxygen_positions;

    for (int i = 0; i < 0; ++i)
    {
        bool valid_position = false;
        float ox, oy;
    
        while (!valid_position)
        {
            ox = dis(gen);
            oy = dis(gen);
            valid_position = true;
    
            for (const auto& pos : oxygen_positions)
            {
                float dx = ox - pos.x;
                float dy = oy - pos.y;
                float distance = std::sqrt(dx * dx + dy * dy);
                if (distance < 25.0f) 
                {
                    valid_position = false;
                    break;
                }
            }
        }
    
        oxygen_positions.push_back({ox, oy});
    
        size_t o_idx = universe.createAtom({ox, oy}, {0.f, 0.f}, 8, 8, 8);
    
        float angle = -45.0f * RADIAN; 
        float bond_length = 8.0f;
    
        size_t h1_idx = universe.createAtom({ox + bond_length * cos(angle), oy + bond_length * sin(angle)}, {1.f, 0.f}, 1);
        size_t h2_idx = universe.createAtom({ox + bond_length * cos(-angle), oy + bond_length * sin(-angle)}, {0.f, 1.f}, 1);
    
        universe.createSubset(o_idx, -1, -1, {h1_idx, h2_idx}, {sim::fun::BondType::SINGLE, sim::fun::BondType::SINGLE});
    }
    
    sim::fun::molecule_structure benzene{};
    benzene.atoms = {
        {6, 6, 0}, // C1
        {6, 6, 0}, // C2
        {6, 6, 0}, // C3
        {6, 6, 0}, // C4
        {6, 6, 0}, // C5
        {6, 6, 0}, // C6
        {1, 1, 0}, // H1
        {1, 1, 0}, // H2
        {1, 1, 0}, // H3
        {1, 1, 0}, // H4
        {1, 1, 0}, // H5
        {1, 1, 0}  // H6
    };
    benzene.bonds = {
        {6, 0},  // H1 to C1
        {7, 1},  // H2 to C2
        {8, 2},  // H3 to C3
        {9, 3},  // H4 to C4
        {10, 4}, // H5 to C5
        {11, 5}, // H6 to C6
        // C-C bonds (adjacent carbons)
        {1, 0},  // C2 to C1
        {2, 1},  // C3 to C2
        {3, 2},  // C4 to C3
        {4, 3},  // C5 to C4
        {5, 4},  // C6 to C5
        {0, 5}   // C1 to C6 (closes the ring)
    };
    benzene.bondTypes = 
    {
        sim::fun::BondType::SINGLE,
        sim::fun::BondType::SINGLE,
        sim::fun::BondType::SINGLE,
        sim::fun::BondType::SINGLE,
        sim::fun::BondType::SINGLE,
        sim::fun::BondType::SINGLE,
        
        sim::fun::BondType::SINGLE,
        sim::fun::BondType::DOUBLE,
        sim::fun::BondType::SINGLE,
        sim::fun::BondType::DOUBLE,
        sim::fun::BondType::SINGLE,
        sim::fun::BondType::DOUBLE,
    };
    benzene.subsets = 
    {
        {0, 1, 5, {6}}, // C1 subset
        {1, 2, 0, {7}}, // C2 subset
        {2, 3, 1, {8}}, // C3 subset
        {3, 4, 2, {9}}, // C4 subset
        {4, 5, 3, {10}}, // C5 subset
        {5, 0, 4, {11}}  // C6 subset
    };

    universe.createMolecule(benzene, {200.5f, 20.5f});

    float targetTemp = 100.0f;

    bool test = false;
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
            //std::cout << delta_time_ms << "ms" << std::endl;
        }
        else
        {
            if (!test)
            {
                //universe.createMolecule(benzene, {100.5f, 20.5f});
                test = true;
            }
        }
        
        window.clear();
        universe.drawDebug(window);
        universe.draw(window, true);
        displayUI(window, universe);
        window.display();
    }
}