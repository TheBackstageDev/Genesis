#include "core/window.hpp"
#include "simulation/fundamental_structures.hpp"

#include <random>

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
    size_t universeSize = 200.f;

    sim::fun::universe universe(universeSize);
    
    setupUI(window);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis(10.0f, universeSize - 10.f); 
    std::uniform_real_distribution<> vel(-5.f, 5.f); 

    float start_x = 100.f;
    float start_y = universeSize / 2.f;
    
    std::vector<sf::Vector2f> oxygen_positions;
    
    for (int i = 0; i < 30; ++i)
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
                if (distance < 20.0f) 
                {
                    valid_position = false;
                    break;
                }
            }
        }
    
        oxygen_positions.push_back({ox, oy});
    
        size_t o_idx = universe.numAtoms();
        universe.createAtom({ox, oy}, {static_cast<float>(vel(gen)), static_cast<float>(vel(gen))}, 8, 8);
    
        float angle = -45.0f * RADIAN; 
        float bond_length = 8.0f;
    
        size_t h1_idx = universe.numAtoms();
        universe.createAtom({ox + bond_length * cos(angle), oy + bond_length * sin(angle)}, {static_cast<float>(vel(gen)), static_cast<float>(vel(gen))}, 1);
    
        size_t h2_idx = universe.numAtoms();
        universe.createAtom({ox + bond_length * cos(-angle), oy + bond_length * sin(-angle)}, {static_cast<float>(vel(gen)), static_cast<float>(vel(gen))}, 1);
    
        universe.createBond(o_idx, h1_idx);
        universe.createBond(o_idx, h2_idx);
    }

    size_t carbonIndex = universe.numAtoms();
    universe.createAtom({20.f, 100.f}, {1.f, 2.f}, 6, 6);
    universe.createAtom({10.f, 100.f}, {0.f, 3.f}, 8, 8);
    universe.createAtom({30.f, 100.f}, {0.f, -1.f}, 8, 8);
    universe.createBond(universe.numAtoms() - 1, carbonIndex, sim::fun::BondType::DOUBLE);
    universe.createBond(universe.numAtoms() - 2, carbonIndex, sim::fun::BondType::DOUBLE);

    float targetTemp = 100.0f;
    
    while (window.isOpen())
    {
        window.pollEvents();
        window.refresh();
        
        if (!window.isPaused())
        universe.update(targetTemp);
        
        window.clear();
        universe.drawDebug(window);
        universe.draw(window, true);
        displayUI(window, universe);
        window.display();
    }
}