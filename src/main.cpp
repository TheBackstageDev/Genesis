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
    size_t universeSize = 1000.f;

    sim::fun::universe universe(universeSize);
    
    setupUI(window);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis(10.0f, universeSize - 10.f); 
    std::uniform_real_distribution<> vel(-3.f, 3.f); 

    std::vector<sf::Vector2f> oxygen_positions;

    for (int i = 0; i < 10; ++i)
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
                if (distance < 30.0f) 
                {
                    valid_position = false;
                    break;
                }
            }
        }

        oxygen_positions.push_back({ox, oy});

        size_t o_idx = universe.numAtoms();
        universe.createAtom({ox, oy}, {static_cast<float>(vel(gen)), static_cast<float>(vel(gen))}, 8);

        float angle = -45.0f * 3.1415926535f / 180.0f; 
        float bond_length = 10.0f;

        size_t h1_idx = universe.numAtoms();
        universe.createAtom({ox + bond_length * cos(angle), oy + bond_length * sin(angle)}, {static_cast<float>(vel(gen)), static_cast<float>(vel(gen))}, 1);

        size_t h2_idx = universe.numAtoms();
        universe.createAtom({ox + bond_length * cos(-angle), oy + bond_length * sin(-angle)}, {static_cast<float>(vel(gen)), static_cast<float>(vel(gen))}, 1);

        universe.createBond(o_idx, h1_idx);
        universe.createBond(o_idx, h2_idx);
    }

    float start_x = 300.f;
    float start_y = universeSize / 2.f;

    std::vector<size_t> carbon_indices;
    size_t numCarbons = 3;
    for (int i = 0; i < numCarbons; ++i)
    {
        float x = start_x + i * 12.0f;
        if (x > universeSize) break;
        size_t c_idx = universe.numAtoms();
        universe.createAtom({x, start_y}, {static_cast<float>(vel(gen)), static_cast<float>(vel(gen))}, 6); 
        carbon_indices.push_back(c_idx);

        float angle1 = -45.0f * 3.1415926535f / 180.0f; 
        float angle2 = 45.0f * 3.1415926535f / 180.0f; 
        float bond_length = 20.0f;

        size_t h1_idx = universe.numAtoms();
        universe.createAtom({x + bond_length * cos(angle1), start_y + bond_length * sin(angle1)}, {0.f, 0.0f}, 1);

        size_t h2_idx = universe.numAtoms();
        universe.createAtom({x + bond_length * cos(angle2), start_y + bond_length * sin(angle2)}, {0.0f, 0.f}, 1);

        universe.createBond(c_idx, h1_idx);
        universe.createBond(c_idx, h2_idx);
        
        if (i == 0)
        {
            size_t h3_idx = universe.numAtoms();
            universe.createAtom({x - bond_length, start_y + bond_length * sin(0.f)}, {0.0f, 0.0f}, 1);
            universe.createBond(c_idx, h3_idx);
        }
        if (i == numCarbons - 1)
        {
            size_t h3_idx = universe.numAtoms();
            universe.createAtom({x + bond_length, start_y - bond_length * sin(0.f)}, {0.0f, 0.0f}, 1);
            universe.createBond(c_idx, h3_idx);
        }
    }

    if (numCarbons > 0)
        for (size_t i = 0; i < carbon_indices.size() - 1; ++i)
        {
            universe.createBond(carbon_indices[i], carbon_indices[i + 1]);
        }

    float targetTemp = 70.f;

    while (window.isOpen())
    {
        window.pollEvents();
        window.refresh();

        if (!window.isPaused())
            universe.update(targetTemp);

        window.clear();
        universe.drawDebug(window);
        universe.draw(window, false);
        displayUI(window, universe);
        window.display();
    }
}