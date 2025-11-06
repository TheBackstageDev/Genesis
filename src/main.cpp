#include "core/window.hpp"
#include "simulation/fundamental_structures.hpp"
#include "simulation/smiles_parser.hpp"

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
    size_t universeSize = 500.f;

    sim::fun::universe universe(universeSize);
    
    setupUI(window);

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dis(10.0f, universeSize - 10.f); 
    std::uniform_real_distribution<> vel(-5.f, 5.f); 

    float targetTemp = 100.0f;

    auto cyclopropane = sim::parseSMILES("CC1(C2CCC1(C=C2)C)C");
    universe.createMolecule(cyclopropane, sf::Vector3f(250.f , 250.f, 0.f));

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
            universe.update(targetTemp);
        
        window.clear();
        universe.drawDebug(window);
        universe.draw(window, true);
        displayUI(window, universe);
        window.display();
    }
}