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
    sim::fun::universe universe(251.f);
    
    setupUI(window);
    
    universe.createAtom({50.0f, 50.0f}, {0.0f, 0.0f}, 8); 
    universe.createAtom({50.0f + LJ_SIGMA_H * 7.f, 50.0f}, {0.0f, 0.0f}, 1); 
    universe.createAtom({50.0f, 50.0f + LJ_SIGMA_H * 7.f}, {0.0f, 0.0f}, 1); 

    universe.createBond(0, 1);
    universe.createBond(0, 2); 

    float targetTemp = 14.f;

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