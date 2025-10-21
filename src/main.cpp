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
    sim::fun::universe universe(201.f);
    
    setupUI(window);
    
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> posDist(1.0f, 200.0f); // Å
    std::uniform_real_distribution<float> velDist(-1.f, 1.f);

    std::vector<sf::Vector2f> positions;
    const int maxAttempts = 100;
    while (positions.size() < 200)
    {
        sf::Vector2f pos(posDist(rng), posDist(rng));
        bool tooClose = false;
        for (const auto &existingPos : positions)
        {
            sf::Vector2f r_vec = pos - existingPos;
            float r = std::sqrt(r_vec.x * r_vec.x + r_vec.y * r_vec.y);
            if (r < 5.f)
            {
                tooClose = true;
                break;
            }
        }
        if (!tooClose)
        {
            positions.push_back(pos);
            universe.createAtom(pos, {velDist(rng), velDist(rng)});
        }
    }

    for (uint32_t Z = 1; Z <= 1; ++Z)
    {
        sf::Vector2f pos(posDist(rng), posDist(rng));
        universe.createAtom(pos, {velDist(rng), velDist(rng)}, 8);
    }

    // universe.createAtom({posDist(rng), posDist(rng)}, {0.f, 0.f}, 2);

    while (window.isOpen())
    {
        window.pollEvents();
        window.refresh();

        if (!window.isPaused())
            universe.update(0.1f);

        window.clear();
        universe.drawDebug(window);
        universe.draw(window, true);
        displayUI(window, universe);
        window.display();
    }
}