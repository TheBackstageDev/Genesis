#include "core/window.hpp"
#include "simulation/fundamental_structures.hpp"

#include <random>

int main()
{
    core::window_t window(500, 500, "Genesis Engine");
    sim::fun::universe universe(20.f);

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<float> posDist(2.0f, 15.0f); // Å

    std::vector<sf::Vector2f> positions;
    const int maxAttempts = 100;
    while (positions.size() < 20)
    {
        sf::Vector2f pos(posDist(rng), posDist(rng));
        bool tooClose = false;
        for (const auto& existingPos : positions)
        {
            sf::Vector2f r_vec = pos - existingPos;
            float r = std::sqrt(r_vec.x * r_vec.x + r_vec.y * r_vec.y);
            if (r < 2.f)
            {
                tooClose = true;
                break;
            }
        }
        if (!tooClose)
        {
            positions.push_back(pos);
            universe.createAtom(pos, {0.f, 0.f});
        }
    }

    while (window.isOpen())
    {
        window.pollEvents();
        window.refresh();

        universe.update();

        window.clear();
        universe.drawDebug(window);
        universe.draw(window);
        window.display();
    }
}