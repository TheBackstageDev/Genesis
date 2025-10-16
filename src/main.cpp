#include "core/window.hpp"

int main()
{
    core::window_t window(500, 500, "Genesis Engine");

    sf::CircleShape shape(0.5f);
    shape.setFillColor(sf::Color::Green);
    shape.setPosition({5.0f, 5.0f});

    while (window.isOpen())
    {
        window.pollEvents();
    
        window.clear();
        window.draw(shape);
        window.display();
    }
}