#include "window.hpp"
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

namespace core
{
    window_t::window_t(uint32_t width, uint32_t height, const std::string& title, float boxSize, uint32_t framerate)
        : width(width), height(height), boxSize(boxSize), window(sf::VideoMode({width, height}), title)
    {
        window.setFramerateLimit(framerate);
        view.setCenter(sf::Vector2f(boxSize / 2.0f, boxSize / 2.0f));
        view.setSize(sf::Vector2f(boxSize, boxSize));
        window.setView(view);
    }

    window_t::~window_t()
    {
    
    }

    bool window_t::pollEvents()
    {
        while (const std::optional event = window.pollEvent())
        {
            if (event->is<sf::Event::Closed>())
            {
                window.close();
                return false;
            }
            if (event->is<sf::Event::Resized>())
            {
                const auto& resizedEvent = event->getIf<sf::Event::Resized>()->size;
                width = resizedEvent.x;
                height = resizedEvent.y;

                view.setSize(sf::Vector2f(boxSize, boxSize));
                float windowAspect = static_cast<float>(width) / height;
                float simAspect = 1.0f;
                if (windowAspect > simAspect)
                {
                    float viewportWidth = simAspect / windowAspect;
                    sf::FloatRect rect(sf::Vector2((1.0f - viewportWidth) / 2.0f, 0.0f), sf::Vector2(viewportWidth, 1.0f));
                    view.setViewport(rect);
                }
                else
                {
                    float viewportHeight = windowAspect / simAspect;
                    sf::FloatRect rect(sf::Vector2(0.0f, (1.0f - viewportHeight) / 2.0f), sf::Vector2(1.0f, viewportHeight));
                    view.setViewport(rect);
                }

                window.setView(view);
            }
        }
        return true;
    }
} // namespace core
