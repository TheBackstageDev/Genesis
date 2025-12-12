#include "window.hpp"
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include <ostream>
#include <iostream>

namespace core
{
    window_t::window_t(uint32_t width, uint32_t height, const std::string& title, float boxSize, uint32_t framerate)
        : width(width), height(height), boxSize(boxSize), window(sf::VideoMode({width, height}), title)
    {
        window.setFramerateLimit(framerate);
        view.setCenter(sf::Vector2f(boxSize / 2.0f, boxSize / 2.0f));
        view.setSize(sf::Vector2f(boxSize, boxSize));
        window.setView(view);

        uiView.setSize({static_cast<float>(width), static_cast<float>(height)});
        uiView.setCenter({width / 2.0f, height / 2.0f});
        uiView.setViewport(sf::FloatRect({0.f, 0.f}, {1.f, 1.f}));

        if (!arial.openFromFile("src/resource/ARIAL.ttf"))
        {
            std::cerr << "Failed to load font. UI and letter mode may not display correctly.\n";
        }
    }

    window_t::~window_t()
    {
    
    }

    bool window_t::pollEvents()
    {
        step = false;

        for (std::optional<sf::Event> optEvent = window.pollEvent(); optEvent.has_value(); optEvent = window.pollEvent())
        {
            const sf::Event& event = *optEvent;
            ImGui::SFML::ProcessEvent(window, event);
            if (event.is<sf::Event::Closed>())
            {
                window.close();
                return false;
            }
            if (event.is<sf::Event::Resized>())
            {
                const auto& resizedEvent = event.getIf<sf::Event::Resized>()->size;
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

                float boxSizeAspect = simAspect * boxSize;
                view.setCenter(sf::Vector2f(boxSizeAspect / 2.0f, boxSizeAspect / 2.0f));
                view.setSize(sf::Vector2f(boxSizeAspect, boxSizeAspect));
                window.setView(view);
            }
            if (event.is<sf::Event::KeyPressed>())
            {
                const auto& key = event.getIf<sf::Event::KeyPressed>()->code;
                if (key == sf::Keyboard::Key::Space)
                    paused = !paused;
                if (key == sf::Keyboard::Key::F)
                    step = true;
            }
            if (event.is<sf::Event::MouseWheelScrolled>())
            {
                float zoomFactor = 1.0f - event.getIf<sf::Event::MouseWheelScrolled>()->delta * ZOOM_SPEED;
                sf::Vector2f currentSize = view.getSize();
                sf::Vector2f newSize = currentSize * zoomFactor;
    
                float minSize = MIN_ZOOM * window.getSize().x;
                float maxSize = MAX_ZOOM * window.getSize().x;
                if (newSize.x >= minSize && newSize.x <= maxSize)
                {
                    view.zoom(zoomFactor);
                    window.setView(view);
                }
            }
        }

        return true;
    }

    void window_t::handleCameraInput()
    {
        sf::Vector2f pan(0.f, 0.f);

        float zoom = view.getSize().length();
        float zoom_speed = zoom / 100.f;

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::W))
            pan.y -= PAN_SPEED * zoom_speed; // Move up
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::S))
            pan.y += PAN_SPEED * zoom_speed; // Move down
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::A))
            pan.x -= PAN_SPEED * zoom_speed; // Move left
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::D))
            pan.x += PAN_SPEED * zoom_speed; // Move right

        view.move(pan);
        window.setView(view);
    }
} // namespace core
