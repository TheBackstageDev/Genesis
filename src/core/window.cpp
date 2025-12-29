#include "window.hpp"
#include <SFML/Window.hpp>

#include <ostream>
#include <iostream>

namespace core
{
    window_t::window_t(uint32_t width, uint32_t height, const std::string& title, float boxSize, uint32_t framerate)
        : width(width), height(height), boxSize(boxSize)
    {
        sf::ContextSettings settings;
        settings.depthBits = 24;
        settings.stencilBits = 8;
        settings.antiAliasingLevel = 2;
        settings.majorVersion = 4;
        settings.minorVersion = 4;

        window.create(sf::VideoMode({width, height}), title, sf::Style::Default, sf::State::Windowed, settings);
        window.setActive(true);
        
        if (!gladLoadGLLoader((GLADloadproc)sf::Context::getFunction)) 
        {
            throw std::exception("[WINDOW]: glad Failed to load OpenGL!");
        }
        
        glEnable(GL_TEXTURE_2D);
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
        glEnable(GL_BLEND);
        glEnable(GL_CULL_FACE);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

        window.setFramerateLimit(framerate);
        view.setCenter(sf::Vector2f(boxSize / 2.0f, boxSize / 2.0f));
        view.setSize(sf::Vector2f(boxSize, boxSize));
        window.setView(view);

        uiView.setSize({static_cast<float>(width), static_cast<float>(height)});
        uiView.setCenter({width / 2.0f, height / 2.0f});
        uiView.setViewport(sf::FloatRect({0.f, 0.f}, {1.f, 1.f}));

        if (!arial.openFromFile("src/resource/fonts/Orbitron-Regular.ttf"))
        {
            std::cerr << "[WINDOW]: Failed to load font. letter mode may not display correctly.\n";
        }
    }

    window_t::~window_t()
    {
        window.setActive(false);
    }

    bool window_t::pollEvents()
    {
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

                glViewport(0, 0, width, height);

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
                window.setSize(sf::Vector2u(width, height));
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
