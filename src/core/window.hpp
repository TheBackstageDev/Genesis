#pragma once

#include <SFML/Graphics.hpp>
#include <functional>

#include <imgui/imgui.h>
#include "imgui/imgui-SFML.h"

#define MIN_ZOOM 0.001f
#define MAX_ZOOM 100.f
#define ZOOM_SPEED 0.1f

#define PAN_SPEED 1.f

namespace core
{
    class window_t
    {
    public:
        window_t(uint32_t width, uint32_t height, const std::string& title, float boxSize = 50.f, uint32_t framerate = 60);
        ~window_t();

        bool pollEvents(); 

        bool isOpen() const { return window.isOpen(); }

        void refresh()
        {
            window.setView(view);
        }

        void draw(sf::Drawable& drawable) { window.draw(drawable); }

        void clear(sf::Color color = sf::Color::Black) { window.clear(sf::Color::Black); }
        void display() { window.display(); }
        void handleCameraInput();

        sf::View& getuiview() { return uiView; };
        sf::Font& getFont() { return arial; }
        
        sf::RenderWindow& getWindow() { return window; }

    private:
        uint32_t height; uint32_t width;
        float boxSize;
        
        sf::RenderWindow window;
        sf::View view;
        sf::View uiView;

        sf::Font arial;
    };
} // namespace core
