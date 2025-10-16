#pragma once

#include <SFML/Graphics.hpp>

namespace core
{
    class window_t
    {
    public:
        window_t(uint32_t width, uint32_t height, const std::string& title, float boxSize = 10.f, uint32_t framerate = 60);
        ~window_t();

        bool pollEvents(); 

        bool isOpen() const
        {
            return window.isOpen();
        }

        void draw(sf::Drawable& drawable)
        {
            window.setView(view);
            window.draw(drawable);
        }

        void clear(sf::Color color = sf::Color::Black)
        {
            window.clear(sf::Color::Black); 
        }

        void display()
        {
            window.display();
        }

        sf::RenderWindow& getWindow() { return window; }

    private:
        uint32_t height; uint32_t width;
        
        sf::RenderWindow window;
        sf::View view;

        float boxSize;
    };
} // namespace core
