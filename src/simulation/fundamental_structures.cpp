#include "fundamental_structures.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <unordered_set>

#include "constants.hpp"

namespace sim
{
    sf::Color lerpColor(sf::Color a, sf::Color b, float t)
    {
        return sf::Color(
            static_cast<uint8_t>(a.r + (b.r - a.r) * sin(t)),
            static_cast<uint8_t>(a.g + (b.g - a.g) * sin(t)),
            static_cast<uint8_t>(a.b + (b.b - a.b) * sin(t)),
            255);
    }

    namespace fun
    {
        void atom::draw(const sf::Vector2f pos, const float camDistance, const float q, core::window_t &window, sf::RenderTarget &target)
        {
            sf::Color atomColor = constants::getElementColor(ZIndex);

            float rad = constants::VDW_RADII[ZIndex] * 3.f / camDistance;

            auto &font = window.getFont();
            std::string name = constants::getAtomLetter(ZIndex);
            sf::Text name_text(font, name, 30.f);
            name_text.setScale({0.05f * rad, 0.05f * rad});
            name_text.setFillColor({255,
                                    255,
                                    255,
                                    255});

            sf::FloatRect b = name_text.getLocalBounds();
            name_text.setOrigin(b.getCenter());
            name_text.setPosition(pos);
            target.draw(name_text);

            // Charge
            if (std::abs(q) >= 1.f)
            {
                sf::Text ion_text(font, q > 0.0f ? "+" : "-", 15.f);

                ion_text.setScale({0.03f * rad, 0.03f * rad});
                ion_text.setFillColor({255, 255, 255, 255});

                sf::FloatRect ionBounds = ion_text.getLocalBounds();
                ion_text.setOrigin(ionBounds.getCenter());

                float offsetScale = rad * 0.5f;
                ion_text.setPosition(pos + sf::Vector2f(offsetScale, -offsetScale));

                target.draw(ion_text);
            }
        }
    } // namespace fun
} // namespace sim
