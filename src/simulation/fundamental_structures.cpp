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
        void atom::draw(float temperature, sf::Vector2f &pos, float camDistance, float q, core::window_t &window, sf::RenderTarget& target, const rendering_info info)
        {
            sf::Color atomColor = constants::getElementColor(ZIndex) + info.color_addition;

            float vdw_r = radius;

            if (info.spaceFilling)
            {
                vdw_r = constants::VDW_RADII[ZIndex] * 3;
            }

            float minDist = 20.0f;
            float maxDist = 300.0f;
            float alpha = info.opacity * 255.0f;

            if (camDistance > minDist)
            {
                float t = (camDistance - minDist) / (maxDist - minDist);
                t = std::clamp(t, 0.0f, 1.0f);
                alpha = info.opacity * 255.0f * (1.0f - t);
            }

            float rad = vdw_r / camDistance * 10.f;
            if (info.lennardBall)
            {
                sf::CircleShape cloud(rad);
                cloud.setOrigin({rad, rad});
                cloud.setPosition(pos);
                cloud.setFillColor(sf::Color(atomColor.r, atomColor.g, atomColor.b, alpha));
                target.draw(cloud);
            }

            if (info.letter)
            {
                auto &font = window.getFont();
                std::string name = constants::getAtomLetter(ZIndex);
                sf::Text name_text(font, name, 25.f);
                name_text.setScale({0.02f * rad, 0.02f * rad});
                name_text.setFillColor({255,
                                        255,
                                        255,
                                        static_cast<uint8_t>(alpha)});

                sf::FloatRect b = name_text.getLocalBounds();
                name_text.setOrigin(b.getCenter());
                name_text.setPosition(pos);
                target.draw(name_text);

                // Charge
                if (std::abs(q) > 1.f)
                {
                    sf::Text ion_text(font, q > 0.0f ? "+" : "-", 15.f);

                    ion_text.setScale({0.018f * rad, 0.018f * rad});
                    ion_text.setFillColor({255, 255, 255, static_cast<uint8_t>(alpha)});

                    sf::FloatRect ionBounds = ion_text.getLocalBounds();
                    ion_text.setOrigin(ionBounds.getCenter());

                    float offsetScale = rad * 0.5f;
                    ion_text.setPosition(pos + sf::Vector2f(offsetScale, -offsetScale));

                    target.draw(ion_text);
                }
            }
        }
    } // namespace fun
} // namespace sim
