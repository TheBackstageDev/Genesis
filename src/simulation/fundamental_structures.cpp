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

    } // namespace fun
} // namespace sim
