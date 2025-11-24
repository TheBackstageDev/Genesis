#pragma once

#include "constants.hpp"
#include "fundamental_structures.hpp"
#include <filesystem>

namespace sim::io
{
    bool loadXYZ(const std::filesystem::path path, std::vector<fun::def_atom>& atoms, std::vector<fun::def_bond>& bonds, std::vector<sf::Vector3f>& positions);
} // namespace sim::io
