#pragma once

#include <glad.h>
#include <filesystem>

namespace core
{
    GLuint loadTexture(std::filesystem::path filename);
} // namespace core
