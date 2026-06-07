#pragma once
#include <string>
#include <filesystem>
#include <vector>
#include <glad/glad.h>

namespace fbutils
{
    std::vector<unsigned char> captureFramebuffer(int width, int height);

    bool savePNG(const std::filesystem::path& filepath,
                 int width, int height,
                 const std::vector<unsigned char>& pixels);

    GLuint createOffscreenFramebuffer(int width, int height, GLuint& colorTex, GLuint& rbo);
}
