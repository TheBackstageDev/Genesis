#include "framebufferUtils.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>
#include <iostream>

namespace fbutils
{
    std::vector<unsigned char> captureFramebuffer(int width, int height)
    {
        std::vector<unsigned char> pixels(width * height * 4);
        glPixelStorei(GL_PACK_ALIGNMENT, 1);
        glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());

        // Flip vertically
        for (int y = 0; y < height / 2; ++y)
        {
            int top = y * width * 4;
            int bottom = (height - 1 - y) * width * 4;
            for (int x = 0; x < width * 4; ++x)
                std::swap(pixels[top + x], pixels[bottom + x]);
        }
        return pixels;
    }

    bool savePNG(const std::filesystem::path& filepath,
                 int width, int height,
                 const std::vector<unsigned char>& pixels)
    {
        if (!stbi_write_png(filepath.string().c_str(), width, height, 4, pixels.data(), width * 4))
        {
            std::cerr << "[FramebufferUtils]: Failed to save PNG " << filepath << std::endl;
            return false;
        }
        return true;
    }

    GLuint createOffscreenFramebuffer(int width, int height, GLuint& colorTex, GLuint& rbo)
    {
        GLuint fbo;
        glGenFramebuffers(1, &fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);

        glGenTextures(1, &colorTex);
        glBindTexture(GL_TEXTURE_2D, colorTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0,
                     GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, colorTex, 0);

        glGenRenderbuffers(1, &rbo);
        glBindRenderbuffer(GL_RENDERBUFFER, rbo);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rbo);

        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        {
            std::cerr << "[FramebufferUtils]: Framebuffer incomplete!" << std::endl;
            return 0;
        }

        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        return fbo;
    }
}
