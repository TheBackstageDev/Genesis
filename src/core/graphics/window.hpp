#pragma once

#include <glad.h>
#include <GLFW/glfw3.h>

#include <vector>

namespace core
{   
    struct extent2D 
    {
        uint32_t width;
        uint32_t height;
    };

    class window_t
    {
    public:
        window_t(int32_t width, int32_t height, const char* title);
        ~window_t();

        window_t(const window_t &) = delete;
        window_t &operator=(const window_t &) = delete;

        const char* get_title() { return title; }

        GLFWwindow* getWindow() { return p_window; }

        extent2D extent() const 
        {
            int32_t width, height;

            glfwGetWindowSize(p_window, &width, &height);

            return extent2D{static_cast<uint32_t>(width), static_cast<uint32_t>(height)}; 
        }

        bool should_close() { return glfwWindowShouldClose(p_window); }
        bool resized() { return _resized; }

        void resetResizedFlag() { _resized = false; }
        void processInput();

        operator GLFWwindow *() const { return p_window; }
    private:
        static void frameBufferResizeCallback(GLFWwindow *window, int width, int height);

        bool _resized = false;
        
        const char* title = nullptr;
        GLFWwindow* p_window = nullptr;
    };
} // namespace core
