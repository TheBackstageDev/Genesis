#include "window.hpp"
#include <iostream>

namespace core
{
    window_t::window_t(int32_t width, int32_t height, const char* title)
    : title(title)
    {
        if (!glfwInit())
        {
            std::cerr << "Failed to initialize GLFW" << std::endl;
            exit(EXIT_FAILURE);
        }

        glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

        p_window = glfwCreateWindow(width, height, title, nullptr, nullptr);
        if (!p_window)
        {
            std::cerr << "Failed to create GLFW window" << std::endl;
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        glfwMakeContextCurrent(p_window);
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) 
            throw std::runtime_error("Failed to initialize GLAD");

        glfwSetWindowUserPointer(p_window, this);
        glfwSetFramebufferSizeCallback(p_window, frameBufferResizeCallback);
    }

    window_t::~window_t()
    {
        glfwDestroyWindow(p_window);
        glfwTerminate();
    }

    void window_t::processInput()
    {

    }

    void window_t::frameBufferResizeCallback(GLFWwindow *pWindow, int width, int height)
    {
        auto win = reinterpret_cast<window_t*>(glfwGetWindowUserPointer(pWindow));
        if (win) win->_resized = true;
        glViewport(0, 0, width, height); 
    }
} // namespace core
