#pragma once

#include <imgui/imgui.h>
#include <imgui/imgui-SFML.h>

#include "core/window.hpp"

#include <memory>

namespace core
{
    class application
    {
    public:
        application(int32_t height, int32_t width, const std::string name);
        ~application();

        void run();
    private:
        std::unique_ptr<core::window_t> window;
    };
} // namespace core
