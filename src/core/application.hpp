#pragma once

#include <imgui/imgui.h>
#include <imgui/imgui-SFML.h>

#include "core/window.hpp"
#include "core/UIHandler.hpp"

#include "json.hpp"

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

        void save();
    
        application_state current_state{application_state::APP_STATE_MENU};
        options app_options{};
    
        std::unique_ptr<core::window_t> window;
        std::unique_ptr<UIHandler> ui;
    };
} // namespace core
