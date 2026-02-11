#pragma once

#include <imgui/imgui.h>
#include <imgui/imgui-SFML.h>

#include "core/window.hpp"
#include "core/UIHandler.hpp"
#include "core/audio.hpp"

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
        bool wasAnyItemHoveredLastFrame = false;
        bool wasMouseActiveLastFrame = false;

        void save();
        void initSounds();
    
        application_state current_state{application_state::APP_STATE_MENU};
        options app_options{};
    
        window_t window;
        AudioEngine audio{};
        UIHandler ui;
    };
} // namespace core
