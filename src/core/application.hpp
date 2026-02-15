#pragma once

#include <imgui/imgui.h>
#include <imgui/imgui-SFML.h>

#include "window.hpp"
#include "UIHandler.hpp"
#include "audio.hpp"

#include "json.hpp"

#include <cstdlib>
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

        const char* appdata = std::getenv("APPDATA");
        const std::filesystem::path savePath = std::filesystem::path(appdata) / "Genesis Molecular Dynamics Engine";

        void save();
        void load();
        void initSounds();
    
        application_state current_state{application_state::APP_STATE_MENU};
        options app_options{};
    
        window_t window;
        AudioEngine audio{};
        UIHandler ui;
    };
} // namespace core
