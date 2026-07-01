#pragma once

#include <imgui.h>

#include "graphics/window.hpp"
#include "simulation/core/simulation_inspector.hpp"
#include "simulation/physics/fields/parameter_table.hpp"
#include "ui/UIHandler.hpp"
#include "core/audio/audio.hpp"

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
        void initIMGUI();

        bool wasAnyItemHoveredLastFrame = false;
        bool wasMouseActiveLastFrame = false;

        const char* appdata = std::getenv("APPDATA");
        const std::filesystem::path savePath = std::filesystem::path(appdata) / "Genesis Molecular Dynamics Engine";

        void save();
        void load();
        void initSounds();

        void runUniverse();

        std::thread simulationThread;
        std::atomic<bool> simulationRunning{false};
        std::mutex simMutex;
    
        application_state current_state{application_state::APP_STATE_MENU};
        options app_options{};
    
        window_t window;
        AudioEngine audio{};
        UIHandler ui;
        
        sim::parameter_table m_parameterTable;
        sim::simulation_inspector m_siminspector{};
    };
} // namespace core
