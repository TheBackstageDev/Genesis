#include "application.hpp"

#include <imgui-SFML.h>

namespace core
{
    application::application(int32_t height, int32_t width, const std::string name)
    {
        window = std::make_unique<core::window_t>(width, height, name, 50.f, app_options.sim_options.target_fps);

        if (!ImGui::SFML::Init(window->getWindow(), false))
            throw std::runtime_error("Error! failed to init Imgui");

        ui.set_language(app_options.lang);
        ui.setApplicationStateCallback([&](application_state newState)
        {
            current_state = newState;
        });
        
        ImGuiIO& io = ImGui::GetIO();
        io.Fonts = new ImFontAtlas();

        io.Fonts->AddFontDefault();
        constexpr float size = 18.5f;

        ImFont* regular = io.Fonts->AddFontFromFileTTF("src/resource/fonts/Orbitron-Regular.ttf", size);
        if (regular == nullptr) throw std::runtime_error("Failed to load Orbitron-Regular.ttf");

        ImFont* medium = io.Fonts->AddFontFromFileTTF("src/resource/fonts/Orbitron-Medium.ttf", size);
        ImFont* bold   = io.Fonts->AddFontFromFileTTF("src/resource/fonts/Orbitron-Bold.ttf", size);
        ImFont* black  = io.Fonts->AddFontFromFileTTF("src/resource/fonts/Orbitron-Black.ttf", size);

        ui.set_regular_font(regular);
        ui.set_medium_font(medium);
        ui.set_bold_font(bold);
        ui.set_black_font(black);

        ImGui::SFML::UpdateFontTexture();
    }

    application::~application()
    {
        ImGui::SFML::Shutdown(window->getWindow());
    }

    void application::run()
    {
        sf::Clock deltaClock;
        while (window->isOpen())
        {
            window->pollEvents(); 
            ImGui::SFML::Update(window->getWindow(), deltaClock.restart());

            window->getWindow().setFramerateLimit(app_options.target_fps);
            window->refresh();
            window->clear();
            
            if (current_state == application_state::APP_STATE_MENU)
                ui.drawMenu();
            if (current_state == application_state::APP_STATE_SIMULATION)
                ui.drawUniverse(*window.get());
            
            ImGui::SFML::Render(window->getWindow());
            window->display();
        }
    }
} // namespace core
