#include "application.hpp"

#include <imgui-SFML.h>
#include <iostream>
#include <fstream>

namespace core
{
    application::application(int32_t height, int32_t width, const std::string name)
    {
        app_options.sim_options.target_fps = 0;
        window = std::make_unique<core::window_t>(width, height, name, 50.f, app_options.sim_options.target_fps);

        if (!ImGui::SFML::Init(window->getWindow(), false))
            throw std::runtime_error("Error! failed to init Imgui");

        ui = std::make_unique<UIHandler>(app_options, *window.get());

        ui->set_language(app_options.lang);
        ui->setApplicationStateCallback([&](application_state newState)
        {
            current_state = newState;
        });
        
        ImGuiIO& io = ImGui::GetIO();
        io.Fonts = new ImFontAtlas();

        io.Fonts->AddFontDefault();
        constexpr float size = 18.5f;

        ImFontGlyphRangesBuilder builder; 

        builder.AddRanges(io.Fonts->GetGlyphRangesDefault());
        ImWchar subscript_range[] = { 0x2080, 0x209F, 0 };
        builder.AddRanges(subscript_range);
        builder.AddText("H₂O CO₂ Na⁺ Cl⁻ → Δ µ ° ± ∞, ™®©℠ sus™ amongAI®");
        builder.AddChar(0x2122);  // ™
        builder.AddChar(0x00AE);  // ®
        builder.AddChar(0x2120);  // ℠
        builder.AddChar(0x00A9);  // © 

        ImVector<ImWchar> ranges;
        builder.BuildRanges(&ranges);

        ImFont* regular = io.Fonts->AddFontFromFileTTF("src/resource/fonts/Orbitron-Regular.ttf", size, nullptr, ranges.Data);
        if (regular == nullptr) throw std::runtime_error("Failed to load Orbitron-Regular.ttf");

        ImFont* medium = io.Fonts->AddFontFromFileTTF("src/resource/fonts/Orbitron-Medium.ttf", size, nullptr, ranges.Data);
        ImFont* bold   = io.Fonts->AddFontFromFileTTF("src/resource/fonts/Orbitron-Bold.ttf", size, nullptr, ranges.Data);
        ImFont* black  = io.Fonts->AddFontFromFileTTF("src/resource/fonts/Orbitron-Black.ttf", size, nullptr, ranges.Data);

        ui->set_regular_font(regular);
        ui->set_medium_font(medium);
        ui->set_bold_font(bold);
        ui->set_black_font(black);

        if (!ImGui::SFML::UpdateFontTexture())
        {
            throw std::exception("[APPLICATION]: IMGUI SFML Couldn't update it's fonts");
        }
    }

    application::~application()
    {
        save();
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
                ui->drawMenu();
            if (current_state == application_state::APP_STATE_SIMULATION)
                ui->drawUniverse(*window.get());
            
            ImGui::SFML::Render(window->getWindow());
            window->display();
        }
    }

    void application::save()
    {
        nlohmann::json app_save{};

        std::filesystem::path path{"src/resource"};

        try
        {
            std::fstream file{path};

            file << app_save.dump();
        }
        catch(const std::exception& e)
        {
            std::cerr << "[APPLICATION]: " << e.what() << '\n';
        }
    }
} // namespace core
