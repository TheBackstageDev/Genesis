#include "application.hpp"

#include <imgui-SFML.h>
#include <implot.h>
#include <iostream>
#include <fstream>
#include <stdexcept>

namespace core
{
    application::application(int32_t height, int32_t width, const std::string name)
        : window(width, height, name, 50.f, 0), ui(app_options, window)
    {
        if (!ImGui::SFML::Init(window.getWindow(), false))
            throw std::runtime_error("Error! failed to init Imgui");

        sf::Image icon;
        if (!icon.loadFromFile("resource/images/icons/window.png"))
        {
            std::cerr << "Error loading icon.png" << std::endl;
        }

        window.getWindow().setIcon(sf::Vector2u(icon.getSize().x, icon.getSize().y), icon.getPixelsPtr());

        ImPlot::CreateContext();

        ui.set_language(app_options.lang);
        ui.setApplicationStateCallback([&](application_state newState)
        {
            current_state = newState;
        });

        ui.setGetApplicationStateCallback([&]()
        {
            return current_state;
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

        ImFont* regular = io.Fonts->AddFontFromFileTTF("resource/fonts/Orbitron-Regular.ttf", size, nullptr, ranges.Data);
        ImFont* regular_small = io.Fonts->AddFontFromFileTTF("resource/fonts/Orbitron-Regular.ttf", size * 0.5f, nullptr, ranges.Data);
        ImFont* regular_big = io.Fonts->AddFontFromFileTTF("resource/fonts/Orbitron-Regular.ttf", size * 2.f, nullptr, ranges.Data);
        if (regular == nullptr) throw std::runtime_error("Failed to load Orbitron-Regular.ttf");

        ImFont* medium = io.Fonts->AddFontFromFileTTF("resource/fonts/Orbitron-Medium.ttf", size, nullptr, ranges.Data);
        ImFont* bold   = io.Fonts->AddFontFromFileTTF("resource/fonts/Orbitron-Bold.ttf", size, nullptr, ranges.Data);
        ImFont* bold_big  = io.Fonts->AddFontFromFileTTF("resource/fonts/Orbitron-Bold.ttf", size * 2.f, nullptr, ranges.Data);
        ImFont* black  = io.Fonts->AddFontFromFileTTF("resource/fonts/Orbitron-Black.ttf", size, nullptr, ranges.Data);

        ui.set_regular_font(regular);
        ui.set_regular_font_small(regular_small);
        ui.set_regular_big_font(regular_big);
        ui.set_medium_font(medium);
        ui.set_bold_font(bold);
        ui.set_bold_big_font(bold_big);
        ui.set_black_font(black);

        if (!ImGui::SFML::UpdateFontTexture())
        {
            throw std::runtime_error("[APPLICATION]: IMGUI SFML Couldn't update it's fonts");
        }
    }

    application::~application()
    {
        save();
        ImPlot::DestroyContext();
        ImGui::SFML::Shutdown(window.getWindow());
    }

    void application::run()
    {
        sf::Clock deltaClock;
        while (window.isOpen())
        {
            window.pollEvents(); 
            
            ui.setDeltaTime(deltaClock.getElapsedTime().asSeconds());
            ImGui::SFML::Update(window.getWindow(), deltaClock.restart());

            window.getWindow().setFramerateLimit(app_options.target_fps);
            window.refresh();
            window.clear();

            auto start = std::chrono::high_resolution_clock::now();

            if (current_state == application_state::APP_STATE_MENU)
                ui.drawMenu();
            if (current_state == application_state::APP_STATE_SIMULATION)
                ui.drawUniverse();

            auto end = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double, std::milli> duration = end - start;

            //std::cout << "UI execution time: " << duration.count() << " milliseconds" << std::endl;
                
            ImGui::SFML::Render(window.getWindow());
            window.display();
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
