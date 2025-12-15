#include "UIHandler.hpp"

#include <fstream>
#include <iostream>

namespace core
{
    UIHandler::UIHandler(options& app_options)
        : app_options(app_options)
    {
        write_localization_json(lang);
    }

    // Helper

    std::filesystem::path getLocalizationFile(localization lang)
    {
        std::filesystem::path localization_path = "src/resource/localizations";
        switch (lang)
        {
            case localization::EN_US: return localization_path / "enus.json";
            case localization::PT_BR: return localization_path / "ptbr.json";
            default: return localization_path / "enus.json";
        }
    }

    void UIHandler::write_localization_json(localization lang)
    {
        std::filesystem::path path = getLocalizationFile(lang);
        std::ifstream file(path);
        if (!file.is_open())
        {
            std::cerr << "[UI Handler] Cannot open: " << path << '\n';
            return;
        }

        try
        {
            file >> localization_json;
        }
        catch (const nlohmann::json::parse_error& e)
        {
            std::cerr << "[UI Handler] JSON parse error: " << e.what() << '\n';
            return;
        }
    }

    // Menu

    void UIHandler::drawMenu()
    {
        ImGuiIO& io = ImGui::GetIO();
        float padding = 50.0f;
        float buttonWidth = 260.0f;
        float buttonHeight = 45.0f;

        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.05f, 0.05f, 0.12f, 0.95f));
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 1.0f, 1.0f, 1.0f));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(15, 15));
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 20));

        {
            std::string menu_title = localization_json["Menu"]["title"].get<std::string>().c_str();
            ImVec2 titleSize = ImGui::CalcTextSize(menu_title.c_str());
                    
            float windowWidth = ImGui::GetMainViewport()->GetCenter().x - titleSize.x;
            ImVec2 titlePos(windowWidth, padding);
            ImGui::SetNextWindowPos(titlePos, ImGuiCond_Always);
            ImGui::Begin("TitleWindow", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar);

            ImGui::SetWindowFontScale(3.f);
            ImGui::SetCursorPosX((ImGui::GetWindowWidth() - titleSize.x) * 0.25f);
            ImGui::Text(menu_title.c_str());
            ImGui::SetWindowFontScale(1.2f);
            std::string menu_subtitle = localization_json["Menu"]["subtitle"].get<std::string>().c_str();

            ImGui::SetCursorPosX((ImGui::GetWindowWidth() - ImGui::CalcTextSize(menu_subtitle.c_str()).x) * 0.5f);
            ImGui::Text(menu_subtitle.c_str());

            ImGui::End();
        }

        {
            ImVec2 panelPos(padding, io.DisplaySize.y - padding);
            ImGui::SetNextWindowPos(panelPos, ImGuiCond_Always, ImVec2(0.0f, 1.0f));
            ImGui::SetNextWindowSize(ImVec2(buttonWidth + 40, 0));
            ImGui::Begin("NavigationPanel", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);

            if (ImGui::Button(localization_json["Menu"]["button_tutorials"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight))) tutorialSelectionOpen = !tutorialSelectionOpen;
            if (ImGui::Button(localization_json["Menu"]["button_scenarios"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight))) sceneSelectionOpen = !sceneSelectionOpen;
            if (ImGui::Button(localization_json["Menu"]["button_sandbox"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight))) sandboxSelectionOpen = !sandboxSelectionOpen;
            if (ImGui::Button(localization_json["Menu"]["button_challenges"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight))) challengeSelectionOpen = !challengeSelectionOpen;
            if (ImGui::Button(localization_json["Menu"]["button_achievements"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight))) challengeViewOpen = !challengeViewOpen;
            if (ImGui::Button(localization_json["Menu"]["button_options"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight))) optionsOpen = !optionsOpen;
            if (ImGui::Button(localization_json["Menu"]["button_quit"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight))) { std::exit(EXIT_SUCCESS); }

            ImGui::End();
        }

        {
            ImVec2 infoPos(io.DisplaySize.x - padding, io.DisplaySize.y - padding);
            ImGui::SetNextWindowPos(infoPos, ImGuiCond_Always, ImVec2(1.0f, 1.0f));
            ImGui::Begin("InfoPanel", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize);
            ImGui::Text(localization_json["Menu"]["info_version"].get<std::string>().c_str());
            ImGui::Text(localization_json["Menu"]["info_author"].get<std::string>().c_str());
            ImGui::End();
        }
        
        if (sceneSelectionOpen) drawSceneSelection();
        if (sandboxSelectionOpen) drawSandboxCreation();
        if (optionsOpen) drawOptions();
        
        ImGui::PopStyleVar(2);
        ImGui::PopStyleColor(2);
    }

    void UIHandler::drawBackgroundDisplay()
    {
        
    }
    
    sim::fun::universe_create_info sandbox_info{};

    void UIHandler::drawSandboxCreation()
    {
        ImGui::SetNextWindowSize(ImVec2(700, 550), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Always, ImVec2(0.5f, 0.5f));

        auto& sandbox_creation = localization_json["Menu"]["Sandbox_Creation"];

        ImGui::Begin(sandbox_creation["title"].get<std::string>().c_str(), nullptr, ImGuiWindowFlags_NoMove);
        if (ImGui::CollapsingHeader(sandbox_creation["header_physics"].get<std::string>().c_str(), ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Checkbox(sandbox_creation["gravity_enabled"].get<std::string>().c_str(), &sandbox_info.has_gravity);
            if (sandbox_info.has_gravity)
            {
                ImGui::Indent();
                ImGui::DragFloat(sandbox_creation["gravity_magnitude"].get<std::string>().c_str(), 
                                &sandbox_info.mag_gravity, 0.1f, 0.0f, 50.0f, "%.2f m/s²");
                ImGui::Unindent();
            }

            ImGui::Checkbox(sandbox_creation["wall_collision"].get<std::string>().c_str(), &sandbox_info.wall_collision);

            ImGui::BeginDisabled();
            ImGui::Checkbox(sandbox_creation["isothermal"].get<std::string>().c_str(), &sandbox_info.isothermal);
            ImGui::Checkbox(sandbox_creation["reactive"].get<std::string>().c_str(), &sandbox_info.reactive);
            ImGui::EndDisabled();
        }

        if (ImGui::CollapsingHeader(sandbox_creation["header_box"].get<std::string>().c_str(), ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::Text(sandbox_creation["box_dimensions"].get<std::string>().c_str());
            ImGui::SliderFloat("X##box", &sandbox_info.box.x, 1.0f, 100.0f, "%.2f");
            ImGui::SliderFloat("Y##box", &sandbox_info.box.y, 1.0f, 100.0f, "%.2f");
            ImGui::SliderFloat("Z##box", &sandbox_info.box.z, 1.0f, 100.0f, "%.2f");

            if (ImGui::Button(sandbox_creation["box_cubic"].get<std::string>().c_str()))
            {
                float size = (sandbox_info.box.x + sandbox_info.box.y + sandbox_info.box.z) / 3.0f;
                sandbox_info.box = sf::Vector3f(size, size, size);
            }
        }

        if (ImGui::CollapsingHeader(sandbox_creation["header_visual"].get<std::string>().c_str()))
        {
            ImGui::Checkbox(sandbox_creation["render_water"].get<std::string>().c_str(), &sandbox_info.render_water);
        }

        if (ImGui::CollapsingHeader(sandbox_creation["header_logging"].get<std::string>().c_str()))
        {
            ImGui::BeginDisabled();
            ImGui::Checkbox(sandbox_creation["log_reactions"].get<std::string>().c_str(), &sandbox_info.log_flags.log_reactions);
            ImGui::EndDisabled();
        }

        ImGui::Separator();
        ImGui::Dummy(ImVec2(0.0f, 20.0f));

        if (ImGui::Button(sandbox_creation["button_create"].get<std::string>().c_str(), ImVec2(200, 50)))
        {
            simulation_universe = std::make_unique<sim::fun::universe>(sandbox_info);
            sandboxSelectionOpen = false;

            setState(application_state::APP_STATE_SIMULATION);
        }

        if (ImGui::Button(sandbox_creation["button_cancel"].get<std::string>().c_str(), ImVec2(200, 50)))
        {
            sandboxSelectionOpen = false;
        }

        ImGui::End();
    }

    void UIHandler::drawSceneFrame(scenario_info& info, int32_t id)
    {
        ImGui::BeginChild(ImGui::GetID(id), ImVec2());
        ImGui::EndChild();
    }

    void UIHandler::drawSceneSelection()
    {
        ImGui::Begin(localization_json["Menu"]["Scene_Selection"]["title"].get<std::string>().c_str(), nullptr, ImGuiWindowFlags_NoMove);

        ImGui::End();
    }

    void UIHandler::drawOptions()
    {
        auto& options = localization_json["Menu"]["Options"];

        ImVec2 viewport_size = ImGui::GetMainViewport()->Size;
        
        ImGui::SetNextWindowSizeConstraints(ImVec2(500, 300), // Min
        ImVec2(800, 600));  // Max
        ImGui::SetNextWindowSize(ImVec2(0, 0), ImGuiCond_Always);
        ImGui::SetNextWindowPos(ImVec2(viewport_size.x / 2.f - 200, viewport_size.y / 3.f), ImGuiCond_Appearing);
        
        ImGui::Begin("Options", &optionsOpen, ImGuiWindowFlags_NoTitleBar);
        ImGui::BeginTabBar(options["title"].get<std::string>().c_str(), ImGuiTabBarFlags_NoCloseWithMiddleMouseButton);
        if (ImGui::BeginTabItem(options["tab_graphics"].get<std::string>().c_str()))
        {
            ImGui::Text(options["tab_graphics"].get<std::string>().c_str());

            static std::string ball_and_stick = "";
            static std::string letter_and_stick = "";
            static std::string space_filling = "";

            ball_and_stick = options["render_modes"]["ball_and_stick"].get<std::string>();
            letter_and_stick = options["render_modes"]["letter_and_stick"].get<std::string>();
            space_filling = options["render_modes"]["space_filling"].get<std::string>();

            const char* modes[] = 
            {
                ball_and_stick.c_str(),
                letter_and_stick.c_str(),
                space_filling.c_str()
            };

            static int32_t current_mode = static_cast<int32_t>(app_options.sim_options.render_mode);
            if (ImGui::Combo("##render_mode", &current_mode, modes, IM_ARRAYSIZE(modes)))
                app_options.sim_options.render_mode = static_cast<simulation_render_mode>(current_mode);

            ImGui::SliderInt(options["tab_target_fps"].get<std::string>().c_str(), &app_options.sim_options.target_fps, 24, 144);

            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem(options["tab_audio"].get<std::string>().c_str()))
        {
            ImGui::SliderFloat(options["tab_audio_master"].get<std::string>().c_str(), &app_options.master_volume, 0.0f, 1.0f, "%.1f");
            ImGui::Checkbox(options["tab_audio_effects"].get<std::string>().c_str(), &app_options.sound_effects);
            ImGui::Checkbox(options["tab_audio_background"].get<std::string>().c_str(), &app_options.background_music);

            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem(options["tab_language"].get<std::string>().c_str()))
        {
            static const char* languages[] = {"English (US)", "Português (BR)"};
            static int32_t selected_language = static_cast<int32_t>(app_options.lang);

            if (ImGui::Combo("##language", &selected_language, languages, IM_ARRAYSIZE(languages)))
            {
                localization new_lang = static_cast<localization>(selected_language);
                
                if (new_lang != app_options.lang)
                {
                    app_options.lang = new_lang;
                    write_localization_json(app_options.lang);

                    ImGui::EndTabItem();
                    ImGui::EndTabBar();
                    ImGui::End();

                    return;
                }
            }
                
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();

        ImGui::Dummy(ImVec2(0.0f, 20.0f));
        if (ImGui::Button(options["tab_exit"].get<std::string>().c_str())) { optionsOpen = false; }

        ImGui::End();
    }

    // Universe Simulation

    void UIHandler::drawUniverseUI()
    {
    
    }

    void UIHandler::drawUniverse(window_t& window)
    {
        simulation_render_mode mode = app_options.sim_options.render_mode;

        bool letter = false;
        bool ball = false;

        if (mode == simulation_render_mode::BALL_AND_STICK)
        {
            ball = true;
            letter = true;
        } 
        else if (mode == simulation_render_mode::LETTER_AND_STICK)
        {
            ball = false;
            letter = true;
        }

        simulation_universe->handleCamera();
        simulation_universe->draw(window, letter, ball);

        if (!pauseMenuOpen)
            drawUniverseUI();
        else
            pauseMenu();

        if (ImGui::IsKeyPressed(ImGuiKey_Escape))
            pauseMenuOpen = !pauseMenuOpen;
    }

    void UIHandler::pauseMenu()
    {
        auto& pause_menu = localization_json["Simulation"]["pause_menu"];

        ImGui::SetNextWindowSize(ImVec2(300, 400), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSizeConstraints(ImVec2(350, 450), ImVec2(450, 600));
        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Always, ImVec2(0.5f, 0.5f));
        
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.05f, 0.05f, 0.12f, 0.95f));
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 1.0f, 1.0f, 1.0f));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(30, 30));
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 20));

        ImGui::Begin("pause menu", &pauseMenuOpen, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize);

        ImGui::SetWindowFontScale(2.8f);
        const std::string title = pause_menu["title"].get<std::string>();
        ImVec2 titleSize = ImGui::CalcTextSize(title.c_str());
        ImGui::SetCursorPosX((ImGui::GetWindowWidth() - titleSize.x) * 0.5f);
        ImGui::Text("%s", title.c_str());

        ImGui::SetWindowFontScale(1.0f);

        ImGui::Dummy(ImVec2(0.0f, 40.0f));

        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 8.0f);
        float buttonWidth = 270.0f;
        ImVec2 buttonSize(buttonWidth, 50.0f);

        if (ImGui::Button(pause_menu["button_resume"].get<std::string>().c_str(), buttonSize))
        {
            pauseMenuOpen = false;
        }

        if (ImGui::Button(pause_menu["button_restart"].get<std::string>().c_str(), buttonSize))
        {
            simulation_universe = std::make_unique<sim::fun::universe>(sandbox_info);
            pauseMenuOpen = false;
        }

        if (ImGui::Button(pause_menu["button_options"].get<std::string>().c_str(), buttonSize))
        {
            optionsOpen = !optionsOpen;
        }

        if (ImGui::Button(pause_menu["button_exit_menu"].get<std::string>().c_str(), buttonSize))
        {
            ImGui::OpenPopup("ConfirmExitToMenu");
        }

        if (ImGui::Button(pause_menu["button_exit_desktop"].get<std::string>().c_str(), buttonSize))
        {
            ImGui::OpenPopup("ConfirmExitToMenu");
            exitDesktop = true;
        }

        if (ImGui::BeginPopupModal("ConfirmExitToMenu", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar))
        {
            float button_width = 130.0f;
            
            ImVec2 center = ImGui::GetMainViewport()->GetCenter();
            ImGui::SetWindowPos({center.x - button_width * 1.5f, center.y - button_width});
            ImGui::TextWrapped(pause_menu["warning"].get<std::string>().c_str());

            ImGui::Dummy(ImVec2(0.0f, 20.0f));

            ImGui::SetCursorPosX((ImGui::GetWindowWidth() - (button_width * 2 + ImGui::GetStyle().ItemSpacing.x)) * 0.5f);

            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.6f, 0.1f, 0.f, 1.f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.9f, 0.1f, 0.f, 1.f));
            if (ImGui::Button(pause_menu["true_exit"].get<std::string>().c_str(), ImVec2(button_width, 30)))
            {
                pauseMenuOpen = false;
                setState(application_state::APP_STATE_MENU);
                
                if (exitDesktop)
                    std::exit(EXIT_SUCCESS);
                ImGui::CloseCurrentPopup();
            }
            ImGui::PopStyleColor(2);

            ImGui::SameLine();

            if (ImGui::Button(pause_menu["cancel"].get<std::string>().c_str(), ImVec2(button_width, 30)))
            {
                exitDesktop = false;
                ImGui::CloseCurrentPopup();
            }

            ImGui::EndPopup();
        }
        
        ImGui::PopStyleVar();

        if (optionsOpen)
        {
            drawOptions();
        }

        ImGui::Dummy(ImVec2(0.0f, 20.0f));
        ImGui::TextWrapped(pause_menu["hint_esc"].get<std::string>().c_str());

        ImGui::End();

        ImGui::PopStyleVar(2);
        ImGui::PopStyleColor(2);
    }

    // Video 

    void UIHandler::drawVideoControls()
    {
        ImGui::SetNextWindowPos(ImVec2(100, 100), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(400, 220), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowBgAlpha(0.85f);

        if (!ImGui::Begin("Video Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
        {
            ImGui::End();
            return;
        }
    }
} // namespace core
