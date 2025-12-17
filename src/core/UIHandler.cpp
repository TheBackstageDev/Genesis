#include "UIHandler.hpp"

#include <fstream>
#include <iostream>

#include "simulation/smiles_parser.hpp"
#include "simulation/format_loader.hpp"

namespace core
{
    UIHandler::UIHandler(options &app_options, window_t& window)
        : app_options(app_options)
    {
        write_localization_json(lang);
        initCompoundPresets();

        initImages(window);
    }

    void UIHandler::initImages(window_t& window)
    {
        // PlaceHolder

        sf::Image placeholder_image;
        placeholder_image.createMaskFromColor(sf::Color(80, 80, 90), 255); // Dark gray
        placeholder_texture_id = placeholder_texture.getNativeHandle();

        initCompoundPresetsImages(window);
    }

    void UIHandler::initCompoundPresetsImages(window_t& window)
    {
        sim::fun::universe_create_info display_info{};
        display_info.box.x = 100.f;
        display_info.box.y = 100.f;
        display_info.box.z = 100.f;
        display_info.wall_collision = false;
        
        display_universe = std::make_unique<sim::fun::universe>(display_info);
        display_universe->camera().target = {0.f, 0.f, 0.f};
        display_universe->camera().distance = 0.f;
        display_universe->camera().azimuth = 200.f;
        display_universe->camera().elevation = 25.f;
        display_universe->camera().fov = 60.f;
        display_universe->camera().nearPlane = 0.1f;
        display_universe->camera().farPlane = 500.f;

        for (auto &compound : compound_presets)
        {
            sf::RenderTexture thumbnail_renderer;
            display_universe->clear();
            display_universe->createMolecule(compound.structure, {0.f, 0.f, 0.f});

            display_universe->camera().distance = 1.f * compound.structure.atoms.size();

            window.clear();
            display_universe->draw(window, window.getWindow(), true, true, true, false);

            sf::Texture thumb_texture;
            if (!thumb_texture.resize(window.getWindow().getSize()))
                continue;

            thumb_texture.update(window.getWindow());
            thumb_textures.emplace_back(std::move(thumb_texture));
        }
    }

    void UIHandler::initCompoundPresets()
    {
        auto &compounds = localization_json["Compounds"];
        auto &compound_names = compounds["compound_names"];
        auto &compound_descriptions = compounds["compound_descriptions"];

        compound_presets.clear();

        constexpr int32_t num_small_molecules = 8;
        std::array<std::string, num_small_molecules> smilesToParse{"O", "O=C=O", "C", "CC", "CCO", "c1ccccc1", "N#N", "C12=C3C4=C5C6=C3C7=C1C8=C9C2=C4C1=C5C2=C6C7=C8C2=C91"};
        std::array<std::string, num_small_molecules> formulas{"H2O", "CO2", "CH4", "C2H6", "C2H5OH", "C6H6", "N2", "C20"};

        using type = compound_type;
        std::array<type, num_small_molecules> types{type::INORGANIC, type::INORGANIC, type::ORGANIC, type::ORGANIC, type::ORGANIC, type::ORGANIC, type::INORGANIC, type::NANOMATERIAL};

        for (int32_t i = 0; i < num_small_molecules; ++i)
        {
            compound_preset_info nInfo{};
            nInfo.name = compound_names[i];
            nInfo.id = i;
            nInfo.type = types[i];
            nInfo.formula = formulas[i];
            nInfo.SMILES = smilesToParse[i];
            nInfo.structure = sim::parseSMILES(nInfo.SMILES);

            for (int32_t i = 0; i < nInfo.structure.atoms.size(); ++i)
                nInfo.molecular_weight += nInfo.structure.atoms[i].NIndex * MASS_NEUTRON + nInfo.structure.atoms[i].ZIndex * MASS_PROTON;

            compound_presets.emplace_back(std::move(nInfo));
        }
    }

    // Helper

    std::filesystem::path getLocalizationFile(localization lang)
    {
        std::filesystem::path localization_path = "src/resource/localizations";
        switch (lang)
        {
        case localization::EN_US:
            return localization_path / "enus.json";
        case localization::PT_BR:
            return localization_path / "ptbr.json";
        default:
            return localization_path / "enus.json";
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
        catch (const nlohmann::json::parse_error &e)
        {
            std::cerr << "[UI Handler] JSON parse error: " << e.what() << '\n';
            return;
        }
    }

    // Menu

    void UIHandler::drawMenu()
    {
        ImGuiIO &io = ImGui::GetIO();
        float padding = 100.0f;
        float buttonWidth = 260.0f;
        float buttonHeight = 45.0f;

        ImGui::PushFont(regular);
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.05f, 0.05f, 0.12f, 0.95f));
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 1.0f, 1.0f, 1.0f));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(15, 15));
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 20));

        {
            std::string menu_title = localization_json["Menu"]["title"].get<std::string>().c_str();
            float titleSize = ImGui::CalcTextSize(menu_title.c_str()).x;

            float windowWidth = ImGui::GetMainViewport()->GetCenter().x - titleSize;
            ImVec2 titlePos(windowWidth, padding);
            ImGui::SetNextWindowPos(titlePos, ImGuiCond_Always);
            ImGui::Begin("TitleWindow", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar);

            ImGui::SetWindowFontScale(2.f);
            ImGui::SetCursorPosX((ImGui::GetWindowWidth() - titleSize) * 0.3f);
            ImGui::Text(menu_title.c_str());
            ImGui::SetWindowFontScale(1.f);

            ImGui::SetCursorPosX((ImGui::GetWindowWidth()));
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

            if (ImGui::Button(localization_json["Menu"]["button_tutorials"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
                tutorialSelectionOpen = !tutorialSelectionOpen;
            if (ImGui::Button(localization_json["Menu"]["button_scenarios"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
                sceneSelectionOpen = !sceneSelectionOpen;
            if (ImGui::Button(localization_json["Menu"]["button_sandbox"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
                sandboxSelectionOpen = !sandboxSelectionOpen;
            if (ImGui::Button(localization_json["Menu"]["button_challenges"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
                challengeSelectionOpen = !challengeSelectionOpen;
            if (ImGui::Button(localization_json["Menu"]["button_achievements"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
                challengeViewOpen = !challengeViewOpen;
            if (ImGui::Button(localization_json["Menu"]["button_options"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
                optionsOpen = !optionsOpen;
            if (ImGui::Button(localization_json["Menu"]["button_quit"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
            {
                std::exit(EXIT_SUCCESS);
            }

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

        if (sceneSelectionOpen)
            drawSceneSelection();
        if (sandboxSelectionOpen)
            drawSandboxCreation();
        if (optionsOpen)
            drawOptions();

        ImGui::PopStyleVar(2);
        ImGui::PopStyleColor(2);
        ImGui::PopFont();
    }

    void UIHandler::drawBackgroundDisplay()
    {
    }

    sim::fun::universe_create_info sandbox_info{};

    void UIHandler::drawSandboxCreation()
    {
        ImGui::SetNextWindowSize(ImVec2(700, 550), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Always, ImVec2(0.5f, 0.5f));

        auto &sandbox_creation = localization_json["Menu"]["Sandbox_Creation"];

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
            ImGui::SliderFloat("X##box", &sandbox_info.box.x, 1.0f, 200.0f, "%.2f");
            ImGui::SliderFloat("Y##box", &sandbox_info.box.y, 1.0f, 200.0f, "%.2f");
            ImGui::SliderFloat("Z##box", &sandbox_info.box.z, 1.0f, 200.0f, "%.2f");

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

        if (ImGui::Button(sandbox_creation["button_create"].get<std::string>().c_str(), ImVec2(300, 50)))
        {
            simulation_universe = std::make_unique<sim::fun::universe>(sandbox_info);
            simulation_universe->camera().target = {sandbox_info.box.x / 2.f, sandbox_info.box.y / 2.f, sandbox_info.box.z / 2.f};

            sandboxSelectionOpen = false;
            paused_simulation = false;

            setState(application_state::APP_STATE_SIMULATION);

            /* sim::fun::molecule_structure structure{};
            sim::io::loadXYZ("src/resource/molecules/presets/dna.xyz", structure.atoms, structure.bonds, structure.positions);
            sim::organizeSubsets(structure.subsets, structure.atoms, structure.bonds);

            simulation_universe->createMolecule(structure, {30, 30, 30}); */
        }

        if (ImGui::Button(sandbox_creation["button_cancel"].get<std::string>().c_str(), ImVec2(200, 50)))
        {
            sandboxSelectionOpen = false;
        }

        ImGui::End();
    }

    void UIHandler::drawSceneFrame(scenario_info &info, int32_t id)
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
        auto &options = localization_json["Menu"]["Options"];

        ImVec2 viewport_size = ImGui::GetMainViewport()->Size;

        ImGui::SetNextWindowSizeConstraints(ImVec2(500, 300),  // Min
                                            ImVec2(800, 600)); // Max
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

            const char *modes[] =
                {
                    ball_and_stick.c_str(),
                    letter_and_stick.c_str(),
                    space_filling.c_str()};

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
            static const char *languages[] = {"English (US)", "Português (BR)"};
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

                    compound_presets.clear();
                    initCompoundPresets();

                    return;
                }
            }

            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();

        ImGui::Dummy(ImVec2(0.0f, 20.0f));
        if (ImGui::Button(options["tab_exit"].get<std::string>().c_str()))
        {
            optionsOpen = false;
        }

        ImGui::End();
    }

    // Universe Simulation

    void UIHandler::drawUniverseUI()
    {
        auto &sim_ui = localization_json["Simulation"]["universe_ui"];

        ImGuiIO &io = ImGui::GetIO();
        const float panel_height = 180.0f;

        ImGui::PushFont(regular);
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.08f, 0.08f, 0.15f, 0.8f));
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.20f, 0.25f, 0.35f, 0.60f));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(25, 20));
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(20, 12));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 12.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 1.0f);

        ImVec2 panel_pos(0, io.DisplaySize.y);
        ImGui::SetNextWindowPos(panel_pos, ImGuiCond_Always, ImVec2(0.0f, 1.0f));
        ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x, panel_height), ImGuiCond_Always);

        ImGuiWindowFlags hud_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                                     ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar |
                                     ImGuiWindowFlags_NoBringToFrontOnFocus;

        ImGui::Begin("HUD", &showHUD, hud_flags);
        ImGui::Columns(2, "HUDColumns", false);
        ImGui::SetColumnWidth(0, io.DisplaySize.x * 0.4f);

        ImGui::TextColored(ImVec4(0.7f, 0.8f, 1.0f, 1.0f), "%s", sim_ui["stats_title"].get<std::string>().c_str());
        ImGui::Separator();
        ImGui::Text("%s %.2f K", sim_ui["temperature"].get<std::string>().c_str(), simulation_universe->temperature());
        ImGui::Text("%s %.2f bar", sim_ui["pressure"].get<std::string>().c_str(), simulation_universe->pressure());
        ImGui::Text("%s %.2f ps", sim_ui["time"].get<std::string>().c_str(), simulation_universe->timestep() * DT);

        ImGui::NextColumn();

        float button_width = 120.0f;
        float button_height = 35.0f;

        std::string timeControlName = paused_simulation ? "Resume" : "Pause";

        if (ImGui::Button(timeControlName.c_str(), ImVec2(button_width, button_height)))
        {
            paused_simulation = !paused_simulation;
        }
        if (ImGui::IsItemHovered())
        {
            std::string tooltip = paused_simulation ? sim_ui["tooltip_resume"] : sim_ui["tooltip_pause"];
            ImGui::SetTooltip(tooltip.c_str());
        }
        ImGui::SameLine();
        if (ImGui::Button("Screenshot", ImVec2(button_width, button_height)))
        {
            // TODO: Implement screenshot capture
        }
        if (ImGui::IsItemHovered())
        {
            ImGui::SetTooltip(sim_ui["tooltip_screenshot"].get<std::string>().c_str());
        }

        ImGui::SameLine();
        if (ImGui::Button("Add Compound", ImVec2(button_width, button_height)))
        {
            compoundSelector = !compoundSelector;
        }
        if (ImGui::IsItemHovered())
        {
            ImGui::SetTooltip(sim_ui["tooltip_compound_selector"].get<std::string>().c_str());
        }

        ImGui::Columns(1);
        ImGui::End();

        ImGui::PopStyleVar(4);
        ImGui::PopStyleColor(2);

        if (compoundSelector)
            drawCompoundSelector();

        ImGui::PopFont();
    }

    constexpr float padding = 20.0f;
    constexpr float row_height = 200.f;
    constexpr float row_width = 200.f;

    void UIHandler::drawCompoundView(const compound_preset_info &compound)
    {
        auto &comp_sel = localization_json["Simulation"]["universe_ui"]["compound_selector"];
        auto &compounds = localization_json["Compounds"];
            
        const ImVec2 card_size(-FLT_MIN, 0.0f);

        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(10, 8));
        
        float avail_width = ImGui::GetContentRegionAvail().x;
        float offset_x = (avail_width - image_size) * 0.5f;
        if (offset_x > 0) ImGui::SetCursorPosX(ImGui::GetCursorPosX() + offset_x);

        ImGui::Image(thumb_textures[compound.id], ImVec2(image_size, image_size));
        ImGui::PopStyleVar();

        ImGui::Dummy(ImVec2(0.0f, 4.0f));

        {
            ImGui::Text(compound.name.c_str());
            ImGui::Text(compound.formula.c_str());

            std::ostringstream info_stream;
            info_stream << comp_sel["molar_weight"].get<std::string>()
                        << " " << std::fixed << std::setprecision(2) << compound.molecular_weight
                        << " g/mol";

            std::string info_line = info_stream.str();

            ImVec2 info_size = ImGui::CalcTextSize(info_line.c_str());
            float info_x = (avail_width - info_size.x) * 0.5f;
            if (info_x > 0) ImGui::SetCursorPosX(ImGui::GetCursorPosX() + info_x);

            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.7f, 0.7f, 0.8f, 1.0f));
            ImGui::TextUnformatted(info_line.c_str());
            ImGui::PopStyleColor();
        }

        ImGui::Dummy(ImVec2(0.0f, 12.0f));
        
        ImGui::PushID(ImGui::GetID(compound.name.c_str()));

        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 6.0f);
        if (ImGui::Button(comp_sel["button_add"].get<std::string>().c_str(), ImVec2(-1.0f, 40.0f)))
        {
            simulation_universe->createMolecule(compound.structure, simulation_universe->camera().target);
        }
        ImGui::PopStyleVar();

        ImGui::PopID();
    }

    void UIHandler::drawCompoundFulLView()
    {
        auto &comp_sel = localization_json["Simulation"]["universe_ui"]["compound_selector"];
        auto &compounds = localization_json["Compounds"];

        compound_preset_info& current_info = compound_presets[selectedCompound];
    }

    void UIHandler::drawCompoundSelector()
    {
        auto &comp_sel = localization_json["Simulation"]["universe_ui"]["compound_selector"];
        auto &compounds = localization_json["Compounds"];

        ImGui::SetNextWindowSize(ImVec2(1000, 650), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

        if (!ImGui::Begin(comp_sel["title"].get<std::string>().c_str(), &compoundSelector))
        {
            ImGui::End();
            return;
        }

        if (ImGui::BeginTabBar("CompoundTabs"))
        {
            for (int32_t i = 0; i < static_cast<int32_t>(compound_type::COUNT); ++i)
            {
                compound_type current_type = static_cast<compound_type>(i);

                if (ImGui::BeginTabItem(compounds["compound_types"][i].get<std::string>().c_str()))
                {

                    int columns = static_cast<int>((ImGui::GetContentRegionAvail().x - padding) / (image_size + padding));
                    columns = std::max(1, columns);

                    if (ImGui::BeginTable("CompoundsGrid", columns, ImGuiTableFlags_ScrollY | ImGuiTableFlags_BordersInnerV))
                    {
                        for (const auto &compound : compound_presets)
                        {
                            if (compound.type != current_type)
                                continue;

                            ImGui::TableNextColumn();

                            ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(15, 15));
                            ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, 8.0f);
                            ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.12f, 0.12f, 0.18f, 0.9f));

                            ImGui::BeginChild(ImGui::GetID(compound.name.c_str()), 
                                            ImVec2(0, 0), ImGuiChildFlags_Border | ImGuiChildFlags_AutoResizeY);
                            drawCompoundView(compound);
                            ImGui::EndChild();

                            ImGui::PopStyleVar(2);
                            ImGui::PopStyleColor();
                        }
                        ImGui::EndTable();
                    }
                    ImGui::EndTabItem();
                }
            }
            ImGui::EndTabBar();
        }
        ImGui::End();
    }

    void UIHandler::runUniverse()
    {
        if (!paused_simulation)
            simulation_universe->update(target_temperature, target_pressure);
    }

    void UIHandler::drawUniverse(window_t &window)
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
        simulation_universe->draw(window, window.getWindow(), letter, ball);

        if (!pauseMenuOpen || !showHUD)
            drawUniverseUI();
        else if (pauseMenuOpen)
            pauseMenu();

        runUniverse();

        if (ImGui::IsKeyPressed(ImGuiKey_Space))
        {
            paused_simulation = !paused_simulation;
        }
        if (ImGui::IsKeyPressed(ImGuiKey_Escape))
        {
            pauseMenuOpen = !pauseMenuOpen;
            paused_simulation = true;
        }
    }

    void UIHandler::pauseMenu()
    {
        auto &pause_menu = localization_json["Simulation"]["pause_menu"];

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

        if (ImGui::Button(pause_menu["button_resume"].get<std::string>().c_str(), buttonSize) || ImGui::IsKeyPressed(ImGuiKey_Space))
        {
            pauseMenuOpen = false;
        }

        if (ImGui::Button(pause_menu["button_restart"].get<std::string>().c_str(), buttonSize))
        {
            simulation_universe = std::make_unique<sim::fun::universe>(sandbox_info);

            paused_simulation = false;
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
        
        ImGui::Dummy(ImVec2(0.0f, 20.0f));
        ImGui::TextWrapped(pause_menu["hint_esc"].get<std::string>().c_str());

        if (optionsOpen)
        {
            drawOptions();
        }

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
