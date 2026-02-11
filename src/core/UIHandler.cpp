#include "UIHandler.hpp"

#include <fstream>
#include <iostream>
#include <random>

#include <implot.h>

#include "simulation/smiles_parser.hpp"
#include "simulation/format_loader.hpp"

namespace core
{
    // Periodic Table

    ImVec4 getElementColor(char type)
    {
        switch (type)
        {
        case 'K':
            return ImVec4(1.00f, 0.40f, 0.40f, 1.00f);

        case 'T':
            return ImVec4(1.00f, 0.65f, 0.00f, 1.00f);

        case 'Z':
            return ImVec4(0.40f, 0.70f, 1.00f, 1.00f);

        case 'M':
            return ImVec4(0.70f, 0.70f, 0.80f, 1.00f);

        case 'N':
            return ImVec4(0.40f, 0.90f, 0.40f, 1.00f);

        case 'O':
            return ImVec4(0.80f, 1.00f, 1.00f, 1.00f);

        case 'L':
            return ImVec4(0.80f, 0.40f, 0.90f, 1.00f);

        case 'A':
            return ImVec4(0.60f, 0.20f, 0.80f, 1.00f);

        default:
            return ImVec4(0.50f, 0.50f, 0.50f, 1.00f);
        }

        return ImVec4(0.5f, 0.5f, 0.5f, 1.0f);
    }

    // Helper

    size_t count_all_items(std::filesystem::path path)
    {
        try
        {
            return static_cast<std::size_t>(std::distance(std::filesystem::directory_iterator{path}, std::filesystem::directory_iterator{}));
        }
        catch (std::filesystem::filesystem_error const &ex)
        {
            std::cerr << "Error accessing directory: " << ex.what() << std::endl;
            return 0;
        }
    }

    std::filesystem::path sandboxSave = std::filesystem::path("scenes/saved");

    std::filesystem::path getLocalizationFile(localization lang)
    {
        std::filesystem::path localization_path = "resource/localizations";
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

    UIHandler::UIHandler(options &app_options, window_t &window, AudioEngine& engine)
        : app_options(app_options), m_window(window), m_rendering_eng(window), m_audio_eng(engine),
          m_scenarioHandler(compound_presets, m_simpacker)
    {
        write_localization_json(lang);
        std::filesystem::path path = getLocalizationFile(localization::EN_US);
        std::ifstream file(path);
        if (!file.is_open())
        {
            std::cerr << "[UI Handler] Cannot open: " << path << '\n';
            return;
        }

        try
        {
            file >> default_json;
        }
        catch (const nlohmann::json::parse_error &e)
        {
            std::cerr << "[UI Handler] JSON parse error: " << e.what() << '\n';
            return;
        }

        initCompoundXYZ();
        initImages();
        initSavedData();
    }

    void UIHandler::initSavedData()
    {
        std::filesystem::path saved_sandbox_dir = "scenes/saved";
        std::filesystem::path scenarios_dir = "scenes/scenarios";
        std::filesystem::path saved_menu_displays_dir = "scenes/menu display";

        loadScenariosFromFolder(saved_sandbox_dir);
        loadScenariosFromFolder(scenarios_dir);
        loadScenariosFromFolder(saved_menu_displays_dir, true);
    }

    void UIHandler::loadScenariosFromFolder(const std::filesystem::path path, bool background)
    {
        if (!std::filesystem::exists(path))
        {
            std::filesystem::create_directories(path);
            return;
        }

        for (const auto &entry : std::filesystem::directory_iterator(path))
        {
            if (entry.is_regular_file() && entry.path().extension() == ".json" && entry.path().filename().string().find("video_") != 0)
            {
                std::ifstream file(entry.path());
                if (file.is_open())
                {
                    try
                    {
                        std::filesystem::path videoPath = path / ("video_" + entry.path().filename().string());

                        if (!background)
                        {
                            nlohmann::json scene_json = nlohmann::json::parse(file);

                            scenario_info nInfo{};
                            nInfo.file = entry;
                            nInfo.title = entry.path().filename().replace_extension("").string();

                            if (std::filesystem::directory_entry(videoPath).exists())
                                nInfo.video = videoPath;
                            else
                                videoPath = "";

                            nInfo.is_sandbox = true;

                            if (path != "scenes/scenarios")
                                m_savedSandbox.emplace_back(std::move(nInfo));
                        }
                        else
                        {
                            m_backgroundUniverses.emplace_back(std::make_shared<sim::fun::universe>(entry.path(), m_rendering_eng));

                            if (std::filesystem::directory_entry(videoPath).exists())
                                m_backgroundUniverses.back()->loadFrames(videoPath);

                            ++m_backgroundDisplays;
                        }

                        std::cout << "Loaded saved scene: " << entry.path() << std::endl;
                    }
                    catch (const nlohmann::json::parse_error &e)
                    {
                        std::cerr << "JSON parse error in " << entry.path() << ": " << e.what() << std::endl;
                    }
                }
            }
        }
    }

    void resize_texture(sf::Texture &tex, sf::Vector2u size)
    {
        if (!tex.resize(size))
        {
            std::cerr << "[UI HANDLER]: icon couldn't resize! \n";
        }
    }

    void load_texture(sf::Texture &tex, std::filesystem::path path)
    {
        if (!tex.loadFromFile(path))
        {
            std::cerr << "[UI HANDLER]: icon couldn't load! \n";
        }
    }

    void UIHandler::initImages()
    {
        initCompoundPresetsImages();

        sf::Image placeholder_image;
        placeholder_image.createMaskFromColor(sf::Color(80, 80, 90), 255); // Dark gray
        placeholder_texture_id = placeholder_texture.getNativeHandle();

        std::filesystem::path icons = "resource/images/icons";
        std::filesystem::path magnifying_glass = icons / "magnifying_glass.png";
        std::filesystem::path pause = icons / "pause.png";
        std::filesystem::path resume = icons / "resume.png";
        std::filesystem::path right_arrow = icons / "right_arrow.png";
        std::filesystem::path left_arrow = icons / "left_arrow.png";
        std::filesystem::path plus = icons / "plus.png";
        std::filesystem::path stats = icons / "stats.png";
        std::filesystem::path genesis_icon = icons / "Genesis.png";

        sf::Texture magnifying_texture{};
        sf::Texture pause_texture{};
        sf::Texture resume_texture{};
        sf::Texture right_arrow_texture{};
        sf::Texture left_arrow_texture{};
        sf::Texture plus_texture{};
        sf::Texture stats_texture{};
        sf::Texture genesis_icon_texture{};

        resize_texture(magnifying_texture, {64, 64});
        resize_texture(pause_texture, {32, 32});
        resize_texture(resume_texture, {32, 32});
        resize_texture(right_arrow_texture, {32, 32});
        resize_texture(left_arrow_texture, {32, 32});
        resize_texture(plus_texture, {32, 32});
        load_texture(magnifying_texture, magnifying_glass);
        load_texture(genesis_icon_texture, genesis_icon);
        load_texture(pause_texture, pause);
        load_texture(resume_texture, resume);
        load_texture(right_arrow_texture, right_arrow);
        load_texture(left_arrow_texture, left_arrow);
        load_texture(plus_texture, plus);
        load_texture(stats_texture, stats);

        textures.emplace("magnifying_texture", std::move(magnifying_texture));
        textures.emplace("genesis_icon", std::move(genesis_icon_texture));
        textures.emplace("pause_icon", std::move(pause_texture));
        textures.emplace("resume_icon", std::move(resume_texture));
        textures.emplace("right_arrow_icon", std::move(right_arrow_texture));
        textures.emplace("left_arrow_icon", std::move(left_arrow_texture));
        textures.emplace("plus_icon", std::move(plus_texture));
        textures.emplace("stats_icon", std::move(stats_texture));

        const std::filesystem::path saved_sandbox_dir = "scenes/saved";
        const std::filesystem::path scenarios_dir = "scenes/scenarios";
        if (!std::filesystem::exists(saved_sandbox_dir))
        {
            std::filesystem::create_directories(saved_sandbox_dir);
            return;
        }

        for (const auto &entry : std::filesystem::directory_iterator(saved_sandbox_dir))
        {
            if (entry.is_regular_file() && entry.path().extension() == ".png")
            {
                sf::Texture scenario_texture{};
                if (!scenario_texture.loadFromFile(entry.path()))
                {
                    std::cerr << "[UI HANDLER]: Couldn't load sandbox image!" << entry.path() << "\n";
                    continue;
                }

                textures.emplace(entry.path().filename().string(), std::move(scenario_texture));
            }
        }

        if (!std::filesystem::exists(scenarios_dir))
        {
            std::filesystem::create_directories(scenarios_dir);
            return;
        }

        for (const auto &entry : std::filesystem::directory_iterator(scenarios_dir))
        {
            if (entry.is_regular_file() && entry.path().extension() == ".png")
            {
                sf::Texture scenario_texture{};
                if (!scenario_texture.loadFromFile(entry.path()))
                {
                    std::cerr << "[UI HANDLER]: Couldn't load scenario image!" << entry.path() << "\n";
                    continue;
                }

                textures.emplace(entry.path().filename().string(), std::move(scenario_texture));
            }
        }
    }

    void UIHandler::initCompoundPresetsImages()
    {
        sim::fun::universe_create_info display_info{};
        display_info.box.x = 100.f;
        display_info.box.y = 100.f;
        display_info.box.z = 100.f;
        display_info.wall_collision = false;

        display_universe = std::make_unique<sim::fun::universe>(display_info, m_rendering_eng);
        auto &cam = m_rendering_eng.camera();
        cam.target = {0.f, 0.f, 0.f};
        cam.distance = 0.f;
        cam.azimuth = 200.f;
        cam.elevation = 25.f;
        cam.fov = 60.f;
        cam.nearPlane = 0.1f;
        cam.farPlane = 500.f;

        for (auto &[name, compound] : compound_presets)
        {
            sf::RenderTexture thumbnail_renderer;
            display_universe->clear();
            display_universe->createMolecule(compound.structure, {0.f, 0.f, 0.f});

            float max_radius = 0.f;
            for (const auto &pos : compound.structure.positions)
            {
                float r = pos.length();
                if (r > max_radius)
                    max_radius = r;
            }

            float molecule_radius = max_radius;
            cam.distance = molecule_radius;
            cam.distance = std::max(cam.distance, 10.f);

            rendering_info info = getSimulationRenderingInfo(simulation_render_mode::SPACE_FILLING);
            info.universeBox = false;

            m_window.clear();
            display_universe->draw(m_window.getWindow(), info);

            sf::Texture thumb_texture;
            if (!thumb_texture.resize(m_window.getWindow().getSize()))
                continue;

            thumb_texture.update(m_window.getWindow());
            thumb_textures.emplace(name, std::move(thumb_texture));
        }

        display_universe->clear();
    }

    void UIHandler::initCompoundXYZ()
    {
        std::filesystem::path folder = "resource/molecules/compounds";
        if (!std::filesystem::exists(folder) || !std::filesystem::is_directory(folder))
        {
            std::cerr << "[Compounds] Folder not found: " << folder << "\n";
            return;
        }

        size_t id_counter = 0;
        for (const auto &entry : std::filesystem::directory_iterator(folder))
        {
            if (!entry.is_regular_file() || entry.path().extension() != ".xyz")
                continue;

            std::string filename = entry.path().stem().string();

            sim::fun::molecule_structure structure{};
            if (!sim::io::loadXYZ(entry.path().string(), structure.atoms, structure.bonds, structure.positions))
            {
                std::cerr << "[Compounds] Failed to load XYZ: " << entry.path() << "\n";
                continue;
            }

            std::string displayName = filename;
            std::string formula = "Unknown";
            compound_type type = compound_type::ORGANIC;

            std::filesystem::path metaPath = entry.path();
            metaPath.replace_extension(".json");
            if (std::filesystem::exists(metaPath))
            {
                std::ifstream metaFile(metaPath);
                nlohmann::json metaJson;
                if (metaFile >> metaJson)
                {
                    displayName = metaJson.value("name", filename);
                    formula = metaJson.value("formula", formula);
                    type = static_cast<compound_type>(metaJson.value("type", 0));
                }
            }

            float molWeight = 0.0f;
            for (const auto &atom : structure.atoms)
            {
                molWeight += atom.ZIndex * MASS_PROTON + atom.NIndex * MASS_NEUTRON;
            }

            sim::organizeSubsets(structure.subsets, structure.atoms, structure.bonds);
            sim::organizeAngles(structure.subsets, structure.atoms, structure.bonds, structure.dihedral_angles, structure.improper_angles, structure.angles, displayName == "Carborane Acid");

            sim::fun::compound_preset_info info{};
            info.id = id_counter++;
            info.name = displayName;
            info.formula = formula;
            info.structure = std::move(structure);
            info.molecular_weight = molWeight;
            info.type = type;

            std::cout << "[Compounds] Loaded XYZ: " << displayName << " (" << formula << ", "
                      << info.structure.atoms.size() << " atoms)\n";

            compound_presets.emplace(displayName, std::move(info));
        }

        std::cout << "[Compounds] Loaded " << compound_presets.size() << " XYZ compounds\n";
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
        float buttonWidth = 200.0f;
        float buttonHeight = 45.0f;

        ImGui::PushFont(regular);
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.04f, 0.04f, 0.10f, 0.92f)); // Darker, more elegant
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.95f, 0.96f, 1.0f, 1.0f));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(25, 25));
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0, 10));

        {
            ImVec2 title_pos(0, (io.DisplaySize.y - 600) * 0.5f);
            ImGui::SetNextWindowPos(title_pos, ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(600, 600));
            ImGui::Begin("TitleOverlay", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                                                ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBackground | ImGuiWindowFlags_NoScrollbar);

            ImVec2 window_size = ImGui::GetWindowSize();

            float aspect_ratio = 5.f / 2.f;

            float image_width = window_size.x * 0.8f;
            float image_height = image_width / aspect_ratio;

            if (image_height > window_size.y * 0.8f) 
            {
                image_height = window_size.y * 0.8f;
                image_width = image_height * aspect_ratio;
            }

            ImGui::Image(textures["genesis_icon"], ImVec2(image_width, image_height));

            ImGui::End();
        }

        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.04f, 0.18f, 0.28f, 0.88f));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.08f, 0.32f, 0.48f, 1.00f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.12f, 0.42f, 0.62f, 1.00f));
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.3f, 0.7f, 1.0f, 0.5f));
        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 12.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 2.0f);

        {
            float total_button_height = 6 * (buttonHeight + 28.0f);
            ImVec2 panel_pos(100, (io.DisplaySize.y - total_button_height) * 0.6f + 80.0f);
            ImGui::SetNextWindowPos(panel_pos, ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(buttonWidth + 60, 0));

            ImGui::Begin("NavigationPanel", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                                                    ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBackground);

            // Buttons with consistent style
            if (ImGui::Button(localization_json["Menu"]["button_tutorials"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
                tutorialSelectionOpen = !tutorialSelectionOpen;

            if (ImGui::Button(localization_json["Menu"]["button_scenarios"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
                sceneSelectionOpen = !sceneSelectionOpen;

            if (ImGui::Button(localization_json["Menu"]["button_load"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
                savesSelectionOpen = !savesSelectionOpen;

            if (ImGui::Button(localization_json["Menu"]["button_sandbox"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
                sandboxSelectionOpen = !sandboxSelectionOpen;

            if (ImGui::Button(localization_json["Menu"]["button_options"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
                optionsOpen = !optionsOpen;

            if (ImGui::Button(localization_json["Menu"]["button_quit"].get<std::string>().c_str(), ImVec2(buttonWidth, buttonHeight)))
            {
                std::exit(EXIT_SUCCESS);
            }
            ImGui::End();
        }

        {
            ImVec2 infoPos(io.DisplaySize.x - padding - 20, io.DisplaySize.y - padding - 20);
            ImGui::SetNextWindowPos(infoPos, ImGuiCond_Always, ImVec2(1.0f, 1.0f));
            ImGui::Begin("InfoPanel", nullptr, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);

            ImGui::TextColored(ImVec4(0.6f, 0.7f, 0.9f, 0.8f), "%s", localization_json["Menu"]["info_version"].get<std::string>().c_str());
            ImGui::TextColored(ImVec4(0.6f, 0.7f, 0.9f, 0.8f), "%s", localization_json["Menu"]["info_author"].get<std::string>().c_str());

            ImGui::End();
        }

        if (tutorialSelectionOpen)
            drawTutorialSelection();
        if (sceneSelectionOpen)
            drawSceneSelection();
        if (sandboxSelectionOpen)
            drawSandboxCreation();
        if (optionsOpen)
            drawOptions();
        if (savesSelectionOpen)
            drawSavedSimulations();

        ImGui::PopStyleVar(4);
        ImGui::PopStyleColor(6);
        ImGui::PopFont();

        if (getState() != application_state::APP_STATE_MENU)
        {
            display_universe->clear();
            m_playingVideo = false;
            m_replaySpeed = 100.0f;

            m_audio_eng.stopSound("MainMenu_Music");
            m_audio_eng.playRandomSong();
        }
        else
        {
            drawMenuBackgroundDisplay();
            m_audio_eng.stopSong();

            if (app_options.background_music)
                m_audio_eng.playPreloaded("MainMenu_Music");
        }
    }

    constexpr float padding = 20.0f;
    constexpr float row_height = 200.f;
    constexpr float row_width = 200.f;

    void UIHandler::chooseNewDisplayScenario()
    {
        constexpr uint32_t maxTries = 5;
        uint32_t newDisplay = 0;
        uint32_t currentTries = 0;

        do
        {
            ++currentTries;

            srand(time(0));
            newDisplay = rand() % m_backgroundDisplays;
        } while (newDisplay == m_currentDisplay && currentTries < maxTries);

        m_currentDisplay = newDisplay;
    }

    void UIHandler::drawMenuBackgroundDisplay()
    {
        auto &camera = m_rendering_eng.camera();

        m_currentDisplayTime += m_deltaTime;
        if (m_currentDisplayTime >= m_displayMaxTime)
        {
            m_currentDisplayTime = 0.f;
            chooseNewDisplayScenario();

            const glm::vec3 boxSizes = m_backgroundUniverses[m_currentDisplay]->boxSizes();
            camera.target = boxSizes * 0.5f;

            float diagonal = glm::length(boxSizes) * 0.8f;
            camera.distance = diagonal;
            camera.azimuth = 45.f;
            camera.elevation = 25.f;
            m_playingVideo = true;
            m_replaySpeed = 100.0f;
            m_currentFrame = 0;
        }

        auto &current_universe = m_backgroundUniverses[m_currentDisplay];
        current_universe->draw(m_window.getWindow(), getSimulationRenderingInfo(app_options.sim_options.render_mode));
        camera.rotateAroundTarget(10.f * m_deltaTime, 0.f);

        if (current_universe->numFrames() > 0)
            playFramesUniverse(*current_universe.get());
    }

    static std::filesystem::path selected_for_delete;

    void UIHandler::drawSandboxSave(scenario_info &info)
    {
        const float image_save_size = image_size;
        const sf::Texture *thumbnail = &placeholder_texture;

        std::string thumb_key = info.title + ".png";
        if (textures.count(thumb_key))
            thumbnail = &textures.at(thumb_key);

        float avail_width = ImGui::GetContentRegionAvail().x;
        float offset_x = (avail_width - image_save_size) * 0.5f;
        if (offset_x > 0)
            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + offset_x);

        ImGui::Image(*thumbnail, ImVec2(image_save_size, image_save_size));

        ImGui::Dummy(ImVec2(0.0f, 10.0f));

        ImGui::PushFont(bold);
        float title_width = ImGui::CalcTextSize(info.title.c_str()).x;
        ImGui::SetCursorPosX((avail_width - title_width) * 0.7f);
        ImGui::Text("%s", info.title.c_str());
        ImGui::PopFont();

        auto &sim_loading = localization_json["Menu"]["Simulation_Loading"];
        float button_width = avail_width - 40.0f;

        ImGui::SetCursorPosX((avail_width - button_width) * 0.7f);
        if (ImGui::Button(sim_loading["button_load"].get<std::string>().c_str(), ImVec2(button_width, 40.0f)))
        {
            simulation_universe = std::make_unique<sim::fun::universe>(info.file, m_rendering_eng);
            savesSelectionOpen = false;
            pauseMenuOpen = false;
            simulation_universe->unpause();
            setState(application_state::APP_STATE_SIMULATION);

            resetVideoData();
        }

        ImGui::SetCursorPosX((avail_width - button_width) * 0.7f);
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.6f, 0.1f, 0.1f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.85f, 0.15f, 0.15f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.9f, 0.2f, 0.2f, 1.0f));
        if (ImGui::Button(sim_loading["button_delete"].get<std::string>().c_str(), ImVec2(button_width, 40.0f)))
        {
            selected_for_delete = info.file;
            ImGui::OpenPopup("ConfirmDelete");
        }
        ImGui::PopStyleColor(3);

        float button_width_popup = 130.0f;

        ImGui::SetNextWindowPos(ImVec2{ImGui::GetMainViewport()->GetCenter().x / 2.f, ImGui::GetMainViewport()->GetCenter().y});
        if (ImGui::BeginPopupModal("ConfirmDelete", nullptr,
                                   ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoTitleBar))
        {
            ImGui::TextWrapped(sim_loading["warning"].get<std::string>().c_str());
            ImGui::Dummy(ImVec2(0.0f, 20.0f));

            float button_width = 120.0f;
            float total_width = button_width * 2 + ImGui::GetStyle().ItemSpacing.x;
            ImGui::SetCursorPosX((ImGui::GetWindowWidth() - total_width) * 0.5f);

            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.65f, 0.15f, 0.15f, 1.0f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.85f, 0.2f, 0.2f, 1.0f));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.95f, 0.25f, 0.25f, 1.0f));
            if (ImGui::Button(sim_loading["true_delete"].get<std::string>().c_str(), ImVec2(button_width, 40)))
            {
                auto it = std::find_if(m_savedSandbox.begin(), m_savedSandbox.end(),
                                       [&](const scenario_info &a)
                                       { return a.file == selected_for_delete; });

                if (it != m_savedSandbox.end())
                {
                    std::filesystem::remove(selected_for_delete);
                    auto thumb_path = selected_for_delete;
                    thumb_path.replace_extension(".png");

                    std::filesystem::remove(thumb_path);

                    textures.erase(it->title + ".png");
                    m_savedSandbox.erase(it);
                }

                selected_for_delete.clear();
                ImGui::CloseCurrentPopup();
            }
            ImGui::PopStyleColor(3);

            ImGui::SameLine();

            if (ImGui::Button(sim_loading["cancel"].get<std::string>().c_str(), ImVec2(button_width, 40)))
            {
                selected_for_delete.clear();
                ImGui::CloseCurrentPopup();
            }

            ImGui::EndPopup();
        }
    }

    void UIHandler::drawSavedSimulations()
    {
        ImGui::SetNextWindowSize(ImVec2(800, 500), ImGuiCond_Appearing);
        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Always, ImVec2(0.5f, 0.5f));

        auto &sim_loading = localization_json["Menu"]["Simulation_Loading"];
        ImGui::Begin(sim_loading["title"].get<std::string>().c_str(), &savesSelectionOpen, ImGuiWindowFlags_NoMove);

        ImGui::Separator();

        const float cell_size = image_size;

        int32_t columns = static_cast<int>((ImGui::GetContentRegionAvail().x - padding) / (image_size + padding));
        columns = std::max(1, columns);
        if (ImGui::BeginTable("SavesGrid", columns, ImGuiTableFlags_BordersInnerV | ImGuiTableFlags_ScrollY))
        {
            for (auto &sandbox : m_savedSandbox)
            {
                std::string lower_title = sandbox.title;
                std::transform(lower_title.begin(), lower_title.end(), lower_title.begin(), ::tolower);

                ImGui::TableNextColumn();

                ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(15, 15));
                ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, 8.0f);
                ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.12f, 0.12f, 0.18f, 0.9f));

                ImGui::BeginChild(ImGui::GetID(sandbox.title.c_str()),
                                  ImVec2(0, 0), ImGuiChildFlags_Border | ImGuiChildFlags_AutoResizeY);
                drawSandboxSave(sandbox);

                ImGui::PopStyleColor();
                ImGui::PopStyleVar(2);

                ImGui::EndChild();
            }
            ImGui::EndTable();
        }

        ImGui::End();
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
                                 &sandbox_info.mag_gravity, 0.1f, 0.0f, 981.0f, "%.2f m/s²");
                ImGui::Unindent();
            }

            ImGui::Checkbox(sandbox_creation["wall_collision"].get<std::string>().c_str(), &sandbox_info.wall_collision);
            ImGui::Checkbox(sandbox_creation["roof_floor_collision"].get<std::string>().c_str(), &sandbox_info.roof_floor_collision);

            ImGui::BeginDisabled();
            ImGui::Checkbox(sandbox_creation["isothermal"].get<std::string>().c_str(), &sandbox_info.isothermal);
            ImGui::Checkbox(sandbox_creation["reactive"].get<std::string>().c_str(), &m_reactive);
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
                sandbox_info.box = glm::vec3(size, size, size);
            }
        }

        if (ImGui::CollapsingHeader(sandbox_creation["header_visual"].get<std::string>().c_str()))
        {
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
            simulation_universe = std::make_unique<sim::fun::universe>(sandbox_info, m_rendering_eng);
            m_rendering_eng.camera().target = sandbox_info.box / 2.f;
            m_rendering_eng.camera().distance = glm::length(sandbox_info.box * 1.2f);

            sandboxSelectionOpen = false;
            simulation_universe->unpause();

            setState(application_state::APP_STATE_SIMULATION);

            target_pressure = 0.f;
            target_temperature = 300.f;

            /* sim::fun::molecule_structure structure{};
            sim::io::loadXYZ("resource/molecules/ice.xyz", structure.atoms, structure.bonds, structure.positions);
            sim::organizeSubsets(structure.subsets, structure.atoms, structure.bonds);
            sim::organizeAngles(structure.subsets, structure.atoms, structure.bonds, structure.dihedral_angles, structure.improper_angles, structure.angles);

            simulation_universe->createMolecule(structure, {10, 10, 10}); */
        }

        if (ImGui::Button(sandbox_creation["button_cancel"].get<std::string>().c_str(), ImVec2(200, 50)))
        {
            sandboxSelectionOpen = false;
        }

        ImGui::End();
    }

    void UIHandler::drawTutorialSelection()
    {
        auto &tutorial_json = localization_json["Menu"]["Tutorial"];
        auto &tutorial_json_default = default_json["Menu"]["Tutorial"];
        auto &descriptions_json = localization_json["Scenarios"];

        ImGui::SetNextWindowSize(ImVec2(1000, 650), ImGuiCond_Appearing);
        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Always, ImVec2(0.5f, 0.5f));
        if (!ImGui::Begin(tutorial_json["title"].get<std::string>().c_str(), &tutorialSelectionOpen,
                          ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse))
        {
            ImGui::End();
            return;
        }

        if (ImGui::BeginTabBar("TutorialTabs", ImGuiTabBarFlags_DrawSelectedOverline))
        {
            for (int32_t i = 0; i < static_cast<int32_t>(tutorial_type::COUNT); ++i)
            {
                std::string category = tutorial_json_default["tutorial_types"][i].get<std::string>();

                if (ImGui::BeginTabItem(category.c_str()))
                {
                    for (auto &[key, scenario] : m_scenarioHandler.getScenarios())
                    {
                        if (key.find("TUTORIAL_" + category + "_") != 0)
                            continue;

                        ImGui::BeginChild(("tab_content_" + category).c_str(), ImVec2(0, 0), true);

                        scenario_info info{};
                        info.title = tutorial_json["titles"][category][key];
                        info.description = descriptions_json[key]["description"];
                        info.is_sandbox = false;
                        info.is_locked = false;

                        drawSceneFrame(info, key);

                        ImGui::Dummy(ImVec2(0, 8));

                        ImGui::EndChild();
                    }
                    ImGui::EndTabItem();
                }
            }
            ImGui::EndTabBar();
        }

        ImGui::End();
    }

    void UIHandler::drawSceneFrame(scenario_info &info, const std::string &id)
    {
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(12, 12));
        ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, 12.0f);
        ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.10f, 0.12f, 0.18f, 0.97f));
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.35f, 0.45f, 0.80f, 0.55f));

        std::string childId = "scene_frame_" + id;
        ImGui::BeginChild(childId.c_str(), ImVec2(-1, 250), true, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);

        ImGui::Columns(2, nullptr, false);
        ImGui::SetColumnWidth(0, 230.0f); 

        const sf::Texture* thumb = &placeholder_texture;
        std::string thumbKey = id + ".png";
        if (textures.count(thumbKey))
            thumb = &textures.at(thumbKey);

        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.3f, 0.7f, 1.0f, 0.6f));
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(4, 4));
        ImGui::Image(*thumb, ImVec2(220, 220));
        ImGui::PopStyleVar();
        ImGui::PopStyleColor();

        ImGui::NextColumn();
        
        ImGui::PushFont(bold);
        ImGui::TextColored(ImVec4(0.85f, 0.95f, 1.0f, 1.0f), "%s", info.title.c_str());
        ImGui::PopFont();

        ImGui::Spacing();
        ImGui::Spacing();

        ImGui::PushTextWrapPos(ImGui::GetContentRegionAvail().x);

        std::string desc = info.description.empty()
                            ? "Explore this simulation in real time. Discover the forces, behaviors and patterns that emerge."
                            : info.description;

        ImGui::BeginChild("desc_scroll", ImVec2(0, -60), false, ImGuiWindowFlags_AlwaysVerticalScrollbar);
        ImGui::TextWrapped("%s", desc.c_str());
        ImGui::EndChild();

        ImGui::PopTextWrapPos();

        ImGui::Spacing();
        ImGui::Spacing();

        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.25f, 0.55f, 0.95f, 0.9f));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.35f, 0.65f, 1.0f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.20f, 0.45f, 0.85f, 1.0f));

        float buttonWidth = ImGui::CalcTextSize("Start Simulation").x + 40.0f;
        ImGui::SetCursorPosX(ImGui::GetColumnWidth() - buttonWidth - 10.0f);

        if (ImGui::Button("Start Simulation", ImVec2(buttonWidth, 38)))
        {
            tutorialSelectionOpen = false;
            sceneSelectionOpen = false;
            startScenario(id);
        }

        ImGui::PopStyleColor(3);

        ImGui::Columns(1);
        ImGui::EndChild();

        ImGui::PopStyleColor(2);
        ImGui::PopStyleVar(2);
    }

    void UIHandler::drawLoadingScreen()
    {
    }

    void UIHandler::drawSceneSelection()
    {
        auto &scene_json = localization_json["Menu"]["Scene_Selection"];
        auto &descriptions_json = localization_json["Scenarios"];
        auto &scene_json_default = default_json["Menu"]["Scene_Selection"];

        ImGui::SetNextWindowSize(ImVec2(1000, 650), ImGuiCond_Appearing);
        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Always, ImVec2(0.5f, 0.5f));
        if (!ImGui::Begin(scene_json["title"].get<std::string>().c_str(), &sceneSelectionOpen,
                          ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse))
        {
            ImGui::End();
            return;
        }

        if (ImGui::BeginTabBar("SceneTabs", ImGuiTabBarFlags_DrawSelectedOverline))
        {
            for (int32_t i = 0; i < static_cast<int32_t>(scene_type::COUNT); ++i)
            {
                std::string category = scene_json_default["scene_types"][i].get<std::string>();

                if (ImGui::BeginTabItem(category.c_str()))
                {
                    for (auto &[key, scenario] : m_scenarioHandler.getScenarios())
                    {
                        if (key.find("SCENARIO_" + category + "_") != 0)
                            continue;

                        ImGui::BeginChild(("tab_content_" + category).c_str(), ImVec2(0, 0), true);

                        scenario_info info{};
                        info.title = scene_json["titles"][category][key];
                        info.is_sandbox = false;
                        info.is_locked = false;
                        info.description = descriptions_json[key]["description"];

                        drawSceneFrame(info, key);

                        ImGui::Dummy(ImVec2(0, 8));

                        ImGui::EndChild();
                    }

                    ImGui::EndTabItem();
                }
            }
            ImGui::EndTabBar();
        }

        ImGui::End();
    }

    void UIHandler::startScenario(const std::string &scenario)
    {
        auto &chosen_scenario = m_scenarioHandler.getScenarios().at(scenario);

        simulation_universe = std::make_unique<sim::fun::universe>(chosen_scenario.file, m_rendering_eng);

        if (!chosen_scenario.file.empty())
            simulation_universe->loadScene(chosen_scenario.file);
        
        if (!chosen_scenario.video.empty())
            simulation_universe->loadFrames(chosen_scenario.video);

        m_scenarioHandler.setCurrentUniverse(simulation_universe.get());
        m_scenarioHandler.chooseScenario(scenario);
        m_scenarioHandler.startScenario();

        auto &cam = m_rendering_eng.camera();
        cam.target = simulation_universe->boxSizes() * 0.5f;
        cam.distance = glm::length(simulation_universe->boxSizes() * 1.2f);

        setState(core::application_state::APP_STATE_SIMULATION);
    }

    void UIHandler::drawOptions()
    {
        auto &options = localization_json["Menu"]["Options"];

        ImVec2 viewport_size = ImGui::GetMainViewport()->Size;

        ImGui::SetNextWindowSizeConstraints(ImVec2(500, 300),  // Min
                                            ImVec2(800, 600)); // Max
        ImGui::SetNextWindowSize(ImVec2(0, 0), ImGuiCond_Always);
        ImGui::SetNextWindowPos(ImVec2(viewport_size.x / 2.f - 200, viewport_size.y / 3.f), ImGuiCond_Appearing);

        ImGui::Begin("Options", &optionsOpen, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse);
        ImGui::BeginTabBar(options["title"].get<std::string>().c_str(), ImGuiTabBarFlags_NoCloseWithMiddleMouseButton);
        if (ImGui::BeginTabItem(options["tab_graphics"].get<std::string>().c_str()))
        {
            ImGui::Text(options["tab_graphics"].get<std::string>().c_str());

            static std::string ball_and_stick = "";
            static std::string licorice = "";
            static std::string letter_and_stick = "";
            static std::string hyper_balls = "";
            static std::string space_filling = "";

            ball_and_stick = options["render_modes"]["ball_and_stick"].get<std::string>();
            letter_and_stick = options["render_modes"]["letter_and_stick"].get<std::string>();
            space_filling = options["render_modes"]["space_filling"].get<std::string>();
            hyper_balls = options["render_modes"]["hyper_balls"].get<std::string>();
            licorice = options["render_modes"]["licorice"].get<std::string>();

            const char *modes[] =
                {
                    ball_and_stick.c_str(),
                    licorice.c_str(),
                    space_filling.c_str()};

            static int32_t current_mode = static_cast<int32_t>(app_options.sim_options.render_mode);
            if (ImGui::Combo("##render_mode", &current_mode, modes, IM_ARRAYSIZE(modes)))
                app_options.sim_options.render_mode = static_cast<simulation_render_mode>(current_mode);

            ImGui::SliderInt(options["tab_target_fps"].get<std::string>().c_str(), &app_options.sim_options.target_fps, 24, 144);

            ImGui::EndTabItem();
        }
        if (ImGui::BeginTabItem(options["tab_audio"].get<std::string>().c_str()))
        {
            ImGui::SliderFloat(options["tab_audio_master"].get<std::string>().c_str(), &app_options.master_volume, 0.0f, 1.0f, "%.2f");
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

    void UIHandler::drawStatsWindow()
    {
        auto &sim_ui = localization_json["Simulation"]["universe_ui"]["stats"];

        ImGuiIO& io = ImGui::GetIO();

        ImVec2 size(600, 600);

        ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f), ImGuiCond_Appearing);
        ImGui::SetNextWindowSize(size, ImGuiCond_Always);

        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 8.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.05f, 0.05f, 0.12f, 0.95f));
        ImGui::PushStyleColor(ImGuiCol_TitleBgActive, ImVec4(0.15f, 0.15f, 0.25f, 1.0f));

        ImGuiWindowFlags flags = ImGuiWindowFlags_NoResize |
                                ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar;

        std::string window_title = sim_ui["title"];
        if (ImGui::Begin(window_title.c_str(), &statsOpen, flags))
        {
            ImGui::PushFont(bold);
            ImGui::TextColored(ImVec4(0.8f, 0.9f, 1.0f, 1.0f), window_title.c_str());
            ImGui::PopFont();
            ImGui::Separator();

            auto &core_json = sim_ui["core_stats"];
            auto &energy_json = sim_ui["energy_stats"];
            auto &graph_json = sim_ui["graphs_stats"];

            std::string time_string = core_json["time"].get<std::string>().c_str();
            std::string temperature_string = core_json["temperature"].get<std::string>().c_str();

            // Basic stats
            if (ImGui::CollapsingHeader(core_json["title"].get<std::string>().c_str(), ImGuiTreeNodeFlags_DefaultOpen))
            {

                ImGui::Text("%s: %.2f K", temperature_string.c_str(), simulation_universe->temperature());
                ImGui::Text("%s:    %.2f Kpa", core_json["pressure"].get<std::string>().c_str(), simulation_universe->pressure());

                ImGui::SetNextItemWidth(100.f);
                ImGui::DragFloat(core_json["slider_temperature"].get<std::string>().c_str(), &target_temperature, 1.0f, 0.0f, 20000.f);
                ImGui::SetNextItemWidth(100.f);
                ImGui::DragFloat(core_json["slider_pressure"].get<std::string>().c_str(), &target_pressure, 1.0f, 0.0f, 1000.f);
                ImGui::Text("%s:        %.2f ps", time_string.c_str(), simulation_universe->getAccumulatedTime());
                ImGui::Text("%s:   %zu", core_json["particles"].get<std::string>().c_str(), simulation_universe->numAtoms());
                ImGui::Text("%s:   %zu", core_json["molecules"].get<std::string>().c_str(), simulation_universe->numMolecules());
            }

            // Energy stats
            if (ImGui::CollapsingHeader(energy_json["title"].get<std::string>().c_str()))
            {
                float ke = simulation_universe->calculateKineticEnergy() / 1000.f;
                ImGui::Text("%s: %.2e kJ", energy_json["kinetic_energy"].get<std::string>().c_str(), ke);
                //ImGui::Text("%s: N/A", energy_json["potential_energy"].get<std::string>().c_str());
                //ImGui::Text("%s:   %.2e kJ", energy_json["total_energy"].get<std::string>().c_str(), ke);
            }

            if (ImGui::CollapsingHeader(graph_json["title"].get<std::string>().c_str())) {
                ImGui::Text(graph_json["warning"].get<std::string>().c_str());

                const float* time = reinterpret_cast<const float*>(time_log.data());
                const float* temperature = reinterpret_cast<const float*>(temperature_log.data());
                int count = static_cast<int>(std::min(time_log.size(), temperature_log.size()));

                if (ImPlot::BeginPlot(graph_json["temperature_v_time"].get<std::string>().c_str())) 
                {
                    std::string x_label = std::format("{} (fs)", time_string); 
                    std::string y_label = std::format("{} (°K)", temperature_string); 
                    ImPlot::SetupAxes(x_label.c_str(), y_label.c_str()); 

                    ImPlot::SetupAxisLimits(ImAxis_Y1, 0.0f, 1000.f, ImPlotCond_Always);

                    ImPlot::SetupLegend(ImPlotLocation_NorthEast, ImPlotLegendFlags_None);

                    ImPlot::PlotLine(time_string.c_str(), time, temperature, count);
                    ImPlot::PlotScatter(temperature_string.c_str(), time, temperature, count);

                    ImPlot::EndPlot(); 
                }
            }

            ImGui::End();
        }

        ImGui::PopStyleColor(2);
        ImGui::PopStyleVar(2);
    }

    void UIHandler::drawTimeControl()
    {
        auto &sim_ui = localization_json["Simulation"]["universe_ui"];

        ImGuiIO &io = ImGui::GetIO();

        bool sim_paused = simulation_universe->isPaused();
        const std::string mode = sim_paused ? "resume_icon" : "pause_icon";

        ImVec2 timecontrol_size(panel_height * 6.f, panel_height);

        ImGui::Dummy(ImVec2(20, 0));

        ImGui::BeginChild("TimeControl", timecontrol_size);

        if (ImGui::ImageButton("##timebutton", textures[mode], ImVec2(32.f, 32.f)))
        {
            sim_paused ? simulation_universe->unpause() : simulation_universe->pause();
        }

        std::string resume_tooltip = sim_ui.value("tooltip_resume", "Resume Simulation");
        std::string pause_tooltip  = sim_ui.value("tooltip_pause", "Pause Simulation");

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip(sim_paused ? resume_tooltip.c_str() : pause_tooltip.c_str());

        ImGui::SameLine();

        ImGui::SetNextItemWidth(150.f);
        ImGui::Text("   %.2f ps   %.1f fs/s   ", simulation_universe->getAccumulatedTime(), simulation_universe->getTimescale());

        ImGui::SameLine();

        if (ImGui::ImageButton("##slowdownbutton", textures["left_arrow_icon"], ImVec2(32.f, 32.f)))
        {
            simulation_universe->setTimescale(simulation_universe->getTimescale() * 0.9f);
        }
        
        if (ImGui::IsItemActive() && ImGui::IsMouseDown(ImGuiMouseButton_Left))
        {
            simulation_universe->setTimescale(simulation_universe->getTimescale() * 0.99f);
        }
        
        std::string slowdown_tooltip  = sim_ui.value("tooltip_slowdown", "Decrease simulation speed");

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip(slowdown_tooltip.c_str());

        ImGui::SameLine();

        if (ImGui::ImageButton("##speedupbutton", textures["right_arrow_icon"], ImVec2(32.f, 32.f)))
        {
            simulation_universe->setTimescale(simulation_universe->getTimescale() * 1.1f);
        }

        if (ImGui::IsItemActive() && ImGui::IsMouseDown(ImGuiMouseButton_Left))
        {
            simulation_universe->setTimescale(simulation_universe->getTimescale() * 1.01f);
        }

        std::string speedup_tooltip  = sim_ui.value("tooltip_speedup", "Increase simulation speed");

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip(speedup_tooltip.c_str());

        ImGui::EndChild();
    }

    void UIHandler::drawHUD()
    {
        auto &sim_ui = localization_json["Simulation"]["universe_ui"];
        ImGuiIO &io = ImGui::GetIO();

        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.20f, 0.25f, 0.35f, 0.60f));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(20, 10));
        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(10, 5));
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 12.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 1.0f);

        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.04f, 0.18f, 0.28f, 0.88f));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.08f, 0.32f, 0.48f, 1.00f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.12f, 0.42f, 0.62f, 1.00f));
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.3f, 0.7f, 1.0f, 0.5f));
        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 12.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 2.0f);

        ImVec2 panel_pos(0, io.DisplaySize.y);
        ImGui::SetNextWindowPos(panel_pos, ImGuiCond_Always, ImVec2(0.0f, 1.0f));
        ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x, panel_height), ImGuiCond_Always);

        ImGuiWindowFlags hud_flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize |
                                    ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoScrollbar |
                                    ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoScrollWithMouse;

        ImGui::Begin("HUD", &showHUD, hud_flags);

        ImGui::SetCursorPosX((io.DisplaySize.x) * 0.5f);

        drawTimeControl();

        ImGui::SameLine();

        ImGui::BeginGroup();

        ImGui::SetCursorPosX((io.DisplaySize.x - 80.0f) * 0.5f);
        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 1.0f);
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.3f, 0.3f, 0.7f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.4f, 0.4f, 0.6f, 1.0f));

        if (ImGui::ImageButton("##compoundsbutton", textures["plus_icon"], ImVec2(32.f, 32.f)))
        {
            compoundSelector = !compoundSelector;
        }

        ImGui::PopStyleColor(2);
        ImGui::PopStyleVar();

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip(sim_ui.value("tooltip_compound_selector", "Add Compound").c_str());

        ImGui::SameLine();
        ImGui::SetCursorPosX((io.DisplaySize.x - 64.0f));

        if (ImGui::ImageButton("##statsbutton", textures["stats_icon"], ImVec2(32.f, 32.f)))
        {
            statsOpen = !statsOpen;
        }

        ImGui::PopStyleColor(4);
        ImGui::PopStyleVar(2);

        ImGui::EndGroup();

        ImGui::End();

        ImGui::PopStyleVar(4);
        ImGui::PopStyleColor();
    }

    void UIHandler::drawUniverseUI()
    {
        ImGui::PushFont(regular);
        drawHUD();

        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.04f, 0.18f, 0.28f, 0.88f));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.08f, 0.32f, 0.48f, 1.00f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.12f, 0.42f, 0.62f, 1.00f));
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.3f, 0.7f, 1.0f, 0.5f));
        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 12.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 2.0f);

        if (statsOpen)
            drawStatsWindow();

        if (videoPlayerOpen)
            drawVideoControls();

        if (compoundSelector)
            drawCompoundSelector();

        ImGui::PopStyleColor(4);
        ImGui::PopStyleVar(2);

        ImGui::PopFont();
    }

    void UIHandler::insertGhost()
    {
        auto &compound = compound_presets[selectedCompound];

        display_universe->clear();
        display_universe->createMolecule(compound.structure, sf::Vector3f(m_rendering_eng.camera().target.x, m_rendering_eng.camera().target.y, m_rendering_eng.camera().target.z));

        sf::Vector3f Centroid = sf::Vector3f{0.f, 0.f, 0.f};
        for (const auto &p : compound.structure.positions)
        {
            Centroid += p;
        }
        Centroid /= static_cast<float>(compound.structure.positions.size());

        for (auto &p : compound.structure.positions)
        {
            p -= Centroid;
        }

        m_currentSelectedCompound = compound;

        ghostDisplay = true;
    }

    void UIHandler::insertGhostElement(std::string symbol)
    {
        molecule_structure element = sim::parseSMILES(symbol, false);

        display_universe->clear();
        display_universe->createMolecule(element, sf::Vector3f(m_rendering_eng.camera().target.x, m_rendering_eng.camera().target.y, m_rendering_eng.camera().target.z));

        ghostDisplay = true;

        sf::Vector3f Centroid = sf::Vector3f{0.f, 0.f, 0.f};
        for (const auto &p : element.positions)
        {
            Centroid += p;
        }
        Centroid /= static_cast<float>(element.positions.size());

        for (auto &p : element.positions)
        {
            p -= Centroid;
        }

        compound_preset_info nInfo{};
        nInfo.structure = element;
        m_currentSelectedCompound = std::move(nInfo);
    }

    void UIHandler::handleGhost()
    {
        if (!ghostDisplay)
            return;

        ghostColliding = false;

        auto &sim_ui = localization_json["Simulation"]["universe_ui"];
        auto &cam = m_rendering_eng.camera();

        ImVec2 mousePos = ImGui::GetMousePos();

        const auto &base_positions = m_currentSelectedCompound.structure.positions;
        for (size_t i = 0; i < base_positions.size(); ++i)
        {
            sf::Vector3f pos = base_positions[i] + sf::Vector3f(cam.target.x, cam.target.y, cam.target.z);
            display_universe->setPosition(i, glm::vec3(pos.x, pos.y, pos.z));
        }

        constexpr float minDistance = 1.6f;
        for (auto &pos : display_universe->positions())
        {
            for (auto &other_pos : simulation_universe->positions())
            {
                glm::vec3 r = pos - other_pos;
                if (simulation_universe->minImageVec(sf::Vector3f(r.x, r.y, r.z)).length() <= minDistance)
                {
                    ghostColliding = true;
                }
            }
        }

        ImGui::BeginTooltip();
        if (ghostColliding)
        {
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.3f, 0.3f, 1.0f));
            ImGui::Text(sim_ui["tooltip_collision_bad"].get<std::string>().c_str());
            ImGui::PopStyleColor();
        }
        else
        {
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.3f, 1.0f, 0.3f, 1.0f));
            ImGui::Text(sim_ui["tooltip_collision_good"].get<std::string>().c_str());
            ImGui::PopStyleColor();
        }

        ImGui::Text(sim_ui["tooltip_emplace_molecule"].get<std::string>().c_str());
        ImGui::EndTooltip();

        if (ImGui::IsKeyPressed(ImGuiKey_G) && !ghostColliding)
        {
            simulation_universe->createMolecule(m_currentSelectedCompound.structure, sf::Vector3f(cam.target.x, cam.target.y, cam.target.z));

            if (app_options.sound_effects)
                m_audio_eng.playPreloaded("Place_Effect", 0.1f);

            if (!ImGui::IsKeyDown(ImGuiKey_LeftShift))
            {
                selectedCompound = "";
                ghostDisplay = false;
            }
        }

        if (ImGui::IsKeyPressed(ImGuiKey_MouseRight))
        {
            ghostDisplay = false;
            display_universe->clear();
            selectedCompound = "";
        }
    }

    void UIHandler::drawCompoundView(const compound_preset_info &compound)
    {
        auto &comp_sel = localization_json["Simulation"]["universe_ui"]["compound_selector"];
        auto &compounds = localization_json["Compounds"];

        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(10, 8));

        float avail_width = ImGui::GetContentRegionAvail().x;
        float offset_x = (avail_width - image_size) * 0.5f;
        if (offset_x > 0)
            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + offset_x);

        ImGui::Image(thumb_textures[compound.name], ImVec2(image_size, image_size));
        ImGui::PopStyleVar();

        ImGui::SameLine(11.f);

        ImGui::PushStyleColor(ImGuiCol_FrameBg, ImVec4(0.8f, 0.8f, 0.8f, 0.0f));
        ImGui::PushStyleColor(ImGuiCol_FrameBgHovered, ImVec4(1.0f, 1.0f, 1.0f, 0.1f));
        ImGui::PushStyleColor(ImGuiCol_FrameBgActive, ImVec4(1.0f, 1.0f, 1.0f, 0.2f));
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.5f, 0.5f, 0.5f, 1.0f));
        ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 1.0f, 1.0f, 1.0f));

        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 1.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 8.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(4, 4));
        if (ImGui::ImageButton("##fullviewbutton", textures["magnifying_texture"], ImVec2(30.f, 30.f), sf::Color::White))
        {
            selectedCompound = compound.name;
            compoundFullView = true;
            newCompoundClicked = true;
        }
        ImGui::PopStyleColor(5);
        ImGui::PopStyleVar(3);

        if (ImGui::IsItemHovered())
        {
            ImGui::SetTooltip(comp_sel["tooltip_fullview"].get<std::string>().c_str());
        }

        ImGui::Dummy(ImVec2(0.0f, 4.0f));

        std::string name = app_options.lang == localization::EN_US ? compound.name : compounds["compound_names"][compound.name].get<std::string>();

        {
            ImVec2 name_size = ImGui::CalcTextSize(name.c_str());
            ImVec2 formula_size = ImGui::CalcTextSize(compound.formula.c_str());

            float content_width = ImGui::GetContentRegionAvail().x;

            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + (content_width - name_size.x) * 0.5f);
            ImGui::Text("%s", name.c_str());

            ImGui::SetCursorPosX(ImGui::GetCursorPosX() + (content_width - formula_size.x) * 0.5f);
            ImGui::Text("%s", compound.formula.c_str());

            std::ostringstream info_stream;
            info_stream << comp_sel["molar_weight"].get<std::string>()
                        << " " << std::fixed << std::setprecision(2) << compound.molecular_weight
                        << " g/mol";

            std::string info_line = info_stream.str();

            ImVec2 info_size = ImGui::CalcTextSize(info_line.c_str());
            float info_x = (avail_width - info_size.x) * 0.5f;
            if (info_x > 0)
                ImGui::SetCursorPosX(ImGui::GetCursorPosX() + info_x);

            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.7f, 0.7f, 0.8f, 1.0f));
            ImGui::TextUnformatted(info_line.c_str());
            ImGui::PopStyleColor();
        }

        ImGui::Dummy(ImVec2(0.0f, 12.0f));
        ImGui::PushID(ImGui::GetID(name.c_str()));

        ImGui::SetCursorPosX(row_width * 0.5f);
        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 6.0f);
        if (ImGui::ImageButton("##addcompound", textures["plus_icon"].getNativeHandle(), ImVec2(32.0f, 32.0f)))
        {
            compoundSelector = false;
            selectedCompound = compound.name;
            insertGhost();
        }

        ImGui::PopStyleVar();

        ImGui::PopID();
    }

    void UIHandler::drawCompoundFullView()
    {
        auto &comp_sel = localization_json["Simulation"]["universe_ui"]["compound_selector"];
        auto &full_view = comp_sel["full_view"];
        auto &compounds = localization_json["Compounds"];

        // Do elements later
        bool isElement = m_selectedElement != UINT32_MAX;

        const sim::fun::compound_preset_info &compound = compound_presets[selectedCompound];

        ImGui::PushStyleColor(ImGuiCol_ModalWindowDimBg, ImVec4(0.0f, 0.0f, 0.0f, 0.7f));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, getMoleculeTypeColor(compound.type));
        ImGui::PushStyleColor(ImGuiCol_Border, getMoleculeTypeColorBorder(compound.type));

        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 4.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 8.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(4, 4));
        ImGui::PushStyleVar(ImGuiStyleVar_ScrollbarSize, 4.0f);

        ImGui::SetNextWindowSize(ImVec2(400, 600), ImGuiCond_Appearing);
        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        if (newCompoundClicked)
        {
            ImGui::SetNextWindowFocus();
            newCompoundClicked = false;
        }

        if (!ImGui::Begin("##FullView", &compoundFullView, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize))
        {
            ImGui::End();
            return;
        }

        ImGui::PushFont(bold);
        ImGui::TextColored(ImVec4(0.8f, 0.9f, 1.0f, 1.0f), "%s", full_view["preview"].get<std::string>().c_str());
        ImGui::PopFont();

        ImGui::SameLine(ImGui::GetWindowWidth() - 100);
        if (ImGui::Button(full_view["close"].get<std::string>().c_str()))
        {
            compoundFullView = false;
        }

        std::string name = app_options.lang == localization::EN_US ? compound.name : compounds["compound_names"][compound.name].get<std::string>();

        ImGui::SetWindowFontScale(1.2f);
        ImGui::Text(name.c_str());
        ImGui::SetWindowFontScale(1.f);

        ImGui::NewLine();
        ImGui::BeginChild("##3DPreview", ImVec2(400, 410), true, ImGuiWindowFlags_HorizontalScrollbar);
        {
            ImGui::SetCursorPosX(15);
            ImGui::Image(thumb_textures[compound.name],
                         ImVec2(360, 360));
        }
        ImGui::EndChild();

        ImGui::BeginChild("InfoPanel", ImVec2(0, 390), true);
        {
            ImGui::PushStyleColor(ImGuiCol_TextDisabled, ImVec4(0.8f, 0.8f, 0.8f, 1.0f));
            ImGui::TextDisabled(comp_sel["formula"].get<std::string>().c_str());
            ImGui::SameLine();
            ImGui::Text("%s", compound.formula.c_str());

            ImGui::TextDisabled(comp_sel["molar_weight"].get<std::string>().c_str());
            ImGui::SameLine();
            ImGui::Text("%.3f g/mol", compound.molecular_weight);

            ImGui::Dummy(ImVec2(0, 15));

            ImGui::TextDisabled(full_view["description"].get<std::string>().c_str());

            // std::string description = localization_json["Compounds"].count(compound.name) == 0 ? compounds["compound_descriptions"][compound.name] : localization_json[default_json["Compounds"]["compound_names"][compound.id]];
            std::string description = compounds["compound_descriptions"][compound.name];
            ImGui::TextWrapped(description.c_str());

            ImGui::Dummy(ImVec2(0, 20));

            ImGui::SetCursorPosY(330);
            if (ImGui::Button(full_view["add_simulation"].get<std::string>().c_str(), ImVec2(-1, 50)))
            {
                compoundSelector = false;
                compoundFullView = false;
                insertGhost();
            }
            ImGui::PopStyleColor();
        }
        ImGui::EndChild();
        ImGui::PopStyleColor(3);
        ImGui::PopStyleVar(4);

        ImGui::End();
    }

    void UIHandler::drawCompoundSelector()
    {
        auto &comp_sel = localization_json["Simulation"]["universe_ui"]["compound_selector"];
        auto &compounds = localization_json["Compounds"];

        ImGui::SetNextWindowSize(ImVec2(1000, 650), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));

        if (!ImGui::Begin(comp_sel["title"].get<std::string>().c_str(), &compoundSelector, ImGuiWindowFlags_NoCollapse))
        {
            ImGui::End();
            return;
        }

        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.04f, 0.18f, 0.28f, 0.88f));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.08f, 0.32f, 0.48f, 1.00f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.12f, 0.42f, 0.62f, 1.00f));
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.3f, 0.7f, 1.0f, 0.5f));
        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 12.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 2.0f);

        if (ImGui::BeginTabBar("##CompoundTabs", ImGuiTabBarFlags_DrawSelectedOverline))
        {
            for (int32_t i = 0; i < static_cast<int32_t>(compound_type::COUNT) - 1; ++i)
            {
                compound_type current_type = static_cast<compound_type>(i);

                if (ImGui::BeginTabItem(compounds["compound_types"][i].get<std::string>().c_str()))
                {
                    int32_t columns = static_cast<int>((ImGui::GetContentRegionAvail().x - padding) / (image_size + padding));
                    columns = std::max(1, columns);

                    if (ImGui::BeginTable("CompoundsGrid", columns, ImGuiTableFlags_ScrollY | ImGuiTableFlags_BordersInnerV))
                    {
                        for (const auto &[name, compound] : compound_presets)
                        {
                            if (compound.type != current_type)
                                continue;

                            ImGui::TableNextColumn();

                            ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(15, 15));
                            ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, 8.0f);
                            ImGui::PushStyleColor(ImGuiCol_ChildBg, ImVec4(0.12f, 0.12f, 0.18f, 0.9f));

                            ImGui::BeginChild(ImGui::GetID(name.c_str()),
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

            if (ImGui::BeginTabItem(compounds["compound_types"][static_cast<uint8_t>(compound_type::ELEMENTS)].get<std::string>().c_str()))
            {
                drawPeriodicTable();
                ImGui::EndTabItem();
            }

            ImGui::EndTabBar();
        }

        if (compoundFullView)
            drawCompoundFullView();

        ImGui::PopStyleColor(4);
        ImGui::PopStyleVar(2);

        ImGui::End();
    }

    void UIHandler::drawPeriodicTable()
    {
        ImGui::BeginChild("##PeriodicTableScroll", ImVec2(0, 0), true);

        constexpr int32_t cols = 18;
        constexpr int32_t rows = 9;

        constexpr uint8_t Zindices[rows][cols] =
            {
                {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2},
                {3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 6, 7, 8, 9, 10},
                {11, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 13, 14, 15, 16, 17, 18},
                {19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36},
                {37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54},
                {55, 56, 57, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86},
                {87, 88, 89, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118},
                {0, 0, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 0, 0},
                {0, 0, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 0, 0}};

        constexpr float atomicMass[9][18] =
            {
                {1.008f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 4.003f},
                {6.941f, 9.012f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 10.811f, 12.011f, 14.007f, 15.999f, 18.998f, 20.180f},
                {22.990f, 24.305f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 26.982f, 28.085f, 30.974f, 32.06f, 35.45f, 39.948f},
                {39.098f, 40.078f, 44.956f, 47.867f, 50.942f, 51.996f, 54.938f, 55.845f, 58.933f, 58.693f, 63.546f, 65.38f, 69.723f, 72.64f, 74.922f, 78.972f, 79.904f, 83.798f},
                {85.468f, 87.62f, 88.906f, 91.224f, 92.906f, 95.95f, 98.0f, 101.07f, 102.906f, 105.905f, 107.868f, 112.414f, 114.818f, 118.711f, 121.760f, 127.60f, 126.904f, 131.293f},
                {132.905f, 137.33f, 138.905f, 178.49f, 180.948f, 183.84f, 186.207f, 190.23f, 192.217f, 195.085f, 196.967f, 200.592f, 204.38f, 207.2f, 208.980f, 209.0f, 210.0f, 222.0f},
                {223.0f, 226.0f, 227.0f, 232.038f, 231.036f, 238.029f, 237.0f, 244.0f, 243.0f, 247.0f, 247.0f, 251.0f, 252.0f, 257.0f, 258.0f, 259.0f, 262.0f, 267.0f},
                {0.0f, 0.0f, 140.116f, 140.908f, 144.24f, 145.0f, 150.36f, 151.964f, 157.25f, 158.925f, 162.500f, 164.930f, 167.259f, 168.934f, 173.045f, 174.967f, 0.0f, 0.0f},
                {0.0f, 0.0f, 227.028f, 232.038f, 231.036f, 238.029f, 237.0f, 244.0f, 243.0f, 247.0f, 247.0f, 251.0f, 252.0f, 257.0f, 258.0f, 259.0f, 262.0f, 267.0f}};

        constexpr char ElementType[rows][cols] =
            {
                //  1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18
                {'N', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'O'},
                {'K', 'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'M', 'N', 'N', 'N', 'N', 'O'},
                {'K', 'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'M', 'M', 'N', 'N', 'N', 'O'},
                {'K', 'T', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'M', 'M', 'M', 'N', 'N', 'O'},
                {'K', 'T', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'M', 'M', 'M', 'M', 'N', 'O'},
                {'K', 'T', 'L', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'M', 'M', 'M', 'M', 'N', 'O'},
                {'K', 'T', 'A', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'Z', 'M', 'M', 'M', 'M', 'N', 'O'},
                {' ', ' ', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', ' ', ' '},
                {' ', ' ', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', ' ', ' '}};

        constexpr char ElementBlock[rows][cols] =
            {
                //  1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18
                {'s', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 's'},
                {'s', 's', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'p', 'p', 'p', 'p', 'p', 'p'},
                {'s', 's', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'p', 'p', 'p', 'p', 'p', 'p'},
                {'s', 's', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'p', 'p', 'p', 'p', 'p', 'p'},
                {'s', 's', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'p', 'p', 'p', 'p', 'p', 'p'},
                {'s', 's', 'f', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'p', 'p', 'p', 'p', 'p', 'p'},
                {'s', 's', 'f', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'd', 'p', 'p', 'p', 'p', 'p', 'p'},
                {' ', ' ', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', ' ', ' '},
                {' ', ' ', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', ' ', ' '}};

        auto &comp_sel = localization_json["Simulation"]["universe_ui"]["compound_selector"];

        constexpr float box_size = 54.f;

        ImGui::PushFont(bold_big);
        std::string title = comp_sel["periodic_title"].get<std::string>();
        std::string subtitle = comp_sel["periodic_subtitle"].get<std::string>();

        ImVec2 titleSize = ImGui::CalcTextSize(title.c_str());
        ImGui::SetCursorPosX(ImGui::GetContentRegionAvail().x * 0.5f - titleSize.x * 0.5f);

        ImGui::BeginGroup();

        ImGui::Text(title.c_str());
        ImGui::PopFont();

        ImVec2 subtitleSize = ImGui::CalcTextSize(subtitle.c_str());

        ImGui::SetCursorPosX(ImGui::GetCursorPosX() + subtitleSize.x * 0.5f);
        ImGui::Text(subtitle.c_str());
        ImGui::EndGroup();

        ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(4, 4));

        ImDrawList *draw_list = ImGui::GetWindowDrawList();

        for (int32_t row = 0; row < rows; ++row)
        {
            for (int32_t col = 0; col < cols; ++col)
            {
                uint8_t Z = Zindices[row][col];
                if (Z == 0)
                {
                    ImGui::Dummy(ImVec2(box_size, box_size));
                    if (col < cols - 1)
                        ImGui::SameLine();
                    continue;
                }

                ImGui::PushID(Z);

                std::string symbol = constants::getAtomLetter(Z);
                std::string name = comp_sel["periodic_names"][row][col];
                char block = ElementBlock[row][col];
                float mass = atomicMass[row][col];

                ImVec4 color = getElementColor(ElementType[row][col]);

                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(color.x * 0.8f, color.y * 0.8f, color.z * 0.8f, 1.0f));
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(color.x * 1.15f, color.y * 1.15f, color.z * 1.15f, 1.0f));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(color.x * 0.8f, color.y * 0.8f, color.z * 0.8f, 1.0f));

                if (ImGui::Button(std::string("##" + symbol).c_str(), ImVec2(box_size, box_size)))
                {
                    m_selectedElement = Z;
                    compoundFullView = true;

                    // upper part is for when implement fullview for elements
                    insertGhostElement(symbol);
                    compoundSelector = false;
                    compoundFullView = false;
                }

                ImVec2 buttonMin = ImGui::GetItemRectMin();
                ImVec2 buttonMax = ImGui::GetItemRectMax();
                ImVec2 buttonCenter = ImVec2((buttonMin.x + buttonMax.x) * 0.5f, (buttonMin.y + buttonMax.y) * 0.5f);

                if (ImGui::IsItemHovered())
                {
                    ImGui::BeginTooltip();
                    ImGui::Text("%s (%s)", name.c_str(), symbol.c_str());
                    ImGui::Text("Z = %d\nMass = %.3f u", Z, mass);
                    ImGui::Text("Block: %c", block);
                    ImGui::EndTooltip();
                }

                ImU32 textColor = ImGui::GetColorU32(ImVec4(1, 1, 1, 1));
                ImU32 smallTextColor = ImGui::GetColorU32(ImVec4(1.0f, 1.0f, 1.0f, 1.0f));
                ImU32 nameColor = ImGui::GetColorU32(ImVec4(1.0f, 1.0f, 1.0f, 1.0f));

                ImGui::PushFont(regular_small);
                draw_list->AddText(ImVec2(buttonMin.x + 4.0f, buttonMin.y + 4.0f),
                                   smallTextColor,
                                   std::to_string(Z).c_str());

                draw_list->AddText(ImVec2(buttonCenter.x + (mass > 100 ? 1.f : 6.f), buttonMin.y + 4.0f),
                                   smallTextColor,
                                   std::format("{:.1f}", mass).c_str());
                ImGui::PopFont();

                ImGui::PushFont(bold);
                ImVec2 symbolSize = ImGui::CalcTextSize(symbol.c_str());
                draw_list->AddText(ImVec2(buttonCenter.x - symbolSize.x * 0.5f, buttonCenter.y - symbolSize.y * 0.5f - 4.0f),
                                   textColor,
                                   symbol.c_str());
                ImGui::PopFont();

                ImGui::PushFont(regular_small);
                ImVec2 nameSize = ImGui::CalcTextSize(name.c_str());
                draw_list->AddText(ImVec2(buttonCenter.x - nameSize.x * 0.5f, buttonMax.y - nameSize.y - 4.0f),
                                   nameColor,
                                   name.c_str());
                ImGui::PopFont();

                ImGui::PopStyleColor(3);

                if (col < cols - 1)
                    ImGui::SameLine();

                ImGui::PopID();
            }
        }

        ImGui::PopStyleVar();
        ImGui::EndChild();
    }

    rendering_info UIHandler::getSimulationRenderingInfo(simulation_render_mode mode)
    {
        rendering_info info{};
        info.opacity = 1.0f;

        if (mode == simulation_render_mode::BALL_AND_STICK)
        {
            info.lennardBall = true;
            info.licorice = false;
        }
        else if (mode == simulation_render_mode::LICORICE)
        {
            info.lennardBall = info.licorice = true;
        }
        else if (mode == simulation_render_mode::HYPER_BALLS)
        {
            info.lennardBall = info.licorice = info.spaceFilling = false;
            info.hyperBalls = true;
        }
        else if (mode == simulation_render_mode::SPACE_FILLING)
        {
            info.lennardBall = info.spaceFilling = true;
            info.licorice = false;
        }

        return info;
    }

    void UIHandler::runUniverse()
    {
        if (simulation_universe->isPaused())
            return;

        auto start = std::chrono::high_resolution_clock::now();

        simulation_universe->update(target_temperature, target_pressure);

        if (m_reactive) m_reaction_eng.update(*simulation_universe.get());

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start;

        //std::cout << "Simulation execution time: " << duration.count() << " milliseconds" << std::endl;

        if (m_recordingFrames && simulation_universe->timestep())
        {
            if (simulation_universe->timestep() % 100 == 0)
            {
                temperature_log.emplace_back(simulation_universe->temperature());
                time_log.emplace_back(simulation_universe->timestep());
            }

            if (simulation_universe->timestep() % 5 == 0)
            {
                simulation_universe->saveFrame();
            }

            if (m_autoFrame)
                m_currentFrame = simulation_universe->numFrames();
        }
    }

    void UIHandler::drawUniverse()
    {
        simulation_render_mode mode = app_options.sim_options.render_mode;

        rendering_info info = getSimulationRenderingInfo(app_options.sim_options.render_mode);
        info.universeBox = !(screenshotToggle && savedSimulation);

        simulation_universe->draw(m_window.getWindow(), info);
        m_rendering_eng.handleCamera();

        if (screenshotToggle)
        {
            std::filesystem::path path = "src/resource/screenshots";
            std::string title = "";

            if (savedSimulation)
            {
                path = sandboxSave;
                title = m_savedSandbox.at(m_savedSandbox.size() - 1).title;
            }

            screenshotWindow(path, title);
            screenshotToggle = false;
        }

        info.universeBox = false;
        info.opacity = 0.7f;
        info.color_addition = ghostColliding ? ImVec4(0.8f, 0.f, 0.f, 0.f) : ImVec4(0.f, 0.f, 0.f, 0.f);
        if (ghostDisplay)
        {
            handleGhost();
            display_universe->draw(m_window.getWindow(), info);
        }

        if (!pauseMenuOpen || !showHUD)
            drawUniverseUI();
        else if (pauseMenuOpen)
            pauseMenu();

        runUniverse();

        playFramesUniverse(*simulation_universe.get());

        if (m_scenarioHandler.inScenario())
        {
            m_scenarioHandler.draw(localization_json, m_deltaTime);
            target_temperature = m_scenarioHandler.getWantedTemperature();
            m_playingVideo = m_scenarioHandler.getWantedVideoPlay();
            m_replaySpeed = m_scenarioHandler.getWantedVideoSpeed();
            
            if (m_scenarioHandler.exit())
            {
                pauseMenuOpen = false;
                compoundFullView = false;
                compoundSelector = false;
                savesSelectionOpen = false;
                optionsOpen = false;

                setState(application_state::APP_STATE_MENU);

                m_currentDisplayTime = m_displayMaxTime + 1.f;

                resetVideoData();

                m_scenarioHandler.clear();
            }
        }

        if (ImGui::IsKeyPressed(ImGuiKey_Space))
        {
            simulation_universe->isPaused() ? simulation_universe->unpause() : simulation_universe->pause();
        }
        if (ImGui::IsKeyPressed(ImGuiKey_Escape))
        {
            pauseMenuOpen = !pauseMenuOpen;
            optionsOpen = false;
            savesSelectionOpen = false;
            compoundSelector = false;
            simulation_universe->pause();
        }
    }

    std::string formatTime()
    {
        const auto now = std::chrono::system_clock::now();

        std::string time_str = std::format("{:%Y_%m_%d%H_%M}", now);
        return time_str;
    }

    void UIHandler::screenshotWindow(std::filesystem::path path, std::string name)
    {
        sf::Texture screenshot{};
        if (!screenshot.resize(m_window.getWindow().getSize()))
        {
            std::cerr << "[UI HANDLER]: Failed to Screenshot!";
            return;
        }

        screenshot.update(m_window.getWindow());
        sf::Image image = screenshot.copyToImage();

        std::filesystem::path filepath = name.empty() ? path / (formatTime() + ".png") : path / (name + ".png");
        if (!image.saveToFile(filepath))
        {
            std::cerr << "[UI HANDLER]: Failed to save Screenshot!";
        }

        if (!name.empty())
        {
            textures.emplace(name + ".png", screenshot);
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

        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.04f, 0.18f, 0.28f, 0.88f));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.08f, 0.32f, 0.48f, 1.00f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.12f, 0.42f, 0.62f, 1.00f));
        ImGui::PushStyleColor(ImGuiCol_Border, ImVec4(0.3f, 0.7f, 1.0f, 0.5f));
        ImGui::PushStyleVar(ImGuiStyleVar_FrameRounding, 12.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 2.0f);

        ImGui::Begin("##pause menu", &pauseMenuOpen, ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize);

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
            simulation_universe = std::make_unique<sim::fun::universe>(sandbox_info, m_rendering_eng);
            m_rendering_eng.camera().target = {sandbox_info.box.x / 2.f, sandbox_info.box.y / 2.f, sandbox_info.box.z / 2.f};

            target_pressure = 0.f;
            target_temperature = 300.f;
            simulation_universe->unpause();
            pauseMenuOpen = false;

            m_scenarioHandler.restart();
            resetVideoData();
        }

        if (ImGui::Button(pause_menu["button_save"].get<std::string>().c_str(), buttonSize))
        {
            ImGui::OpenPopup("SavePopup");
        }

        if (ImGui::Button(pause_menu["button_load"].get<std::string>().c_str(), buttonSize))
        {
            savesSelectionOpen = !savesSelectionOpen;
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

        float button_width_popup = 130.0f;

        if (ImGui::BeginPopupModal("SavePopup", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar))
        {
            static char input[256];

            ImGui::Text(pause_menu["simulation_name"].get<std::string>().c_str());
            ImGui::InputText("##InputName", input, 256 * sizeof(char), ImGuiInputTextFlags_EnterReturnsTrue);

            ImGui::SetNextItemWidth(150.f);
            if (ImGui::Button(pause_menu["confirm_save"].get<std::string>().c_str(), ImVec2(-1.f, 30)))
            {
                scenario_info saved{};
                saved.title = input;
                saved.is_sandbox = true;
                saved.file = sandboxSave / (std::string(input) + ".json");

                if (std::find(m_savedSandbox.begin(), m_savedSandbox.end(), saved) == m_savedSandbox.end())
                    m_savedSandbox.emplace_back(std::move(saved));

                simulation_universe->saveScene(sandboxSave, std::string(input));
                (void)simulation_universe->saveAsVideo(sandboxSave, std::string(input));

                savedSimulation = true;
                screenshotToggle = true;
            }

            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.6f, 0.1f, 0.f, 1.f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.9f, 0.1f, 0.f, 1.f));
            if (ImGui::Button(pause_menu["cancel"].get<std::string>().c_str(), ImVec2(-1.f, 30)))
            {
                ImGui::CloseCurrentPopup();
            }
            ImGui::PopStyleColor(2);

            ImGui::EndPopup();
        }

        if (ImGui::BeginPopupModal("ConfirmExitToMenu", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar))
        {
            ImVec2 center = ImGui::GetMainViewport()->GetCenter();
            ImGui::SetWindowPos({center.x - button_width_popup * 1.5f, center.y - button_width_popup});
            ImGui::TextWrapped(pause_menu["warning"].get<std::string>().c_str());

            ImGui::Dummy(ImVec2(0.0f, 20.0f));

            ImGui::SetCursorPosX((ImGui::GetWindowWidth() - (button_width_popup * 2 + ImGui::GetStyle().ItemSpacing.x)) * 0.5f);

            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.6f, 0.1f, 0.f, 1.f));
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.9f, 0.1f, 0.f, 1.f));
            if (ImGui::Button(pause_menu["true_exit"].get<std::string>().c_str(), ImVec2(button_width_popup, 30)))
            {
                pauseMenuOpen = false;
                compoundFullView = false;
                compoundSelector = false;
                savesSelectionOpen = false;
                optionsOpen = false;
                ghostDisplay = false;
                statsOpen = false;

                setState(application_state::APP_STATE_MENU);

                m_currentDisplayTime = m_displayMaxTime + 1.f;

                resetVideoData();

                m_scenarioHandler.clear();

                if (exitDesktop)
                    std::exit(EXIT_SUCCESS);
                ImGui::CloseCurrentPopup();
            }
            ImGui::PopStyleColor(2);

            ImGui::SameLine();

            if (ImGui::Button(pause_menu["cancel"].get<std::string>().c_str(), ImVec2(button_width_popup, 30)))
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

        if (savesSelectionOpen)
        {
            drawSavedSimulations();
        }

        ImGui::End();

        ImGui::PopStyleVar(4);
        ImGui::PopStyleColor(6);
    }

    // Video

    void UIHandler::playFramesUniverse(sim::fun::universe &u)
    {
        if (m_playingVideo && u.numFrames() > 0)
        {
            const float timePerFrame = 1.0f / m_replaySpeed;

            m_frameAccumulator += m_deltaTime;

            if (m_frameAccumulator > 2.f)
                m_frameAccumulator = 1.f; // Prevent too big
            while (m_frameAccumulator >= timePerFrame)
            {
                m_frameAccumulator -= timePerFrame;
                ++m_currentFrame;

                m_currentFrame = m_currentFrame % u.numFrames();
            }

            u.setDisplayPositions(u.getFrame(m_currentFrame).positions);
        }
    }

    void UIHandler::drawVideoControls()
    {
        auto &video_controller = localization_json["Simulation"]["universe_ui"]["video_controller"];

        ImGui::SetNextWindowPos(ImVec2(0, ImGui::GetMainViewport()->Size.y - panel_height * 1.8f), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(ImGui::GetMainViewport()->Size.x, 50), ImGuiCond_Always);
        ImGui::SetNextWindowBgAlpha(0.85f);

        if (!ImGui::Begin("Video Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar))
        {
            ImGui::End();
            return;
        }

        if (ImGui::Button("Play"))
        {
            m_playingVideo = !m_playingVideo;
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip(m_playingVideo ? video_controller["pause_tooltip"].get<std::string>().c_str() : video_controller["play_tooltip"].get<std::string>().c_str());

        ImGui::SameLine();

        if (ImGui::Button("Restart"))
        {
            m_currentFrame = 0;
        }

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip(video_controller["restart_tooltip"].get<std::string>().c_str());

        ImGui::SameLine();

        if (ImGui::Button("Record"))
        {
            m_recordingFrames = !m_recordingFrames;

            if (!m_recordingFrames)
                simulation_universe->clearDisplayPositions();
        }

        ImGui::SameLine();

        if (ImGui::IsItemHovered())
            ImGui::SetTooltip(m_recordingFrames ? video_controller["record_stop_tooltip"].get<std::string>().c_str() : video_controller["record_start_tooltip"].get<std::string>().c_str());

        ImGui::DragFloat(video_controller["play_speed"].get<std::string>().c_str(), &m_replaySpeed, 0.5f, 1.0f, 1000.f, "%.1f");

        ImGui::SameLine();

        ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x - 50);
        if (ImGui::SliderInt("##Frame", &m_currentFrame, 0, simulation_universe->numFrames(),
                             "Frame %d", ImGuiSliderFlags_AlwaysClamp))
        {
            m_playingVideo = false;

            if (m_currentFrame != simulation_universe->numFrames())
            {
                simulation_universe->setDisplayPositions(simulation_universe->getFrame(m_currentFrame).positions);
                m_autoFrame = false;
            }
            else
            {
                m_autoFrame = true;
                simulation_universe->clearDisplayPositions();
            }
        }

        ImGui::End();
    }

    // Saving

    nlohmann::json UIHandler::save()
    {
        nlohmann::json UI_json{};

        return UI_json;
    }
} // namespace core
