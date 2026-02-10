#pragma once

#include <SFML/Graphics.hpp>
#include <imgui/imgui.h>

#include <vector>
#include <map>
#include <memory>

#include <functional>

#include "json.hpp"

#include "core/ScenarioHandler.hpp"
#include "simulation/reaction_engine.hpp"
#include "simulation/universe.hpp"
#include "simulation/simulation_packer.hpp"

namespace core
{
    enum class application_state
    {
        APP_STATE_MENU,
        APP_STATE_SIMULATION
    };

    enum class localization
    {
        EN_US,
        PT_BR,
        COUNT
    };

    struct scenario_info
    {
        std::string title;
        std::string description;

        std::filesystem::path file;
        std::filesystem::path video;

        bool is_locked = true;
        bool is_sandbox = false;

        bool operator==(const scenario_info& other) const
        {
            return other.file == file && other.title == title;
        }
    };

    enum class tutorial_type : uint32_t
    {
        ENGINE,
        CLASSICAL,
        QUANTUM,
        CHEMISTRY,
        COUNT
    };

    enum class scene_type : uint32_t
    {
        DYNAMICS,
        MOLECULE,
        REACTIONS,
        COUNT
    };

    inline ImVec4 getMoleculeTypeColor(sim::fun::compound_type _type)
    {
        using type = sim::fun::compound_type;

        switch(_type)
        {
        case type::ORGANIC:      return ImVec4(0.080f, 0.380f, 0.120f, 1.0f);  // Dark Green      #146824
        case type::INORGANIC:    return ImVec4(0.040f, 0.220f, 0.450f, 1.0f);  // Deep Blue       #0A3780
        case type::BIOMOLECULE:  return ImVec4(0.650f, 0.300f, 0.000f, 1.0f);  // Burnt Orange    #A64D00
        case type::ION:          return ImVec4(0.380f, 0.080f, 0.450f, 1.0f);  // Deep Purple     #611475
        case type::NANOMATERIAL: return ImVec4(0.300f, 0.300f, 0.300f, 1.0f);  // Dark Gray       #4D4D4D
        default:                 return ImVec4(0.200f, 0.200f, 0.200f, 1.0f);  // Neutral dark
        }
    }

    inline ImVec4 getMoleculeTypeColorBorder(sim::fun::compound_type _type)
    {
        using type = sim::fun::compound_type;
        switch(_type)
        {
        case type::ORGANIC:      return ImVec4(0.294f, 0.686f, 0.314f, 1.0f);  // #4CAF50 Vibrant Green
        case type::INORGANIC:    return ImVec4(0.129f, 0.588f, 0.953f, 1.0f);  // #2196F3 Cool Blue
        case type::BIOMOLECULE:  return ImVec4(1.000f, 0.596f, 0.000f, 1.0f);  // #FF9800 Warm Orange
        case type::ION:          return ImVec4(0.612f, 0.153f, 0.690f, 1.0f);  // #9C27B0 Electric Purple
        case type::NANOMATERIAL: return ImVec4(0.620f, 0.620f, 0.620f, 1.0f);  // #9E9E9E Sleek Gray
        default:                 return ImVec4(0.600f, 0.600f, 0.600f, 1.0f);
        }
    }

    // All options other than BALL_AND_STICK and Letter Mode are WIP
    enum class simulation_render_mode
    {
        BALL_AND_STICK,
        LICORICE,
        SPACE_FILLING,
        HYPER_BALLS,
        COUNT
    };

    struct simulation_options
    {
        int32_t target_fps{60};
        simulation_render_mode render_mode{simulation_render_mode::BALL_AND_STICK};
    };

    struct options
    {
        localization lang = localization::EN_US;

        // Visual
        simulation_options sim_options{};
        bool vsync{true};
        bool fullscreen{false};
        int32_t target_fps{60};

        // Audio
        float master_volume{0.8f};
        bool sound_effects{true};
        bool background_music{true};
    };

    class UIHandler
    {
    public:
        UIHandler(options& app_options, window_t& window);

        void set_language(localization new_lang) { lang = new_lang; write_localization_json(lang); }
        void setDeltaTime(float dt) { m_deltaTime = dt; }

        // Menu
        
        void drawMenu();

        // Universe

        void drawUniverse();
        void screenshotWindow(std::filesystem::path path, std::string name = "");

        // Callback

        void setApplicationStateCallback(std::function<void(application_state newState)> callback)
        {
            setState = callback;
        }

        void setGetApplicationStateCallback(std::function<application_state()> callback)
        {
            getState = callback;
        }

        // Fonts

        void set_regular_font (ImFont* regular_font) { regular = regular_font; }
        void set_regular_font_small (ImFont* regular_font_small) { regular_small = regular_font_small; }
        void set_regular_big_font (ImFont* regular_big_font) { regular_big = regular_big_font; }
        void set_medium_font (ImFont* medium_font) { medium = medium_font; }
        void set_bold_font (ImFont* bold_font) { bold = bold_font; }
        void set_bold_big_font (ImFont* bold_big_font) { bold_big = bold_big_font; }
        void set_black_font (ImFont* black_font) { black = black_font; }

        ImFont* regular_font () { return regular; }
        ImFont* regular_font_small () { return regular_small; }
        ImFont* regular_big_font () { return regular_big; }
        ImFont* medium_font() { return medium; }
        ImFont* bold_font() { return bold; }
        ImFont* bold_big_font() { return bold_big; }
        ImFont* black_font() { return black; }

        nlohmann::json save();

    private:
        localization lang{localization::EN_US};
        options& app_options;
        nlohmann::json localization_json{};
        nlohmann::json default_json{};

        window_t& m_window;

        void write_localization_json(localization lang);
        void initImages();
        void initSavedData();
        void loadScenariosFromFolder(const std::filesystem::path path, bool background = false);

        // Menu

        sf::Texture placeholder_texture;
        uint32_t placeholder_texture_id = 0;

        std::unordered_map<std::string, sf::Texture> thumb_textures;
        std::unordered_map<std::string, sf::Texture> textures;

        const float image_size = 200.f;

        const float m_displayMaxTime = 7.f;
        float m_currentDisplayTime = m_displayMaxTime + 1.f; // trigger initial setup on menu

        size_t m_backgroundDisplays = 0;
        void chooseNewDisplayScenario();
        void drawMenuBackgroundDisplay();
        void drawLoadingScreen();

        uint32_t m_currentDisplay = UINT32_MAX;
        std::vector<std::shared_ptr<sim::fun::universe>> m_backgroundUniverses;

        std::vector<scenario_info> m_savedSandbox;
        std::vector<sim::fun::video> m_SimulationVideos;

        ScenarioHandler m_scenarioHandler;
        
        void drawTutorialSelection();
        void drawSceneSelection();
        void drawSceneFrame(scenario_info& info, const std::string& id);
        void drawScene(scenario_info& info);
        void drawSandboxCreation();

        void startScenario(const std::string& scenario);

        void drawSandboxSave(scenario_info& info);
        void drawSavedSimulations();
        void drawOptions();

        bool sceneSelectionOpen = false;
        bool savesSelectionOpen = false;
        bool tutorialSelectionOpen = false;
        bool sandboxSelectionOpen = false;
        bool challengeSelectionOpen = false;
        bool challengeViewOpen = false;
        bool optionsOpen = false;
        bool statsOpen = false;
        bool videoPlayerOpen = true;

        bool savedSimulation = false;

        float m_deltaTime = 0.0f;

        const float panel_height = 64.0f;

        // Universe

        float target_temperature = 300.f;
        float target_pressure = 0.f;
        
        void runUniverse();
        void drawUniverseUI();
        void drawStatsWindow();
        void drawHUD();
        void drawTimeControl();

        sim::fun::rendering_info getSimulationRenderingInfo(simulation_render_mode mode);

        bool screenshotToggle = false;

        sim::fun::compound_preset_info m_currentSelectedCompound{};
        std::string selectedCompound = "";
        uint8_t m_selectedElement = UINT32_MAX;
        bool compoundFullView = false;
        bool newCompoundClicked = false;
        bool ghostDisplay = false;
        bool ghostColliding = false;

        void insertGhost();
        void insertGhostElement(std::string symbol = "H");
        void handleGhost();
        void drawCompoundSelector();
        void drawCompoundView(const sim::fun::compound_preset_info& compound);
        void drawCompoundFullView();
        void drawPeriodicTable();

        void pauseMenu();

        void initCompoundPresetsImages();
        void initCompoundXYZ();
        std::unordered_map<std::string, sim::fun::compound_preset_info> compound_presets{};

        bool compoundSelector = false;
        bool showHUD = true;

        bool pauseMenuOpen = false;
        bool exitDesktop = false;

        bool m_reactive = false;

        std::unique_ptr<sim::fun::universe> simulation_universe;
        std::unique_ptr<sim::fun::universe> display_universe;

        sim::rendering_engine m_rendering_eng;
        sim::reaction_engine m_reaction_eng;
        sim::simulation_packer m_simpacker{};

        // Callbacks
        
        std::function<void(application_state newState)> setState;
        std::function<application_state()> getState;

        // Video

        std::vector<float> temperature_log{};
        std::vector<int32_t> time_log{};

        void resetVideoData() 
        { 
            m_autoFrame = true;
            m_recordingFrames = false;
            m_rewinding = false;
            m_playingVideo = false;

            m_currentFrame = 0;

            temperature_log.clear();

            simulation_universe->clearDisplayPositions();
            display_universe->clearDisplayPositions();

            simulation_universe->clearFrames();
            display_universe->clearFrames();
        }

        void playFramesUniverse(sim::fun::universe& u);
        void drawVideoControls();
        
        float m_replaySpeed = 60.f; // in FPS
        float m_frameAccumulator  = 0.0f; // seconds
        bool m_autoFrame = true; // lock it when at last frame;
        bool m_recordingFrames = false;
        bool m_rewinding = false;
        bool m_playingVideo = false;

        int32_t m_currentFrame = 0;

        // Fonts
        ImFont* regular = nullptr;
        ImFont* regular_small = nullptr;
        ImFont* regular_big = nullptr;
        ImFont* medium = nullptr;
        ImFont* bold = nullptr;
        ImFont* bold_big = nullptr;
        ImFont* black = nullptr;
    };
} // namespace core
