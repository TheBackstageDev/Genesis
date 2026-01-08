#pragma once

#include <SFML/Graphics.hpp>
#include <imgui/imgui.h>

#include <vector>
#include <map>
#include <memory>

#include <functional>

#include "json.hpp"

#include "core/ScenarioHandler.hpp"
#include "simulation/universe.hpp"

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

    enum class compound_type : uint32_t
    {
        ORGANIC,
        BIOMOLECULE,
        INORGANIC,
        ION,
        NANOMATERIAL,
        POLYMER,
        COUNT
    };

    enum class tutorial_type : uint32_t
    {
        ENGINE,
        CLASSICAL,
        QUANTUM,
        CHEMISTRY,
        COUNT
    };

    enum class scenario_type : uint32_t
    {
        DYNAMICS,
        MOLECULE,
        REACTIONS,
        COUNT
    };

    inline ImVec4 getMoleculeTypeColor(compound_type _type)
    {
        using type = compound_type;

        switch(_type)
        {
        case type::ORGANIC:      return ImVec4(0.080f, 0.380f, 0.120f, 1.0f);  // Dark Green      #146824
        case type::INORGANIC:    return ImVec4(0.040f, 0.220f, 0.450f, 1.0f);  // Deep Blue       #0A3780
        case type::BIOMOLECULE:  return ImVec4(0.650f, 0.300f, 0.000f, 1.0f);  // Burnt Orange    #A64D00
        case type::ION:          return ImVec4(0.380f, 0.080f, 0.450f, 1.0f);  // Deep Purple     #611475
        case type::NANOMATERIAL: return ImVec4(0.300f, 0.300f, 0.300f, 1.0f);  // Dark Gray       #4D4D4D
        case type::POLYMER:      return ImVec4(0.650f, 0.100f, 0.080f, 1.0f);  // Dark Red        #A61A14
        default:                 return ImVec4(0.200f, 0.200f, 0.200f, 1.0f);  // Neutral dark
        }
    }

    inline ImVec4 getMoleculeTypeColorBorder(compound_type _type)
    {
        using type = compound_type;
        switch(_type)
        {
        case type::ORGANIC:      return ImVec4(0.294f, 0.686f, 0.314f, 1.0f);  // #4CAF50 Vibrant Green
        case type::INORGANIC:    return ImVec4(0.129f, 0.588f, 0.953f, 1.0f);  // #2196F3 Cool Blue
        case type::BIOMOLECULE:  return ImVec4(1.000f, 0.596f, 0.000f, 1.0f);  // #FF9800 Warm Orange
        case type::ION:          return ImVec4(0.612f, 0.153f, 0.690f, 1.0f);  // #9C27B0 Electric Purple
        case type::NANOMATERIAL: return ImVec4(0.620f, 0.620f, 0.620f, 1.0f);  // #9E9E9E Sleek Gray
        case type::POLYMER:      return ImVec4(0.957f, 0.263f, 0.212f, 1.0f);  // #F44336 Deep Red
        default:                 return ImVec4(0.600f, 0.600f, 0.600f, 1.0f);
        }
    }

    struct compound_preset_info
    {
        std::string name;
        std::string SMILES;
        std::string formula;
        float molecular_weight = 0.0f;

        compound_type type;
        sim::fun::molecule_structure structure;
        uint32_t id = 0;
    };

    // All options other than BALL_AND_STICK and Letter Mode are WIP
    enum class simulation_render_mode
    {
        BALL_AND_STICK,
        LETTER_AND_STICK,
        SPACE_FILLING,
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
        void set_medium_font (ImFont* medium_font) { medium = medium_font; }
        void set_bold_font (ImFont* bold_font) { bold = bold_font; }
        void set_black_font (ImFont* black_font) { black = black_font; }

        ImFont* regular_font () { return regular; }
        ImFont* medium_font() { return medium; }
        ImFont* bold_font() { return bold; }
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

        std::vector<sf::Texture> thumb_textures;
        std::unordered_map<std::string, sf::Texture> textures;

        const float image_size = 200.f;

        const float m_displayMaxTime = 15.f;
        float m_currentDisplayTime = m_displayMaxTime + 1.f; // trigger initial setup on menu
        void chooseNewDisplayScenario();
        void drawMenuBackgroundDisplay();

        uint32_t m_currentDisplay = 0;
        std::vector<scenario_info> m_backgroundDisplays;
        std::vector<scenario_info> m_savedSandbox;
        std::vector<sim::fun::video> m_SimulationVideos;

        void drawScene(scenario_info& info);
        void drawSceneSelection();
        void drawSceneFrame(scenario_info& info, int32_t id);
        void drawSandboxCreation();

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

        bool savedSimulation = false;

        float m_deltaTime = 0.0f;

        // Universe

        float target_temperature = 300.f;
        float target_pressure = 0.f;
        
        void runUniverse();
        void drawUniverseUI();
        sim::fun::rendering_info getSimulationRenderingInfo(simulation_render_mode mode);

        bool screenshotToggle = false;

        uint32_t selectedCompound = UINT32_MAX;
        bool compoundFullView = false;
        bool newCompoundClicked = false;
        bool ghostDisplay = false;
        bool ghostColliding = false;

        void insertGhost();
        void handleGhost();
        void drawCompoundSelector();
        void drawCompoundView(const compound_preset_info& compound);
        void drawCompoundFulLView();

        void pauseMenu();

        void initCompoundPresetsImages();
        void initCompoundPresets();
        std::vector<compound_preset_info> compound_presets{};

        bool compoundSelector = false;
        bool showHUD = true;

        bool pauseMenuOpen = false;
        bool exitDesktop = false;

        bool paused_simulation = false;

        std::unique_ptr<sim::fun::universe> simulation_universe;
        std::unique_ptr<sim::fun::universe> display_universe;

        sim::rendering_engine rendering_eng;

        // Callbacks
        
        std::function<void(application_state newState)> setState;
        std::function<application_state()> getState;

        // Video

        void drawVideoControls();

        sim::fun::video* m_currentVideo = nullptr;
        bool m_rewinding = false;
        size_t m_currentFrame = 0;

        // Fonts
        ImFont* regular = nullptr;
        ImFont* medium = nullptr;
        ImFont* bold = nullptr;
        ImFont* black = nullptr;
    };
} // namespace core
