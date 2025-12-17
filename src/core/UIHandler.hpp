#pragma once

#include <SFML/Graphics.hpp>
#include <imgui/imgui.h>

#include <vector>
#include <map>
#include <memory>

#include <functional>

#include "json.hpp"
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

    struct videoMetaData
    {
        std::string title;
        std::string description;

        sf::Vector3f box;
        size_t num_atoms;
        size_t num_frames;
    };

    struct frame
    {
        std::vector<sf::Vector3f> positions;
        std::map<size_t, float> temperatures;

        float global_temperature;
    };

    struct video
    {
        std::vector<frame> frames; // index is defined by order
        std::vector<size_t> keyFrames; // for events

        std::vector<std::string> text; // changes every keyframe (or continue);

        videoMetaData metadata;
    };

    struct scenario_info
    {
        std::string title;
        std::string description;

        std::filesystem::path snapshot;
        std::filesystem::path file;
        std::filesystem::path video;

        bool has_video = true;
        bool is_locked = true;
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

        // Video

        void drawVideoControls();

        // Menu
        
        void drawMenu();
        bool isrewinding() { return rewinding; }

        // Universe

        void drawUniverse(window_t& window);

        // Callback

        void setApplicationStateCallback(std::function<void(application_state newState)> callback)
        {
            setState = callback;
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

    private:
        localization lang{localization::EN_US};
        options& app_options;
        nlohmann::json localization_json{};

        void write_localization_json(localization lang);
        void initImages(window_t& window);

        // Menu

        sf::Texture placeholder_texture;
        uint32_t placeholder_texture_id = 0;

        std::vector<sf::Texture> thumb_textures;

        const float image_size = 200.f;

        int32_t currentDisplayScenario = 0;
        float currentDisplayTime = 0.0f;
        const float displayMaxTime = 10.f;
        void chooseNewDisplayScenario();
        void drawBackgroundDisplay();

        std::vector<scenario_info> backgroundDisplays;

        void drawSceneSelection();
        void drawSceneFrame(scenario_info& info, int32_t id);
        void drawSandboxCreation();
        void drawOptions();

        bool sceneSelectionOpen = false;
        bool tutorialSelectionOpen = false;
        bool sandboxSelectionOpen = false;
        bool challengeSelectionOpen = false;
        bool challengeViewOpen = false;
        bool optionsOpen = false;

        // Universe

        float target_temperature = 300.f;
        float target_pressure = 0.f;
        
        void runUniverse();
        void drawUniverseUI();
        
        uint32_t selectedCompound = UINT32_MAX;
        void drawCompoundSelector();
        void drawCompoundView(const compound_preset_info& compound);
        void drawCompoundFulLView();

        void pauseMenu();

        void initCompoundPresetsImages(window_t& window);
        void initCompoundPresets();
        std::vector<compound_preset_info> compound_presets{};

        bool compoundSelector = false;
        bool showHUD = true;

        bool pauseMenuOpen = false;
        bool exitDesktop = false;

        bool paused_simulation = false;

        std::unique_ptr<sim::fun::universe> simulation_universe;
        std::unique_ptr<sim::fun::universe> display_universe;

        // Callbacks
        
        std::function<void(application_state newState)> setState;

        // Video

        video* currentVideo = nullptr;
        bool rewinding = false;

        // Fonts
        ImFont* regular = nullptr;
        ImFont* medium = nullptr;
        ImFont* bold = nullptr;
        ImFont* black = nullptr;
    };
} // namespace core
