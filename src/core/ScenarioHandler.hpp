#pragma once

#include "simulation/universe.hpp"
#include <functional>

#include <filesystem>
#include <chrono>

namespace core
{
    struct ScenarioAction
    {
        enum class specialScenarioAction
        {
            NONE,
            PAUSE_SIMULATION,
            RESUME_SIMULATION,
            DRAW_ARROW,
            HIGHLIGHT
        };

        specialScenarioAction special = specialScenarioAction::NONE;

        // Note: Text will be displayed on the position of an atom, if text is present and the indices are set, to either be between them or 
        // To be on one of them;

        bool specialActionProcessed = false;
        std::string text = "";
        uint32_t fromAtom = UINT32_MAX, toAtom = UINT32_MAX; // for atom position arrows
        uint32_t fromMolecule = UINT32_MAX, toMolecule = UINT32_MAX; // for Molecule position arrows
        std::vector<uint32_t> atomIndices; // For Highlight Atoms

        float duration_s = 0.0f; // 0.0f means it'll last until something happens or it is clicked away
        int32_t simulateSteps = 0;                // if > 0, run sim N steps during this action
        std::function<void(sim::fun::universe& universe)> customAction;
    };

    struct ScenarioStep
    {
        std::vector<ScenarioAction> actions;

        bool autoAdvanceAfterNarration = false;
        float minDisplayTime_s = 1.8f;              // prevent too-fast skipping

        bool simulateDuringStep = false;
        float maxStepDuration_s = 0.0f;

        std::function<bool(sim::fun::universe& universe)> advanceWhen = nullptr;
        std::function<void(sim::fun::universe& universe, std::vector<sim::fun::compound_preset_info>& compounds)> onEnter = nullptr;
        std::function<void(sim::fun::universe& universe)> onExit = nullptr;
    };

    struct Scenario
    {
        std::vector<ScenarioStep> steps;

        std::filesystem::path file = "";
        sim::fun::rendering_info m_render_info{};

        bool allow_previous = false;
        bool allow_compoundSelector = false;
        bool close_not_exit = false;
    };
        
    class ScenarioHandler
    {
    public:
        ScenarioHandler(std::vector<sim::fun::compound_preset_info>& compounds);
        ~ScenarioHandler();

        void chooseScenario(const std::string& scenario)
        {
            m_currentScenarioKey = scenario;
        }

        void update(float dt, nlohmann::json& localization);
        
        void setCurrentUniverse(sim::fun::universe* newUniverse) 
        {
            m_universe = newUniverse;
            m_currentStepIdx = 0;
            m_stepStartTime = 0.f;
            m_narrationFinished = false;
        }

        void loadVideo(const std::filesystem::path path);
        void loadVideo(sim::fun::video* newVideo) { m_video = newVideo; }

        void restart()
        {
            enterStep(0);
        }

        void clear()
        {
            m_universe = nullptr;
            m_video = nullptr;
            m_currentScenarioKey.clear();
            m_currentStepIdx = 0;
            m_stepStartTime = 0.f;
            m_exitFlag = false;
            m_narrationFinished = false;
        }

        bool isActive() const 
        {
            auto it = m_Scenarios.find(m_currentScenarioKey);
            return it != m_Scenarios.end() && m_currentStepIdx >= 0 && m_currentStepIdx < it->second.steps.size();
        }

        bool inScenario() { return !m_currentScenarioKey.empty(); }
        void draw(nlohmann::json& localization, float dt);
        void startScenario();

        void nextStep();
        void lastStep();
        void skipToStep(size_t index) { enterStep(index); }

        const std::map<std::string, Scenario>& getScenarios() const { return m_Scenarios; }

        bool exit() { return m_exitFlag; }
        bool allowedCompoundSelector() 
        { 
            if (!m_currentScenarioKey.empty())
                return m_Scenarios[m_currentScenarioKey].allow_compoundSelector;

            return true;
        }

        // Gets

        float getWantedTemperature() const { return m_wantedTemperature; }
    private:
        friend ScenarioStep;

        void initScenarios();
        void initTutorials();

        void enterStep(size_t idx);
        void executeActions(std::vector<ScenarioAction>& actions);
        void executeSpecialAction(ScenarioAction& action);

        void drawNarration(const std::string& narration, nlohmann::json& localization);
        void updateNarration(float dt, const std::string& narration);
        void updateVisuals();

        void resetStepActions();

        const ScenarioStep& getCurrentStep() const
        {
            auto it = m_Scenarios.find(m_currentScenarioKey);
            if (it == m_Scenarios.end() || m_currentStepIdx >= it->second.steps.size())
            {
                static const ScenarioStep dummy{};
                return dummy;
            }
            return it->second.steps[m_currentStepIdx];
        }

        ScenarioStep& getCurrentStep()
        {
            return const_cast<ScenarioStep&>(std::as_const(*this).getCurrentStep());
        }
        
        size_t getCurrentScenarioStepCount() const
        {
            auto it = m_Scenarios.find(m_currentScenarioKey);
            return it != m_Scenarios.end() ? it->second.steps.size() : 0;
        }

        float getCurrentTime() { return std::chrono::duration<float>(std::chrono::system_clock::now().time_since_epoch()).count(); }
        float getTimeWhen(float future = 0.0f) { return std::chrono::duration<float>(std::chrono::system_clock::now().time_since_epoch()).count() + future; }

        std::map<std::string, Scenario> m_Scenarios;

        std::string m_currentScenarioKey = "";

        struct ActiveVisual 
        {
            ScenarioAction action;
            float endTime = 0.0f; // m_stepStartTime + LifeTime, that is when will it'll stop appearing, unless 0.0f
        };
        std::vector<ActiveVisual> m_activeVisuals;

        bool m_narrationFinished = false;
        bool m_exitFlag = false;
        size_t m_currentStepIdx = 0;
        
        float m_timeSinceScenarioStart = 0.0;
        float m_stepStartTime = 0.f;
        float m_narrationProgress = 0.f;

        const float m_narration_interval = 0.01f; // 1.0 = second, per letter

        sim::fun::video* m_video = nullptr;        
        sim::fun::universe* m_universe = nullptr;

        std::vector<sim::fun::compound_preset_info>& m_compounds;
        
        // "Flags"
        float m_wantedTemperature = 300.f;
    };
} // namespace core
