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

        std::string text = "";
        uint32_t fromAtom = UINT32_MAX, toAtom = UINT32_MAX; // for atom position arrows
        uint32_t fromSubset = UINT32_MAX, toSubset = UINT32_MAX; // for subset position arrows
        uint32_t fromMolecule = UINT32_MAX, toMolecule = UINT32_MAX; // for Molecule position arrows
        std::vector<uint32_t> atomIndices; // For Highlight Atoms

        float duration_s = 0.0f; // 0.0f means it'll last until something happens or it is clicked away
        int32_t simulateSteps = 0;                // if > 0, run sim N steps during this action
        std::function<void(sim::fun::universe& universe)> customAction;
    };

    struct ScenarioStep
    {
        std::string narration; // Main Text
        std::vector<ScenarioAction> actions;

        bool autoAdvanceAfterNarration = false;
        float minDisplayTime_s = 1.8f;              // prevent too-fast skipping

        bool simulateDuringStep = false;
        float maxStepDuration_s = 0.0f;

        std::function<bool(sim::fun::universe& universe)> advanceWhen = nullptr;
        std::function<void(sim::fun::universe& universe)> onEnter = nullptr;
        std::function<void(sim::fun::universe& universe)> onExit = nullptr;
    };
        
    class ScenarioHandler
    {
    public:
        void update(float dt);
        
        void setCurrentUniverse(sim::fun::universe* newUniverse) 
        {
            m_universe = newUniverse;
            m_currentStepIdx = 0;
            m_stepStartTime = 0.f;
            m_narrationFinished = false;
        }

        void loadVideo(const std::filesystem::path path);
        void loadVideo(sim::fun::video* newVideo) { m_video = newVideo; }

        void clear()
        {
            m_universe = nullptr;
            m_scenario.clear();
            m_currentStepIdx = 0;
            m_stepStartTime = 0.f;
            m_narrationFinished = false;
        }

        bool isActive() const { return m_currentStepIdx > 0 && m_currentStepIdx < m_scenario.size(); }

        void nextStep();
        void lastStep();
        void skipToStep(size_t index) { enterStep(index); }

    private:
        void enterStep(size_t idx);
        void executeActions(const std::vector<ScenarioAction>& actions);
        void updateNarration(float dt);
        void updateVisuals();

        float getCurrentTime() { return std::chrono::duration<float>(std::chrono::system_clock::now().time_since_epoch()).count(); }

        using Scenario = std::vector<ScenarioStep>;
        Scenario m_scenario;

        struct ActiveVisual 
        {
            ScenarioAction action;
            float endTime = 0.0f; // m_stepStartTime + LifeTime, that is when will it'll stop appearing, unless 0.0f
        };
        std::vector<ActiveVisual> m_activeVisuals;

        bool m_narrationFinished = false;
        size_t m_currentStepIdx = 0;
        float m_stepStartTime = 0.f;
        float m_narrationProgress = 0.f;

        const float m_narration_interval = 0.1f; // 1.0 = second, per letter

        sim::fun::video* m_video = nullptr;        
        sim::fun::universe* m_universe = nullptr;
    };
} // namespace core
