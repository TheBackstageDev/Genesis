#include "ScenarioHandler.hpp"

namespace core
{
    void ScenarioHandler::nextStep() 
    { 
        if (!isActive()) return;

        auto& current = m_scenario[m_currentStepIdx];
        if (current.onExit)
            current.onExit(*m_universe);

        ++m_currentStepIdx;

        if (isActive())
            enterStep(m_currentStepIdx); 
    }

    void ScenarioHandler::lastStep()
    { 
        if (m_currentStepIdx == 0) return;

        auto& current = m_scenario[m_currentStepIdx];
        if (current.onExit)
            current.onExit(*m_universe);

        --m_currentStepIdx;

        if (isActive())
            enterStep(m_currentStepIdx);
    }

    void ScenarioHandler::enterStep(size_t idx)
    {
        m_stepStartTime = 0.0f;

        m_narrationFinished = false;
        m_narrationProgress = 0.f;

        auto& step = m_scenario[idx];
        if (step.onEnter)
            step.onEnter(*m_universe);
    }

    void ScenarioHandler::updateVisuals()
    {
        float now = getCurrentTime();

        std::erase_if(m_activeVisuals, [&](const ActiveVisual& v)
        {
            return v.endTime > 0.0f && now >= v.endTime;
        });
    }

    void ScenarioHandler::update(float dt)
    {
        if (!m_narrationFinished) updateNarration(dt);
    }

    void ScenarioHandler::executeActions(const std::vector<ScenarioAction>& actions)
    {
        for (const ScenarioAction& action : actions)
        {
            if (action.customAction) action.customAction(*m_universe);
        }
    }

    void ScenarioHandler::updateNarration(float dt)
    {
        auto& currentStep = m_scenario[m_currentStepIdx];
        const std::string& current_narration = currentStep.narration;
    }

    void ScenarioHandler::loadVideo(const std::filesystem::path path)
    {
    
    }
} // namespace core
