#include "ScenarioHandler.hpp"

namespace core
{
    ScenarioHandler::ScenarioHandler(std::vector<sim::fun::compound_preset_info>& compounds)
        : m_compounds(compounds)
    {
        initScenarios();
    }

    ScenarioHandler::~ScenarioHandler()
    {
        clear();
    }

    void ScenarioHandler::initScenarios()
    {
        const std::filesystem::path scenario_path = "src/scenes/scenarios";
        Scenario scenario_brownian_motion{};
        scenario_brownian_motion.file = scenario_path / "SCENARIO_DYNAMICS_brownian_motion.json";
        scenario_brownian_motion.steps = 
        {
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 5.0f,
            },
            {
                .autoAdvanceAfterNarration = false,
                .minDisplayTime_s = 5.0f,
            }
        };
        
        m_Scenarios.emplace("SCENARIO_DYNAMICS_browinian_motion", std::move(scenario_brownian_motion));

        initTutorials();
    }

    void ScenarioHandler::initTutorials()
    {
        Scenario tutorial_engine_controls{};
        tutorial_engine_controls.steps = 
        {
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 3.0f,
            },  
            {
                .autoAdvanceAfterNarration = false,
                .minDisplayTime_s = 3.0f,
                .advanceWhen = [](sim::fun::universe& u) 
                { 
                    auto& cam = u.getRenderingCamera();

                    static glm::vec3 initialPosition = cam.eye();
                    static glm::vec3 initialTarget   = cam.target;
                    static float initialAzimuth      = cam.azimuth;
                    static float initialElevation    = cam.elevation;

                    static bool firstCall = true;
                    if (firstCall)
                    {
                        initialPosition = cam.eye();
                        initialTarget   = cam.target;
                        initialAzimuth  = cam.azimuth;
                        initialElevation = cam.elevation;
                        firstCall = false;
                    }

                    bool moved = glm::distance(cam.eye(), initialPosition) > 2.f &&
                                glm::distance(cam.target,   initialTarget)   > 2.f &&
                                std::abs(cam.azimuth   - initialAzimuth)   > 10.f &&
                                std::abs(cam.elevation - initialElevation) > 10.f;

                    if (moved)
                    {
                        firstCall = true;
                        return true;
                    }

                    return false;
                },
            },
            {
                .actions = 
                {
                    { .special = ScenarioAction::specialScenarioAction::PAUSE_SIMULATION },
                },
                .autoAdvanceAfterNarration = false,
                .minDisplayTime_s = 3.0f,
                .advanceWhen = [](sim::fun::universe& u) { return !u.isPaused(); },
                .onEnter = [](sim::fun::universe& u, std::vector<sim::fun::compound_preset_info>& compounds)
                {
                    u.createMolecule(compounds.at(6).structure, sf::Vector3f(u.boxSizes().x * 0.5f, u.boxSizes().y * 0.5f, u.boxSizes().z * 0.5f));
                }
            }, 
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 5.0f,
                .onEnter = [](sim::fun::universe& u, std::vector<sim::fun::compound_preset_info>& compounds)
                {
                    u.unpause();
                }
            },
            {
                .autoAdvanceAfterNarration = false,
            }
        };

        m_Scenarios.emplace("TUTORIAL_ENGINE_Controls", std::move(tutorial_engine_controls));
    }

    void ScenarioHandler::startScenario()
    {
        m_currentStepIdx = 0;
        m_stepStartTime = 0.0f;
        m_narrationFinished = false;
        m_activeVisuals.clear();

        enterStep(0);
    }

    void ScenarioHandler::nextStep() 
    { 
        if (!isActive()) return;

        auto& current = getCurrentStep();
        if (current.onExit)
            current.onExit(*m_universe);

        resetStepActions();

        ++m_currentStepIdx;

        if (isActive())
            enterStep(m_currentStepIdx); 
    }

    void ScenarioHandler::lastStep()
    { 
        if (m_currentStepIdx == 0) return;

        auto& current = getCurrentStep();
        if (current.onExit)
            current.onExit(*m_universe);

        resetStepActions();

        --m_currentStepIdx;

        if (isActive())
            enterStep(m_currentStepIdx);
    }

    void ScenarioHandler::enterStep(size_t idx)
    {
        m_stepStartTime = m_timeSinceScenarioStart;

        m_narrationFinished = false;
        m_narrationProgress = 0.f;

        auto& step = getCurrentStep();
        if (step.onEnter)
            step.onEnter(*m_universe, m_compounds);
    }

    void ScenarioHandler::resetStepActions()
    {
        auto& step = getCurrentStep();
        
        for (ScenarioAction& action : step.actions)
        {
            action.specialActionProcessed = false;
        }
    }

    void ScenarioHandler::updateVisuals()
    {
        float now = m_timeSinceScenarioStart;

        std::erase_if(m_activeVisuals, [&](const ActiveVisual& v)
        {
            return v.endTime > 0.0f && now >= v.endTime;
        });
    }

    void ScenarioHandler::update(float dt, nlohmann::json& localization)
    {
        m_timeSinceScenarioStart += dt;

        updateNarration(dt, localization["narrations"][m_currentStepIdx].get<std::string>());

        float timeInStep = m_timeSinceScenarioStart - m_stepStartTime;
        auto& current_step = getCurrentStep();
        
        if (current_step.advanceWhen && current_step.minDisplayTime_s <= timeInStep)
        {
            if (current_step.advanceWhen(*m_universe) && m_currentStepIdx < getCurrentScenarioStepCount() - 1) 
            {
                nextStep();
                return;
            }
        }

        executeActions(current_step.actions);
    }

    void ScenarioHandler::executeActions(std::vector<ScenarioAction>& actions)
    {
        for (ScenarioAction& action : actions)
        {
            if (action.customAction) action.customAction(*m_universe);

            if (!action.specialActionProcessed)
            {
                executeSpecialAction(action);
                action.specialActionProcessed = true;
            }
        }
    }

    void ScenarioHandler::executeSpecialAction(ScenarioAction& action)
    {
        using type = ScenarioAction::specialScenarioAction;

        switch (action.special)
        {
        case type::PAUSE_SIMULATION:
            m_universe->pause(); break;
        case type::RESUME_SIMULATION:
            m_universe->unpause(); break;
        case type::HIGHLIGHT:
        {
            for (uint32_t& index : action.atomIndices)
                m_universe->highlightAtom(index);
         
            if (action.toAtom && action.fromAtom)
                m_universe->highlightBond(action.fromAtom, action.toAtom);
                
            break;
        }
        case type::DRAW_ARROW:
        {
            if (action.toAtom && action.fromAtom)
                m_universe->createArrow(action.fromAtom, action.toAtom);

            break;
        }
        default:
            break;
        }
    }

    void ScenarioHandler::updateNarration(float dt, const std::string& narration)
    {
        auto& step = getCurrentStep();
        const std::string& fullText = narration;

        if (fullText.empty())
        {
            m_narrationFinished = true;
            return;
        }

        m_narrationProgress += dt / m_narration_interval;
        size_t visibleChars = static_cast<size_t>(m_narrationProgress);

        if (visibleChars >= fullText.size())
        {
            m_narrationFinished = true;
            visibleChars = fullText.size();
        }
        
        float timeInStep = m_timeSinceScenarioStart - m_stepStartTime;
        if (m_narrationFinished && step.autoAdvanceAfterNarration && 
            timeInStep >= step.minDisplayTime_s)
        {
            nextStep();
        }
    }

    void ScenarioHandler::draw(nlohmann::json& localization, float dt)
    {   
        auto& current_json = localization["Scenarios"][m_currentScenarioKey]; 
        
        drawNarration(current_json["narrations"][m_currentStepIdx].get<std::string>(), localization);
        update(dt, current_json);
    }

    void ScenarioHandler::drawNarration(const std::string& narration, nlohmann::json& localization)
    {
        if (!m_universe || !isActive())
            return;

        auto& step = getCurrentStep();
        if (narration.empty() && step.actions.empty())
            return;

            ImGuiIO& io = ImGui::GetIO();

        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(),
                                ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x * 0.25f, 400), ImGuiCond_Always);
        ImGui::SetNextWindowBgAlpha(0.80f);

        ImGuiWindowFlags flags = ImGuiWindowFlags_NoTitleBar
                            | ImGuiWindowFlags_NoResize
                            | ImGuiWindowFlags_NoScrollbar
                            | ImGuiWindowFlags_AlwaysAutoResize
                            | ImGuiWindowFlags_NoSavedSettings;

        auto& current_json = localization["Simulation"]["narration_menu"]; 

        if (ImGui::Begin("ScenarioNarration", nullptr, flags))
        {
            ImGui::BeginChild("NarrationContent", ImVec2(0, -80), true);

            float progress = static_cast<float>(m_currentStepIdx + 1) / getCurrentScenarioStepCount();
            ImGui::ProgressBar(progress, ImVec2(-1, 8));
            ImGui::Dummy(ImVec2(0, 8));

            ImGui::PushTextWrapPos(ImGui::GetContentRegionAvail().x);

            std::string visibleText = narration.substr(0, static_cast<size_t>(m_narrationProgress));
            ImGui::TextWrapped("%s", visibleText.c_str());

            ImGui::PopTextWrapPos();
            ImGui::EndChild();

            ImGui::Separator();
            ImGui::TextColored(ImVec4(0.7f, 0.9f, 1.0f, 1.0f),
                            "Step %zu / %zu",
                            m_currentStepIdx + 1,
                            getCurrentScenarioStepCount());

            if (step.autoAdvanceAfterNarration && m_narrationFinished)
            {
                float remaining = step.minDisplayTime_s - (getCurrentTime() - m_stepStartTime);
                if (remaining > 0)
                    ImGui::TextDisabled("Auto-advances in %.1fs", remaining);
            }

            ImGui::Dummy(ImVec2(0, 12));

            float buttonWidth = ImGui::GetContentRegionAvail().x * 0.48f;

            bool allow_previous = m_Scenarios[m_currentScenarioKey].allow_previous;
            if (m_currentStepIdx == 0 && !m_Scenarios[m_currentScenarioKey].allow_previous)
            {
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.3f, 0.3f, 0.3f, 0.8f));  // dark gray
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.3f, 0.3f, 0.3f, 0.8f));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.3f, 0.3f, 0.3f, 0.8f));
            }
            else
            {
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.2f, 0.4f, 0.6f, 1.0f));   // blue
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.3f, 0.5f, 0.7f, 1.0f));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.2f, 0.4f, 0.6f, 1.0f));   // blue
            }

            ImGui::BeginDisabled(!allow_previous);

            if (ImGui::Button(current_json["button_previous"].get<std::string>().c_str(), ImVec2(buttonWidth, 45)))
                if (m_currentStepIdx > 0) lastStep();

            ImGui::EndDisabled();

            ImGui::PopStyleColor(3);

            ImGui::SameLine();

            std::string button_next = "";

            if (m_currentStepIdx >= getCurrentScenarioStepCount() - 1)
            {
                button_next = current_json["button_exit"].get<std::string>().c_str();
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.7f, 0.2f, 0.2f, 0.8f));  // red
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.8f, 0.3f, 0.3f, 0.9f));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.7f, 0.2f, 0.2f, 0.8f));
            }
            else
            {
                button_next = current_json["button_next"].get<std::string>();
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.3f, 0.6f, 0.3f, 1.0f));   // green
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.3f, 0.6f, 0.3f, 1.0f));   // green
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.4f, 0.7f, 0.4f, 1.0f));   // green
            }

            if (ImGui::Button(button_next.c_str(), ImVec2(buttonWidth, 45)))
            {
                if (m_currentStepIdx < getCurrentScenarioStepCount() - 1) nextStep();
                else
                {
                    clear();
                    m_exitFlag = true;
                }
            }

            ImGui::PopStyleColor(3);
            ImGui::End();
        }
    }

    void ScenarioHandler::loadVideo(const std::filesystem::path path)
    {
    
    }
} // namespace core
