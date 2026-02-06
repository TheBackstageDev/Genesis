#include "ScenarioHandler.hpp"

#include <random>

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
        const std::filesystem::path scenario_path = "scenes/scenarios";
        Scenario SCENARIO_DYNAMICS_BrownianMotion{};
        SCENARIO_DYNAMICS_BrownianMotion.close_not_exit = true;
        SCENARIO_DYNAMICS_BrownianMotion.file = scenario_path / "SCENARIO_DYNAMICS_browinian_motion.json";
        SCENARIO_DYNAMICS_BrownianMotion.steps = 
        {
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 5.0f,
                .onEnter = [&](sim::fun::universe& u, std::vector<sim::fun::compound_preset_info>& compounds)
                {
                    u.setTimescale(3.f);
                    m_wantedTemperature = 30.f;
                },
            },
            {
                .autoAdvanceAfterNarration = false,
                .minDisplayTime_s = 5.0f,
            }
        };

        Scenario SCENARIO_DYNAMICS_CrystalNucleation{};
        SCENARIO_DYNAMICS_CrystalNucleation.close_not_exit = true;
        SCENARIO_DYNAMICS_CrystalNucleation.file = scenario_path / "SCENARIO_DYNAMICS_CrystalNucleation.json";
        SCENARIO_DYNAMICS_CrystalNucleation.steps =
        {
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 6.0f,
                .onEnter = [&](sim::fun::universe& u, std::vector<sim::fun::compound_preset_info>& compounds)
                {
                    u.setTimescale(200.f);
                    m_wantedTemperature = 200.f;
                }
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
                .onEnter = [&](sim::fun::universe& u, std::vector<sim::fun::compound_preset_info>& compounds)
                {
                    u.setTimescale(500.f);
                    m_wantedTemperature = 100.f;
                }
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .autoAdvanceAfterNarration = false,
                .minDisplayTime_s = 10.0f,  
            }
        };
        
        m_Scenarios.emplace("SCENARIO_DYNAMICS_browinian_motion", std::move(SCENARIO_DYNAMICS_BrownianMotion));
        m_Scenarios.emplace("SCENARIO_DYNAMICS_CrystalNucleation", std::move(SCENARIO_DYNAMICS_CrystalNucleation));

        initTutorials();
    }

    void ScenarioHandler::initTutorials()
    {
        const std::filesystem::path scenario_path = "scenes/scenarios";

        Scenario TUTORIAL_ENGINE_Controls{};
        TUTORIAL_ENGINE_Controls.steps = 
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
                .autoAdvanceAfterNarration = false,
                .minDisplayTime_s = 3.0f,
                .advanceWhen = [](sim::fun::universe& u) { return !u.isPaused(); },
                .onEnter = [](sim::fun::universe& u, std::vector<sim::fun::compound_preset_info>& compounds)
                {
                    u.pause();
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

        Scenario TUTORIAL_CLASSICAL_VanDerWaals{};
        TUTORIAL_CLASSICAL_VanDerWaals.file = scenario_path / "TUTORIAL_CLASSICAL_VDW.json";
        TUTORIAL_CLASSICAL_VanDerWaals.steps = 
        {
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 5.0f,
                .onEnter = [&](sim::fun::universe& u, std::vector<sim::fun::compound_preset_info>& compounds)
                {
                    u.pause();
                    u.setTimescale(500.f);
                    m_wantedTemperature = 2.0f;
                }
            },
            {
                .actions = 
                {
                    /* { 
                        .special = ScenarioAction::specialScenarioAction::DRAW_ARROW, 
                        .fromAtom = 0, .toAtom = 1
                    } */
                },
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 5.0f,
                .onExit =  [&](sim::fun::universe& u)
                {
                    u.unpause();
                }
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 20.0f,
            },
            {
                .actions =
                {
                    {
                        .customAction = [&](sim::fun::universe& u)
                        {
                            static int32_t numAdded = 0;
                            const int32_t targetCount = 50;

                            if (numAdded >= targetCount) return;

                            const glm::vec3 center = (u.getPosition(0) + u.getPosition(1)) * 0.5f;

                            std::random_device rd;
                            std::mt19937 gen(rd());
                            std::uniform_real_distribution<float> distPos(-12.0f, 12.0f);
                            std::uniform_real_distribution<float> distVel(-0.1f, 0.1f);

                            auto positions = u.positions();

                            for (int32_t attempt = 0; numAdded < targetCount && attempt < 300; ++attempt)
                            {
                                glm::vec3 candidate = center + glm::vec3(distPos(gen), distPos(gen), distPos(gen));

                                bool tooClose = false;
                                float minDist = 9999.0f;
                                for (const auto& p : positions)
                                {
                                    float d = glm::length(p - candidate);
                                    if (d < 2.8f)
                                    {
                                        tooClose = true;
                                        break;
                                    }
                                    minDist = std::min(minDist, d);
                                }

                                if (!tooClose)
                                {
                                    u.createAtom(
                                        candidate,
                                        glm::vec3(distVel(gen), distVel(gen), distVel(gen)),
                                        18,
                                        18,
                                        18,
                                        0
                                    );

                                    ++numAdded;
                                    positions = u.positions();
                                }
                            }
                        }
                    }
                },
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
                .onEnter = [&](sim::fun::universe& u, std::vector<sim::fun::compound_preset_info>& compounds)
                {
                    u.setTimescale(500.f);
                }
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
                .onEnter = [&](sim::fun::universe& u, std::vector<sim::fun::compound_preset_info>& compounds)
                {
                    u.setTimescale(30.f);
                    m_wantedTemperature = 120.f;
                }
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
                .onEnter = [&](sim::fun::universe& u, std::vector<sim::fun::compound_preset_info>& compounds)
                {
                    u.setTimescale(500.f);
                    m_wantedTemperature = 2.0f;
                }
            },
            {
                .autoAdvanceAfterNarration = false,
            },
        };

        Scenario TUTORIAL_CLASSICAL_CoulombForce{};
        TUTORIAL_CLASSICAL_CoulombForce.steps = 
        {
        };

        m_Scenarios.emplace("TUTORIAL_ENGINE_Controls", std::move(TUTORIAL_ENGINE_Controls));
        m_Scenarios.emplace("TUTORIAL_CLASSICAL_VanDerWaals", std::move(TUTORIAL_CLASSICAL_VanDerWaals));
        m_Scenarios.emplace("TUTORIAL_CLASSICAL_CoulombForce", std::move(TUTORIAL_CLASSICAL_CoulombForce));
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
            if (action.toAtom != UINT32_MAX && action.fromAtom != UINT32_MAX)
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

            bool close = m_Scenarios[m_currentScenarioKey].close_not_exit;
            if (m_currentStepIdx >= getCurrentScenarioStepCount() - 1)
            {
                std::string exit_str = close ? current_json["button_close"] : current_json["button_exit"];
                button_next = exit_str.c_str();
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
                    if (!close)
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
