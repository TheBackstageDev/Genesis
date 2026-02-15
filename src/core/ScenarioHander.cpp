#include "ScenarioHandler.hpp"

#include "simulation/smiles_parser.hpp"
#include <random>

namespace core
{
    ScenarioHandler::ScenarioHandler(std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds, sim::simulation_packer& simpacker)
        : m_compounds(compounds), m_simpacker(simpacker)
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
        SCENARIO_DYNAMICS_BrownianMotion.file = scenario_path / "SCENARIO_DYNAMICS_browinian_motion.json";
        SCENARIO_DYNAMICS_BrownianMotion.steps =
            {
                {
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 5.0f,
                    .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                    {
                        u.setTimescale(10.f);
                        m_wantedTemperature = 3000.f;
                    },
                },
                {
                    .autoAdvanceAfterNarration = false,
                    .minDisplayTime_s = 5.0f,
                }};

        Scenario SCENARIO_DYNAMICS_CrystalNucleation{};
        SCENARIO_DYNAMICS_CrystalNucleation.file = scenario_path / "SCENARIO_DYNAMICS_CrystalNucleation.json";
        SCENARIO_DYNAMICS_CrystalNucleation.steps =
            {
                {.autoAdvanceAfterNarration = true,
                 .minDisplayTime_s = 6.0f,
                 .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                 {
                     u.setTimescale(200.f);
                     m_wantedTemperature = 200.f;
                 }},
                {
                    .autoAdvanceAfterNarration = true, 
                    .minDisplayTime_s = 15.0f, 
                    .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                    {
                        u.setTimescale(400.f);
                        m_wantedTemperature = 90.f;
                    }
                },
                {
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 15.0f,
                },
                {
                    .autoAdvanceAfterNarration = false,
                    .minDisplayTime_s = 10.0f,
                }};

        Scenario SCENARIO_DYNAMICS_DropletFormation{};
        SCENARIO_DYNAMICS_DropletFormation.file = scenario_path / "SCENARIO_DYNAMICS_DropletFormation.json";
        SCENARIO_DYNAMICS_DropletFormation.steps = 
        {
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 12.0f,
                .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                {
                    u.setTimescale(200.f);
                    m_wantedTemperature = 1000.f;
                },
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 25.0f,
                .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                {
                    u.setTimescale(750.f);
                    m_wantedTemperature = 40.f;
                },
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 25.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 20.0f,
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
        };

        Scenario SCENARIO_DYNAMICS_AtmosphereLayers{};
        SCENARIO_DYNAMICS_AtmosphereLayers.file = scenario_path / "SCENARIO_DYNAMICS_AtmosphereLayers.json";
        SCENARIO_DYNAMICS_AtmosphereLayers.steps = 
        {
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
                .onEnter = [&](sim::fun::universe& u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                {
                    u.pause();
                    u.setTimescale(1.f);

                    auto argon = sim::parseSMILES("Ar", false);

                    std::vector<sim::fun::molecule_structure> molecules =
                    {
                        compounds["Nitrogen"].structure,
                        compounds["Oxygen"].structure,
                        argon,
                        compounds["Carbon Dioxide"].structure,
                        compounds["Methane"].structure
                    };

                    std::vector<float> chances =
                    {
                        0.78f,
                        0.21f,
                        0.0093f,
                        0.004f,
                        0.0001f
                    };

                    m_simpacker.pack(u, molecules, chances, u.boxSizes() * 0.5f, u.boxSizes());
                }
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
                .onEnter = [&](sim::fun::universe& u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                {
                    u.setTimescale(2.f);
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
                .autoAdvanceAfterNarration = false,
                .minDisplayTime_s = 15.0f,
            }
        };

        m_Scenarios.emplace("SCENARIO_DYNAMICS_browinian_motion", std::move(SCENARIO_DYNAMICS_BrownianMotion));
        m_Scenarios.emplace("SCENARIO_DYNAMICS_DropletFormation", std::move(SCENARIO_DYNAMICS_DropletFormation));
        m_Scenarios.emplace("SCENARIO_DYNAMICS_CrystalNucleation", std::move(SCENARIO_DYNAMICS_CrystalNucleation));
        m_Scenarios.emplace("SCENARIO_DYNAMICS_AtmosphereLayers", std::move(SCENARIO_DYNAMICS_AtmosphereLayers));

        initMoleculeScenarios();
        initTutorials();
    }

    void ScenarioHandler::initMoleculeScenarios()
    {
        Scenario SCENARIO_MOLECULE_Water{};
        SCENARIO_MOLECULE_Water.steps = 
        {
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,

                .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                {
                    const sf::Vector3f center{u.boxSizes().x * 0.5f, u.boxSizes().y * 0.5f, u.boxSizes().z * 0.5f};

                    u.setTimescale(1.5f);
                    u.createMolecule(compounds["Water"].structure, center);
                },
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .actions = 
                {
                    {
                        .action = [&](sim::fun::universe& u)
                        {
                            u.getRenderingEngine().drawAngle(ImGui::GetBackgroundDrawList(), u.getPosition(0), u.getPosition(1), u.getPosition(2));
                        }
                    }
                },

                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .actions = 
                {
                    {
                        .action = [&](sim::fun::universe& u)
                        {
                            u.getRenderingEngine().drawAngle(ImGui::GetBackgroundDrawList(), u.getPosition(0), u.getPosition(1), u.getPosition(2));
                        }
                    }
                },
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .actions = 
                {
                    {
                        .action = [&](sim::fun::universe& u)
                        {
                            u.getRenderingEngine().drawAngle(ImGui::GetBackgroundDrawList(), u.getPosition(0), u.getPosition(1), u.getPosition(2));
                        }
                    }
                },

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
                .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                {   
                    m_simpacker.pack(u, compounds["Water"].structure, u.boxSizes() * 0.5f, u.boxSizes() * 0.3f, 5);
                },
            },
            {
                .autoAdvanceAfterNarration = false,
                .minDisplayTime_s = 15.0f,
            },
        };

        Scenario SCENARIO_MOLECULE_Benzene{};
        SCENARIO_MOLECULE_Benzene.steps = 
        {
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,

                .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                {
                    const sf::Vector3f center{u.boxSizes().x * 0.5f, u.boxSizes().y * 0.5f, u.boxSizes().z * 0.5f};

                    u.setTimescale(1.f);
                    u.createMolecule(compounds.at("Benzene").structure, center);
                },
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .actions = 
                {
                    {
                        .action = [&](sim::fun::universe& u)
                        {
                            u.getRenderingEngine().drawAngle(ImGui::GetBackgroundDrawList(), u.getPosition(0), u.getPosition(1), u.getPosition(2));
                        }
                    }
                },
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .actions = 
                {
                    {
                        .action = [&](sim::fun::universe& u)
                        {
                            u.getRenderingEngine().drawAngle(ImGui::GetBackgroundDrawList(), u.getPosition(0), u.getPosition(1), u.getPosition(2));
                        }
                    }
                },
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .actions = 
                {
                    {
                        .action = [&](sim::fun::universe& u)
                        {
                            u.getRenderingEngine().drawAngle(ImGui::GetBackgroundDrawList(), u.getPosition(0), u.getPosition(1), u.getPosition(2));
                        }
                    }
                },
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .actions = 
                {
                    {
                        .action = [&](sim::fun::universe& u)
                        {
                            u.getRenderingEngine().drawAngle(ImGui::GetBackgroundDrawList(), u.getPosition(0), u.getPosition(1), u.getPosition(2));
                        }
                    }
                },
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .autoAdvanceAfterNarration = false,
                .minDisplayTime_s = 15.0f,
            },
        };

        Scenario SCENARIO_MOLECULE_CarbonDioxide{};
        SCENARIO_MOLECULE_CarbonDioxide.steps = 
        {
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,

                .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                {
                    const sf::Vector3f center{u.boxSizes().x * 0.5f, u.boxSizes().y * 0.5f, u.boxSizes().z * 0.5f};

                    u.setTimescale(1.f);
                    u.createMolecule(compounds.at("Carbon Dioxide").structure, center);
                },
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .actions = 
                {
                    {
                        .action = [&](sim::fun::universe& u)
                        {
                            u.getRenderingEngine().drawAngle(ImGui::GetBackgroundDrawList(), u.getPosition(1), u.getPosition(0), u.getPosition(2));
                        }
                    }
                },
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .actions = 
                {
                    {
                        .action = [&](sim::fun::universe& u)
                        {
                            u.getRenderingEngine().drawAngle(ImGui::GetBackgroundDrawList(), u.getPosition(1), u.getPosition(0), u.getPosition(2));
                        }
                    }
                },
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .actions = 
                {
                    {
                        .action = [&](sim::fun::universe& u)
                        {
                            u.getRenderingEngine().drawAngle(ImGui::GetBackgroundDrawList(), u.getPosition(1), u.getPosition(0), u.getPosition(2));
                        }
                    }
                },
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .actions = 
                {
                    {
                        .action = [&](sim::fun::universe& u)
                        {
                            u.getRenderingEngine().drawAngle(ImGui::GetBackgroundDrawList(), u.getPosition(1), u.getPosition(0), u.getPosition(2));
                        }
                    }
                },
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 15.0f,
            },
            {
                .autoAdvanceAfterNarration = false,
                .minDisplayTime_s = 15.0f,
            },
        };

        m_Scenarios.emplace("SCENARIO_MOLECULE_Water", std::move(SCENARIO_MOLECULE_Water));
        m_Scenarios.emplace("SCENARIO_MOLECULE_CarbonDioxide", std::move(SCENARIO_MOLECULE_CarbonDioxide));
        m_Scenarios.emplace("SCENARIO_MOLECULE_Benzene", std::move(SCENARIO_MOLECULE_Benzene));
    }

    void ScenarioHandler::initTutorials()
    {
        const std::filesystem::path scenario_path = "scenes/scenarios";

        Scenario TUTORIAL_ENGINE_Controls{};
        TUTORIAL_ENGINE_Controls.steps =
            {
                {
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 10.0f,
                },
                {
                    .autoAdvanceAfterNarration = false,
                    .minDisplayTime_s = 3.0f,
                    .advanceWhen = [&](sim::fun::universe &u)
                    {
                        auto &cam = u.getRenderingCamera();

                        static glm::vec3 initialPosition = cam.eye();
                        static glm::vec3 initialTarget = cam.target;
                        static float initialAzimuth = cam.azimuth;
                        static float initialElevation = cam.elevation;

                        static bool firstCall = true;
                        if (firstCall)
                        {
                            initialPosition = cam.eye();
                            initialTarget = cam.target;
                            initialAzimuth = cam.azimuth;
                            initialElevation = cam.elevation;
                            firstCall = false;
                        }

                        bool moved = glm::distance(cam.eye(), initialPosition) > 2.f &&
                                     glm::distance(cam.target, initialTarget) > 2.f &&
                                     std::abs(cam.azimuth - initialAzimuth) > 10.f &&
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
                    .minDisplayTime_s = 1.0f,
                    .advanceWhen = [&](sim::fun::universe& u)
                    {
                        return u.numAtoms() > 0;
                    }
                },
                {
                    .autoAdvanceAfterNarration = false,
                    .minDisplayTime_s = 1.0f,
                    .advanceWhen = [&](sim::fun::universe& u)
                    {
                        return u.isPaused();
                    }
                },
                {
                    .autoAdvanceAfterNarration = false,
                    .minDisplayTime_s = 1.0f,
                    .advanceWhen = [&](sim::fun::universe& u)
                    {
                        return !u.isPaused();
                    }
                },
                {
                    .autoAdvanceAfterNarration = false,
                    .minDisplayTime_s = 1.0f,
                    .advanceWhen = [&](sim::fun::universe& u)
                    {
                        return u.getTimescale() > 1.5f || u.getTimescale() < 0.5f;
                    }
                },
                {
                    .autoAdvanceAfterNarration = false,
                    .minDisplayTime_s = 1.0f,
                }
            };

        Scenario TUTORIAL_CLASSICAL_VanDerWaals{};
        TUTORIAL_CLASSICAL_VanDerWaals.file = scenario_path / "TUTORIAL_CLASSICAL_VDW.json";
        TUTORIAL_CLASSICAL_VanDerWaals.steps =
            {
                {
                 .autoAdvanceAfterNarration = true,
                 .minDisplayTime_s = 5.0f,
                 .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                 {
                    u.pause();
                    u.setTimescale(500.f);
                    m_wantedTemperature = 2.0f;
                 }
                },
                {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 5.0f,
                .onExit = [&](sim::fun::universe &u)
                {
                    u.unpause();
                }},
                {
                    .actions =
                    {
                        {
                            .action = [&](sim::fun::universe& u)
                            {
                                u.clearArrows();

                                for (size_t i = 0; i < u.numAtoms(); ++i) 
                                {
                                    glm::vec3 pos = u.getPosition(i);
                                    glm::vec3 force = u.getForce(i);

                                    float scale = 1.5f;  

                                    u.createArrow(pos, pos + scale * force);
                                }
                            }
                        }
                    },
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 15.0f,
                    .onExit = [&](sim::fun::universe& u)
                    {
                        u.clearArrows();
                    }
                },
                {
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 20.0f,
                },
                {
                 .autoAdvanceAfterNarration = true,
                 .minDisplayTime_s = 10.0f,
                 .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                 {
                    u.setTimescale(500.f);
                    auto argon = sim::parseSMILES("Ar", false);

                    m_simpacker.pack(u, argon, u.boxSizes() * 0.5f, u.boxSizes(), 50);
                 }},
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
                {.autoAdvanceAfterNarration = true, .minDisplayTime_s = 10.0f, .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                                                                               {
                                                                                   u.setTimescale(30.f);
                                                                                   m_wantedTemperature = 120.f;
                                                                               }},
                {.autoAdvanceAfterNarration = true, .minDisplayTime_s = 10.0f, .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                                                                               {
                                                                                   u.setTimescale(500.f);
                                                                                   m_wantedTemperature = 2.0f;
                                                                               }},
                {
                    .autoAdvanceAfterNarration = false,
                },
            };

        Scenario TUTORIAL_CLASSICAL_CoulombForce{};
        TUTORIAL_CLASSICAL_CoulombForce.file = scenario_path / "TUTORIAL_CLASSICAL_COULOMB.json";
        TUTORIAL_CLASSICAL_CoulombForce.steps =
            {
                {.autoAdvanceAfterNarration = true,
                 .minDisplayTime_s = 15.0f,
                 .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                 {
                     u.pause();
                     u.setTimescale(10.f);
                     m_wantedTemperature = 250.f;
                 }},
                {
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 10.0f,
                },
                {
                    .actions =
                    {
                        {
                            .action = [&](sim::fun::universe& u)
                            {
                                u.clearArrows();

                                for (size_t i = 0; i < u.numAtoms(); ++i) 
                                {
                                    glm::vec3 pos = u.getPosition(i);
                                    glm::vec3 force = u.getForce(i);

                                    float scale = 0.005f;  

                                    u.createArrow(pos, pos + scale * force);
                                }
                            }
                        }
                    },
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 15.0f,
                    .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                    { 
                        u.unpause(); 
                    },
                    .onExit = [&](sim::fun::universe &u)
                    {
                        u.pause();
                        u.clear(); 
                    },
                },
                {
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 5.0f,
                    .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                    {
                        const glm::vec3 center = u.boxSizes() * 0.5f;

                        u.createAtom(center - glm::vec3(1.f, 0.f, 0.f), glm::vec3(1.f, 0.f, 0.f), 11, 11, 10, 0);
                        u.createAtom(center + glm::vec3(1.f, 0.f, 0.f), glm::vec3(-1.f, 0.f, 0.f), 11, 11, 10, 0);
                    },
                },
                {
                    .actions =
                    {
                        {
                            .action = [&](sim::fun::universe& u)
                            {
                                u.clearArrows();

                                for (size_t i = 0; i < u.numAtoms(); ++i) 
                                {
                                    glm::vec3 pos = u.getPosition(i);
                                    glm::vec3 force = u.getForce(i);

                                    float scale = 0.005f;  

                                    u.createArrow(pos, pos + scale * force);
                                }
                            }
                        }
                    },
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 25.0f,
                    .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                    {
                        u.unpause();
                        u.clearArrows();
                    },
                },
                {
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 25.0f,
                    .onExit = [&](sim::fun::universe &u)
                    {
                        u.clear();
                    },
                },
                {
                    .actions =
                    {
                        {
                            .action = [&](sim::fun::universe& u)
                            {
                                u.clearArrows();

                                for (size_t i = 0; i < u.numAtoms(); ++i) 
                                {
                                    glm::vec3 pos = u.getPosition(i);
                                    glm::vec3 force = u.getForce(i);

                                    float scale = 0.02f;  

                                    u.createArrow(pos, pos + scale * force);
                                }
                            }
                        }
                    },
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 25.0f,
                    .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                    {
                        u.setTimescale(1.f);
                        m_wantedTemperature = 100.f;

                        const glm::vec3 center = u.boxSizes() * 0.5f;

                        const float spacing = 3.f;
                        const int nx = 2;
                        const int ny = 2;
                        const int nz = 2;

                        for (int ix = 0; ix <= nx; ++ix)
                        {
                            for (int iy = 0; iy <= ny; ++iy)
                            {
                                for (int iz = 0; iz <= nz; ++iz)
                                {
                                    glm::vec3 offset(
                                        (ix - nx / 2.0f) * spacing,
                                        (iy - ny / 2.0f) * spacing,
                                        (iz - nz / 2.0f) * spacing);

                                    glm::vec3 pos = center + offset;

                                    bool is_sodium = ((ix + iy + iz) % 2 == 0);

                                    if (is_sodium)
                                    {
                                        u.createAtom(pos, glm::vec3(0), 11, 11, 10, 0);
                                    }
                                    else
                                    {
                                        u.createAtom(pos, glm::vec3(0), 17, 18, 18, 0);
                                    }
                                }
                            }
                        }

                        m_wantedTemperature = 200.f;
                    },
                },
                {
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 25.0f,
                    .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                    {
                        m_wantedTemperature = 350.f;

                        u.clearArrows();
                        u.setTimescale(3.f);
                        m_simpacker.pack(u, compounds["Water"].structure, u.boxSizes() * 0.5f, u.boxSizes() * 0.3f, 50);
                    }
                },
                {
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 30.0f,
                },
                {
                    .autoAdvanceAfterNarration = true,
                    .minDisplayTime_s = 30.0f,
                },
                {
                    .autoAdvanceAfterNarration = false,
                    .minDisplayTime_s = 10.0f,
                }
            };

        Scenario TUTORIAL_CLASSICAL_Phases{};
        TUTORIAL_CLASSICAL_Phases.file = scenario_path / "TUTORIAL_CLASSICAL_Phases.json";
        TUTORIAL_CLASSICAL_Phases.video = scenario_path / "video_TUTORIAL_CLASSICAL_Phases.json";
        TUTORIAL_CLASSICAL_Phases.steps =
        {
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
                .onEnter = [&](sim::fun::universe &u, std::unordered_map<std::string, sim::fun::compound_preset_info> &compounds)
                {
                    u.pause();
                    m_wantedVideoPlay = true;
                    m_wantedFramesPerSecond = 30.f;
                }
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
            },
            {
                .autoAdvanceAfterNarration = true,
                .minDisplayTime_s = 10.0f,
            },
            {
                .autoAdvanceAfterNarration = false,
                .minDisplayTime_s = 10.0f,
            },
        };

        m_Scenarios.emplace("TUTORIAL_ENGINE_Controls", std::move(TUTORIAL_ENGINE_Controls));
        m_Scenarios.emplace("TUTORIAL_CLASSICAL_VanDerWaals", std::move(TUTORIAL_CLASSICAL_VanDerWaals));
        m_Scenarios.emplace("TUTORIAL_CLASSICAL_CoulombForce", std::move(TUTORIAL_CLASSICAL_CoulombForce));
        m_Scenarios.emplace("TUTORIAL_CLASSICAL_Phases", std::move(TUTORIAL_CLASSICAL_Phases));
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
        if (!isActive())
            return;

        auto &current = getCurrentStep();
        if (current.onExit)
            current.onExit(*m_universe);

        resetStepActions();

        ++m_currentStepIdx;

        if (isActive())
            enterStep(m_currentStepIdx);
    }

    void ScenarioHandler::lastStep()
    {
        if (m_currentStepIdx == 0)
            return;

        auto &current = getCurrentStep();
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

        auto &step = getCurrentStep();
        if (step.onEnter)
            step.onEnter(*m_universe, m_compounds);
    }

    void ScenarioHandler::resetStepActions()
    {
        auto &step = getCurrentStep();
    }

    void ScenarioHandler::updateVisuals()
    {
        float now = m_timeSinceScenarioStart;

        std::erase_if(m_activeVisuals, [&](const ActiveVisual &v)
                      { return v.endTime > 0.0f && now >= v.endTime; });
    }

    void ScenarioHandler::update(float dt, nlohmann::json &localization)
    {
        m_timeSinceScenarioStart += dt;

        updateNarration(dt, localization["narrations"][m_currentStepIdx].get<std::string>());

        float timeInStep = m_timeSinceScenarioStart - m_stepStartTime;
        auto &current_step = getCurrentStep();

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

    void ScenarioHandler::executeActions(std::vector<ScenarioAction> &actions)
    {
        for (ScenarioAction &action : actions)
        {
            if (action.action)
                action.action(*m_universe);
        }
    }

    void ScenarioHandler::updateNarration(float dt, const std::string &narration)
    {
        auto &step = getCurrentStep();
        const std::string &fullText = narration;

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

    void ScenarioHandler::draw(nlohmann::json &localization, float dt)
    {
        auto &current_json = localization["Scenarios"][m_currentScenarioKey];

        drawNarration(current_json["narrations"][m_currentStepIdx].get<std::string>(), localization);
        update(dt, current_json);
    }

    void ScenarioHandler::drawNarration(const std::string &narration, nlohmann::json &localization)
    {
        if (!m_universe || !isActive())
            return;

        auto &step = getCurrentStep();
        if (narration.empty() && step.actions.empty())
            return;

        ImGuiIO &io = ImGui::GetIO();

        ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(),
                                ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x * 0.25f, 400), ImGuiCond_Always);
        ImGui::SetNextWindowBgAlpha(0.80f);

        ImGuiWindowFlags flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings;

        auto &current_json = localization["Simulation"]["narration_menu"];

        if (ImGui::Begin("ScenarioNarration", &m_exitFlag, flags))
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
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.3f, 0.3f, 0.3f, 0.8f)); // dark gray
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.3f, 0.3f, 0.3f, 0.8f));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.3f, 0.3f, 0.3f, 0.8f));
            }
            else
            {
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.2f, 0.4f, 0.6f, 1.0f)); // blue
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.3f, 0.5f, 0.7f, 1.0f));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.2f, 0.4f, 0.6f, 1.0f)); // blue
            }

            ImGui::BeginDisabled(!allow_previous);

            if (ImGui::Button(current_json["button_previous"].get<std::string>().c_str(), ImVec2(buttonWidth, 45)))
                if (m_currentStepIdx > 0)
                    lastStep();

            ImGui::EndDisabled();

            ImGui::PopStyleColor(3);

            ImGui::SameLine();

            std::string button_next = "";

            if (m_currentStepIdx >= getCurrentScenarioStepCount() - 1)
            {
                std::string exit_str = current_json["button_close"];
                button_next = exit_str.c_str();
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.7f, 0.2f, 0.2f, 0.8f)); // red
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.8f, 0.3f, 0.3f, 0.9f));
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.7f, 0.2f, 0.2f, 0.8f));
            }
            else
            {
                button_next = current_json["button_next"].get<std::string>();
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.3f, 0.6f, 0.3f, 1.0f));        // green
                ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.3f, 0.6f, 0.3f, 1.0f));  // green
                ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.4f, 0.7f, 0.4f, 1.0f)); // green
            }

            if (ImGui::Button(button_next.c_str(), ImVec2(buttonWidth, 45)))
            {
                if (m_currentStepIdx < getCurrentScenarioStepCount() - 1)
                    nextStep();
                else
                {
                    clear();
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
