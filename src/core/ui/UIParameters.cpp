#include "UIParameters.hpp"

namespace ui
{
    UIParameters::UIParameters(sim::parameter_table& paramTable)
        : m_paramTable(paramTable)
    {
    }
    
    void UIParameters::draw(sim::sim_dynamics& dynamics, sim::fun::universe& universe)
    {
        if (!m_active) return;

        if (!ImGui::Begin("##Parameters", &m_active))
        {
            ImGui::End();
            return;
        }

        if (ImGui::BeginTabBar("ParamsTabs"))
        {
            if (ImGui::BeginTabItem("Thermodynamics"))
            {
                ImGui::Text("Temperature");
                ImGui::SliderFloat("Target (K)", &m_targetTemp, 10.0f, 2000.0f);
                if (ImGui::Button("Apply Temperature"))
                    dynamics.setTargetTemperature(m_targetTemp);

                ImGui::Separator();
                ImGui::Text("Pressure");
                ImGui::SliderFloat("Target (kPa)", &m_targetPres, 0.0f, 10000.0f);
                if (ImGui::Button("Apply Pressure"))
                    dynamics.setTargetPressure(m_targetPres);

                ImGui::EndTabItem();
            }

            if (ImGui::BeginTabItem("Physics"))
            {
                ImGui::SliderFloat("Timestep (fs)", &m_targetTimeStep, 0.1f, 5.0f);
                dynamics.setTimescale(m_targetTimeStep);

                ImGui::EndTabItem();
            }

            if (ImGui::BeginTabItem("Interactions"))
            {
                if (ImGui::BeginTabBar("InteractionsSubTabs"))
                {
                    if (ImGui::BeginTabItem("Lennard-Jones"))
                    {
                        

                        ImGui::EndTabItem();
                    }

                    if (ImGui::BeginTabItem("Morse (Bonds)"))
                    {


                        ImGui::EndTabItem();
                    }

                    ImGui::EndTabBar();
                }

                ImGui::EndTabItem();
            }

            ImGui::EndTabBar();
        }

        ImGui::End();
    }
} // namespace ui
