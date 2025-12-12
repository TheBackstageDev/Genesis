#include "UIHandler.hpp"

namespace core
{
    void UIHandler::drawVideoControls()
    {
        ImGui::SetNextWindowPos(ImVec2(100, 100), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(400, 220), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowBgAlpha(0.85f);

        if (!ImGui::Begin("Video Controls", nullptr, ImGuiWindowFlags_AlwaysAutoResize))
        {
            ImGui::End();
            return;
        }
    }
} // namespace core
