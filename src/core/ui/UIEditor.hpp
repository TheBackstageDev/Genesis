#pragma once
#include <memory>
#include <imgui.h>

#include "simulation/core/universe.hpp"
#include "simulation/physics/dynamics.hpp"

#include "UIHierarchy.hpp"
#include "UIParameters.hpp"
#include "UICompoundBrowser.hpp"

namespace ui
{
    class UIEditor
    {
    public:
        UIEditor(sim::parameter_table& paramTable);
        ~UIEditor();

        void draw(sim::sim_dynamics& dynamics, sim::fun::universe& universe);
        bool isVisible() const { return m_visible; }
        void setVisible(bool visible) { m_visible = visible; }

    private:

        bool m_visible = true;

        std::unique_ptr<UIHierarchy> m_hierarchy;
        std::unique_ptr<UIParameters> m_parameters;
        std::unique_ptr<UICompoundBrowser> m_compoundBrowser;

        void drawLeftPanel();      // Hierarchy
        void drawRightPanel();     // Parameters + Inspector
        void drawBottomToolbar();
    };
} // namespace ui