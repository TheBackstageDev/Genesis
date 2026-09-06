#include "UIEditor.hpp"

namespace ui
{
    UIEditor::UIEditor(sim::parameter_table& paramTable)
    {
        m_hierarchy = std::make_unique<UIHierarchy>();
        m_parameters = std::make_unique<UIParameters>(paramTable);
        m_compoundBrowser = std::make_unique<UICompoundBrowser>();
    }

    UIEditor::~UIEditor() = default;

    void UIEditor::draw(sim::sim_dynamics& dynamics, sim::fun::universe& universe)
    {
        if (!m_visible)
            return;


        m_parameters->draw(dynamics, universe);
    }
} // namespace ui