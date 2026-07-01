#pragma once

#include "simulation/physics/fields/parameter_table.hpp"
#include "simulation/physics/dynamics.hpp"

#include <imgui.h>

namespace ui
{
    class UIParameters
    {
    public:
        UIParameters(sim::parameter_table& paramTable);
        
        void draw(sim::sim_dynamics& dynamics, sim::fun::universe& universe);

        void open() { m_active = true; }
        void close() { m_active = false; }

        float targetTemp() { return m_targetTemp; }
        float targetPressure() { return m_targetPres; }

        void setTargetTemp(float nTemp) { m_targetTemp = nTemp; }
        void setTargetPres(float nPres) { m_targetPres = nPres; }
    private:
        bool m_active = false;

        float m_targetTemp = 297.3f;
        float m_targetPres = 100.f;
        float m_targetTimeStep = 1.0f;

        sim::parameter_table& m_paramTable;
    };
    
} // namespace ui
