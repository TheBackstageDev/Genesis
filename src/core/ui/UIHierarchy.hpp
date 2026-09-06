#pragma once

#include "simulation/core/universe.hpp"
#include <imgui.h>

namespace ui
{
    class UIHierarchy
    {
    public:
        UIHierarchy();

        void draw(sim::fun::universe& universe);
    private:
    };
} // namespace ui
