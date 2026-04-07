#pragma once

#include "fundamental_structures.hpp"
#include "universe.hpp"

#include "core/verletlist.hpp"
#include "core/spatialgrid.hpp"

#include <vector>
#include <functional>

namespace sim
{
    class reaction_engine
    {
    public:
        reaction_engine();
        ~reaction_engine() = default;

        void update(fun::universe& u, core::verlet_list& verlet_list);
    private:
        fun::universe* m_universe = nullptr;

        void redistributeCharge();

        // Bond-related
        void createBond(const uint32_t a, const uint32_t b);
        void breakBond(const uint32_t a, const uint32_t b);
    };
} // namespace sim
