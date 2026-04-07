#include "reaction_engine.hpp"
#include <iostream>
#include <numeric>

namespace sim
{
    reaction_engine::reaction_engine()
    {
    }

    void reaction_engine::update(fun::universe& u, core::verlet_list& verlet_list)
    {
        
    }
    
    // Bond Related

    void reaction_engine::createBond(const uint32_t a, const uint32_t b)
    {
        if (!m_universe) return;

        m_universe->createBond(a, b);
        m_universe->updateAngles(a);
    }

    void reaction_engine::breakBond(const uint32_t a, const uint32_t b)
    {
        if (!m_universe) return;
        
        m_universe->breakBond(a, b);
        m_universe->updateAngles(a);
        m_universe->updateAngles(b);
    }
} // namespace sim
