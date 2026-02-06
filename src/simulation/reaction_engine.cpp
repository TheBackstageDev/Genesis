#include "reaction_engine.hpp"
#include <iostream>

namespace sim
{
    reaction_engine::reaction_engine()
    {
        setupRules();
    }

    void reaction_engine::setupRules()
    {
    }

    void reaction_engine::update(fun::universe& u)
    {

    }
    
    // Bond Related

    void reaction_engine::createBond(fun::subset& incoming, fun::subset& target)
    {
        if (!m_universe) return;

        if (incoming.mainAtomIdx == UINT32_MAX || incoming.mainAtomIdx >= m_universe->getAtoms().size())
            return;
        
        uint32_t incomingAtom = incoming.mainAtomIdx;
        uint32_t targetAtom = target.mainAtomIdx;

        if (m_universe->areBonded(incomingAtom, targetAtom))
            return;

        m_universe->createBond(incomingAtom, targetAtom);
    }

    void reaction_engine::breakBond(fun::subset& a, fun::subset& b)
    {
        if (!m_universe) return;
        
    }

    uint32_t reaction_engine::getClosestAtom(const fun::subset& target, const fun::subset& reference)
    {
        if (!m_universe) return UINT32_MAX;

        glm::vec3 referencePos = m_universe->getPosition(reference.mainAtomIdx);

        uint32_t closestAtom = target.mainAtomIdx;
        float smallest_distance = glm::length(m_universe->minImageVec(m_universe->getPosition(target.mainAtomIdx) - referencePos));

        for (int32_t i = 0; i < target.hydrogenCount; ++i)
        {
            uint32_t currentIndex = target.hydrogenBegin + i;
            float hydrogen_distance = glm::length(m_universe->minImageVec(m_universe->getPosition(currentIndex) - referencePos));

            if (hydrogen_distance < smallest_distance)
            {
                closestAtom = currentIndex;
                smallest_distance = hydrogen_distance;
            }
        }

        return closestAtom;
    }

    bool reaction_engine::isRadical(const fun::subset& sub, fun::universe& u)
    {
        const auto& atoms = u.getAtoms();
        int32_t unpaired = 0;

        if (sub.hydrogenCount == 0 && atoms[sub.mainAtomIdx].ZIndex == 1) // Hydrogen H* Radical
            return true;

        if (sub.mainAtomIdx < atoms.size())
        {
            const fun::atom& main = atoms[sub.mainAtomIdx];
            int expected = constants::getUsualBonds(main.ZIndex);
            unpaired += std::max(0, expected - main.bondCount);
        }

        if (sub.connectedCount != UINT32_MAX && sub.connectedBegin != UINT32_MAX)
        {
            for (uint32_t k = 0; k < sub.connectedCount; ++k)
            {
                uint32_t idx = sub.connectedBegin + k;
                if (idx >= atoms.size()) continue;

                const fun::atom& a = atoms[idx];
                int expected = constants::getUsualBonds(a.ZIndex);
                unpaired += std::max(0, expected - a.bondCount);
            }
        }

        if (sub.hydrogenCount != UINT32_MAX && sub.hydrogenBegin != UINT32_MAX)
        {
            for (uint32_t k = 0; k < sub.hydrogenCount; ++k)
            {
                uint32_t idx = sub.hydrogenBegin + k;
                if (idx >= atoms.size()) continue;

                const fun::atom& h = atoms[idx];
                unpaired += std::max(0, 1 - h.bondCount); 
            }
        }

        return unpaired > 1;
    }
} // namespace sim
