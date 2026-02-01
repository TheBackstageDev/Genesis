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
        reaction_rule radical_radical{};
        radical_radical.subsetA = fun::functionalGroup::RADICAL;
        radical_radical.subsetB = fun::functionalGroup::RADICAL;
        radical_radical.priority = 1.0f;
        radical_radical.action = 
            [](fun::subset& a, fun::subset& b, reaction_engine& engine)
            {
                engine.createBond(a, b);
            };

        register_rule(radical_radical);
    }

    void reaction_engine::update(fun::universe& u)
    {
        m_universe = &u;

        auto& subsets = u.getSubsets();
        for (auto& sub : subsets)
        {
            bool was_radical = (sub.group == fun::functionalGroup::RADICAL);

            if (isRadical(sub, u))
            {
                if (!was_radical)
                {
                    sub.group = fun::functionalGroup::RADICAL;
                    std::cout << "[Radical detected] Subset centered on atom " 
                              << sub.mainAtomIdx << " (Z=" 
                              << u.getAtoms()[sub.mainAtomIdx].ZIndex << ")\n";
                }
            }
            else
            {
                sub.group = fun::functionalGroup::ANY;
            }
        }

        for (size_t i = 0; i < subsets.size(); ++i)
        {
            auto& subA = subsets[i];
            glm::vec3 posA = u.getPosition(subA.mainAtomIdx);

            for (size_t j = i + 1; j < subsets.size(); ++j)
            {
                auto& subB = subsets[j];
                glm::vec3 posB = u.getPosition(subB.mainAtomIdx);

                float dist = glm::length(u.minImageVec(posB - posA));

                if (dist > 2.f) continue;

                for (auto& rule : m_rules)
                {
                    if (rule.match(subA.group, subA.hydrogenCount,
                                   subB.group, subB.hydrogenCount))
                    {
                        std::cout << "[Reaction check] Rule matched between subsets " 
                                  << i << " and " << j << " (dist = " << dist << " Å)\n";

                        rule.action(subA, subB, *this);
                    }
                }
            }
        }
        m_universe = nullptr;
    }
    
    // Bond Related

    void reaction_engine::createBond(fun::subset& incoming, fun::subset& target)
    {
        if (!m_universe) return;

        if (incoming.mainAtomIdx == UINT32_MAX || incoming.mainAtomIdx >= m_universe->getAtoms().size())
            return;
        
        uint32_t incomingAtom = incoming.mainAtomIdx;
        uint32_t bestAtom = getClosestAtom(target, incoming);
        
        if (m_universe->areBonded(incomingAtom, bestAtom))
            return;

        m_universe->createBond(incomingAtom, bestAtom);
    }

    void reaction_engine::breakBond(fun::subset& a, fun::subset& b)
    {
        if (!m_universe) return;
        
    }

    void reaction_engine::transferHydrogen(fun::subset& donor, fun::subset& acceptor)
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
