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
        for (int32_t i = 0; i < subsets.size(); ++i)
        {
            auto& subset_i = subsets[i];
            if (subset_i.group != fun::functionalGroup::RADICAL)
            {
                if (isRadical(subset_i, u))
                {
                    subset_i.group = fun::functionalGroup::RADICAL;
                    std::cout << "subset " << subset_i.mainAtomIdx << " is a radical! \n";
                }
                else
                {
                    subset_i.group = fun::functionalGroup::ANY;
                }
            }

            for (int32_t j = i + 1; j < subsets.size(); ++j)
            {
                auto& subset_j = subsets[j];

                float dr = glm::length(u.minImageVec(u.getPosition(subset_i.mainAtomIdx) - u.getPosition(subset_j.mainAtomIdx)));
                if (dr < 2.5f)
                for (auto& rule : m_rules)
                {
                    if (rule.match(subset_i.group, subset_i.hydrogenCount, subset_j.group, subset_j.hydrogenCount))
                        rule.action(subset_i, subset_j, *this);
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

        int32_t unpairedElectrons = 0;

        const fun::atom& mainAtom = atoms[sub.mainAtomIdx]; 
        int32_t expectedMainCount = constants::getUsualBonds(mainAtom.ZIndex);
        unpairedElectrons += std::max(0, expectedMainCount - mainAtom.bondCount);

        return unpairedElectrons > 0;
    }
} // namespace sim
