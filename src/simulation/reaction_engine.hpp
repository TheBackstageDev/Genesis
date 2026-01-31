#pragma once

#include "fundamental_structures.hpp"
#include "universe.hpp"
#include <vector>
#include <functional>

namespace sim
{
    enum class reaction_settings
    {
        REACTION_RULES,
        REAXFF,
        COUNT
    };

    class reaction_engine;

    struct reaction_rule
    {
        fun::functionalGroup subsetA;
        uint32_t numHydrogenA = 0;

        fun::functionalGroup subsetB;
        uint32_t numHydrogenB = 0;

        std::function<void(fun::subset&, fun::subset&, reaction_engine&)> action;
        float priority = 1.0f;

        bool match(fun::functionalGroup groupA, uint32_t numhydrogenA, fun::functionalGroup groupB, uint32_t numhydrogenB)
        {
            bool matchesA = groupA == subsetA && numhydrogenA == numHydrogenA, matchesB = groupB == subsetB && numhydrogenB == numHydrogenB;
            bool anyA = groupA == fun::functionalGroup::ANY, anyB = groupB == fun::functionalGroup::ANY;
            
            return matchesA && matchesB || (anyA && anyB)
                    || anyA && matchesB || anyB && matchesA || groupA == fun::functionalGroup::RADICAL && groupB == fun::functionalGroup::RADICAL;
        }
    };

    class reaction_engine
    {
    public:
        reaction_engine();
        ~reaction_engine() = default;

        void register_rule(reaction_rule r) { m_rules.emplace_back(std::move(r)); }

        void update(fun::universe& u);
    private:
        friend reaction_rule;

        fun::universe* m_universe = nullptr;
        std::vector<reaction_rule> m_rules;

        void setupRules();

        bool isRadical(const fun::subset& sub, fun::universe& u);

        uint32_t getClosestAtom(const fun::subset& target, const fun::subset& reference);

        // Bond-related
        void createBond(fun::subset& incoming, fun::subset& target);
        void breakBond(fun::subset& a, fun::subset& b);
        void transferHydrogen(fun::subset& donor, fun::subset& acceptor);
        void breakWeakestBond(const fun::subset& root);
        void adjustBondOrders(std::vector<fun::subset>& subsets, std::vector<fun::bond>& bonds);
    };
} // namespace sim
