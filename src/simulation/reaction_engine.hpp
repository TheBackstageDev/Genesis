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

    struct reaction_rule
    {

    };

    class reaction_engine
    {
    public:
        reaction_engine();
        ~reaction_engine() = default;

        void register_rule(reaction_rule r) { m_rules.emplace_back(std::move(r)); }

        void update(fun::universe& u);
    private:
        fun::universe* m_universe = nullptr;
        std::vector<reaction_rule> m_rules;

        void setupRules();

        bool isRadical(const fun::subset& sub, fun::universe& u);

        uint32_t getClosestAtom(const fun::subset& target, const fun::subset& reference);

        // Bond-related
        void createBond(fun::subset& incoming, fun::subset& target);
        void breakBond(fun::subset& a, fun::subset& b);
    };
} // namespace sim
