#pragma once

#include <random>
#include <glm/glm.hpp>
#include "fundamental_structures.hpp"
#include "universe.hpp"

namespace sim
{
    class simulation_packer
    {
    public:
        simulation_packer();

        void pack(fun::universe& u, const fun::molecule_structure& molecule, const glm::vec3 center, const float radius);
    private:
        void optimize(fun::universe& u, const fun::molecule_structure& molecule, const glm::vec3 molecule_position);
    };
}