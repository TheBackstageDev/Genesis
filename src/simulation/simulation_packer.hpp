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

        void pack(fun::universe& u, const fun::molecule_structure& molecule, const glm::vec3 center, const glm::vec3 box, const uint32_t targetAmmount = 0);
        void pack(fun::universe& u, const std::vector<fun::molecule_structure>& molecules, const std::vector<float>& chances, const glm::vec3 center, const glm::vec3 box, const uint32_t targetAmmount = 0);
    private:
        void optimize(fun::universe& u, const fun::molecule_structure& molecule, const glm::vec3 molecule_position);

        std::vector<float> normalizeChances(const std::vector<float>& chances);

        uint32_t getSpawnAmmount(const glm::vec3 box, const std::vector<fun::molecule_structure>& molecules, const std::vector<float>& chances);
        uint32_t getSpawnAmmount(const glm::vec3 box, const fun::molecule_structure& molecule);

        float getVDWVolume(const fun::molecule_structure& molecule);

        float randFloat(float min, float max);
        glm::quat randomRotation();
    };
}