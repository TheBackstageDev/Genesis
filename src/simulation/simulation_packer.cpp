#include "simulation_packer.hpp"
#include <cassert>

#include <glm/gtc/quaternion.hpp>
#include <iostream>

namespace sim
{
    simulation_packer::simulation_packer()
    {
    }

    void simulation_packer::pack(fun::universe &u, const fun::molecule_structure &molecule, const glm::vec3 center, const glm::vec3 box)
    {
        pack(u, {molecule}, {1.0f}, center, box);
    }

    void simulation_packer::pack(fun::universe &u, const std::vector<fun::molecule_structure> &molecules, const std::vector<float> &chances, const glm::vec3 center, const glm::vec3 box)
    {
        assert(molecules.size() == chances.size() && "Molecules size and Chances must match!");
        assert(!molecules.empty() && !chances.empty() && "Empty molecules or Chances!");

        uint32_t target_ammount = getSpawnAmmount(box, molecules, chances);

        std::vector<float> normalized = normalizeChances(chances);
        std::vector<float> cumulative(normalized.size());

        cumulative[0] = normalized[0];
        for (size_t i = 1; i < normalized.size(); ++i) 
        {
            cumulative[i] = cumulative[i - 1] + normalized[i];
        }

        uint32_t created = 0;
        uint32_t attempts = 0;

        constexpr uint32_t maxAttempts = 1000;

        float maxExtent = 0.0f;
        for (const auto& mol : molecules) 
        {
            for (const auto& pos : mol.positions) 
            {
                float r = pos.length();
                if (r > maxExtent) maxExtent = r;
            }
        }

        glm::vec3 halfBox = box * 0.5f;

        glm::vec3 effectiveHalf = halfBox - maxExtent * 1.2f;

        if (effectiveHalf.x <= 0.0f || effectiveHalf.y <= 0.0f || effectiveHalf.z <= 0.0f)
        {
            std::cerr << "[Simulation Packer]: Box too small for molecules!\n";
            return;
        }

        auto positions = u.positions();
        positions.reserve(positions.size() + target_ammount * molecules.size());

        while (created < target_ammount && attempts < maxAttempts)
        {
            float r = randFloat(0.0f, 1.0f);
            size_t typeIndex = 0;
            for (; typeIndex < cumulative.size(); ++typeIndex) 
                if (r <= cumulative[typeIndex]) break;
            if (typeIndex >= molecules.size()) typeIndex = molecules.size() - 1;

            auto chosenMol = molecules[typeIndex];

            glm::vec3 offset(
                randFloat(-effectiveHalf.x, effectiveHalf.x),
                randFloat(-effectiveHalf.y, effectiveHalf.y),
                randFloat(-effectiveHalf.z, effectiveHalf.z)
            );
            glm::vec3 spawnCenter = center + offset;

            bool tooClose = false;
            for (const auto& prev : positions) 
            {
                if (glm::length(spawnCenter - prev) < 2.8f) 
                {
                    tooClose = true;
                    break;
                }
            }
            if (tooClose) 
            {
                ++attempts;
                continue;
            }

            glm::quat rotation = randomRotation();
            glm::mat4 rotationMatrix = glm::mat4_cast(rotation);

            for (auto& pos : chosenMol.positions)
            {
                glm::vec3 local_pos = glm::vec3(pos.x, pos.y, pos.z);
                glm::vec3 rotatedPos = glm::vec3(rotationMatrix * glm::vec4(local_pos, 1.0f));
                glm::vec3 finalPos = spawnCenter + rotatedPos;

                pos = sf::Vector3f(finalPos.x, finalPos.y, finalPos.z);

                positions.emplace_back(finalPos);
            }

            u.createMolecule(chosenMol, sf::Vector3f(spawnCenter.x, spawnCenter.y, spawnCenter.z), sf::Vector3f(randFloat(-0.1, 0.1), randFloat(-0.1, 0.1), randFloat(-0.1, 0.1)));
            
            ++created;
            attempts = 0;
        }
    }

    uint32_t simulation_packer::getSpawnAmmount(const glm::vec3 box, const std::vector<fun::molecule_structure> &molecules, const std::vector<float> &chances)
    {
        float total_volume = box.x * box.y * box.z;

        constexpr float packing_density = 0.55f;

        float avg_volume = 0.0f;
        float sum_chances = 0.0f;

        for (size_t i = 0; i < molecules.size(); ++i)
        {
            float v = getVDWVolume(molecules[i]);
            avg_volume += v * chances[i];
            sum_chances += chances[i];
        }

        if (sum_chances <= 0.0f)
            return 0;
        avg_volume /= sum_chances;

        uint32_t estimated = static_cast<uint32_t>(total_volume / avg_volume * packing_density);

        return estimated;
    }

    uint32_t simulation_packer::getSpawnAmmount(const glm::vec3 box, const fun::molecule_structure &molecule)
    {
        return getSpawnAmmount(box, {molecule}, {1.0f});
    }

    float simulation_packer::getVDWVolume(const fun::molecule_structure &molecule)
    {
        float total_radii = 0.0f;

        for (auto &atom : molecule.atoms)
            total_radii += constants::VDW_RADII[atom.ZIndex];

        return (4.f / 3.f) * M_PI * std::powf(total_radii / molecule.atoms.size(), 3);
    }

    std::vector<float> simulation_packer::normalizeChances(const std::vector<float> &chances)
    {
        float total_chance = std::accumulate(chances.begin(), chances.end(), 0.f, std::plus<float>());

        if (total_chance <= 0.0f)
        {
            std::cerr << "[SIMULATION PACKER]: All chances are zero or negative. Nothing will be spawned.\n";
            return {};
        }

        std::vector<float> normalizedChances(chances.size(), 0.0f);

        for (int32_t i = 0; i < normalizedChances.size(); ++i)
            normalizedChances[i] = chances[i] / total_chance;

        return normalizedChances;
    }

    float simulation_packer::randFloat(float min, float max)
    {
        static std::mt19937 gen(std::random_device{}());
        static std::uniform_real_distribution<float> dist;
        dist.param(std::uniform_real_distribution<float>::param_type(min, max));
        return dist(gen);
    }

    glm::quat simulation_packer::randomRotation()
    {
        float u1 = randFloat(0.0f, 1.0f);
        float u2 = randFloat(0.0f, glm::two_pi<float>());
        float u3 = randFloat(0.0f, glm::two_pi<float>());

        float sqrt1 = std::sqrt(1.0f - u1);
        float sqrt2 = std::sqrt(u1);

        return glm::quat(
            std::cos(u3) * sqrt2,
            std::sin(u2) * sqrt1,
            std::cos(u2) * sqrt1,
            std::sin(u3) * sqrt2);
    }
}