#include "simulation_packer.hpp"
#include <cassert>

#include <glm/gtc/quaternion.hpp>
#include <iostream>

namespace sim
{
    simulation_packer::simulation_packer()
    {
    }

    double gcm3_to_da3(double density_g_cm3) 
    {
        constexpr double CM3_TO_A3 = 1e24;
        return density_g_cm3 * (AVOGADRO / CM3_TO_A3);
    }

    void simulation_packer::pack(fun::universe &u, const fun::molecule_structure &molecule, 
                const glm::vec3 center, const glm::vec3 box, const uint32_t targetAmmount)
    {
        pack(u, {molecule}, {1.0f}, center, box, targetAmmount);
    }

    void packDensity(fun::universe& u, const fun::molecule_structure& molecule, const glm::vec3 center, const glm::vec3 box, const float density)
    {
        packDensity(u, {molecule}, center, box, density);
    }

    void simulation_packer::pack(fun::universe &u, const std::vector<fun::molecule_structure> &molecules, 
                const std::vector<float> &chances, const glm::vec3 center, const glm::vec3 box, const uint32_t targetAmmount)
    {
        assert(molecules.size() == chances.size() && "Molecules size and Chances must match!");
        
        if (molecules.empty() || chances.empty()) return;

        uint32_t target_ammount = targetAmmount == 0 ? getSpawnAmmount(box, molecules, chances) : targetAmmount;

        std::vector<float> normalized = normalizeChances(chances);
        std::vector<float> cumulative(normalized.size());

        cumulative[0] = normalized[0];
        for (size_t i = 1; i < normalized.size(); ++i) 
        {
            cumulative[i] = cumulative[i - 1] + normalized[i];
        }

        std::vector<float> mol_radii = getMoleculesRadii(molecules);
        std::vector<uint32_t> mol_ammounts(molecules.size());
        
        for (int32_t i = 0; i < target_ammount; ++i)
        {
            float r = randFloat(0.0f, 1.0f);
            size_t typeIndex = 0;
            for (; typeIndex < cumulative.size(); ++typeIndex) 
                if (r <= cumulative[typeIndex]) break;
            if (typeIndex >= molecules.size()) typeIndex = molecules.size() - 1;

            ++mol_ammounts[typeIndex];
        }

        std::vector<uint32_t> order(molecules.size()); // the order that the molecules will be placed, from biggest to smallest.
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](const uint32_t a, const uint32_t b){ return mol_radii[a] > mol_radii[b]; });

        glm::vec3 halfBox = box * 0.5f;

        positionMolecules(u, molecules, order, mol_radii, mol_ammounts, halfBox, center);
    }

    void simulation_packer::packDensity(fun::universe &u, const std::vector<fun::molecule_structure> &molecules, 
                const glm::vec3 center, const glm::vec3 box, const float density)
    {
        /* assert(molecules.size() == chances.size() && "Molecules size and Chances must match!");
        
        if (molecules.empty() || chances.empty()) return;

        uint32_t target_ammount = targetAmmount == 0 ? getSpawnAmmount(box, molecules, chances) : targetAmmount;

        std::vector<float> normalized = normalizeChances(chances);
        std::vector<float> cumulative(normalized.size());

        cumulative[0] = normalized[0];
        for (size_t i = 1; i < normalized.size(); ++i) 
        {
            cumulative[i] = cumulative[i - 1] + normalized[i];
        }

        std::vector<float> mol_radii = getMoleculesRadii(molecules);
        std::vector<uint32_t> mol_ammounts(molecules.size());
        
        for (int32_t i = 0; i < target_ammount; ++i)
        {
            float r = randFloat(0.0f, 1.0f);
            size_t typeIndex = 0;
            for (; typeIndex < cumulative.size(); ++typeIndex) 
                if (r <= cumulative[typeIndex]) break;
            if (typeIndex >= molecules.size()) typeIndex = molecules.size() - 1;

            ++mol_ammounts[typeIndex];
        }

        std::vector<uint32_t> order(molecules.size()); // the order that the molecules will be placed, from biggest to smallest.
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](const uint32_t a, const uint32_t b){ return mol_radii[a] > mol_radii[b]; });

        glm::vec3 halfBox = box * 0.5f;

        positionMolecules(u, molecules, order, mol_radii, mol_ammounts, halfBox, center); */
    }
    
    void simulation_packer::positionMolecules(fun::universe& u, const std::vector<fun::molecule_structure>& molecules, const std::vector<uint32_t> mol_order, 
                                                const std::vector<float> mol_radii, const std::vector<uint32_t> mol_ammounts, const glm::vec3 halfBox, const glm::vec3 center)
    {
        constexpr uint32_t maxAttempts = 200;

        for (uint32_t ord : mol_order)
        {
            uint32_t to_place = mol_ammounts[ord];
            if (to_place == 0) continue;

            auto& mol = molecules[ord];
            float this_radius = mol_radii[ord];

            uint32_t created_this = 0;
            uint32_t attempts = 0;

            glm::vec3 effectiveBox = halfBox - this_radius * 0.5f;

            if (effectiveBox.x <= 0.f || effectiveBox.y <= 0.f || effectiveBox.z <= 0.f) continue;

            while (created_this < to_place && attempts < maxAttempts)
            {
                ++attempts;

                glm::vec3 offset(
                    randFloat(-effectiveBox.x, effectiveBox.x),
                    randFloat(-effectiveBox.y, effectiveBox.y),
                    randFloat(-effectiveBox.z, effectiveBox.z)
                );
                glm::vec3 spawn_center = center + offset;

                bool quick_reject = false;
                for (const auto& existing : u.positions()) 
                {
                    if (glm::length(spawn_center - existing) < 2.2f) 
                    {
                        quick_reject = true;
                        break;
                    }
                }
                if (quick_reject) continue;

                glm::quat rot = randomRotation();
                glm::mat4 rot_mat = glm::mat4_cast(rot);

                bool collision = false;
                std::vector<glm::vec3> temp_positions;

                for (const auto& local_pos : mol.positions)
                {
                    glm::vec3 rotated = glm::vec3(rot_mat * glm::vec4(local_pos.x, local_pos.y, local_pos.z, 1.0f));
                    glm::vec3 world_pos = spawn_center + rotated;
                    
                    for (const auto& existing : u.positions())
                    {
                        glm::vec3 dr = u.minImageVec(world_pos - existing);
                        if (glm::length(dr) < 2.f) 
                        {
                            collision = true;
                            break;
                        }
                    }
                    if (collision) break;

                    temp_positions.push_back(world_pos);
                }

                if (collision) continue;

                u.createMolecule(mol, spawn_center,
                                glm::vec3(randFloat(-0.1f, 0.1f), randFloat(-0.1f, 0.1f), randFloat(-0.1f, 0.1f)));

                auto& data = u.getData();

                size_t start_idx = u.numAtoms() - mol.positions.size();
                for (size_t k = start_idx; k < u.numAtoms(); ++k) {
                    glm::vec3 local = data.positions[k] - spawn_center;
                    glm::vec3 rotated = glm::vec3(rot_mat * glm::vec4(local, 1.0f));
                    data.positions[k] = spawn_center + rotated;

                    data.positions[k] += glm::vec3(
                        randFloat(-0.05f, 0.05f),
                        randFloat(-0.05f, 0.05f),
                        randFloat(-0.05f, 0.05f)
                    );
                }

                ++created_this;
                attempts = 0;
            }

            if (created_this < to_place) 
            {
                std::cerr << "[Packer] Only placed " << created_this << "/" << to_place 
                        << " molecules of type " << ord << "\n";
            }
        }
    }

    uint32_t simulation_packer::getSpawnAmmount(const glm::vec3 box, const std::vector<fun::molecule_structure> &molecules, const std::vector<float> &chances)
    {
        float total_volume = box.x * box.y * box.z;

        constexpr float packing_density = 0.5f;

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

    std::vector<float> simulation_packer::getMoleculesRadii(const std::vector<fun::molecule_structure> &molecules)
    {
        std::vector<float> mol_radii(molecules.size());

        for (size_t m = 0; m < molecules.size(); ++m) 
        {
            for (const auto& p : molecules[m].positions) 
            {
                float r = glm::length(glm::vec3(p.x, p.y, p.z));
                if (r > mol_radii[m]) mol_radii[m] = r;
            }
            mol_radii[m] += 2.2f;
        }

        return mol_radii;
    }

    std::vector<float> getMoleculesMass(const std::vector<fun::molecule_structure>& molecules)
    {
        std::vector<float> mol_mass(molecules.size());

        for (size_t m = 0; m < molecules.size(); ++m) 
        {
            float mass = 0.f;
            for (const auto& atom : molecules[m].atoms)
                mass += MASS_NEUTRON * atom.NIndex + MASS_PROTON * atom.ZIndex;

            mol_mass[m] += mass;
        }

        return mol_mass;
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
            return chances;
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