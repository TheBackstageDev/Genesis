#pragma once

#include "core/shader.hpp"
#include "core/buffer.hpp"

#include "core/spatialgrid.hpp"
#include "core/verletlist.hpp"
#include "universe.hpp"

#include <glad/glad.h>
#include <glm/glm.hpp>

namespace sim
{
    class sim_dynamics
    {
    public:
        sim_dynamics(fun::universe &u);
        ~sim_dynamics();

        void step(float target_temp = 300.0f, float target_pressure = 1.01325f);

        void pause() { m_paused = true; }
        void unpause() { m_paused = false; }

        size_t step_count() const { return m_step_count; }
        bool isPaused() { return m_paused; }

        void activateGPU() { m_GPU = true; }
        void deactivateGPU() { m_GPU = false; }

        void setTimescale(float timescale = 1.0f) { m_timescale = timescale; }

        size_t timestep() { return m_step_count; }
        float timescale() { return m_timescale; }
        float accumulated_time() { return m_accumulated_time; }
        float temperature() { return m_temperature; }
        float pressure() { return m_pressure; }

        core::verlet_list& getVerlet() { return universe_verlet; }
        fun::universe &getUniverse() { return m_universe; }
        
        void destroySSBOs();
    private:
        fun::universe &m_universe;

        core::glProgram unbonded_program;
        core::glProgram bonded_program;
        core::glProgram angles_program;
        core::glProgram integrate_program;

        bool m_builtSSBOs = false;
        core::glBuffer ssbo_pos, ssbo_vel, ssbo_force, ssbo_lj, ssbo_q, ssbo_bonds, ssbo_angles, ssbo_bondedPairs;

        struct bufferBond
        {
            uint32_t A, B;
            float K, r0;
        };

        struct bufferAngle
        {
            uint32_t A, B, C;
            float k, theta;
        };

        void createComputeShaders();
        void updateSSBOs();

        std::vector<glm::vec3> processCellUnbonded(int32_t ix, int32_t iy, int32_t iz, int32_t atom_start = 0, int32_t atom_end = 0);
        void computeUnbondedCPU();
        void computeBondedCPU();
        void computeExternalForces(uint32_t i);

        glm::vec3 computeLJforce(uint32_t i, uint32_t j, glm::vec3 &dr_vec);
        glm::vec3 computeCoulombForce(uint32_t i, uint32_t j, glm::vec3 &dr_vec);
        float computeDihedral(const glm::vec3 &pa, const glm::vec3 &pb, const glm::vec3 &pc, const glm::vec3 &pd);
        float computePressure();
        void COMDrift();

        void computeUnbondedGPU();
        void computeBondedGPU();
        void computeAnglesGPU();
        void syncBuffers();

        float m_dt = FEMTOSECOND;
        float m_timescale = 1.0f;
        float m_accumulated_time = 0.0f;
        uint64_t m_step_count = 0;

        std::atomic<float> total_virial{0.0f};
        float m_temperature = 0.0f;
        float m_pressure = 0.0f;

        size_t m_lastBondCount = 0;
        size_t m_lastAngleCount = 0;

        float gauss_random()
        {
            static std::random_device rd;
            static std::mt19937 gen(rd());
            static std::normal_distribution<float> dist(0.0f, 1.0f);

            return dist(gen);
        }

        const std::array<glm::vec3, 6> wall_directions = 
                {glm::vec3(1.f, 0.f, 0.f), glm::vec3(-1.f, 0.f, 0.f), 
                 glm::vec3(0.f, 1.f, 0.f), glm::vec3(0.f, -1.f, 0.f),
                 glm::vec3(0.f, 0.f, 1.f), glm::vec3(0.f, 0.f, -1.f)};

        bool m_GPU = true;
        bool m_paused = false;
        core::SpatialGrid universe_grid{};
        core::verlet_list universe_verlet{};

        void setPressure(float bar = 100);
        void setTemperature(float kelvin = 0.f);

        void integrate();
        void computeForces();
        void zeroForces();
    };
} // namespace sim
