#pragma once

#include "core/shader.hpp"
#include "core/buffer.hpp"

#include "core/spatialgrid.hpp"
#include "universe.hpp"

#include <glad/glad.h>
#include <glm/glm.hpp>

namespace sim
{
    class sim_dynamics
    {
    public:
        sim_dynamics(fun::universe &u);

        void step(float target_temp = 300.0f, float target_pressure = 1.01325f);

        void pause() { m_paused = true; }
        void unpause() { m_paused = false; }

        size_t step_count() const { return m_step_count; }
        bool isPaused() { return m_paused; }

        fun::universe& getUniverse() { return m_universe; }

        glm::vec3 minImageVec(glm::vec3 dr)
        {
            bool wall_col = m_universe.wallcollision();
            bool rf_col = m_universe.rooffloorcollision();

            if (wall_col && rf_col) return dr;

            glm::vec3 box_sizes = m_universe.boxSizes();

            if (!rf_col)
                dr.z -= box_sizes.z * std::round(dr.z / box_sizes.z);

            if (!wall_col)
            {
                dr.x -= box_sizes.x * std::round(dr.x / box_sizes.x);
                dr.y -= box_sizes.y * std::round(dr.y / box_sizes.y);
            }

            return dr;
        }
    private:
        fun::universe &m_universe;

        core::glProgram unbonded_program;
        core::glProgram bonded_program;
        core::glProgram integrate_program;

        core::glBuffer ssbo_pos, ssbo_vel, ssbo_force, ssbo_lj, ssbo_q;

        void createComputeShaders();
        void updateSSBOs();

        std::vector<glm::vec3> processCellUnbonded(int32_t ix, int32_t iy, int32_t iz, int32_t atom_start = 0, int32_t atom_end = 0);
        void computeUnbondedCPU();
        void computeBondedCPU();

        glm::vec3 computeLJforce(uint32_t i, uint32_t j, glm::vec3 &dr_vec);
        glm::vec3 computeCoulombForce(uint32_t i, uint32_t j, glm::vec3 &dr_vec);
        float computeDihedral(const glm::vec3& pa, const glm::vec3& pb, const glm::vec3& pc, const glm::vec3& pd);
        float computePressure();

        void computeUnbondedGPU();
        void computeBondedGPU();
        
        float m_dt = FEMTOSECOND;
        float m_Timescale = 1.0f;
        uint64_t m_step_count = 0;
        
        bool m_paused = false;
        core::SpatialGrid universe_grid{};
        
        void setPressure(float bar = 100);
        void setTemperature(float kelvin = 0.f);
        
        void computeForces();
        void zeroForces();
        void syncBuffers();
        void integrate();
    };
} // namespace sim
