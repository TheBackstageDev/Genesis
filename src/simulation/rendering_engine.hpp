#pragma once

#include <glad/glad.h>

#include "core/camera.hpp"
#include "core/window.hpp"
#include "simulation/fundamental_structures.hpp"

#define GLM_FORCE_RADIANS
#define GLM_FORCE_INTRINSICS
#define GLM_FORCE_ALIGNED_GENTYPES

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>

#include <functional>
#include <SFML/Graphics.hpp>

namespace sim
{
    class rendering_engine
    {
    public:
        rendering_engine(core::window_t &window);
        ~rendering_engine() = default;

        void draw(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);

        void drawDebug();
        void handleCamera();
        core::camera_t &camera() { return cam; }
        
    private:
        core::window_t &window;

        GLuint vao, vbo;
        GLuint atom_program;

        struct AtomInstance 
        {
            alignas(16) glm::vec3 position;
            alignas(16) float radius;
            alignas(16) glm::vec4 color;
        };

        void drawHydrogenBond(sf::RenderTarget &target, int32_t H, const fun::rendering_simulation_info &sim_info);
        void drawBonds(sf::RenderTarget &target, const fun::rendering_simulation_info &sim_info);
        void drawRings(sf::RenderTarget &target, const fun::rendering_simulation_info &sim_info);
        void drawChargeField(sf::RenderTarget &target, const fun::rendering_simulation_info &sim_info);
        void drawBox(const glm::vec3 &box);
        
        // Callbacks

        std::function<uint32_t(int32_t ix, int32_t iy, int32_t iz)> getCellID;
        std::function<glm::vec3(glm::vec3 dir)> minImageVec;

        // Camera
        ImVec2 lastMouse;
        glm::vec2 project(const glm::vec3 &p) const;
        core::camera_t cam;
    };
} // namespace sim
