#pragma once

#include <glad/glad.h>

#include "core/camera.hpp"
#include "core/window.hpp"
#include "core/shader.hpp"
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
        void drawAngle(ImDrawList* draw_list, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c);

        void drawDebug();
        void handleCamera();
        core::camera_t &camera() { return cam; }
        
    private:
        core::window_t &window;

        GLuint box_vao, box_vbo, color_vao, color_vbo, bond_vao, bond_vbo, arrow_vao, arrow_vbo;
        std::unordered_map<std::string, core::glProgram> programs;

        struct AtomInstance 
        {
            alignas(16) glm::vec3 position;
            alignas(16) float radius;
            alignas(16) glm::vec4 color;
        };

        struct BondInstance
        {
            glm::vec4 posA{0.0f};
            glm::vec4 posB{0.0f};
            glm::vec4 colorA{1.0f};
            glm::vec4 colorB{1.0f};
            alignas(16) float radius = 0.0f;
            alignas(16) float radiusA = 0.5f;
            alignas(16) float radiusB = 0.5f;
        };

        struct ArrowInstance
        {
            glm::vec4 posA{0.0f};
            glm::vec4 posB{0.0f};
            glm::vec4 color{1.0f};
            alignas(16) float radius = 0.0f;
        };

        const float licorice_radius = 0.2f;

        void drawHydrogenBond(sf::RenderTarget &target, int32_t H, const fun::rendering_simulation_info &sim_info);
        void drawChargeField(sf::RenderTarget &target, const fun::rendering_simulation_info &sim_info);
        void drawBox(const glm::vec3 &box, sf::RenderTarget &target);

        // Init

        void loadProgram(const std::string key, const std::string& vert, const std::string& frag);

        void bindColor(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);
        void bindBond(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);
        void bindArrow(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);

        // Others

        void drawHighlight(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);
        void drawBondHighlight(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);
        void drawAtomHighlight(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);

        glm::vec4 getAtomColor(const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info, const uint32_t i);

        // Callbacks

        std::function<uint32_t(int32_t ix, int32_t iy, int32_t iz)> getCellID;
        std::function<glm::vec3(glm::vec3 dir)> minImageVec;

        // Camera
        ImVec2 lastMouse;
        glm::vec2 project(const glm::vec3 &p) const;
        core::camera_t cam;
    };
} // namespace sim
