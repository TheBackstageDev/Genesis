#pragma once

#include <glad/glad.h>
#include <imgui.h>

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

namespace sim
{   
    class rendering_engine
    {
    public:
        rendering_engine(core::window_t &window);
        ~rendering_engine() = default;

        void draw(const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);
        void drawAngle(ImDrawList* draw_list, const glm::vec3& a, const glm::vec3& b, const glm::vec3& c);

        void drawDebug();
        void handleCamera();
        core::camera_t &camera() { return cam; }

        int32_t pickAtom(ImVec2 mousePos, const std::vector<sim::fun::atom>& atoms);
        
    private:
        core::window_t &window;

        GLuint box_vao, box_vbo, color_vao, color_vbo, bond_vao, bond_vbo, arrow_vao, arrow_vbo;
        std::unordered_map<std::string, core::glProgram> programs;

        struct AtomInstance 
        {
            alignas(16) glm::vec3 position;
            alignas(16) float radius;
            alignas(16) glm::vec4 color;
            alignas(16) uint32_t id;
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

        struct SelectBuffer 
        {
            GLuint fbo = 0;
            GLuint texture = 0;
            int32_t width = 0;
            int32_t height = 0;
        } selectBuffer;

        std::vector<AtomInstance> m_atomInstances{};

        const float licorice_radius = 0.2f;

        void drawHydrogenBond(int32_t H, const fun::rendering_simulation_info &sim_info);
        void drawChargeField(const fun::rendering_simulation_info &sim_info);
        void drawBox(const glm::vec3 &box);

        // Init

        void loadProgram(const std::string key, const std::string& vert, const std::string& frag);

        void bindColor(const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);
        void bindBond(const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);
        void bindArrow(const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);

        void createSelectBuffer(core::extent2D windowExtent);

        // Others

        void drawHighlight(const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);
        void drawBondHighlight(const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);
        void drawAtomHighlight(const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info);

        glm::vec4 getAtomColor(const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info, const uint32_t i);

        // Callbacks

        std::function<uint32_t(int32_t ix, int32_t iy, int32_t iz)> getCellID;
        std::function<glm::vec3(glm::vec3 dir)> minImageVec;

        // Camera

        glm::mat4 lastProj{1.f};
        glm::mat4 lastView{1.f};

        ImVec2 lastMouse;
        glm::vec2 project(const glm::vec3 &p) const;
        core::camera_t cam;
    };
} // namespace sim
