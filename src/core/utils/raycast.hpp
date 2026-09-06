#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <imgui.h>

#include <optional>
#include "core/graphics/camera.hpp"

namespace core
{
    struct Ray
    {
        glm::vec3 origin{0.0f};
        glm::vec3 direction{0.0f, 0.0f, -1.0f};   // default forward

        Ray() = default;
        Ray(glm::vec3 o, glm::vec3 d) : origin(o), direction(glm::normalize(d)) {}
    };

    struct RayHit
    {
        float distance = 0.0f;
        glm::vec3 point{0.0f};
        glm::vec3 normal{0.0f};
        bool hit = false;
    };

    std::optional<RayHit> raySphereIntersect(const Ray& ray, const glm::vec3& center, float radius);
    Ray screenToRay(ImVec2 mousePos, const core::extent2D& viewportSize, const camera_t& cam);
} // namespace core
