#pragma once
#include <SFML/Graphics.hpp>
#include "simulation/constants.hpp"
#include "core/window.hpp"
#include <algorithm>
#include <cmath>

#define GLM_FORCE_RADIANS
#define GLM_FORCE_INTRINSICS
#define GLM_FORCE_ALIGNED_GENTYPES

#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace core
{
    struct camera_t
    {
        glm::vec3 target{0.0f, 0.0f, 0.0f};  // Point the camera looks at
        float distance = 30.0f;             // Distance from target
        float azimuth = 45.0f;              // Horizontal rotation (degrees)
        float elevation = 25.0f;            // Vertical rotation (degrees), -89 to +89
        float fov = 50.0f;                  // Field of view in degrees
        float nearPlane = 1.0f;
        float farPlane = 1000.0f;

        glm::vec3 eye() const
        {
            const float radAz = glm::radians(azimuth);
            const float radEl = glm::radians(elevation);

            const float cosEl = std::cos(radEl);
            const float sinEl = std::sin(radEl);

            const float x = distance * cosEl * std::sin(radAz);
            const float y = distance * cosEl * std::cos(radAz);
            const float z = distance * sinEl;

            return target + glm::vec3(x, y, z);
        }

        void rotateAzimuth(float deltaDegrees)
        {
            azimuth += deltaDegrees;

            if (azimuth >= 360.0f) azimuth -= 360.0f;
            if (azimuth < 0.0f)    azimuth += 360.0f;
        }

        void rotateElevation(float deltaDegrees)
        {
            elevation += deltaDegrees;
            elevation = std::clamp(elevation, -89.0f, 89.0f);
        }

        void rotateAroundTarget(float deltaAzimuth, float deltaElevation, float sensitivity = 1.0f)
        {
            rotateAzimuth(deltaAzimuth * sensitivity);
            rotateElevation(deltaElevation * sensitivity);
        }

        void setSideView(char axis)
        {
            switch (axis)
            {
                case 'X': azimuth = 0.0f;   elevation = 0.0f;   break;  // +X view
                case 'Y': azimuth = 90.0f;  elevation = 0.0f;   break;  // +Y view
                case 'Z': azimuth = 45.0f;  elevation = 89.9f;  break;  // Top-down (+Z)
                default: break;
            }
        }

        glm::vec2 project(const glm::vec3& worldPos, float viewportWidth, float viewportHeight) const
        {
            glm::mat4 view = getViewMatrix();
            glm::mat4 proj = getProjectionMatrix(viewportWidth, viewportHeight);

            glm::vec4 clip = proj * view * glm::vec4(worldPos, 1.0f);

            if (clip.w <= 0.0f)
                return glm::vec2(-9999.f, -9999.f);

            glm::vec3 ndc = glm::vec3(clip) / clip.w;

            float screenX = (ndc.x * 0.5f + 0.5f) * viewportWidth;
            float screenY = (1.0f - ndc.y * 0.5f - 0.5f) * viewportHeight;

            return glm::vec2(screenX, screenY);
        }

        glm::mat4 getProjectionMatrix(float viewportWidth, float viewportHeight) const 
        {
            float aspect = static_cast<float>(viewportWidth) / static_cast<float>(viewportHeight);

            return glm::perspective(
                glm::radians(fov),
                aspect,
                nearPlane,
                farPlane
            );
        }

        glm::mat4 getViewMatrix() const
        {
            const glm::vec3 e = eye();
            const glm::vec3 f = glm::normalize(target - e);
            const glm::vec3 u = glm::vec3(0.0f, 0.0f, 1.0f);
            const glm::vec3 r = glm::normalize(glm::cross(f, u));
            const glm::vec3 actualUp = glm::cross(r, f);

            glm::mat4 view(1.0f);
            view[0] = glm::vec4(r, 0.0f);
            view[1] = glm::vec4(actualUp, 0.0f);
            view[2] = glm::vec4(-f, 0.0f);
            view[3] = glm::vec4(e, 1.0f);

            return glm::inverse(view);
        }
    };
}