#pragma once
#include <SFML/Graphics.hpp>
#include "simulation/constants.hpp"
#include "core/window.hpp"
#include <cmath>

namespace core
{
    struct camera_t
    {
        sf::Vector3f target{20, 10, 0};      // look-at point
        float distance = 1.f;            // eye distance from target
        float azimuth = 0.f;              // horizontal angle (deg)
        float elevation = 0.f;            // vertical angle (deg)
        float fov = 50.f;                  // field of view (deg)
        float nearPlane = 1.f;
        float farPlane = 1000.f;

        sf::Vector3f eye() const
        {
            float x = distance * cos(elevation * RADIAN) * sin(azimuth * RADIAN);
            float y = distance * cos(elevation * RADIAN) * cos(azimuth * RADIAN);
            float z = distance * sin(elevation * RADIAN);
            return target + sf::Vector3f{x, y, z};
        }

        void setSideView(char axis)
        {
            switch (axis)
            {
                case 'X': azimuth = 0.f;   elevation = 0.f;   break;
                case 'Y': azimuth = 90.f;  elevation = 0.f;   break;
                case 'Z': azimuth = 45.f;  elevation = 90.f;  break;
            }
        }

        sf::Vector2f project(const sf::Vector3f& p, float w, float h) const
        {
            sf::Vector3f eyePos = eye();
            sf::Vector3f view = p - eyePos;
            sf::Vector3f forward = distance == 0.f ? sf::Vector3f(1.f, 0.f, 0.f) : (target - eyePos).normalized();
            sf::Vector3f up{0, 0, 1};
            sf::Vector3f right = forward.cross(up).normalized();
            up = right.cross(forward).normalized();

            float depth = view.dot(forward);
            if (depth <= nearPlane || depth >= farPlane)
                return {-9999.f, -9999.f};

            float f = 1.0f / std::tan(fov * 0.5f * M_PI / 180.0f);
            float aspect = w / h;

            float px = view.dot(right) * f / depth;
            float py = view.dot(up)    * f / depth;

            float x = (px + 1.0f) * 0.5f * w;
            float y = (1.0f - py) * 0.5f * h;

            return {x, y};
        }

        sf::Transform getViewMatrix() const
        {
            sf::Vector3f e  = eye();
            sf::Vector3f f  = (target - e).normalized();
            sf::Vector3f up = sf::Vector3f(0, 0, 1);
            sf::Vector3f r  = f.cross(up).normalized();
            up = r.cross(f);

            sf::Transform rotation(
                r.x, r.y, r.z,
                up.x, up.y, up.z,
                -f.x, -f.y, -f.z
            );

            sf::Transform translation;
            translation.translate({-r.dot(e), -up.dot(e)});

            return translation * rotation;
        }
    };
}