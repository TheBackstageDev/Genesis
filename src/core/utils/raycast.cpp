#include "raycast.hpp"

namespace core
{
    std::optional<RayHit> raySphereIntersect(const Ray& ray, const glm::vec3& center, float radius)
    {
        glm::vec3 oc = ray.origin - center;
        glm::vec3 d = ray.direction;
        
        float a = glm::dot(d, d);
        float b = 2.f * glm::dot(oc, d);
        float c = glm::dot(oc, oc) - radius * radius;

        float discriminant = b*b - 4.f*a*c;

        if (discriminant < 0.f)
            return std::nullopt;

        float sqrtDiscriminant = sqrtf(discriminant);

        float t1 = (-b - sqrtDiscriminant) / (2.f * a);
        float t2 = (-b + sqrtDiscriminant) / (2.f * a);

        if (t1 < 0.f && t2 < 0.f) 
            return std::nullopt;
        
        float t = std::numeric_limits<float>::infinity();

        if (t1 > 0.f) t = t1;
        if (t2 > 0.f && t2 < t) t = t2;

        if (t == std::numeric_limits<float>::infinity())
            return std::nullopt;

        RayHit hit{};
        hit.point = ray.origin + t * ray.direction;
        hit.normal = glm::normalize(hit.point - center);
        hit.distance = t;
        return hit;
    }

    Ray screenToRay(ImVec2 mousePos, const core::extent2D& viewportSize, const camera_t& cam)
    {
        float ndcX = (2.0f * mousePos.x) / viewportSize.width - 1.0f;
        float ndcY = 1.0f - (2.0f * mousePos.y) / viewportSize.height;

        glm::vec4 viewPort(0.0f, 0.0f, viewportSize.width, viewportSize.height);

        glm::vec3 nearPoint = glm::unProject(glm::vec3(ndcX, ndcY, 0.0f),
                                            cam.getViewMatrix(),
                                            cam.getProjectionMatrix(viewportSize),
                                            viewPort);

        glm::vec3 farPoint = glm::unProject(glm::vec3(ndcX, ndcY, 1.0f),
                                            cam.getViewMatrix(),
                                            cam.getProjectionMatrix(viewportSize),
                                            viewPort);

        Ray ray;
        ray.origin = nearPoint;
        ray.direction = glm::normalize(farPoint - nearPoint);
        return ray;
    }
} // namespace core
