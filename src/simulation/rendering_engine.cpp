#include "rendering_engine.hpp"
#include <algorithm>
#include <numeric>

#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <stdexcept>

#include <string>
#include <fstream>
#include <sstream>

namespace sim
{
    std::string loadShaderSource(const std::filesystem::path &path)
    {
        std::ifstream file(path, std::ios::ate | std::ios::binary);
        if (!file.is_open())
        {
            throw std::runtime_error("Cannot open shader: " + path.string());
        }

        auto size = file.tellg();
        std::string content(size, '\0');
        file.seekg(0);
        file.read(content.data(), size);

        if (file.fail())
        {
            throw std::runtime_error("Failed to read shader: " + path.string());
        }

        return content;
    }

    rendering_engine::rendering_engine(core::window_t &window)
        : window(window)
    {
        const std::filesystem::path shader_root = "src/shaders";
        std::string vertex = loadShaderSource(shader_root / "color.vert");
        std::string fragment = loadShaderSource(shader_root / "color.frag");

        GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
        const char *vertex_src = vertex.c_str();
        glShaderSource(vertex_shader, 1, &vertex_src, nullptr);
        glCompileShader(vertex_shader);

        GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
        const char *fragment_src = fragment.c_str();
        glShaderSource(fragment_shader, 1, &fragment_src, nullptr);
        glCompileShader(fragment_shader);

        atom_program = glCreateProgram();
        glAttachShader(atom_program, vertex_shader);
        glAttachShader(atom_program, fragment_shader);
        glLinkProgram(atom_program);

        GLint success = 0;
        glGetProgramiv(atom_program, GL_LINK_STATUS, &success);
        if (!success)
        {
            char info[1024] = {};
            glGetProgramInfoLog(atom_program, 1024, nullptr, info);
            std::cerr << info << "\n";
            throw std::runtime_error(std::string("[Shader link failed]: ") + info);
        }

        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);

        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);

        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomInstance), (void *)offsetof(AtomInstance, position));
        glVertexAttribDivisor(0, 1);

        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(AtomInstance), (void *)offsetof(AtomInstance, radius));
        glVertexAttribDivisor(1, 1);

        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(AtomInstance), (void *)offsetof(AtomInstance, color));
        glVertexAttribDivisor(2, 1);
    }

    void rendering_engine::drawBox(const glm::vec3 &box)
    {
        const std::array<glm::vec3, 8> corners =
            {{{0.f, 0.f, 0.f},
              {box.x, 0.f, 0.f},
              {box.x, box.y, 0.f},
              {0.f, box.y, 0.f},
              {0.f, 0.f, box.z},
              {box.x, 0.f, box.z},
              {box.x, box.y, box.z},
              {0.f, box.y, box.z}}};

        const std::array<std::pair<int32_t, int32_t>, 12> edges =
            {{
                {0, 1}, {1, 2}, {2, 3}, {3, 0}, // bottom face
                {4, 5},
                {5, 6},
                {6, 7},
                {7, 4}, // top face
                {0, 4},
                {1, 5},
                {2, 6},
                {3, 7} // vertical pillars
            }};

        sf::VertexArray lines(sf::PrimitiveType::Lines, edges.size() * 2);
        sf::Color edgeColor(200, 200, 200, 180);

        int32_t idx = 0;
        for (const auto &[i, j] : edges)
        {
            glm::vec2 p1 = project(corners[i]);
            glm::vec2 p2 = project(corners[j]);

            if (p1.x < -1000.f || p2.x < -1000.f)
                continue;

            const glm::vec3 &a = corners[i];
            const glm::vec3 &b = corners[j];

            glm::vec3 cam_eye = cam.eye();
            glm::vec3 camToA = a - glm::vec3(cam_eye.x, cam_eye.y, cam_eye.z);
            glm::vec3 camToB = b - glm::vec3(cam_eye.x, cam_eye.y, cam_eye.z);

            float dot = glm::dot(glm::normalize(camToA), glm::normalize(camToB));
            const float THRESHOLD = 0.3f;
            if (dot < THRESHOLD)
                continue;

            lines[idx++] = sf::Vertex(sf::Vector2f(p1.x, p1.y), edgeColor);
            lines[idx++] = sf::Vertex(sf::Vector2f(p2.x, p2.y), edgeColor);
        }

        window.getWindow().draw(lines);
    }

    void rendering_engine::draw(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info)
    {
        std::vector<int32_t> drawOrder(sim_info.positions.size());
        std::iota(drawOrder.begin(), drawOrder.end(), 0);

        std::sort(drawOrder.begin(), drawOrder.end(),
                  [&](int32_t a, int32_t b)
                  { return (sim_info.positions[a] - cam.eye()).length() > (sim_info.positions[b] - cam.eye()).length(); });

        std::vector<int32_t> no_draw{};

        if (!sim_info.renderWater)
            for (auto &m : sim_info.molecules)
            {
                if (m.exclude)
                {
                    no_draw.emplace_back(m.angleBegin);
                    no_draw.emplace_back(m.angleBegin + 1);
                    no_draw.emplace_back(m.angleBegin + 2);
                }
            }

        std::vector<AtomInstance> instances;
        instances.reserve(sim_info.atoms.size());

        if (info.lennardBall)
        for (int32_t sortedIdx : drawOrder)
        {
            const auto& atom = sim_info.atoms[sortedIdx];
            float radius = info.spaceFilling
                            ? constants::VDW_RADII[atom.ZIndex] * 1.5f
                            : atom.radius / 3.f;

            sf::Color col = constants::getElementColor(atom.ZIndex);
            glm::vec4 color_norm(col.r / 255.f, col.g / 255.f, col.b / 255.f, 1.0f);

            instances.emplace_back(sim_info.positions[sortedIdx], radius, color_norm);
        }

        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     instances.size() * sizeof(AtomInstance),
                     instances.data(),
                     GL_DYNAMIC_DRAW);

        glBindVertexArray(vao);
        glUseProgram(atom_program);

        GLint loc_projection = glGetUniformLocation(atom_program, "u_projection");
        GLint loc_view = glGetUniformLocation(atom_program, "u_view");
        GLint loc_cam_pos = glGetUniformLocation(atom_program, "u_campos");

        if (loc_projection != -1)
            glUniformMatrix4fv(loc_projection, 1, GL_FALSE, glm::value_ptr(cam.getProjectionMatrix(target)));

        if (loc_view != -1)
            glUniformMatrix4fv(loc_view, 1, GL_FALSE, glm::value_ptr(cam.getViewMatrix()));

        if (loc_cam_pos != -1)
            glUniform3fv(loc_cam_pos, 1, glm::value_ptr(cam.eye()));

        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, static_cast<GLsizei>(instances.size()));

        GLenum err = glGetError();
        if (err != GL_NO_ERROR)
        {
            std::cerr << "[RENDERING ENGINE] OpenGL error after draw: 0x" << std::hex << err << "\n";
        }

        glBindVertexArray(0);
        glUseProgram(0);

        target.pushGLStates();

        if (info.universeBox)
            drawBox(sim_info.box);
    
        glm::vec3 cam_eye = cam.eye();

        if (info.letter)
            for (int32_t i = 0; i < sim_info.positions.size(); ++i)
            {
                if (std::find(no_draw.begin(), no_draw.end(), drawOrder[i]) != no_draw.end())
                    continue;

                glm::vec2 p = project(sim_info.positions[drawOrder[i]]);
                if (p.x < -1000)
                    continue;

                glm::vec3 r_vec = sim_info.positions[drawOrder[i]] - cam_eye;

                sim_info.atoms[i].draw(sf::Vector2f(p.x, p.y), glm::length(r_vec), sim_info.q[drawOrder[i]], window, window.getWindow());
            }

        drawBonds(window.getWindow(), sim_info);

        target.popGLStates();
    }

    float smootherstep(float edge0, float edge1, float x)
    {
        x = std::clamp((x - edge0) / (edge1 - edge0), 0.f, 1.f);

        return x * x * (3.0f - 2.0f * x);
    }

    void rendering_engine::drawChargeField(sf::RenderTarget &target, const fun::rendering_simulation_info &sim_info)
    {
    }

    void rendering_engine::drawBonds(sf::RenderTarget &target, const fun::rendering_simulation_info &sim_info)
    {
        sf::Vector2f dimensions = window.getWindow().getView().getSize();

        for (int32_t b = 0; b < sim_info.bonds.size(); ++b)
        {
            const fun::bond &bond = sim_info.bonds[b];

            const glm::vec3 &pCentral = sim_info.positions[bond.centralAtom]; // pB
            const glm::vec3 &pBonded = sim_info.positions[bond.bondedAtom];   // pA

            glm::vec2 s1 = cam.project(pCentral, dimensions.x, dimensions.y); // center
            glm::vec2 s2 = cam.project(pBonded, dimensions.x, dimensions.y);  // end

            if (s1.x <= -9999 || s2.x <= -9999)
                continue;
            if ((pCentral - pBonded).length() > 3.f)
                continue;
            glm::vec3 cam_eye = cam.eye();
            if ((pCentral - cam_eye).length() > 50.f)
                continue;

            glm::vec2 dir = s2 - s1;
            float len = dir.length();
            if (len < 0.5f)
                continue;
            dir /= len;

            glm::vec2 perp{-dir.y, dir.x};

            uint8_t lines = static_cast<uint8_t>(bond.type);
            float distance = (pCentral - cam_eye).length();

            float shrink = 1.5f;  // pixels
            float spacing = 0.5f; // pixels between lines

            std::vector<sf::Vertex> verts;
            verts.reserve(lines * 2);

            for (int32_t i = 0; i < lines; ++i)
            {
                float offset = spacing * (i - (lines - 1) * 0.5f);

                glm::vec2 off = perp * offset;

                glm::vec2 start = s1 + dir * shrink + off;
                glm::vec2 end = s2 - dir * shrink + off;

                verts.emplace_back(sf::Vector2f(start.x, start.y), sf::Color::White);
                verts.emplace_back(sf::Vector2f(end.x, end.y), sf::Color::White);
            }

            target.draw(verts.data(), verts.size(), sf::PrimitiveType::Lines);
        }
    }

    void rendering_engine::drawRings(sf::RenderTarget &target, const fun::rendering_simulation_info &sim_info)
    {
        constexpr float visualRadius = 0.5f;
        constexpr float visualThickness = 0.1f;
        constexpr uint32_t segments = 60;
    }

    void rendering_engine::drawHydrogenBond(sf::RenderTarget &target, int32_t H, const fun::rendering_simulation_info &sim_info)
    {
    }

    void rendering_engine::drawDebug()
    {
    }

    glm::vec2 rendering_engine::project(const glm::vec3 &p) const
    {
        auto size = window.getWindow().getView().getSize();
        return cam.project(p, size.x, size.y);
    }

    void rendering_engine::handleCamera()
    {
        ImGuiIO &io = ImGui::GetIO();

        if (io.WantCaptureMouse || io.WantCaptureKeyboard)
            return;

        if (ImGui::IsAnyItemHovered() || ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow))
            return;

        if (ImGui::IsMouseDragging(ImGuiMouseButton_Left))
        {
            ImVec2 dragDelta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Left, 0.0f);

            cam.azimuth -= dragDelta.x * 0.25f;                                             // yaw
            cam.elevation = std::clamp(cam.elevation - dragDelta.y * 0.25f, -89.0f, 89.0f); // pitch

            ImGui::ResetMouseDragDelta(ImGuiMouseButton_Left);
        }

        glm::vec3 cam_eye = cam.eye();
        glm::vec3 forward = glm::normalize(cam.target - cam_eye);
        glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0, 0, 1)));
        glm::vec3 up = glm::normalize(glm::cross(right, forward));

        if (ImGui::IsMouseDragging(ImGuiMouseButton_Right))
        {
            ImVec2 dragDelta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Right, 0.0f);
            float panSpeed = cam.distance * 0.002f;

            cam.target += right * (-dragDelta.x * panSpeed);
            cam.target += up * (dragDelta.y * panSpeed);

            ImGui::ResetMouseDragDelta(ImGuiMouseButton_Right);
        }

        float wheel = ImGui::GetIO().MouseWheel;
        if (wheel != 0.0f)
        {
            const float zoomSpeed = 6.0f;
            float factor = (wheel > 0) ? (1.0f / 1.1f) : 1.1f;
            cam.distance = std::clamp(cam.distance * std::pow(factor, zoomSpeed * std::abs(wheel)),
                                      5.0f, 1000.0f);
        }

        ImVec2 ioMousePos = ImGui::GetMousePos();
        (void)ioMousePos; // silence unused warning

        float moveSpeed = 0.5f;

        if (ImGui::IsKeyDown(ImGuiKey_W))
            cam.target += forward * moveSpeed;
        if (ImGui::IsKeyDown(ImGuiKey_S))
            cam.target -= forward * moveSpeed;
        if (ImGui::IsKeyDown(ImGuiKey_A))
            cam.target -= right * moveSpeed;
        if (ImGui::IsKeyDown(ImGuiKey_D))
            cam.target += right * moveSpeed;
        if (ImGui::IsKeyDown(ImGuiKey_Q))
            cam.target -= up * moveSpeed;
        if (ImGui::IsKeyDown(ImGuiKey_E))
            cam.target += up * moveSpeed;
    }
} // namespace sim
