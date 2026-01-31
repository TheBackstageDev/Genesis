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
        initShaders();

        glGenVertexArrays(1, &color_vao);
        glGenBuffers(1, &color_vbo);
        glBindVertexArray(color_vao);
        glBindBuffer(GL_ARRAY_BUFFER, color_vbo);

        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(AtomInstance), (void *)offsetof(AtomInstance, position));
        glVertexAttribDivisor(0, 1);

        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(AtomInstance), (void *)offsetof(AtomInstance, radius));
        glVertexAttribDivisor(1, 1);

        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(AtomInstance), (void *)offsetof(AtomInstance, color));
        glVertexAttribDivisor(2, 1);

        glGenVertexArrays(1, &bond_vao);
        glGenBuffers(1, &bond_vbo);
        glBindVertexArray(bond_vao);
        glBindBuffer(GL_ARRAY_BUFFER, bond_vbo);

        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(BondInstance), (void *)offsetof(BondInstance, posA));
        glVertexAttribDivisor(0, 1);

        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(BondInstance), (void *)offsetof(BondInstance, posB));
        glVertexAttribDivisor(1, 1);

        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(BondInstance), (void *)offsetof(BondInstance, colorA));
        glVertexAttribDivisor(2, 1);

        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(BondInstance), (void *)offsetof(BondInstance, colorB));
        glVertexAttribDivisor(3, 1);

        glEnableVertexAttribArray(4);
        glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(BondInstance), (void *)offsetof(BondInstance, radius));
        glVertexAttribDivisor(4, 1);

        glEnableVertexAttribArray(5);
        glVertexAttribPointer(5, 1, GL_FLOAT, GL_FALSE, sizeof(BondInstance), (void *)offsetof(BondInstance, radiusA));
        glVertexAttribDivisor(5, 1);

        glEnableVertexAttribArray(6);
        glVertexAttribPointer(6, 1, GL_FLOAT, GL_FALSE, sizeof(BondInstance), (void *)offsetof(BondInstance, radiusB));
        glVertexAttribDivisor(6, 1);
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

    void rendering_engine::initShaders()
    {
        const std::filesystem::path shader_root = "src/shaders";

        initShaders(shader_root / "color.vert", shader_root / "color.frag", atom_program);
        initShaders(shader_root / "bond.vert", shader_root / "bond.frag", bond_program);
        initShaders(shader_root / "hyper_balls.vert", shader_root / "hyper_balls.frag", hyperballs_program);
    }

    void rendering_engine::initShaders(const std::filesystem::path vert, const std::filesystem::path frag, GLuint& program)
    {
        std::string vertex = loadShaderSource(vert);
        std::string fragment = loadShaderSource(frag);

        GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
        const char *vertex_src = vertex.c_str();
        glShaderSource(vertex_shader, 1, &vertex_src, nullptr);
        glCompileShader(vertex_shader);

        GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
        const char *fragment_src = fragment.c_str();
        glShaderSource(fragment_shader, 1, &fragment_src, nullptr);
        glCompileShader(fragment_shader);

        program = glCreateProgram();
        glAttachShader(program, vertex_shader);
        glAttachShader(program, fragment_shader);
        glLinkProgram(program);

        GLint success = 0;
        glGetProgramiv(program, GL_LINK_STATUS, &success);
        if (!success)
        {
            char info[1024] = {};
            glGetProgramInfoLog(program, 1024, nullptr, info);
            std::cerr << info << "\n";
            throw std::runtime_error(std::string("[Shader link failed]: ") + info);
        }

        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);
    }

    void rendering_engine::bindColor(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info)
    {
        if (info.hyperBalls) return;

        std::vector<AtomInstance> instances;
        instances.reserve(sim_info.atoms.size());

        if (info.lennardBall || info.licorice)
        {
            for (int32_t i = 0; i < sim_info.atoms.size(); ++i)
            {
                if (sim_info.positions.size() < i)
                    break;

                const auto &atom = sim_info.atoms[i];
                float radius = 0.f;

                if (info.licorice)
                    radius = licorice_radius;
                else
                    radius = info.spaceFilling
                                 ? constants::VDW_RADII[atom.ZIndex] * 0.8f
                                 : constants::covalent_radius[atom.ZIndex] * 0.7f;

                sf::Color col = constants::getElementColor(atom.ZIndex) + info.color_addition;
                glm::vec4 color_norm(col.r / 255.f, col.g / 255.f, col.b / 255.f, info.opacity);

                instances.emplace_back(glm::vec3(glm::vec4(sim_info.positions[i], 1.0)), radius, color_norm);
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, color_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     instances.size() * sizeof(AtomInstance),
                     instances.data(),
                     GL_DYNAMIC_DRAW);

        glBindVertexArray(color_vao);
        glUseProgram(atom_program);

        GLint loc_projection = glGetUniformLocation(atom_program, "u_proj");
        GLint loc_view = glGetUniformLocation(atom_program, "u_view");

        if (loc_projection != -1)
            glUniformMatrix4fv(loc_projection, 1, GL_FALSE, glm::value_ptr(cam.getProjectionMatrix(target)));

        if (loc_view != -1)
            glUniformMatrix4fv(loc_view, 1, GL_FALSE, glm::value_ptr(cam.getViewMatrix()));

        glEnable(GL_DEPTH_TEST);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, static_cast<GLsizei>(instances.size()));

        GLenum err = glGetError();
        if (err != GL_NO_ERROR)
        {
            std::cerr << "[RENDERING ENGINE] OpenGL error after draw: 0x" << std::hex << err << "\n";
        }

        drawHighlight(target, info, sim_info);
    }

    void rendering_engine::bindBond(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info)
    {
        if (info.spaceFilling || (!info.lennardBall && !info.licorice && !info.hyperBalls))
            return;

        std::vector<BondInstance> instances;
        instances.reserve(sim_info.bonds.size() * 4);

        for (int32_t i = 0; i < sim_info.bonds.size(); ++i)
        {
            const auto &bond = sim_info.bonds[i];

            sf::Color colA = constants::getElementColor(sim_info.atoms[bond.bondedAtom].ZIndex) + info.color_addition;
            sf::Color colB = constants::getElementColor(sim_info.atoms[bond.centralAtom].ZIndex) + info.color_addition;

            glm::vec4 colorA_norm(colA.r / 255.f, colA.g / 255.f, colA.b / 255.f, info.opacity);
            glm::vec4 colorB_norm(colB.r / 255.f, colB.g / 255.f, colB.b / 255.f, info.opacity);

            int32_t order = static_cast<int32_t>(bond.type);

            glm::vec3 posA = sim_info.positions[bond.bondedAtom];
            glm::vec3 posB = sim_info.positions[bond.centralAtom];

            glm::vec3 r_vec = posA - posB;
            if (glm::length(r_vec) > 5.f)
                continue;

            /* if (info.licorice)
            {
                instances.emplace_back(
                glm::vec4(posA, 1.0f),
                glm::vec4(posB, 1.0f),
                colorA_norm, colorB_norm, licorice_radius);
                
                continue;
            } */

            float radiusA = constants::covalent_radius[sim_info.atoms[bond.bondedAtom].ZIndex];
            float radiusB = constants::covalent_radius[sim_info.atoms[bond.centralAtom].ZIndex];
            float baseRadius = 0.15f;
            float bondR = info.licorice ? licorice_radius / order : baseRadius / order;
            glm::vec3 bondDir = glm::normalize(r_vec);

            glm::vec3 arbitrary = (std::abs(bondDir.x) < 0.9f) ? glm::vec3(1, 0, 0) : glm::vec3(0, 1, 0);
            glm::vec3 perp1 = glm::normalize(glm::cross(bondDir, arbitrary));
            glm::vec3 perp2 = glm::cross(bondDir, perp1);

            const float offsetStep = 0.10f;

            if (order == 1 || info.hyperBalls)
            {
                instances.emplace_back(
                    glm::vec4(posA, 1.0f),
                    glm::vec4(posB, 1.0f),
                    colorA_norm, colorB_norm, bondR, radiusA, radiusB);
            }
            else
            {
                for (int k = 0; k < order; ++k)
                {
                    float angle = (static_cast<float>(k) / order) * glm::two_pi<float>();
                    float offsetAmount = offsetStep * (order == 2 ? 1.0f : 1.0f);

                    glm::vec3 offset = offsetAmount * (std::cos(angle) * perp1 + std::sin(angle) * perp2);

                    glm::vec3 offsetposA = posA + offset;
                    glm::vec3 offsetposB = posB + offset;

                    instances.emplace_back(
                        glm::vec4(offsetposA, 1.0f),
                        glm::vec4(offsetposB, 1.0f),
                        colorA_norm, colorB_norm, bondR, radiusA, radiusB);
                }
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, bond_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     instances.size() * sizeof(BondInstance),
                     instances.data(),
                     GL_DYNAMIC_DRAW);

        glBindVertexArray(bond_vao);

        info.hyperBalls ? glUseProgram(hyperballs_program) : glUseProgram(bond_program);

        GLint loc_projection = glGetUniformLocation(bond_program, "u_proj");
        GLint loc_view = glGetUniformLocation(bond_program, "u_view");
        GLint loc_licorice = glGetUniformLocation(bond_program, "licorice");

        if (loc_projection != -1)
            glUniformMatrix4fv(loc_projection, 1, GL_FALSE, glm::value_ptr(cam.getProjectionMatrix(target)));

        if (loc_view != -1)
            glUniformMatrix4fv(loc_view, 1, GL_FALSE, glm::value_ptr(cam.getViewMatrix()));

        if (!info.hyperBalls) 
        {
            GLint loc_licorice = glGetUniformLocation(bond_program, "licorice");
            if (loc_licorice != -1)
                glUniform1i(loc_licorice, static_cast<uint8_t>(info.licorice));
        }

        glEnable(GL_DEPTH_TEST);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, static_cast<GLsizei>(instances.size()));

        GLenum err = glGetError();
        if (err != GL_NO_ERROR)
        {
            std::cerr << "[RENDERING ENGINE] OpenGL error after draw: 0x" << std::hex << err << "\n";
        }
    }

    void rendering_engine::drawHighlight(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info)
    {
        if (!info.flag_highlights)
            return;
        drawAtomHighlight(target, info, sim_info);
        drawBondHighlight(target, info, sim_info);
    }

    void rendering_engine::drawBondHighlight(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info)
    {
        if (info.spaceFilling || !info.lennardBall)
            return;

        std::vector<BondInstance> instances;
        instances.reserve(info.highlight_bonds.size() * 4);

        for (int32_t i = 0; i < info.highlight_bonds.size(); ++i)
        {
            const auto &predict_bond = info.highlight_bonds[i];
            const auto &bond = *std::find_if(sim_info.bonds.begin(), sim_info.bonds.end(), [&](const sim::fun::bond &b)
                                             { return b.bondedAtom == predict_bond.first, b.centralAtom == predict_bond.second; });
            glm::vec4 col(1.0, 1.0, 0.0, 0.3);

            int32_t order = static_cast<int32_t>(bond.type);

            glm::vec3 posA = sim_info.positions[bond.bondedAtom];
            glm::vec3 posB = sim_info.positions[bond.centralAtom];

            glm::vec3 r_vec = posA - posB;
            if (glm::length(r_vec) > 5.f)
                continue;

            float baseRadius = 0.2f;
            float bondR = baseRadius / order;
            glm::vec3 bondDir = glm::normalize(r_vec);

            glm::vec3 arbitrary = (std::abs(bondDir.x) < 0.9f) ? glm::vec3(1, 0, 0) : glm::vec3(0, 1, 0);
            glm::vec3 perp1 = glm::normalize(glm::cross(bondDir, arbitrary));
            glm::vec3 perp2 = glm::cross(bondDir, perp1);

            const float offsetStep = 0.10f;

            if (order == 1)
            {
                instances.emplace_back(
                    glm::vec4(posA, 1.0f),
                    glm::vec4(posB, 1.0f),
                    col, col, bondR);
            }
            else
            {
                for (int k = 0; k < order; ++k)
                {
                    float angle = (static_cast<float>(k) / order) * glm::two_pi<float>();
                    float offsetAmount = offsetStep * (order == 2 ? 1.0f : 1.0f);

                    glm::vec3 offset = offsetAmount * (std::cos(angle) * perp1 + std::sin(angle) * perp2);

                    glm::vec3 offsetposA = posA + offset;
                    glm::vec3 offsetposB = posB + offset;

                    instances.emplace_back(
                        glm::vec4(offsetposA, 1.0f),
                        glm::vec4(offsetposB, 1.0f),
                        col, col, bondR);
                }
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, bond_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     instances.size() * sizeof(BondInstance),
                     instances.data(),
                     GL_DYNAMIC_DRAW);

        glBindVertexArray(bond_vao);
        glUseProgram(bond_program);

        GLint loc_projection = glGetUniformLocation(bond_program, "u_proj");
        GLint loc_view = glGetUniformLocation(bond_program, "u_view");

        if (loc_projection != -1)
            glUniformMatrix4fv(loc_projection, 1, GL_FALSE, glm::value_ptr(cam.getProjectionMatrix(target)));

        if (loc_view != -1)
            glUniformMatrix4fv(loc_view, 1, GL_FALSE, glm::value_ptr(cam.getViewMatrix()));

        glEnable(GL_DEPTH_TEST);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, static_cast<GLsizei>(instances.size()));

        GLenum err = glGetError();
        if (err != GL_NO_ERROR)
        {
            std::cerr << "[RENDERING ENGINE] OpenGL error after draw: 0x" << std::hex << err << "\n";
        }
    }

    void rendering_engine::drawAtomHighlight(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info)
    {
        std::vector<AtomInstance> instances;
        instances.reserve(info.highlight_indices.size());

        for (int32_t i = 0; i < info.highlight_indices.size(); ++i)
        {
            auto &highlight = info.highlight_indices[i];

            if (sim_info.positions.size() < highlight)
                break;

            const auto &atom = sim_info.atoms[highlight];
            float radius = info.spaceFilling
                               ? constants::VDW_RADII[atom.ZIndex]
                               : constants::covalent_radius[atom.ZIndex];

            glm::vec4 color_norm(1.f, 1.f, 0.f, 0.3f);

            instances.emplace_back(glm::vec3(glm::vec4(sim_info.positions[i], 1.0)), radius, color_norm);
        }

        glBindBuffer(GL_ARRAY_BUFFER, color_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     instances.size() * sizeof(AtomInstance),
                     instances.data(),
                     GL_DYNAMIC_DRAW);

        glBindVertexArray(color_vao);
        glUseProgram(atom_program);

        GLint loc_projection = glGetUniformLocation(atom_program, "u_proj");
        GLint loc_view = glGetUniformLocation(atom_program, "u_view");

        if (loc_projection != -1)
            glUniformMatrix4fv(loc_projection, 1, GL_FALSE, glm::value_ptr(cam.getProjectionMatrix(target)));

        if (loc_view != -1)
            glUniformMatrix4fv(loc_view, 1, GL_FALSE, glm::value_ptr(cam.getViewMatrix()));

        glEnable(GL_DEPTH_TEST);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, static_cast<GLsizei>(instances.size()));

        GLenum err = glGetError();
        if (err != GL_NO_ERROR)
        {
            std::cerr << "[RENDERING ENGINE] OpenGL error after draw: 0x" << std::hex << err << "\n";
        }
    }

    void rendering_engine::draw(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info)
    {
        if (info.universeBox)
            drawBox(sim_info.box);

        bindBond(target, info, sim_info);
        bindColor(target, info, sim_info);

        glBindVertexArray(0);
        glUseProgram(0);
    }

    float smootherstep(float edge0, float edge1, float x)
    {
        x = std::clamp((x - edge0) / (edge1 - edge0), 0.f, 1.f);

        return x * x * (3.0f - 2.0f * x);
    }

    void rendering_engine::drawChargeField(sf::RenderTarget &target, const fun::rendering_simulation_info &sim_info)
    {
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
