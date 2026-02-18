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
    rendering_engine::rendering_engine(core::window_t &window)
        : window(window)
    {
        loadProgram("box",        "box.vert",        "box.frag");
        loadProgram("atom",       "color.vert",      "color.frag");
        loadProgram("bond",       "bond.vert",       "bond.frag");
        loadProgram("hyperballs", "hyper_balls.vert","hyper_balls.frag");
        loadProgram("arrow",      "arrow.vert",      "arrow.frag");

        glGenVertexArrays(1, &box_vao);
        glGenBuffers(1, &box_vbo);
        glBindVertexArray(box_vao);
        glBindBuffer(GL_ARRAY_BUFFER, box_vbo);

        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
        glVertexAttribDivisor(0, 0);

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
        
        glGenVertexArrays(1, &arrow_vao);
        glGenBuffers(1, &arrow_vbo);
        glBindVertexArray(arrow_vao);
        glBindBuffer(GL_ARRAY_BUFFER, arrow_vbo);

        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(ArrowInstance), (void *)offsetof(ArrowInstance, posA));
        glVertexAttribDivisor(0, 1);

        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(ArrowInstance), (void *)offsetof(ArrowInstance, posB));
        glVertexAttribDivisor(1, 1);

        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, sizeof(ArrowInstance), (void *)offsetof(ArrowInstance, color));
        glVertexAttribDivisor(2, 1);

        glEnableVertexAttribArray(3);
        glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(ArrowInstance), (void *)offsetof(ArrowInstance, radius));
        glVertexAttribDivisor(3, 1);
    }

    void rendering_engine::loadProgram(const std::string key, const std::string& vert, const std::string& frag)
    {
        programs.emplace(key,
            core::glProgram{
                core::glShader{GL_VERTEX_SHADER,   "resource/shaders/" + vert},
                core::glShader{GL_FRAGMENT_SHADER, "resource/shaders/" + frag}
            });
    }

    void rendering_engine::drawBox(const glm::vec3 &box, sf::RenderTarget &target)
    {
        const std::array<glm::vec3, 8> corners = {{
            {0.f, 0.f, 0.f},
            {box.x, 0.f, 0.f},
            {box.x, box.y, 0.f},
            {0.f, box.y, 0.f},
            {0.f, 0.f, box.z},
            {box.x, 0.f, box.z},
            {box.x, box.y, box.z},
            {0.f, box.y, box.z}
        }};

        const std::array<std::pair<int32_t, int32_t>, 12> edges = {{
            {0, 1}, {1, 2}, {2, 3}, {3, 0}, // bottom
            {4, 5}, {5, 6}, {6, 7}, {7, 4}, // top
            {0, 4}, {1, 5}, {2, 6}, {3, 7}  // verticals
        }};

        std::vector<glm::vec3> lineVertices;
        lineVertices.reserve(edges.size() * 2);
        for (auto [i, j] : edges) 
        {
            lineVertices.push_back(corners[i]);
            lineVertices.push_back(corners[j]);
        }

        glBindBuffer(GL_ARRAY_BUFFER, box_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                    lineVertices.size() * sizeof(glm::vec3),
                    lineVertices.data(),
                    GL_DYNAMIC_DRAW);

        glBindVertexArray(box_vao);

        auto& box_program = programs["box"];
        box_program.use();

        box_program.setUniform("u_proj", cam.getProjectionMatrix(target.getView().getSize().x, target.getView().getSize().y));
        box_program.setUniform("u_view", cam.getViewMatrix());

        glEnable(GL_DEPTH_TEST);
        glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(lineVertices.size()));

        GLenum err = glGetError();
        if (err != GL_NO_ERROR) 
        {
            std::cerr << "[RENDERING ENGINE] OpenGL error after draw: 0x"
                    << std::hex << err << "\n";
        }
    }

    void rendering_engine::bindColor(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info)
    {
        if (info.hyperBalls || info.wireframe && !(info.lennardBall || info.licorice)) return;

        std::vector<AtomInstance> instances;
        instances.reserve(sim_info.atoms.size());

        for (int32_t i = 0; i < sim_info.atoms.size(); ++i)
        {
            if (sim_info.positions.size() - 1 < i)
                break;

            const auto &atom = sim_info.atoms[i];
            float radius = 0.f;

            if (info.licorice)
                radius = licorice_radius * 1.18f;
            else
                radius = info.spaceFilling
                             ? constants::VDW_RADII[atom.ZIndex]
                              : constants::covalent_radius[atom.ZIndex] * 0.5f;

            glm::vec4 col_addition = glm::vec4(info.color_addition.x, info.color_addition.y, info.color_addition.z, 1.0);
            glm::vec4 col = getAtomColor(info, sim_info, i);

            col += col_addition;
            col.w = info.opacity;

            instances.emplace_back(glm::vec3(glm::vec4(sim_info.positions[i], 1.0)), radius, col);
        }

        glBindBuffer(GL_ARRAY_BUFFER, color_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     instances.size() * sizeof(AtomInstance),
                     instances.data(),
                     GL_DYNAMIC_DRAW);

        glBindVertexArray(color_vao);

        auto& atom_program = programs["atom"];
        atom_program.use();

        atom_program.setUniform("u_proj", cam.getProjectionMatrix(target.getView().getSize().x, target.getView().getSize().y));
        atom_program.setUniform("u_view", cam.getViewMatrix());

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

            if (sim_info.positions.size() - 1 < bond.bondedAtom || sim_info.positions.size() - 1 < bond.centralAtom) break;

            glm::vec4 col_addition = glm::vec4(info.color_addition.x, info.color_addition.y, info.color_addition.z, 1.0);
            glm::vec4 colA = getAtomColor(info, sim_info, bond.bondedAtom);
            glm::vec4 colB = getAtomColor(info, sim_info, bond.centralAtom);

            colA += col_addition;
            colB += col_addition;

            colA.w = info.opacity;
            colB.w = info.opacity;

            int32_t order = static_cast<int32_t>(bond.type);

            glm::vec3 posA = sim_info.positions[bond.bondedAtom];
            glm::vec3 posB = sim_info.positions[bond.centralAtom];

            glm::vec3 r_vec = posA - posB;
            if (glm::length(r_vec) > 5.f)
                continue;

            float radiusA = constants::covalent_radius[sim_info.atoms[bond.bondedAtom].ZIndex];
            float radiusB = constants::covalent_radius[sim_info.atoms[bond.centralAtom].ZIndex];
            float baseRadius = 0.10f;
            float bondR = info.licorice ? licorice_radius / order : baseRadius / order;
            glm::vec3 bondDir = glm::normalize(r_vec);

            glm::vec3 arbitrary = (std::abs(bondDir.x) < 0.9f) ? glm::vec3(1, 0, 0) : glm::vec3(0, 1, 0);
            glm::vec3 perp1 = glm::normalize(glm::cross(bondDir, arbitrary));
            glm::vec3 perp2 = glm::cross(bondDir, perp1);

            const float offsetStep = 0.08f;

            if (order == 1 || info.hyperBalls)
            {
                instances.emplace_back(
                    glm::vec4(posA, 1.0f),
                    glm::vec4(posB, 1.0f),
                    colA, colB, bondR, radiusA, radiusB);
            }
            else
            {
                for (int k = 0; k < order; ++k)
                {
                    float angle = (static_cast<float>(k) / order) * glm::two_pi<float>();
                    float offsetAmount = offsetStep * (order == 2 ? 1.1f : 1.0f);

                    glm::vec3 offset = offsetAmount * (std::cos(angle) * perp1 + std::sin(angle) * perp2);

                    glm::vec3 offsetposA = posA + offset;
                    glm::vec3 offsetposB = posB + offset;

                    instances.emplace_back(
                        glm::vec4(offsetposA, 1.0f),
                        glm::vec4(offsetposB, 1.0f),
                        colA, colB, bondR, radiusA, radiusB);
                }
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, bond_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     instances.size() * sizeof(BondInstance),
                     instances.data(),
                     GL_DYNAMIC_DRAW);

        glBindVertexArray(bond_vao);

        auto& bond_program = programs["bond"];
        info.hyperBalls ? programs["hyper_balls"].use() : bond_program.use();

        bond_program.setUniform("u_proj", cam.getProjectionMatrix(target.getView().getSize().x, target.getView().getSize().y));
        bond_program.setUniform("u_view", cam.getViewMatrix());

        if (!info.hyperBalls) 
            bond_program.setUniform("licorice", static_cast<uint8_t>(info.licorice));

        glEnable(GL_DEPTH_TEST);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, static_cast<GLsizei>(instances.size()));

        GLenum err = glGetError();
        if (err != GL_NO_ERROR)
        {
            std::cerr << "[RENDERING ENGINE] OpenGL error after draw: 0x" << std::hex << err << "\n";
        }
    }

    void rendering_engine::bindArrow(sf::RenderTarget &target, const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info)
    {
        if (!info.flag_arrows) return;

        std::vector<ArrowInstance> instances;
        instances.reserve(info.arrows.size());

        for (int32_t i = 0; i < info.arrows.size(); ++i)
        {
            auto& positions = info.arrows[i];

            ArrowInstance nInstance{};
            nInstance.color = glm::vec4(0.8, 0.8, 0.8, 1.0);
            nInstance.posA = glm::vec4(positions.first, 1.0);
            nInstance.posB = glm::vec4(positions.second, 1.0);
            nInstance.radius = .1f;

            instances.emplace_back(std::move(nInstance));
        }

        glBindBuffer(GL_ARRAY_BUFFER, arrow_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     instances.size() * sizeof(ArrowInstance),
                     instances.data(),
                     GL_DYNAMIC_DRAW);

        glBindVertexArray(arrow_vao);

        auto& arrow_program = programs["arrow"];
        arrow_program.use();

        arrow_program.setUniform("u_proj", cam.getProjectionMatrix(target.getView().getSize().x, target.getView().getSize().y));
        arrow_program.setUniform("u_view", cam.getViewMatrix());

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

        auto& bond_program = programs["bond"];
        bond_program.use();

        bond_program.setUniform("u_proj", cam.getProjectionMatrix(target.getView().getSize().x, target.getView().getSize().y));
        bond_program.setUniform("u_view", cam.getViewMatrix());

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

        auto& atom_program = programs["atom"];
        atom_program.use();

        atom_program.setUniform("u_proj", cam.getProjectionMatrix(target.getView().getSize().x, target.getView().getSize().y));
        atom_program.setUniform("u_view", cam.getViewMatrix());

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
            drawBox(sim_info.box, target);

        bindBond(target, info, sim_info);
        bindColor(target, info, sim_info);
        bindArrow(target, info, sim_info);

        glBindVertexArray(0);
        glUseProgram(0);
    }

    float smootherstep(float edge0, float edge1, float x)
    {
        x = std::clamp((x - edge0) / (edge1 - edge0), 0.f, 1.f);

        return x * x * (3.0f - 2.0f * x);
    }

    void rendering_engine::drawAngle(ImDrawList* draw_list,
                                    const glm::vec3& a,
                                    const glm::vec3& b,
                                    const glm::vec3& c)
    {
        constexpr ImU32 arc_color  = IM_COL32(220, 60, 40, 220);
        constexpr ImU32 text_color = IM_COL32(255, 255, 255, 255);
        constexpr int arc_segments = 32;
        
        glm::vec3 cam_pos = cam.eye();
        glm::vec3& cam_target = cam.target;
        float dist_to_cam = glm::length(b - cam_pos);

        float arc_radius = 450.0f / dist_to_cam;
        float text_offset = 80.0f / dist_to_cam;

        glm::vec2 screen_a = project(a);
        glm::vec2 screen_b = project(b);
        glm::vec2 screen_c = project(c);

        if (glm::length(screen_a) < 1e-6f || glm::length(screen_b) < 1e-6f || glm::length(screen_c) < 1e-6f)
            return;

        glm::vec2 dir_ba = glm::normalize(screen_a - screen_b);
        glm::vec2 dir_bc = glm::normalize(screen_c - screen_b);

        float angle_rad = std::acosf(glm::clamp(glm::dot(dir_ba, dir_bc), -1.0f, 1.0f));

        float cross = dir_ba.x * dir_bc.y - dir_ba.y * dir_bc.x;
        bool clockwise = true;

        float start_angle = atan2f(dir_ba.y, dir_ba.x);
        float delta_angle = clockwise ? -angle_rad : angle_rad;
        float end_angle   = start_angle + delta_angle;

        if (angle_rad > 179.0f || angle_rad < 181.0f)
            delta_angle = -angle_rad;

        draw_list->PathArcTo(
            ImVec2(screen_b.x, screen_b.y),
            arc_radius,
            start_angle,
            end_angle,
            arc_segments
        );
        draw_list->PathStroke(arc_color, false, 2.5f);

        float mid_angle = start_angle + delta_angle * 0.5f;
        glm::vec2 text_pos = screen_b + glm::vec2(cosf(mid_angle), sinf(mid_angle)) * (arc_radius + text_offset);

        glm::vec3 dir_ba_world = glm::normalize(a - b);
        glm::vec3 dir_bc_world = glm::normalize(c - b);

        float angle_rad_world = std::acosf(glm::clamp(glm::dot(dir_ba_world, dir_bc_world), -1.0f, 1.0f));

        char text[16];
        std::snprintf(text, sizeof(text), "%.1f°", glm::degrees(angle_rad_world));

        ImVec2 text_size = ImGui::CalcTextSize(text);
        draw_list->AddText(
            ImVec2(text_pos.x - text_size.x * 0.5f, text_pos.y - text_size.y * 0.5f),
            text_color,
            text
        );
    }

    void rendering_engine::drawChargeField(sf::RenderTarget &target, const fun::rendering_simulation_info &sim_info)
    {
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

    glm::vec3 hsv2rgb(glm::vec3 c) 
    {
        glm::vec3 rgb = glm::clamp(glm::abs(glm::mod(c.x * 6.0f + glm::vec3(0,4,2), 6.0f) - 3.0f) - 1.0f,
                                0.0f, 1.0f);
        return c.z * glm::mix(glm::vec3(0.8f), rgb, c.y);
    }

    glm::vec4 rendering_engine::getAtomColor(const fun::rendering_info &info, const fun::rendering_simulation_info &sim_info, const uint32_t i)
    {
        const auto &atom = sim_info.atoms[i];
        
        sf::Color default_color = constants::getElementColor(atom.ZIndex);
        glm::vec4 default_norm(default_color.r / 255.f, default_color.g / 255.f, default_color.b / 255.f, info.opacity);

        switch(info.color_mode)
        {
            case fun::color_rendering_mode::COLOR:
                return default_norm;
            case fun::color_rendering_mode::CHARGE:
            {
                float t = (sim_info.q[i] + 1.0) * 0.5;
                glm::vec3 color = glm::mix(glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(1.0f, 0.0f, 0.0f), t);

                return glm::vec4(color, 1.0f);
            }
            case fun::color_rendering_mode::VELOCITY:
            {
                constexpr float d_velocity = 40.f;
                float v = glm::length(sim_info.velocities[i]);
                float t = std::clamp(v / d_velocity, 0.0f, 1.0f);

                float hue = glm::mix(0.0f, 0.75f, t); 
                glm::vec3 rgb = hsv2rgb(glm::vec3(hue, 1.0f, 1.0f)); 
                return glm::vec4(rgb, 1.0f);
            }
            default:
                return default_norm;
        }
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
