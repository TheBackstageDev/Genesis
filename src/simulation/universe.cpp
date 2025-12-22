#include "universe.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <queue>
#include <fstream>

#include "constants.hpp"

namespace sim
{
    namespace fun
    {

        universe::universe(const universe_create_info &create_info)
            : box(create_info.box), react(create_info.reactive),
              gravity(create_info.has_gravity), mag_gravity(create_info.mag_gravity),
              wall_collision(create_info.wall_collision), isothermal(create_info.isothermal),
              render_water(create_info.render_water), HMassRepartitioning(create_info.HMassRepartitioning), 
              log_flags(create_info.log_flags)
        {
            if (react)
                initReaxParams();
        }

        universe::universe(const std::filesystem::path path)
        {
            loadScene(path);
        }

        void universe::drawBox(core::window_t &window)
        {
            const std::array<sf::Vector3f, 8> corners =
                {{{0.f, 0.f, 0.f},
                  {box.x, 0.f, 0.f},
                  {box.x, box.y, 0.f},
                  {0.f, box.y, 0.f},
                  {0.f, 0.f, box.z},
                  {box.x, 0.f, box.z},
                  {box.x, box.y, box.z},
                  {0.f, box.y, box.z}}};

            const std::array<std::pair<int, int>, 12> edges =
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
                sf::Vector2f p1 = project(window, corners[i]);
                sf::Vector2f p2 = project(window, corners[j]);

                if (p1.x < -1000.f || p2.x < -1000.f)
                    continue;

                const sf::Vector3f &a = corners[i];
                const sf::Vector3f &b = corners[j];

                glm::vec3 cam_eye = cam.eye();
                sf::Vector3f camToA = a - sf::Vector3f(cam_eye.x, cam_eye.y, cam_eye.z);
                sf::Vector3f camToB = b - sf::Vector3f(cam_eye.x, cam_eye.y, cam_eye.z);

                float dot = camToA.normalized().dot(camToB.normalized());
                const float THRESHOLD = 0.3f;
                if (dot < THRESHOLD)
                    continue;

                lines[idx++] = sf::Vertex(p1, edgeColor);
                lines[idx++] = sf::Vertex(p2, edgeColor);
            }

            window.getWindow().draw(lines);
        }

        void universe::draw(core::window_t &window, sf::RenderTarget& target, const rendering_info info)
        {
            std::vector<int32_t> drawOrder(atoms.size());
            std::iota(drawOrder.begin(), drawOrder.end(), 0);
            sf::Vector3f eye = sf::Vector3f(cam.eye().x, cam.eye().y, cam.eye().z);

            std::sort(drawOrder.begin(), drawOrder.end(),
                      [&](int32_t a, int32_t b)
                      { return (pos(a) - eye).lengthSquared() > (pos(b) - eye).lengthSquared(); });

            if (info.universeBox)
                drawBox(window);

            std::vector<int32_t> no_draw{};

            if (!render_water)
                for (auto &m : molecules)
                {
                    if (m.exclude)
                    {
                        no_draw.emplace_back(m.angleBegin);
                        no_draw.emplace_back(m.angleBegin + 1);
                        no_draw.emplace_back(m.angleBegin + 2);
                    }
                }

            if (!react && !info.spaceFilling)
            {
                drawBonds(window, target, no_draw);
                //drawRings(window, target);
            }

            for (int32_t i = 0; i < atoms.size(); ++i)
            {
                if (std::find(no_draw.begin(), no_draw.end(), drawOrder[i]) != no_draw.end())
                    continue;

                sf::Vector2f p = project(window, pos(drawOrder[i]));
                if (p.x < -1000)
                    continue;

                float temp = calculateAtomTemperature(drawOrder[i]);
                if (atoms[drawOrder[i]].ZIndex == 1)
                    drawHydrogenBond(window, target, drawOrder[i]);

                glm::vec3 cam_eye = cam.eye();
                float camDistance = (pos(drawOrder[i]) - sf::Vector3f(cam_eye.x, cam_eye.y, cam_eye.z)).length();

                float T = calculateAtomTemperature(drawOrder[i]);
                atoms[drawOrder[i]].draw(T, p, camDistance, data.q[drawOrder[i]], window, target, info);
            }
        }

        float smootherstep(float edge0, float edge1, float x)
        {
            x = std::clamp((x - edge0) / (edge1 - edge0), 0.f, 1.f);

            return x * x * (3.0f - 2.0f * x);
        }

        void universe::drawChargeField(core::window_t &window, sf::RenderTarget& target)
        {
            if (!(cx > 0 && cy > 0 && cz > 0))
                return;

            sf::RenderStates states;
            states.blendMode = sf::BlendAlpha;

            constexpr int32_t GRID = 50;
            float step = box.length() / GRID;

            const float isovalue_pos = 1.f;
            const float isovalue_neg = -1.f;
            const float thickness = .7f;

            struct Blob
            {
                sf::Vector2f screenPos;
                float depth;
                float size;
                sf::Color color;
            };

            std::vector<Blob> blobs;
            blobs.reserve(GRID * GRID * GRID);

            for (int iz = 0; iz < GRID; ++iz)
                for (int iy = 0; iy < GRID; ++iy)
                    for (int ix = 0; ix < GRID; ++ix)
                    {
                        sf::Vector3f worldPos(
                            (ix + 0.5f) * step,
                            (iy + 0.5f) * step,
                            (iz + 0.5f) * step);

                        sf::Vector2f screen = project(window, worldPos);
                        if (screen.x < -500)
                            continue;

                        float density = 0.0f;

                        for (int32_t d = 0; d < 14; ++d)
                        {
                            auto &cell = cells[getCellID(ix + offsets[d][0], iy + offsets[d][1], iz + offsets[d][2])];
                            for (int32_t i = 0; i < cell.size(); ++i)
                            {
                                sf::Vector3f dr = minImageVec(worldPos - pos(cell[i]));
                                float r2 = dr.lengthSquared();
                                density += data.q[cell[i]] * std::exp(-r2 / (1.2f * 1.2f));
                            }
                        }

                        if (std::abs(density) < 0.1f)
                            continue;

                        float dist_to_pos = std::abs(density - isovalue_pos);
                        float dist_to_neg = std::abs(density - isovalue_neg);
                        float dist = std::min(dist_to_pos, dist_to_neg);

                        if (dist > thickness)
                            continue;

                        sf::Color col;
                        if (density <= -1.f)
                            col = sf::Color(220, 20, 60, 180); // deep red
                        else if (density <= -0.5f)
                            col = sf::Color(255, 80, 80, 180); // red
                        else if (density <= -0.2f)
                            col = sf::Color(255, 147, 147, 150); // mild red / light pink
                        else if (density < 0.2f)
                            col = sf::Color(120, 240, 150, 200); // green
                        else if (density < 0.5f)
                            col = sf::Color(135, 206, 255, 160); // light sky blue
                        else if (density < 1.f)
                            col = sf::Color(70, 130, 255, 180); // medium blue
                        else
                            col = sf::Color(30, 70, 200, 180); // deep royal blue

                        float alpha = 1.0f - smootherstep(0.0f, thickness, dist);
                        float t = (density - isovalue_neg) / (isovalue_pos - isovalue_neg);
                        t = std::clamp(t, 0.0f, 1.0f);

                        float opacity = 1.0f - smootherstep(0.0f, thickness, dist);

                        col.a = static_cast<uint8_t>(255 * opacity);

                        glm::vec3 cam_eye = cam.eye();
                        float distance = (worldPos - sf::Vector3f(cam_eye.x, cam_eye.y, cam_eye.z)).length();
                        float size = step * GRID / distance;

                        blobs.emplace_back(screen, distance, size, col);
                    }

            std::sort(blobs.begin(), blobs.end(), [=](const Blob &a, const Blob &b)
                      { return a.depth > b.depth; });

            sf::VertexArray billboards(sf::PrimitiveType::Triangles);
            billboards.resize(blobs.size() * 6);

            int32_t idx = 0;
            for (const auto &b : blobs)
            {
                sf::Vector2f p = b.screenPos;
                float s = b.size;
                sf::Color c = b.color;

                billboards[idx + 0] = {p + sf::Vector2f(-s, -s), c};
                billboards[idx + 1] = {p + sf::Vector2f(s, -s), c};
                billboards[idx + 2] = {p + sf::Vector2f(s, s), c};

                billboards[idx + 3] = {p + sf::Vector2f(-s, -s), c};
                billboards[idx + 4] = {p + sf::Vector2f(s, s), c};
                billboards[idx + 5] = {p + sf::Vector2f(-s, s), c};

                idx += 6;
            }
            target.draw(billboards, states);
        }

        void universe::drawBonds(core::window_t &window, sf::RenderTarget& target, const std::vector<int32_t> &no_draw)
        {
            sf::Vector2f dimensions = window.getWindow().getView().getSize();

            for (int32_t b = 0; b < bonds.size(); ++b)
            {
                const bond &bond = bonds[b];
                if (std::find(no_draw.begin(), no_draw.end(), bond.centralAtom) != no_draw.end())
                    continue;

                const sf::Vector3f &pCentral = pos(bond.centralAtom); // pB
                const sf::Vector3f &pBonded = pos(bond.bondedAtom);   // pA

                sf::Vector2f s1 = cam.project(glm::vec3(pCentral.x, pCentral.y, pCentral.z), dimensions.x, dimensions.y); // center
                sf::Vector2f s2 = cam.project(glm::vec3(pBonded.x, pBonded.y, pBonded.z), dimensions.x, dimensions.y);  // end

                if (s1.x <= -9999 || s2.x <= -9999)
                    continue;
                if ((pCentral - pBonded).length() > 5.f)
                    continue;
                glm::vec3 cam_eye = cam.eye();
                if ((pCentral - sf::Vector3f(cam_eye.x, cam_eye.y, cam_eye.z)).length() > 200.f)
                    continue;

                sf::Vector2f dir = s2 - s1;
                float len = dir.length();
                if (len < 0.5f)
                    continue;
                dir /= len;

                sf::Vector2f perp{-dir.y, dir.x};

                uint8_t lines = static_cast<uint8_t>(bond.type);
                float distance = (pCentral - sf::Vector3f(cam_eye.x, cam_eye.y, cam_eye.z)).length();

                float shrink = 5.f / distance;   // pixels
                float spacing = 10.f / distance; // pixels between lines

                std::vector<sf::Vertex> verts;
                verts.reserve(lines * 2);

                for (int32_t i = 0; i < lines; ++i)
                {
                    float offset = spacing * (i - (lines - 1) * 0.5f);

                    sf::Vector2f off = perp * offset;

                    sf::Vector2f start = s1 + dir * shrink + off;
                    sf::Vector2f end = s2 - dir * shrink + off;

                    verts.emplace_back(start, sf::Color::White);
                    verts.emplace_back(end, sf::Color::White);
                }

                target.draw(verts.data(), verts.size(), sf::PrimitiveType::Lines);
            }
        }

        void universe::highlightAtom(core::window_t& window, size_t i)
        {
            const atom &a = atoms[i];

            constexpr uint32_t pointCount = 60;
            constexpr float thickness = 4.0f;

            float ringRadius = a.radius * 1.2f;

            sf::Vector3f worldPos = pos(i);
            sf::Vector2f currentPos = project(window, worldPos);

            glm::vec3 cam_eye = cam.eye();
            float distance = (pos(i) - sf::Vector3f(cam_eye.x, cam_eye.y, cam_eye.z)).length();

            sf::CircleShape circle(ringRadius / distance);
            circle.setPosition(currentPos - sf::Vector2f(ringRadius / distance, ringRadius / distance));
            circle.setFillColor(sf::Color::Transparent);
            circle.setOutlineColor(sf::Color::White);
            circle.setOutlineThickness(thickness);
            circle.setPointCount(pointCount);

            window.getWindow().draw(circle);
        }

        void universe::drawRings(core::window_t &window, sf::RenderTarget& target)
        {
            constexpr float visualRadius = 0.5f;
            constexpr float visualThickness = 0.1f;
            constexpr uint32_t segments = 60;

            for (const auto& ring : rings)
            {
                if (ring.size() < 3) continue;

                sf::Vector3f center(0.0f, 0.0f, 0.0f);
                for (uint32_t idx : ring) center += pos(idx);
                center /= static_cast<float>(ring.size());

                sf::Vector3f normal(0.0f, 0.0f, 0.0f);
                float radius3D = 0.0f;
                sf::Vector3f first = pos(ring[0]);
                for (size_t i = 0; i < ring.size(); ++i)
                {
                    sf::Vector3f p1 = pos(ring[i]);
                    sf::Vector3f p2 = pos(ring[(i + 1) % ring.size()]);
                    radius3D += (p1 - center).length();

                    normal.x += (p1.y - p2.y) * (p1.z + p2.z);
                    normal.y += (p1.z - p2.z) * (p1.x + p2.x);
                    normal.z += (p1.x - p2.x) * (p1.y + p2.y);
                }
                radius3D /= static_cast<float>(ring.size());
                if (normal.lengthSquared() > 0.0f) normal /= normal.length();
                else normal = sf::Vector3f(0.0f, 0.0f, 1.0f); 

                glm::vec3 cam_eye = cam.eye();
                float distance = (center - sf::Vector3f(cam_eye.x, cam_eye.y, cam_eye.z)).length();
                if (distance < 0.1f) continue;

                float scaledRadius = visualRadius * radius3D;
                float scaledHalfThick = visualThickness * 0.5f;

                sf::Vector2f screenCenter = project(window, center);

                if (screenCenter.x < -1000.f) continue;

                sf::Vector3f arb = (std::abs(normal.x) > 0.5f) ? sf::Vector3f(0,1,0) : sf::Vector3f(1,0,0);

                sf::Vector3f u = normal.cross(arb);
                if (u.lengthSquared() > 0.0f) u /= u.length();
                sf::Vector3f v = normal.cross(u);

                sf::VertexArray verts(sf::PrimitiveType::TriangleStrip, 2 * (segments + 1));
                for (uint32_t i = 0; i <= segments; ++i)
                {
                    float angle = 2.0f * static_cast<float>(M_PI) * i / segments;
                    float c = std::cos(angle);
                    float s = std::sin(angle);
                    sf::Vector3f dir = u * c + v * s;

                    // Outer
                    sf::Vector3f outer3D = center + dir * (radius3D + 0.1f * scaledHalfThick);
                    sf::Vector2f outer2D = project(window, outer3D) - screenCenter;
                    verts[2 * i].position = screenCenter + outer2D * (scaledRadius / radius3D);
                    verts[2 * i].color = sf::Color::White;

                    // Inner
                    sf::Vector3f inner3D = center + dir * (radius3D - scaledHalfThick);
                    sf::Vector2f inner2D = project(window, inner3D) - screenCenter;
                    verts[2 * i + 1].position = screenCenter + inner2D * (scaledRadius / radius3D);
                    verts[2 * i + 1].color = sf::Color::White;
                }

                target.draw(verts);
            }
        }
        void universe::drawHydrogenBond(core::window_t &window, sf::RenderTarget& target, int32_t H)
        {
            if (timeStep == 0 || cells.size() == 0)
                return;

            if (data.q[H] > -0.2f && data.q[H] < 0.2f)
                return;

            int32_t D = SIZE_MAX;
            for (const auto &b : bonds)
            {
                if (b.bondedAtom == H)
                {
                    D = b.centralAtom;
                    break;
                }
            }
            if (D == SIZE_MAX)
                return;

            sf::Vector2f pH = project(window, pos(H));

            constexpr float MAX_HA_DISTANCE = 2.8f;
            constexpr float MIN_COS_ANGLE = 0.8f;
            constexpr float STRONG_ENERGY = -700.0f;
            constexpr float WEAK_ENERGY = -5.0f;

            sf::Vector3f p = pos(H);
            p.x -= box.x * std::floor(p.x * (1.0f / box.x));
            p.y -= box.y * std::floor(p.y * (1.0f / box.y));
            p.z -= box.z * std::floor(p.z * (1.0f / box.z));

            int32_t ix = static_cast<int32_t>(p.x * 1 / box.x);
            int32_t iy = static_cast<int32_t>(p.y * 1 / box.y);
            int32_t iz = static_cast<int32_t>(p.z * 1 / box.z);

            ix = ix >= static_cast<int32_t>(cx) ? static_cast<int32_t>(cx) - 1 : ix;
            iy = iy >= static_cast<int32_t>(cy) ? static_cast<int32_t>(cy) - 1 : iy;
            iz = iz >= static_cast<int32_t>(cz) ? static_cast<int32_t>(cz) - 1 : iz;

            for (int32_t dx = -1; dx <= 1; ++dx)
            for (int32_t dy = -1; dy <= 1; ++dy)
            for (int32_t dz = -1; dz <= 1; ++dz)
            {
                int32_t x = ix + dx, y = iy + dy, z = iz + dz;
                uint32_t cell_id = getCellID(x, y, z);

                for (uint32_t A : cells[cell_id])
                {
                    if (A == D || A == H)
                        continue;

                    if (std::abs(data.q[A]) < 0.15f)
                        continue;
                    if (data.q[H] * data.q[A] > 0.0f)
                        continue;

                    sf::Vector3f rHD = minImageVec(pos(D) - pos(H));
                    sf::Vector3f rHA = minImageVec(pos(A) - pos(H));
                    float dHD = rHD.length();
                    float dHA = rHA.length();

                    if (dHA > MAX_HA_DISTANCE || dHA < 0.7f)
                        continue;
                    if (dHD < 0.5f)
                        continue;

                    float cos_angle = rHD.dot(rHA) / (dHD * dHA);
                    if (cos_angle < MIN_COS_ANGLE)
                        continue;

                    float energy_kJ = COULOMB_K * data.q[H] * data.q[A] / dHA;
                    if (energy_kJ > WEAK_ENERGY)
                        continue;

                    float strength = std::clamp((energy_kJ - WEAK_ENERGY) / (STRONG_ENERGY - WEAK_ENERGY), 0.0f, 1.0f);
                    if (strength < 0.15f)
                        continue;

                    uint8_t alpha = static_cast<uint8_t>(255 * strength);

                    sf::Vector2f pA = project(window, pos(A));
                    sf::Vector2f dir = (pA - pH).normalized();
                    float dist2D = (pA - pH).length();
                    if (dist2D < 1e-3f)
                        continue;

                    constexpr int NUM_DASHES = 8;
                    constexpr float DASH_FRACTION = 0.55f;
                    constexpr float GAP_FRACTION = 0.45f;

                    float dash_len = dist2D * DASH_FRACTION / NUM_DASHES;
                    float gap_len = dist2D * GAP_FRACTION / NUM_DASHES;

                    sf::Vertex line[2];
                    line[0].color = line[1].color = sf::Color(255, 255, 255, alpha);

                    for (int i = 0; i < NUM_DASHES; ++i)
                    {
                        float start = (dash_len + gap_len) * i;
                        float end = start + dash_len;

                        if (end > dist2D)
                            end = dist2D;
                        if (start >= dist2D)
                            break;

                        line[0].position = pH + dir * start;
                        line[1].position = pH + dir * end;
                        window.getWindow().draw(line, 2, sf::PrimitiveType::Lines);
                    }
                }
            }
        }

        void universe::drawDebug(core::window_t &window)
        {
            if (cx == 0 || box.x <= 0.f)
                return;

            sf::VertexArray grid(sf::PrimitiveType::Lines, 0);
            grid.resize(24 * cx * cy * cz);

            sf::Color cellColor(100, 149, 237, 60);
            sf::Color hotColor(255, 69, 0, 120);

            int32_t vert_idx = 0;

            for (int32_t ix = 0; ix < cx; ++ix)
                for (int32_t iy = 0; iy < cy; ++iy)
                    for (int32_t iz = 0; iz < cz; ++iz)
                    {
                        int32_t cell_id = getCellID(ix, iy, iz);
                        bool has_atoms = !cells[cell_id].empty();

                        sf::Color color = has_atoms ? hotColor : cellColor;

                        float x0 = ix * CELL_CUTOFF;
                        float y0 = iy * CELL_CUTOFF;
                        float z0 = iz * CELL_CUTOFF;
                        float x1 = x0 + CELL_CUTOFF;
                        float y1 = y0 + CELL_CUTOFF;
                        float z1 = z0 + CELL_CUTOFF;

                        // 8 corners of the cell
                        std::array<sf::Vector3f, 8> c = {{{x0, y0, z0}, {x1, y0, z0}, {x1, y1, z0}, {x0, y1, z0}, {x0, y0, z1}, {x1, y0, z1}, {x1, y1, z1}, {x0, y1, z1}}};

                        // 12 edges
                        const int edges[12][2] = {
                            {0, 1}, {1, 2}, {2, 3}, {3, 0}, // bottom
                            {4, 5},
                            {5, 6},
                            {6, 7},
                            {7, 4}, // top
                            {0, 4},
                            {1, 5},
                            {2, 6},
                            {3, 7} // vertical
                        };

                        for (int e = 0; e < 12; ++e)
                        {
                            sf::Vector2f a = project(window, c[edges[e][0]]);
                            sf::Vector2f b = project(window, c[edges[e][1]]);

                            if (a.x < -1000.f || b.x < -1000.f)
                                continue;

                            glm::vec3 cam_eye = cam.eye();
                            sf::Vector3f camToA = c[edges[e][0]] - sf::Vector3f(cam_eye.x, cam_eye.y, cam_eye.z);
                            sf::Vector3f camToB = c[edges[e][1]] - sf::Vector3f(cam_eye.x, cam_eye.y, cam_eye.z);
                            if (camToA.normalized().dot(camToB.normalized()) < 0.2f)
                                continue;

                            grid[vert_idx++] = sf::Vertex(a, color);
                            grid[vert_idx++] = sf::Vertex(b, color);
                        }
                    }

            window.getWindow().draw(grid);
        }

        int32_t universe::createAtom(sf::Vector3f p, sf::Vector3f v, uint8_t ZIndex, uint8_t numNeutrons, uint8_t numElectron, int32_t chirality)
        {
            atom newAtom{};
            newAtom.ZIndex = ZIndex;

            emplace_vel(v);
            emplace_pos(p);

            std::pair<float, float> constants = constants::getAtomConstants(ZIndex);

            newAtom.sigma = constants.first;
            newAtom.epsilon = constants.second;
            newAtom.radius = constants::VDW_RADII[ZIndex] * 1.5f;
            newAtom.electrons = numElectron;
            newAtom.NCount = numNeutrons;
            newAtom.mass = ZIndex * MASS_PROTON + numNeutrons * MASS_NEUTRON + numElectron * MASS_ELECTRON;
            newAtom.chirality = chirality;
            newAtom.bondCount = 0;

            if (HMassRepartitioning)
            {
                if (ZIndex == 1) newAtom.mass *= 3.f;
                if (ZIndex > 5) newAtom.mass -= 3.f;
            }

            atoms.emplace_back(std::move(newAtom));
            data.q.emplace_back(ZIndex - numElectron);
            data.fx.resize(atoms.size());
            data.fy.resize(atoms.size());
            data.fz.resize(atoms.size());
            data.bond_orders.resize(atoms.size());
            data.temperature.resize(atoms.size());

            frozen_atoms.emplace_back(false);

            return atoms.size() - 1;
        }

        void universe::createBond(int32_t idx1, int32_t idx2, BondType type)
        {
            if (idx1 >= atoms.size() || idx2 >= atoms.size() || idx1 == idx2)
                return;

            if (std::find_if(bonds.begin(), bonds.end(), [&](bond &a)
                             { return a.bondedAtom == idx1 && a.centralAtom == idx2 || a.bondedAtom == idx2 && a.centralAtom == idx1; }) != bonds.end())
                return;

            bond nBond{};
            nBond.bondedAtom = idx1;
            nBond.centralAtom = idx2;
            nBond.type = type;
            nBond.equilibriumLength = constants::getBondLength(atoms[idx1].ZIndex, atoms[idx2].ZIndex, type);
            nBond.k = constants::getBondHarmonicConstantFromEnergy(atoms[idx1].ZIndex, atoms[idx2].ZIndex, type);

            int8_t bondCount = static_cast<int8_t>(type);

            atoms[idx1].bondCount += bondCount;
            atoms[idx2].bondCount += bondCount;

            float EN1 = constants::getElectronegativity(atoms[idx1].ZIndex);
            float EN2 = constants::getElectronegativity(atoms[idx2].ZIndex);
            float deltaEN = std::abs(EN1 - EN2);

            if (deltaEN > 0.1f) // Significant electronegativity difference
            {
                float charge = deltaEN;

                if (EN1 > EN2)
                {
                    data.q[idx1] -= charge;
                    data.q[idx2] += charge;
                }
                else
                {
                    data.q[idx2] -= charge;
                    data.q[idx1] += charge;
                }
            }

            bonds.emplace_back(std::move(nBond));
        }

        int32_t universe::createSubset(const def_subset &nSub, const int32_t baseAtom, const int32_t baseSubset)
        {
            subset nSubset{};
            nSubset.mainAtomIdx = nSub.mainAtomIdx + baseAtom;
            nSubset.bondedSubsetIdx = nSub.bondedSubset == SIZE_MAX ? SIZE_MAX : nSub.bondedSubset + baseSubset;
            nSubset.bondingSubsetIdx = nSub.bondingSubset == SIZE_MAX ? SIZE_MAX : nSub.bondingSubset + baseSubset;

            std::vector<uint8_t> neighbourZs;
            neighbourZs.reserve(nSub.connectedIdx.size() + nSub.hydrogensIdx.size());

            for (int32_t i = 0; i < nSub.connectedIdx.size(); ++i)
            {
                const int32_t bondedAtom = nSub.connectedIdx[i] + baseAtom;
                neighbourZs.emplace_back(atoms[bondedAtom].ZIndex);
            }

            for (int32_t h = 0; h < nSub.hydrogensIdx.size(); ++h)
            {
                neighbourZs.emplace_back(1);
            }

            if (nSub.hydrogensIdx.size() > 0)
            {
                nSubset.hydrogenBegin = nSub.hydrogensIdx[0] + baseAtom;
                nSubset.hydrogenCount = nSub.hydrogensIdx.size();
            }
            if (nSub.connectedIdx.size() > 0)
            {
                nSubset.connectedBegin = nSub.connectedIdx[0] + baseAtom;
                nSubset.connectedCount = nSub.connectedIdx.size();
            }

            subsets.emplace_back(std::move(nSubset));
            return subsets.size() - 1; // Index
        }

        void universe::createMolecule(molecule_structure structure, sf::Vector3f pos, sf::Vector3f vel)
        {
            int32_t baseAtomIndex = atoms.size();

            molecule nMolecule{};

            sf::Vector3f centroid{0.f, 0.f, 0.f};
            for (const auto& p : structure.positions)
            {
                centroid += p;
            }
            centroid /= static_cast<float>(structure.positions.size());

            for (auto& p : structure.positions)
            {
                p -= centroid;
            }

            for (int32_t i = 0; i < structure.atoms.size(); ++i)
            {
                const def_atom &a = structure.atoms[i];

                createAtom(structure.positions[i] + pos, vel, a.ZIndex, a.NIndex, a.ZIndex - a.charge, a.chirality);
                data.q[i] += structure.atoms[i].charge;
            }

            nMolecule.atomBegin = baseAtomIndex;
            nMolecule.atomCount = structure.atoms.size();

            if (!react)
            {
                int32_t baseBondIndex = bonds.size();
                for (int32_t b = 0; b < structure.bonds.size(); ++b)
                {
                    const def_bond &db = structure.bonds[b];
                    int32_t central = baseAtomIndex + db.centralAtomIdx;
                    int32_t bonded = baseAtomIndex + db.bondingAtomIdx;

                    createBond(bonded, central, db.type);
                }

                nMolecule.bondBegin = baseBondIndex;
                nMolecule.bondCount = structure.bonds.size();

                int32_t baseSubset = subsets.size();
                for (int32_t s = 0; s < structure.subsets.size(); ++s)
                    createSubset(structure.subsets[s], baseAtomIndex, baseSubset);

                nMolecule.subsetBegin = baseSubset;
                nMolecule.subsetCount = structure.subsets.size();

                if (nMolecule.subsetCount > 0 && structure.atoms[0].charge == 0)
                    balanceMolecularCharges(subsets[baseSubset]);

                rebuildBondTopology();

                int32_t baseAngle = angles.size();
                for (int32_t a = 0; a < structure.angles.size(); ++a)
                {
                    angle angle = structure.angles[a];
                    angle.A += baseAtomIndex;
                    angle.B += baseAtomIndex;
                    angle.C += baseAtomIndex;

                    angles.emplace_back(angle);
                }

                rings.insert(rings.end(), structure.rings_aromatic.begin(), structure.rings_aromatic.end());

                nMolecule.angleBegin = baseAngle;
                nMolecule.angleCount = structure.angles.size();

                int32_t baseDihedral = dihedral_angles.size();
                for (int32_t a = 0; a < structure.dihedral_angles.size(); ++a)
                {
                    dihedral_angle angle = structure.dihedral_angles[a];
                    angle.A += baseAtomIndex;
                    angle.B += baseAtomIndex;
                    angle.C += baseAtomIndex;
                    angle.D += baseAtomIndex;

                    dihedral_angles.emplace_back(angle);
                }

                nMolecule.dihedralBegin = baseDihedral;
                nMolecule.dihedralCount = structure.dihedral_angles.size();

                if (structure.atoms.size() == 3 && structure.atoms[0].ZIndex == 8 && structure.atoms[1].ZIndex == 1 && structure.atoms[2].ZIndex == 1)
                    nMolecule.exclude = true;
            }

            molecules.emplace_back(std::move(nMolecule));
        }

        void universe::balanceMolecularCharges(subset &mol)
        {
            auto it = std::find_if(subsets.begin(), subsets.end(), [&](const subset &s)
                                   { return &s == &mol; });
            if (it == subsets.end())
                return;

            int32_t startIdx = std::distance(subsets.begin(), it);

            std::vector<int32_t> atomsToBalance;
            std::unordered_set<int32_t> visited;
            std::queue<int32_t> q;
            q.push(startIdx);
            visited.insert(startIdx);

            while (!q.empty())
            {
                int32_t currentIdx = q.front();
                q.pop();
                const auto &currentSubset = subsets[currentIdx];

                atomsToBalance.push_back(currentSubset.mainAtomIdx);

                if (currentSubset.connectedCount != UINT32_MAX)
                    for (int32_t i = 0; i < currentSubset.connectedCount; ++i)
                    {
                        atomsToBalance.emplace_back(currentSubset.connectedBegin + i);
                    }

                if (currentSubset.hydrogenCount != UINT32_MAX)
                    for (int32_t i = 0; i < currentSubset.hydrogenCount; ++i)
                    {
                        atomsToBalance.emplace_back(currentSubset.hydrogenCount + i);
                    }

                int32_t nextIdx = currentSubset.bondedSubsetIdx;
                if (nextIdx < subsets.size() && nextIdx != SIZE_MAX && visited.insert(nextIdx).second)
                {
                    q.push(nextIdx);
                }
            }

            if (atomsToBalance.empty())
                return;

            std::sort(atomsToBalance.begin(), atomsToBalance.end());
            auto last = std::unique(atomsToBalance.begin(), atomsToBalance.end());
            atomsToBalance.erase(last, atomsToBalance.end());

            std::vector<std::vector<int32_t>> neighborLists(atoms.size());
            for (const auto &bond : bonds)
            {
                int32_t a = bond.bondedAtom;
                int32_t b = bond.centralAtom;
                bool aIn = std::binary_search(atomsToBalance.begin(), atomsToBalance.end(), a);
                bool bIn = std::binary_search(atomsToBalance.begin(), atomsToBalance.end(), b);
                if (aIn && bIn)
                {
                    neighborLists[a].push_back(b);
                    neighborLists[b].push_back(a);
                }
            }

            float totalCharge = std::accumulate(data.q.begin(), data.q.end(), 0.f);

            if (std::abs(totalCharge) < 1e-6f)
                return;

            std::vector<float> valenceDeficit(atomsToBalance.size(), 0.0f);
            for (int32_t i = 0; i < atomsToBalance.size(); ++i)
            {
                int32_t idx = atomsToBalance[i];
                if (idx >= atoms.size()) continue;
                
                uint8_t usualValence = constants::getUsualBonds(atoms[idx].ZIndex);
                int deficit = usualValence - static_cast<int>(atoms[idx].bondCount);
                valenceDeficit[i] = static_cast<float>(deficit);
            }

            float chargeToDistribute = -totalCharge;
            float totalDeficit = 0.0f;
            for (float d : valenceDeficit)
            {
                if (d > 0)
                    totalDeficit += d;
            }

            if (totalDeficit > 0.5f)
            {
                for (int32_t i = 0; i < atomsToBalance.size(); ++i)
                {
                    if (valenceDeficit[i] > 0)
                    {
                        float share = chargeToDistribute * (valenceDeficit[i] / totalDeficit);
                        data.q[atomsToBalance[i]] += share;
                    }
                }
            }
            else
            {
                float adjustment = chargeToDistribute / atomsToBalance.size();
                for (int32_t idx : atomsToBalance)
                {
                    data.q[idx] += adjustment;
                }
            }
        }

        void universe::boundCheck(uint32_t i)
        {
            float &x = data.x[i];
            float &y = data.y[i];
            float &z = data.z[i];

            float &vx = data.vx[i];
            float &vy = data.vy[i];
            float &vz = data.vz[i];

            if (wall_collision)
            {
                if (x < 0.0f)
                {
                    x = 0.0f;
                    vx = -vx;
                }
                else if (x > box.x)
                {
                    x = box.x;
                    vx = -vx;
                }

                if (y < 0.0f)
                {
                    y = 0.0f;
                    vy = -vy;
                }
                else if (y > box.y)
                {
                    y = box.y;
                    vy = -vy;
                }

                if (z < 0.0f)
                {
                    z = 0.0f;
                    vz = -vz;
                }
                else if (z > box.z)
                {
                    z = box.z;
                    vz = -vz;
                }

                vx *= 0.99f;
                vy *= 0.99f;
                vz *= 0.99f;
            }
            else
            {
                x = std::fmod(std::fmod(x, box.x) + box.x, box.x);
                y = std::fmod(std::fmod(y, box.y) + box.y, box.y);
                z = std::fmod(std::fmod(z, box.z) + box.z, box.z);
            }
        }

        float universe::ljPot(uint32_t i, uint32_t j)
        {
            float potential = 0.f;

            float dr = minImageVec(pos(i) - pos(j)).length();

            auto [sigma_i, epsilon_i] = constants::getAtomConstants(atoms[i].ZIndex);
            auto [sigma_j, epsilon_j] = constants::getAtomConstants(atoms[j].ZIndex);
            float sigma = (sigma_i + sigma_j) / 2.0f;

            if (dr < sigma * CUTOFF && dr > EPSILON)
            {
                float epsilon = sqrtf(epsilon_i * epsilon_j);
                float r6 = powf((sigma / dr), 6);
                float r12 = r6 * r6;
                potential += 4.0f * epsilon * (r12 - r6);
            }

            return potential;
        }

        sf::Vector3f universe::ljForce(uint32_t i, uint32_t j)
        {
            const atom &a1 = atoms[i];
            const atom &a2 = atoms[j];

            const float sigma_i = a1.sigma;
            const float sigma_j = a2.sigma;
            const float epsilon_i = a1.epsilon;
            const float epsilon_j = a2.epsilon;

            sf::Vector3f dr_vec = minImageVec(pos(i) - pos(j));
            float dr = dr_vec.length();

            const float sigma = sqrtf(sigma_i * sigma_j);

            if (dr < sigma * CUTOFF && dr > EPSILON)
            {
                float sigma6 = powf(sigma, 6);
                float sigma12 = sigma6 * sigma6;

                float epsilon = sqrtf(epsilon_i * epsilon_j);
                float r8 = sigma6 / powf(dr, 8);
                float r14 = 2.f * sigma12 / powf(dr, 14);
                float du_dr = 24.0f * epsilon * (r14 - r8);
                return (du_dr / dr) * dr_vec;
            }

            return sf::Vector3f{0.f, 0.f, 0.f};
        }

        // Wolf method
        float universe::wolfForce(float r, float qi_qj)
        {
            constexpr float alpha = 0.25f; // damping
            constexpr float rc = COULOMB_CUTOFF;
            if (r >= rc)
                return 0.f;

            float ir = 1.f / r;
            float ir_c = 1.f / rc;

            float erfc_ar = std::erfc(alpha * r);
            float exp_term = std::exp(-(alpha * alpha * r * r));

            float force = COULOMB_K * qi_qj * ir * ir * (erfc_ar + 2.f * alpha / std::sqrt(M_PI) * r * exp_term);

            return force;
        }

        sf::Vector3f universe::coulombForce(uint32_t i, uint32_t j, sf::Vector3f &dr_vec)
        {
            float dr = dr_vec.length();

            if (dr < EPSILON || dr > COULOMB_CUTOFF)
                return {0.f, 0.f, 0.f};

            float qq = data.q[i] * data.q[j];
            if (qq == 0.f)
                return {0.f, 0.f, 0.f};

            float forceMag = COULOMB_K * qq / (dr * dr * dr);
            return -forceMag * dr_vec;
        }

        void universe::calcBondForces()
        {
            for (int32_t i = 0; i < bonds.size(); ++i)
            {
                bond &bond = bonds[i];

                int32_t idx1 = bond.bondedAtom;
                int32_t idx2 = bond.centralAtom;
                sf::Vector3f r_vec = minImageVec(pos(idx2) - pos(idx1));
                float dr = r_vec.length();
                if (dr <= EPSILON)
                    continue;

                float delta_r = dr - bond.equilibriumLength;
                float force_magnitude = bond.k * delta_r;

                sf::Vector3f force_dir = r_vec / dr;
                sf::Vector3f force = force_magnitude * force_dir;

                add_force(idx1, force);
                add_force(idx2, -force);
            }
        }

        void universe::calcAngleForces()
        {
            // int32_t count = 0;

            for (const angle &ang : angles)
            {
                int32_t i = ang.A; // left atom
                int32_t j = ang.B; // central atom
                int32_t k = ang.C; // right atom

                sf::Vector3f r_ji = minImageVec(pos(i) - pos(j));
                sf::Vector3f r_jk = minImageVec(pos(k) - pos(j));

                float len_ji = r_ji.length();
                float len_jk = r_jk.length();

                if (len_ji < EPSILON || len_jk < EPSILON)
                    continue;

                sf::Vector3f u_ji = r_ji / len_ji;
                sf::Vector3f u_jk = r_jk / len_jk;

                float cos_theta = std::clamp(u_ji.dot(u_jk), -1.0f, 1.0f);
                float sin_theta = 1.0f - cos_theta * cos_theta;
                sin_theta = std::max(sin_theta, 0.0f);
                sin_theta = std::sqrt(sin_theta);
                if (sin_theta < 1e-6f)
                    sin_theta = 1e-6f;

                float theta = std::acos(cos_theta);
                float delta_theta = theta - ang.rad;

                sf::Vector3f dtheta_dri = (cos_theta * u_ji - u_jk) / (len_ji * sin_theta);
                sf::Vector3f dtheta_drk = (cos_theta * u_jk - u_ji) / (len_jk * sin_theta);

                sf::Vector3f F_i = -ang.K * delta_theta * dtheta_dri;
                sf::Vector3f F_k = -ang.K * delta_theta * dtheta_drk;
                sf::Vector3f F_j = -F_i - F_k; // Newton’s 3rd law

                add_force(i, F_i);
                add_force(j, F_j);
                add_force(k, F_k);

                //++count;
            }

            // std::cout << "angle forces applied: " << count << std::endl;
        }

        void universe::calcLjForces()
        {
            // int32_t count = 0;

            for (int32_t ix = 0; ix < cx; ++ix)
                for (int32_t iy = 0; iy < cy; ++iy)
                    for (int32_t iz = 0; iz < cz; ++iz)
                    {
                        int32_t cell_id = ix + cx * (iy + cy * iz);

                        const std::vector<uint32_t> &cell = cells[cell_id];

                        for (int32_t dz = -1; dz <= 1; ++dz)
                            for (int32_t dy = -1; dy <= 1; ++dy)
                                for (int32_t dx = -1; dx <= 1; ++dx)
                                {
                                    int32_t n_ix = (ix + dx + cx) % cx;
                                    int32_t n_iy = (iy + dy + cy) % cy;
                                    int32_t n_iz = (iz + dz + cz) % cz;
                                    int32_t neighbour_id = getCellID(n_ix, n_iy, n_iz);

                                    const std::vector<uint32_t> &neighbour_cell = cells[neighbour_id];

                                    for (int32_t ii = 0; ii < cell.size(); ++ii)
                                    {
                                        const int32_t &i = cell[ii];
                                        for (int32_t jj = 0; jj < neighbour_cell.size(); ++jj)
                                        {
                                            const int32_t &j = neighbour_cell[jj];

                                            if (j <= i)
                                                continue;

                                            sf::Vector3f f = ljForce(i, j);
                                            add_force(i, f);
                                            add_force(j, -f);
                                        }
                                    }
                                }
                    }

            // std::cout << "count: " << count << std::endl;
        }

        float universe::calculateDihedral(const sf::Vector3f &pa, const sf::Vector3f &pb, const sf::Vector3f &pc, const sf::Vector3f &pd)
        {
            sf::Vector3f v1 = pb - pa;
            sf::Vector3f v2 = pc - pb;
            sf::Vector3f v3 = pd - pc;

            sf::Vector3f n1 = v1.cross(v2);
            sf::Vector3f n2 = v2.cross(v3);

            n1 = n1.length() == 0.f ? sf::Vector3f{0.f, 0.f, 0.f} : n1.normalized();
            n2 = n2.length() == 0.f ? sf::Vector3f{0.f, 0.f, 0.f} : n2.normalized();

            float cosPhi = std::clamp(n1.dot(n2), -1.f, 1.f);
            float phi = std::acos(cosPhi);

            // sign via right-hand rule
            if (v2.dot(n1.cross(n2)) < 0.f)
                phi = -phi;
            if (phi < 0.f)
                phi += 2.f * M_PI;
            return phi;
        }

        void universe::calcDihedralForces()
        {
            for (int32_t d = 0; d < dihedral_angles.size(); ++d)
            {
                const dihedral_angle &d_angle = dihedral_angles[d];

                const sf::Vector3f &pa = pos(d_angle.A);
                const sf::Vector3f &pb = pos(d_angle.B);
                const sf::Vector3f &pc = pos(d_angle.C);
                const sf::Vector3f &pd = pos(d_angle.D);

                float phi = calculateDihedral(pa, pb, pc, pd);

                float target = d_angle.rad;
                if (target == 0.f && d_angle.periodicity == 1)
                {
                    int32_t chi = atoms[d_angle.B].chirality ? atoms[d_angle.B].chirality : atoms[d_angle.C].chirality;
                    if (chi == 1)
                        target = (phi < static_cast<float>(M_PI)) ? 1.047f : 5.236f; // ~60° / ~300°
                    else if (chi == 2)
                        target = (phi < static_cast<float>(M_PI)) ? 5.236f : 1.047f;
                    else
                        continue;
                }

                float diff = phi - target;
                while (diff > M_PI)
                    diff -= 2.f * M_PI;
                while (diff < -M_PI)
                    diff += 2.f * M_PI;

                sf::Vector3f b1 = pb - pa;
                sf::Vector3f b2 = pc - pb;
                sf::Vector3f b3 = pd - pc;

                sf::Vector3f n1 = b1.cross(b2).normalized();
                sf::Vector3f n2 = b2.cross(b3).normalized();
                sf::Vector3f u2 = b2.normalized();

                float sin_term = std::sin(d_angle.periodicity * phi - d_angle.rad);
                float torque_mag = d_angle.K * d_angle.periodicity * sin_term;

                sf::Vector3f axis = (pc - pb).normalized();

                sf::Vector3f rA = (pa - (pb.length() > 0.f ? pb.normalized() : sf::Vector3f{0.f, 0.f, 0.f}));
                sf::Vector3f torqueA = axis.cross(rA).normalized() * torque_mag;

                sf::Vector3f rD = (pd - pc).normalized();
                sf::Vector3f torqueD = axis.cross(rD).normalized() * (-torque_mag);

                sf::Vector3f f1 = (torque_mag / b1.length()) * (n1.cross(u2));
                sf::Vector3f f4 = (torque_mag / b3.length()) * (u2.cross(n2));

                add_force(d_angle.A, -f1);
                add_force(d_angle.D, -f4);
                add_force(d_angle.B, f1);
                add_force(d_angle.C, f4);
                add_force(d_angle.D, -torqueD * 0.5f);
            }
        }

        void universe::calcElectrostaticForces()
        {
            for (int32_t ix = 0; ix < cx; ++ix)
                for (int32_t iy = 0; iy < cy; ++iy)
                    for (int32_t iz = 0; iz < cz; ++iz)
                    {
                        int32_t cell_id = getCellID(ix, iy, iz);

                        const std::vector<uint32_t> &cell = cells[cell_id];

                        for (int32_t ii = 0; ii < cell.size(); ++ii)
                        {
                            uint32_t i = cell[ii];
                            if (data.q[i] == 0.f)
                                continue;
                            for (int32_t jj = ii + 1; jj < cell.size(); ++jj)
                            {
                                const uint32_t &j = cell[jj];

                                if (j <= i)
                                    continue;
                                if (data.q[j] == 0.f)
                                    continue;

                                sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                                sf::Vector3f force = coulombForce(i, j, dr);

                                add_force(i, force);
                                add_force(j, -force);
                            }
                        }

                        for (int32_t d = 1; d < 14; ++d)
                        {
                            int32_t n_ix = ix + offsets[d][0], n_iy = iy + offsets[d][1], n_iz = iz + offsets[d][2];
                            int32_t neighbour_id = getCellID(n_ix, n_iy, n_iz);

                            const std::vector<uint32_t> &neighbour_cell = cells[neighbour_id];

                            for (int32_t ii = 0; ii < cell.size(); ++ii)
                            {
                                const uint32_t &i = cell[ii];
                                if (data.q[i] == 0.f)
                                    continue;
                                for (int32_t jj = 0; jj < neighbour_cell.size(); ++jj)
                                {
                                    const uint32_t &j = neighbour_cell[jj];

                                    if (data.q[j] == 0.f)
                                        continue;

                                    sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                                    sf::Vector3f force = coulombForce(i, j, dr);

                                    add_force(i, force);
                                    add_force(j, -force);
                                }
                            }
                        }
                    }
        }

        void universe::calcBondedForces()
        {
            calcBondForces();
            calcAngleForces();
            calcDihedralForces();
        }

        std::vector<sf::Vector3f> universe::processCellUnbonded(int32_t ix, int32_t iy, int32_t iz)
        {
            float local_virial = 0.f;

            int32_t cellID = getCellID(ix, iy, iz);
            const auto &cell = cells[cellID];
            thread_local std::vector<sf::Vector3f> local_forces(atoms.size(), {0, 0, 0});

            if (local_forces.size() != atoms.size())
                local_forces.assign(atoms.size(), {0.0f, 0.0f, 0.0f});

            std::fill(local_forces.begin(), local_forces.end(), sf::Vector3f{0, 0, 0});

            for (int32_t ii = 0; ii < cell.size(); ++ii)
            {
                const int32_t i = cell[ii];

                for (int32_t jj = ii + 1; jj < cell.size(); ++jj)
                {
                    const int32_t j = cell[jj];
                    sf::Vector3f dr = minImageVec(pos(j) - pos(i));

                    float r2 = dr.lengthSquared();
                    if (r2 > CELL_CUTOFF)
                        continue;

                    if (areBonded(i, j))
                        continue;

                    sf::Vector3f cForce = coulombForce(i, j, dr);
                    sf::Vector3f lForce = ljForce(i, j);

                    sf::Vector3f total_force = cForce + lForce;
                    local_forces[i] += total_force;
                    local_forces[j] -= total_force;

                    local_virial += dr.x * total_force.x +
                                    dr.y * total_force.y +
                                    dr.z * total_force.z;
                }
            }

            for (int32_t dx = 0; dx <= 1; ++dx)
                for (int32_t dy = -1; dy <= 1; ++dy)
                    for (int32_t dz = -1; dz <= 1; ++dz)
                    {
                        if (dx == 0 && dy == 0 && dz == 0)
                            continue;

                        int32_t n_ix = ix + dx;
                        int32_t n_iy = iy + dy;
                        int32_t n_iz = iz + dz;
                        int32_t neighbor_id = getCellID(n_ix, n_iy, n_iz);

                        if (neighbor_id == cellID)
                            continue;

                        const auto &neighbor_cell = cells[neighbor_id];

                        for (int32_t ii = 0; ii < cell.size(); ++ii)
                        {
                            int32_t i = cell[ii];

                            for (int32_t jj = 0; jj < neighbor_cell.size(); ++jj)
                            {
                                int32_t j = neighbor_cell[jj];

                                if (j <= i)
                                    continue;

                                sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                                float r = dr.length();

                                if (r > CELL_CUTOFF)
                                    continue;

                                if (areBonded(i, j))
                                    continue;

                                sf::Vector3f cForce = coulombForce(i, j, dr);
                                sf::Vector3f lForce = ljForce(i, j);

                                sf::Vector3f total_force = cForce + lForce;
                                local_forces[i] += total_force;
                                local_forces[j] -= total_force;

                                local_virial += dr.x * total_force.x +
                                                dr.y * total_force.y +
                                                dr.z * total_force.z;
                            }
                        }
                    }

            total_virial.fetch_add(local_virial, std::memory_order_relaxed);
            return local_forces;
        }

        void universe::calcBondedForcesParallel()
        {
            const uint32_t n_threads = std::max(1u, std::thread::hardware_concurrency());
            std::vector<std::future<std::vector<sf::Vector3f>>> futures;

            auto make_task = [this](auto &&func)
            {
                return [this, func](int32_t start, int32_t end) -> std::vector<sf::Vector3f>
                {
                    thread_local std::vector<sf::Vector3f> local_forces;
                    if (local_forces.size() != atoms.size())
                    {
                        local_forces.assign(atoms.size(), {0.0f, 0.0f, 0.0f});
                    }
                    std::fill(local_forces.begin(), local_forces.end(), sf::Vector3f{0, 0, 0});

                    func(start, end, local_forces);

                    return local_forces;
                };
            };

            auto bond_func = [this](int32_t start, int32_t end, std::vector<sf::Vector3f> &lf)
            {
                for (int32_t i = start; i < end; ++i)
                {
                    const bond &b = bonds[i];
                    int32_t a = b.bondedAtom;
                    int32_t c = b.centralAtom;

                    sf::Vector3f dr = minImageVec(pos(c) - pos(a));
                    float len = dr.length();
                    if (len <= EPSILON)
                        continue;

                    float dl = len - b.equilibriumLength;
                    sf::Vector3f force = (b.k * dl / len) * dr;

                    lf[a] += force;
                    lf[c] -= force;
                }
            };

            for (int32_t t = 0; t < n_threads; ++t)
            {
                int32_t start = bonds.size() * t / n_threads;
                int32_t end = bonds.size() * (t + 1) / n_threads;
                futures.emplace_back(std::async(std::launch::async, make_task(bond_func), start, end));
            }

            auto angle_func = [this](int32_t start, int32_t end, std::vector<sf::Vector3f> &lf)
            {
                for (int32_t a = start; a < end; ++a)
                {
                    const angle &ang = angles[a];
                    int32_t i = ang.A, j = ang.B, k = ang.C;

                    sf::Vector3f r_ji = minImageVec(pos(i) - pos(j));
                    sf::Vector3f r_jk = minImageVec(pos(k) - pos(j));
                    float len_ji = r_ji.length();
                    float len_jk = r_jk.length();
                    if (len_ji < EPSILON || len_jk < EPSILON)
                        continue;

                    sf::Vector3f u_ji = r_ji / len_ji;
                    sf::Vector3f u_jk = r_jk / len_jk;
                    float cos_theta = std::clamp(u_ji.dot(u_jk), -1.0f, 1.0f);
                    float sin_theta = std::sqrt(std::max(1.0f - cos_theta * cos_theta, 0.0f));
                    if (sin_theta < 1e-6f)
                        sin_theta = 1e-6f;

                    float theta = std::acos(cos_theta);
                    float delta_theta = theta - ang.rad;

                    sf::Vector3f dtheta_dri = (cos_theta * u_ji - u_jk) / (len_ji * sin_theta);
                    sf::Vector3f dtheta_drk = (cos_theta * u_jk - u_ji) / (len_jk * sin_theta);

                    sf::Vector3f F_i = -ang.K * delta_theta * dtheta_dri;
                    sf::Vector3f F_k = -ang.K * delta_theta * dtheta_drk;
                    sf::Vector3f F_j = -F_i - F_k;

                    lf[i] += F_i;
                    lf[j] += F_j;
                    lf[k] += F_k;
                }
            };

            for (int32_t t = 0; t < n_threads; ++t)
            {
                int32_t start = angles.size() * t / n_threads;
                int32_t end = angles.size() * (t + 1) / n_threads;
                futures.emplace_back(std::async(std::launch::async, make_task(angle_func), start, end));
            }

            auto dihedral_func = [this](int32_t start, int32_t end, std::vector<sf::Vector3f> &lf)
            {
                for (int32_t d = start; d < end; ++d)
                {
                    const dihedral_angle &da = dihedral_angles[d];
                    const sf::Vector3f &pa = pos(da.A);
                    const sf::Vector3f &pb = pos(da.B);
                    const sf::Vector3f &pc = pos(da.C);
                    const sf::Vector3f &pd = pos(da.D);

                    float phi = calculateDihedral(pa, pb, pc, pd);

                    float target = da.rad;
                    if (target == 0.0f && da.periodicity == 1)
                    {
                        int32_t chi = atoms[da.B].chirality ? atoms[da.B].chirality : atoms[da.C].chirality;
                        if (chi == 1)
                            target = (phi < M_PI) ? 1.047f : 5.236f;
                        if (chi == 2)
                            target = (phi < M_PI) ? 5.236f : 1.047f;
                    }

                    float diff = phi - target;
                    while (diff > M_PI)
                        diff -= 2.0f * M_PI;
                    while (diff < -M_PI)
                        diff += 2.0f * M_PI;

                    sf::Vector3f b1 = pb - pa;
                    sf::Vector3f b2 = pc - pb;
                    sf::Vector3f b3 = pd - pc;

                    sf::Vector3f n1 = b1.cross(b2).normalized();
                    sf::Vector3f n2 = b2.cross(b3).normalized();
                    sf::Vector3f u2 = b2.normalized();

                    float sin_term = std::sin(da.periodicity * phi - da.rad);
                    float torque_mag = da.K * da.periodicity * sin_term;

                    sf::Vector3f axis = (pc - pb).normalized();

                    sf::Vector3f rA = (pa - (pb.length() > 0.f ? pb.normalized() : sf::Vector3f{0.f, 0.f, 0.f}));
                    sf::Vector3f torqueA = axis.cross(rA).normalized() * torque_mag;

                    sf::Vector3f rD = (pd - pc).normalized();
                    sf::Vector3f torqueD = axis.cross(rD).normalized() * (-torque_mag);

                    sf::Vector3f f1 = (torque_mag / b1.length()) * (n1.cross(u2));
                    sf::Vector3f f4 = (torque_mag / b3.length()) * (u2.cross(n2));

                    lf[da.A] += -f1;
                    lf[da.D] += -f4;
                    lf[da.B] +=  f1;
                    lf[da.C] +=  f4;
                    lf[da.D] +=  -torqueD * 0.5f;
                }
            };

            for (int32_t t = 0; t < n_threads; ++t)
            {
                int32_t start = dihedral_angles.size() * t / n_threads;
                int32_t end = dihedral_angles.size() * (t + 1) / n_threads;
                futures.emplace_back(std::async(std::launch::async, make_task(dihedral_func), start, end));
            }

            for (auto &fut : futures)
            {
                auto local_f = fut.get();
                for (int32_t i = 0; i < atoms.size(); ++i)
                    add_force(i, local_f[i]);
            }
        }

        void universe::calcUnbondedForcesParallel()
        {
            // NOTE: eventually distribute cells with more work over threads instead of evenly giving to all threads, threads with no work are gonna be dismissed.
            const int32_t n_threads = std::thread::hardware_concurrency();

            std::vector<std::future<std::vector<sf::Vector3f>>> futures;

            std::vector<uint32_t> work(cells.size());
            for (int32_t c = 0; c < cells.size(); ++c)
            {
                uint32_t n = cells[c].size();
                work[c] = n * n;
            }

            auto worker = [this](int32_t start_flat, int32_t end_flat) -> std::vector<sf::Vector3f>
            {
                std::vector<sf::Vector3f> thread_forces(atoms.size(), {0, 0, 0});

                for (int32_t flat = start_flat; flat < end_flat; ++flat)
                {
                    int32_t iz = flat % cz;
                    int32_t iy = (flat / cz) % cy;
                    int32_t ix = flat / (cz * cy);

                    auto cell_forces = processCellUnbonded(ix, iy, iz);
                    for (int32_t i = 0; i < atoms.size(); ++i)
                        thread_forces[i] += cell_forces[i];
                }
                return thread_forces;
            };

            for (int32_t t = 0; t < n_threads; ++t)
            {
                int32_t start = (cells.size() * t) / n_threads;
                int32_t end = (cells.size() * (t + 1)) / n_threads;

                futures.push_back(std::async(std::launch::async, worker, start, end));
            }

            for (auto &fut : futures)
            {
                std::vector<sf::Vector3f> local_f = fut.get();
                for (int32_t i = 0; i < atoms.size(); ++i)
                    add_force(i, local_f[i]);
            }
        }

        void universe::calcUnbondedForces()
        {
            for (int32_t ix = 0; ix < cx; ++ix)
                for (int32_t iy = 0; iy < cy; ++iy)
                    for (int32_t iz = 0; iz < cz; ++iz)
                    {
                        int32_t cell_id = getCellID(ix, iy, iz);
                        const auto &cell = cells[cell_id];

                        for (int32_t ii = 0; ii < cell.size(); ++ii)
                        {
                            int32_t i = cell[ii];

                            for (int32_t jj = ii + 1; jj < cell.size(); ++jj)
                            {
                                int32_t j = cell[jj];

                                sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                                float r2 = dr.lengthSquared();
                                if (r2 > CELL_CUTOFF)
                                    continue;

                                if (areBonded(i, j))
                                    continue;

                                sf::Vector3f cForce = coulombForce(i, j, dr);
                                sf::Vector3f lForce = ljForce(i, j);

                                sf::Vector3f total_force = cForce + lForce;
                                add_force(i, total_force);
                                add_force(j, -total_force);

                                total_virial += dr.x * total_force.x +
                                                dr.y * total_force.y +
                                                dr.z * total_force.z;
                            }
                        }

                        for (int32_t dx = 0; dx <= 1; ++dx)
                            for (int32_t dy = -1; dy <= 1; ++dy)
                                for (int32_t dz = -1; dz <= 1; ++dz)
                                {
                                    if (dx == 0 && dy == 0 && dz == 0)
                                        continue;

                                    int32_t n_ix = ix + dx;
                                    int32_t n_iy = iy + dy;
                                    int32_t n_iz = iz + dz;
                                    int32_t neighbor_id = getCellID(n_ix, n_iy, n_iz);

                                    if (neighbor_id == cell_id)
                                        continue;

                                    const auto &neighbor_cell = cells[neighbor_id];

                                    for (int32_t ii = 0; ii < cell.size(); ++ii)
                                    {
                                        int32_t i = cell[ii];

                                        for (int32_t jj = 0; jj < neighbor_cell.size(); ++jj)
                                        {
                                            int32_t j = neighbor_cell[jj];

                                            if (j <= i)
                                                continue;

                                            sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                                            float r = dr.length();

                                            if (r > CELL_CUTOFF)
                                                continue;
                                            if (areBonded(i, j))
                                                continue;

                                            sf::Vector3f cForce = coulombForce(i, j, dr);
                                            sf::Vector3f lForce = ljForce(i, j);

                                            sf::Vector3f total_force = cForce + lForce;
                                            add_force(i, total_force);
                                            add_force(j, -total_force);

                                            total_virial += dr.x * total_force.x +
                                                            dr.y * total_force.y +
                                                            dr.z * total_force.z;
                                        }
                                    }
                                }
                    }
        }

        float universe::calculateUncorrectedBondOrder(int32_t i, int32_t j)
        {
            sf::Vector3f dr_vec = minImageVec(pos(i) - pos(j));
            float r = dr_vec.length();
            if (r > REACTION_CUTOFF)
                return 0.0f;

            uint8_t Zi = atoms[i].ZIndex;
            uint8_t Zj = atoms[j].ZIndex;

            auto &pi = constants::getParams(Zi);
            auto &pj = constants::getParams(Zj);

            float r_sigma = (pi.r0_sigma + pj.r0_sigma) / 2.0f;
            float r_pi = (pi.r0_pi + pj.r0_pi) / 2.0f;
            float r_pp = (pi.r0_pp + pj.r0_pp) / 2.0f;

            float bo_sigma = std::exp(pi.p_bo1 * std::pow(r / r_sigma, pi.p_bo2));
            float bo_pi = std::exp(pi.p_bo3 * std::pow(r / r_pi, pi.p_bo4));
            float bo_pp = std::exp(pi.p_bo5 * std::pow(r / r_pp, pi.p_bo6));

            return bo_sigma + bo_pi + bo_pp;
        }

        float universe::calculateBondEnergy(int32_t i, int32_t j, float bo_sigma, float bo_pi, float bo_pp)
        {
            auto &parami = constants::getParams(atoms[i].ZIndex);
            auto &paramj = constants::getParams(atoms[j].ZIndex);

            return -(sqrt(parami.De_sigma * paramj.De_sigma) * bo_sigma + sqrt(parami.De_pi * paramj.De_pi) * bo_pi + sqrt(parami.De_pp * paramj.De_pp) * bo_pp);
        }

        void universe::processReactivePair(int32_t i, int32_t j, float cutoff, float vis_thresh)
        {
        }

        void universe::handleReactiveForces()
        {
        }

        void universe::update(float targetTemperature, float targetPressure)
        {
            int32_t N = atoms.size();

            if (N == 0)
                return;

            data.fx.assign(N, 0.0f);
            data.fy.assign(N, 0.0f);
            data.fz.assign(N, 0.0f);

            total_virial = 0.0f;

            if (!react)
            {
                calcBondedForcesParallel();
                calcUnbondedForcesParallel();
            }
            else
                handleReactiveForces();

            setPressure(targetPressure);
            setTemperature(targetTemperature);

            for (int32_t i = 0; i < N; ++i)
            {
                sf::Vector3f accel = force(i) / atoms[i].mass;
                add_vel(i, accel * (0.5f * DT));

                if (gravity)
                    add_vel(i, sf::Vector3f(0.f, 0.f, -mag_gravity * 0.5f * DT));
            }

            for (int32_t i = 0; i < N; ++i)
            {
                sf::Vector3f dPos = vel(i) * DT;
                add_pos(i, dPos);
                boundCheck(i);
            }

            data.fx.assign(N, 0.0f);
            data.fy.assign(N, 0.0f);
            data.fz.assign(N, 0.0f);

            if (!react)
            {
                calcBondedForcesParallel();
                calcUnbondedForcesParallel();
            }
            else
                handleReactiveForces();

            for (int32_t i = 0; i < N; ++i) {
                sf::Vector3f accel = force(i) / atoms[i].mass;
                add_vel(i, accel * (0.5f * DT));

                if (gravity)
                    add_vel(i, sf::Vector3f(0.f, 0.f, -mag_gravity * 0.5f * DT));
            }

            buildCells();

            /* if (timeStep % 2 == 0)
                saveFrame(); */

            ++timeStep;
        }

        float universe::calculatePressure()
        {
            if (atoms.empty())
                return 0.f;

            float kinetic = calculateKineticEnergy();
            float volume = box.x * box.y * box.z;
            float temperature = (2.0f / 3.0f) * kinetic / atoms.size();

            float P_ideal = atoms.size() * temperature / volume;

            float P_virial = -total_virial / (3.0f * volume);

            float pressure = P_ideal + P_virial;
            return pressure;
        }

        void universe::setPressure(float Target_P_Bar)
        {
            if (Target_P_Bar <= 0.0f)
                return;
            if (timeStep % BAROSTAT_INTERVAL != 0)
                return;

            constexpr float beta_T = 4.5e-5f;
            constexpr float tau_P = 1.0f;

            pres = calculatePressure();

            float delta_P = Target_P_Bar - pres;
            float mu = 1.0f - (2.f / tau_P) * beta_T * delta_P;

            mu = std::clamp(mu, 0.5f, 1.5f);

            float scale = std::cbrt(mu);

            box.z *= scale;

            for (int32_t i = 0; i < atoms.size(); ++i)
            {
                data.z[i] *= scale;
            }
        }

        void universe::setTemperature(float kelvin)
        {
            if (timeStep % THERMOSTAT_INTERVAL != 0)
                return;

            if (atoms.size() == 0)
                return;

            float avg_KE = calculateKineticEnergy() / atoms.size();
            temp = (2.0f / 3.0f) * (avg_KE * KB);
            float lambda = sqrtf(kelvin / (temp + 1e-5f));

            if (timeStep > 100)
            {
                constexpr float tau = 0.6f;
                lambda = 1.0f + (lambda - 1.0f) / tau;
            }

            for (int32_t i = 0; i < data.vx.size(); ++i)
            {
                data.vx[i] *= lambda;
                data.vy[i] *= lambda;
                data.vz[i] *= lambda;
            }
        }

        // Energy Calculation
        float universe::calculateKineticEnergy()
        {
            float kinetic_energy = 0.0f;
            for (int32_t i = 0; i < atoms.size(); ++i)
                kinetic_energy += 0.5f * atoms[i].mass * vel(i).lengthSquared();

            return kinetic_energy;
        }

        void universe::buildCells()
        {
            const float cutoff = react ? REACTION_CUTOFF : CELL_CUTOFF;
            if ((cutoff - box.x) > EPSILON || timeStep == 0)
            {
                cx = static_cast<int32_t>(std::ceil(cutoff / box.x));
                cy = static_cast<int32_t>(std::ceil(cutoff / box.y));
                cz = static_cast<int32_t>(std::ceil(cutoff / box.z));
            }

            int32_t ncells = cx * cy * cz;
            if (cells.size() != ncells)
                cells.assign(ncells, {});
            else
                for (auto &c : cells)
                    c.clear();

            const float inv = 1.0f / box.x;

            for (int32_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector3f p = pos(i);
                p.x -= cutoff * std::floor(p.x / cutoff);
                p.y -= cutoff * std::floor(p.y / cutoff);
                p.z -= cutoff * std::floor(p.z / cutoff);

                int32_t ix = static_cast<int32_t>(p.x * inv);
                int32_t iy = static_cast<int32_t>(p.y * inv);
                int32_t iz = static_cast<int32_t>(p.z * inv);

                ix = std::clamp(ix, 0, (int32_t)cx - 1);
                iy = std::clamp(iy, 0, (int32_t)cy - 1);
                iz = std::clamp(iz, 0, (int32_t)cz - 1);

                int32_t id = static_cast<int32_t>(ix) + cx * (static_cast<int32_t>(iy) + cy * static_cast<int32_t>(iz));
                cells[id].emplace_back(i);
            }
        }

        float universe::calculateAtomTemperature(int32_t i)
        {
            float ke = 0.5f * atoms[i].mass * vel(i).lengthSquared();
            return (2.0f / 3.0f) * ke * KB;
        }

        // Camera
        sf::Vector2f universe::project(core::window_t &window, const sf::Vector3f &p) const
        {
            auto size = window.getWindow().getView().getSize();
            return cam.project(glm::vec3{p.x, p.y, p.z}, size.x, size.y);
        }

        void universe::handleCamera()
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

        // TO DO
        std::string universe::moleculeName(const std::vector<uint32_t> &subsetIdx)
        {
            std::map<uint8_t, int32_t> ZIndices;

            for (int32_t i = 0; i < subsetIdx.size(); ++i)
            {
                auto &s = subsets[subsetIdx[i]];
                int32_t centralAtom = s.mainAtomIdx;
                ++ZIndices[atoms[centralAtom].ZIndex];
                ZIndices[1] += s.connectedCount;
            }

            std::unordered_map<int32_t, std::string> prefixes{
                {1, "mono"},
                {2, "di"},
                {3, "tri"},
                {4, "tetra"},
                {5, "penta"},
                {6, "hexa"},
                {7, "hepta"},
                {8, "octa"},
                {9, "nona"},
                {10, "deca"}};

            std::string name("");

            bool organic = ZIndices[6] > 0;

            if (!organic)
            {
                for (int32_t i = 1; i <= COUNT_ATOMS; ++i)
                {
                    if (ZIndices[i] <= 0)
                        continue;
                    std::string name = constants::getAtomName(i);
                    std::transform(name.begin(), name.end(), name.begin(), [](uint8_t c)
                                   { return std::tolower(c); });

                    name += prefixes[ZIndices[i]] + name + " ";
                }
            }
            else // organic nomeclature
            {
            }

            return name;
        }

        void universe::initReaxParams()
        {

        }

        // loading and saving scenes

        void universe::saveFrame()
        {
            positionxLog.emplace_back(data.x);
            positionyLog.emplace_back(data.y);
            positionzLog.emplace_back(data.z);

            temperatureLog.emplace_back(temp);
        }

        // Helper
        std::string formatTime()
        {
            const auto now = std::chrono::system_clock::now();

            std::string time_str = std::format("{:%Y_%m_%d%H_%M}", now);
            return time_str;
        }

        void universe::saveAsVideo(const std::filesystem::path path, const std::string name)
        {
            nlohmann::json video{};

            for (int32_t i = 0; i < positionxLog.size(); ++i)
            {
                video["posx"].emplace_back(positionxLog[i]);
                video["posy"].emplace_back(positionyLog[i]);
                video["posz"].emplace_back(positionzLog[i]);

                video["temperature"].emplace_back(temperatureLog[i]);
            }

            video["metadata"] =
                {
                    {"title", formatTime()},
                    {"description", "Molecular dynamics trajectory"},
                    {"atoms", atoms.size()},
                    {"frames", positionxLog.size()},
                    {"box", {box.x, box.y, box.z}}};

            try
            {
                if (!std::filesystem::is_directory(path))
                    std::filesystem::create_directory(path);
                
                std::filesystem::path filepath = name.empty() ? path / ("video_" + formatTime() + ".json") : path / ("video_" + name + ".json");
                std::ofstream file(filepath);
                file.flush();
                file << video.dump();
                file.close();

                std::cout << "[Recorder] Saved trajectory with " << video["metadata"]["frames"] << " frames to " << path << '\n';
            }
            catch (std::exception &e)
            {
                std::cout << "[Recorder] Failed to save trajectory: " << e.what() << std::endl;
                return;
            }

            positionxLog.clear();
            positionyLog.clear();
            positionzLog.clear();

            temperatureLog.clear();
            energyLog.clear();
        }

        void universe::saveScene(const std::filesystem::path path, const std::string name)
        {
            if (!std::filesystem::is_directory(path))
                std::filesystem::create_directory(path);

            nlohmann::json scene{};

            for (int32_t x = 0; x < data.x.size(); ++x)
            {
                scene["posx"].emplace_back(data.x[x]);
                scene["posy"].emplace_back(data.y[x]);
                scene["posz"].emplace_back(data.z[x]);

                scene["velx"].emplace_back(data.vx[x]);
                scene["vely"].emplace_back(data.vy[x]);
                scene["velz"].emplace_back(data.vz[x]);

                scene["charge"].emplace_back(data.q[x]);
                scene["temperature"].emplace_back(data.temperature[x]);
                scene["bondOrder"].emplace_back(data.bond_orders[x]);

                scene["ZIndex"].emplace_back(atoms[x].ZIndex);
                scene["neutrons"].emplace_back(atoms[x].NCount);
                scene["electrons"].emplace_back(atoms[x].electrons);
                scene["boundCount"].emplace_back(atoms[x].bondCount);
                scene["chirality"].emplace_back(atoms[x].chirality);
            }

            for (int32_t b = 0; b < bonds.size(); ++b)
            {
                auto &bond = bonds[b];
                nlohmann::json b_json{};
                b_json["central"] = bond.centralAtom;
                b_json["bonded"] = bond.bondedAtom;
                b_json["equilibrium"] = bond.equilibriumLength;
                b_json["type"] = static_cast<uint8_t>(bond.type);
                b_json["k"] = bond.k;

                scene["bonds"].emplace_back(b_json);
            }

            for (int32_t a = 0; a < angles.size(); ++a)
            {
                auto &angle = angles[a];
                nlohmann::json a_json{};
                a_json["A"] = angle.A;
                a_json["B"] = angle.B;
                a_json["C"] = angle.C;
                a_json["K"] = angle.K;
                a_json["rad"] = angle.rad;

                scene["angles"].emplace_back(a_json);
            }

            for (int32_t d = 0; d < dihedral_angles.size(); ++d)
            {
                auto &dihedral = dihedral_angles[d];
                nlohmann::json d_json{};

                d_json["A"] = dihedral.A;
                d_json["B"] = dihedral.B;
                d_json["C"] = dihedral.C;
                d_json["D"] = dihedral.D;
                d_json["K"] = dihedral.K;
                d_json["periodicity"] = dihedral.periodicity;
                d_json["rad"] = dihedral.rad;

                scene["dihedrals"].emplace_back(d_json);
            }

            for (int32_t s = 0; s < subsets.size(); ++s)
            {
                auto &subset = subsets[s];
                nlohmann::json s_json{};

                s_json["mainAtom"] = subset.mainAtomIdx;
                s_json["bonded"] = subset.bondedSubsetIdx;
                s_json["bonding"] = subset.bondingSubsetIdx;
                s_json["connected_begin"] = subset.connectedBegin;
                s_json["connected_count"] = subset.connectedCount;
                s_json["hydrogen_begin"] = subset.hydrogenBegin;
                s_json["hydrogen_count"] = subset.hydrogenCount;

                scene["subsets"].emplace_back(s_json);
            }

            for (int32_t m = 0; m < molecules.size(); ++m)
            {
                auto &molecule = molecules[m];
                nlohmann::json m_json{};

                m_json["atomBegin"] = molecule.atomBegin;
                m_json["atomCount"] = molecule.atomCount;
                m_json["bondBegin"] = molecule.bondBegin;
                m_json["bondCount"] = molecule.bondCount;
                m_json["angleBegin"] = molecule.angleBegin;
                m_json["angleCount"] = molecule.angleCount;
                m_json["dihedralBegin"] = molecule.dihedralBegin;
                m_json["dihedralCount"] = molecule.dihedralCount;
                m_json["subsetBegin"] = molecule.subsetBegin;
                m_json["subsetCount"] = molecule.subsetCount;
                m_json["exclude"] = molecule.exclude;

                scene["molecules"].emplace_back(m_json);
            }

            scene["gravity"] = gravity;
            scene["react"] = react;
            scene["isothermal"] = isothermal;
            scene["render_water"] = render_water;
            scene["wall_collision"] = wall_collision;
            scene["mag_gravity"] = mag_gravity;
            scene["boxx"] = box.x;
            scene["boxy"] = box.y;
            scene["boxz"] = box.z;

            try
            {
                std::filesystem::path filepath = name.empty() ? path / (formatTime() + ".json") : path / (name + ".json");
                std::ofstream file(filepath);
                if (!file.is_open())
                {
                    std::cerr << "[Save] Failed to create file: " << filepath << "\n";
                    return;
                }

                file << scene.dump();
                file.flush();
                file.close();

                std::cout << "[Save] Successfully saved: " << filepath << "\n";
            }
            catch (const std::exception &e)
            {
                std::cerr << "[Save] Exception while saving: " << e.what() << "\n";
            }
        }

        void universe::loadScene(const std::filesystem::path path)
        {
            std::ifstream file(path);
            if (!file.is_open() || path.extension() != ".json")
            {
                std::cerr << "[Simulation] Cannot open: " << path << '\n';
                return;
            }

            nlohmann::json scene;
            try
            {
                file >> scene;
            }
            catch (const nlohmann::json::parse_error &e)
            {
                std::cerr << "[Simulation] JSON parse error: " << e.what() << '\n';
                return;
            }

            atoms.clear();
            bonds.clear();
            angles.clear();
            dihedral_angles.clear();
            subsets.clear();
            molecules.clear();
            reactive_bonds.clear();
            data.x.clear();
            data.y.clear();
            data.z.clear();
            data.vx.clear();
            data.vy.clear();
            data.vz.clear();
            data.q.clear();
            data.temperature.clear();
            data.bond_orders.clear();

            box.x = scene.value("boxx", 50.0f);
            box.y = scene.value("boxy", 50.0f);
            box.z = scene.value("boxz", 50.0f);

            gravity = scene.value("gravity", false);
            react = scene.value("react", false);
            isothermal = scene.value("isothermal", true);
            render_water = scene.value("render_water", true);
            wall_collision = scene.value("wall_collision", false);
            mag_gravity = scene.value("mag_gravity", 9.81f);

            if (react)
                initReaxParams();

            const auto &posx = scene["posx"];
            const auto &posy = scene["posy"];
            const auto &posz = scene["posz"];
            const auto &velx = scene["velx"];
            const auto &vely = scene["vely"];
            const auto &velz = scene["velz"];
            const auto &charges = scene["charge"];
            const auto &ZIndices = scene["ZIndex"];
            const auto &neutrons = scene["neutrons"];
            const auto &electrons = scene["electrons"];
            const auto &bondCounts = scene["boundCount"];
            const auto &chiralities = scene["chirality"];
            const auto &temperatures = scene["temperature"];
            const auto &bondorders = scene["bondOrder"];

            int32_t N = posx.size();

            atoms.reserve(N);
            data.x.reserve(N);
            data.y.reserve(N);
            data.z.reserve(N);
            data.vx.reserve(N);
            data.vy.reserve(N);
            data.vz.reserve(N);
            data.q.reserve(N);
            data.temperature.reserve(N);
            data.bond_orders.reserve(N);

            for (int32_t i = 0; i < N; ++i)
            {
                uint8_t Z = ZIndices[i];
                auto [sigma, epsilon] = constants::getAtomConstants(Z);
                float radius = sigma / 1.3f;
                uint8_t nNeutrons = neutrons[i];

                atom newAtom{};
                newAtom.ZIndex = Z;
                newAtom.sigma = sigma;
                newAtom.epsilon = epsilon;
                newAtom.radius = radius;
                newAtom.electrons = electrons[i];
                newAtom.NCount = nNeutrons;
                newAtom.mass = Z * MASS_PROTON + nNeutrons * MASS_NEUTRON + electrons[i] * MASS_ELECTRON;
                newAtom.chirality = chiralities[i];
                newAtom.bondCount = bondCounts[i];

                atoms.emplace_back(std::move(newAtom));

                emplace_pos(sf::Vector3f(posx[i], posy[i], posz[i]));
                emplace_vel(sf::Vector3f(velx[i], vely[i], velz[i]));

                data.q.emplace_back(charges[i]);
                data.temperature.emplace_back(temperatures[i]);
                data.bond_orders.emplace_back(bondorders[i]);
            }

            data.fx.assign(N, 0.0f);
            data.fy.assign(N, 0.0f);
            data.fz.assign(N, 0.0f);

            if (scene.contains("bonds"))
            {
                for (const auto &b_json : scene["bonds"])
                {
                    bond b{};
                    b.centralAtom = b_json["central"];
                    b.bondedAtom = b_json["bonded"];
                    b.equilibriumLength = b_json["equilibrium"];
                    b.type = static_cast<BondType>(b_json["type"]);
                    b.k = b_json["k"];
                    bonds.emplace_back(b);
                }
            }

            if (scene.contains("angles"))
            {
                for (const auto &a_json : scene["angles"])
                {
                    angle a{};
                    a.A = a_json["A"];
                    a.B = a_json["B"];
                    a.C = a_json["C"];
                    a.K = a_json["K"];
                    a.rad = a_json["rad"];
                    angles.emplace_back(a);
                }
            }

            if (scene.contains("dihedrals"))
            {
                for (const auto &d_json : scene["dihedrals"])
                {
                    dihedral_angle d{};
                    d.A = d_json["A"];
                    d.B = d_json["B"];
                    d.C = d_json["C"];
                    d.D = d_json["D"];
                    d.K = d_json["K"];
                    d.periodicity = d_json["periodicity"];
                    d.rad = d_json["rad"];
                    dihedral_angles.emplace_back(d);
                }
            }

            if (scene.contains("subsets"))
            {
                for (const auto &s_json : scene["subsets"])
                {
                    subset s{};
                    s.mainAtomIdx = s_json["mainAtom"];
                    s.bondedSubsetIdx = s_json.value("bonded", UINT32_MAX);
                    s.bondingSubsetIdx = s_json.value("bonding", UINT32_MAX);

                    s.connectedBegin = s_json.value("connected_begin", UINT32_MAX);
                    s.connectedCount = s_json.value("connected_count", UINT32_MAX);
                    s.hydrogenCount = s_json.value("hydrogen_count", UINT32_MAX);
                    s.hydrogenBegin = s_json.value("hydrogen_begin", UINT32_MAX);

                    subsets.emplace_back(s);
                }
            }

            if (scene.contains("molecules"))
            {
                for (const auto &m_json : scene["molecules"])
                {
                    molecule m{};

                    m.atomBegin = m_json["atomBegin"];
                    m.atomCount = m_json["atomCount"];
                    m.bondBegin = m_json["bondBegin"];
                    m.bondCount = m_json["bondCount"];
                    m.angleBegin = m_json["angleBegin"];
                    m.angleCount = m_json["angleCount"];
                    m.dihedralBegin = m_json["dihedralBegin"];
                    m.dihedralCount = m_json["dihedralCount"];
                    m.subsetBegin = m_json["subsetBegin"];
                    m.subsetCount = m_json["subsetCount"];
                    m.exclude = m_json.value("exclude", false);

                    m.exclude = m_json.value("exclude", false);
                    molecules.emplace_back(std::move(m));
                }
            }

            rebuildBondTopology();

            std::cout << "[Simulation] Successfully loaded scene: " << path.filename()
                      << " (" << atoms.size() << " atoms, "
                      << molecules.size() << " molecules)\n";
        }

    } // namespace fun
} // namespace sim
