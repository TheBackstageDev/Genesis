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

        universe::universe(universe_create_info &create_info)
            : box(create_info.box), react(create_info.reactive),
              gravity(create_info.has_gravity), mag_gravity(create_info.mag_gravity),
              wall_collision(create_info.wall_collision), isothermal(create_info.isothermal),
              render_water(create_info.render_water)
        {
            if (react)
                initReaxParams();
        }

        universe::universe(const std::filesystem::path path)
        {
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

            size_t idx = 0;
            for (const auto &[i, j] : edges)
            {
                sf::Vector2f p1 = project(window, corners[i]);
                sf::Vector2f p2 = project(window, corners[j]);

                if (p1.x < -1000.f || p2.x < -1000.f)
                    continue;

                const sf::Vector3f &a = corners[i];
                const sf::Vector3f &b = corners[j];

                sf::Vector3f camToA = a - cam.eye();
                sf::Vector3f camToB = b - cam.eye();

                float dot = camToA.normalized().dot(camToB.normalized());
                const float THRESHOLD = 0.3f;
                if (dot < THRESHOLD)
                    continue;

                lines[idx++] = sf::Vertex(p1, edgeColor);
                lines[idx++] = sf::Vertex(p2, edgeColor);
            }

            window.getWindow().draw(lines);
        }

        void universe::draw(core::window_t &window, bool letter, bool lennardBall)
        {
            std::vector<size_t> drawOrder(atoms.size());
            std::iota(drawOrder.begin(), drawOrder.end(), 0);
            sf::Vector3f eye = cam.eye();

            std::sort(drawOrder.begin(), drawOrder.end(),
                      [&](size_t a, size_t b)
                      { return (pos(a) - eye).lengthSquared() > (pos(b) - eye).lengthSquared(); });

            drawBox(window);

            std::vector<size_t> no_draw{};

            if (!render_water)
                for (auto &m : molecules)
                {
                    if (m.water)
                    {
                        no_draw.emplace_back(m.atomIdx[0]);
                        no_draw.emplace_back(m.atomIdx[1]);
                        no_draw.emplace_back(m.atomIdx[2]);
                    }
                }

            if (!react)
                drawBonds(window, no_draw);
            else
                drawReactiveBonds(window);

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                if (std::find(no_draw.begin(), no_draw.end(), drawOrder[i]) != no_draw.end())
                    continue;

                sf::Vector2f p = project(window, pos(drawOrder[i]));
                if (p.x < -1000)
                    continue;

                float temp = calculateAtomTemperature(drawOrder[i]);
                if (atoms[drawOrder[i]].ZIndex == 1)
                    drawHydrogenBond(window, drawOrder[i]);

                float camDistance = (pos(drawOrder[i]) - cam.eye()).length();

                float T = calculateAtomTemperature(i);
                atoms[drawOrder[i]].draw(T, p, camDistance, data.q[drawOrder[i]], window, letter, lennardBall);
            }
        }

        float smootherstep(float edge0, float edge1, float x)
        {
            x = std::clamp((x - edge0) / (edge1 - edge0), 0.f, 1.f);

            return x * x * (3.0f - 2.0f * x);
        }

        void universe::drawChargeField(core::window_t &window)
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

                        for (size_t d = 0; d < 14; ++d)
                        {
                            auto &cell = cells[getCellID(ix + offsets[d][0], iy + offsets[d][1], iz + offsets[d][2])];
                            for (size_t i = 0; i < cell.size(); ++i)
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

                        float distance = (worldPos - cam.eye()).length();
                        float size = step * GRID / distance;

                        blobs.emplace_back(screen, distance, size, col);
                    }

            std::sort(blobs.begin(), blobs.end(), [=](const Blob &a, const Blob &b)
                      { return a.depth > b.depth; });

            sf::VertexArray billboards(sf::PrimitiveType::Triangles);
            billboards.resize(blobs.size() * 6);

            size_t idx = 0;
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
            window.getWindow().draw(billboards, states);
        }

        void universe::drawBonds(core::window_t &window, const std::vector<size_t> &no_draw)
        {
            sf::Vector2f dimensions = window.getWindow().getView().getSize();

            for (size_t b = 0; b < bonds.size(); ++b)
            {
                const bond &bond = bonds[b];
                if (std::find(no_draw.begin(), no_draw.end(), bond.centralAtom) != no_draw.end())
                    continue;

                const sf::Vector3f &pCentral = pos(bond.centralAtom); // pB
                const sf::Vector3f &pBonded = pos(bond.bondedAtom);   // pA

                sf::Vector2f s1 = cam.project(pCentral, dimensions.x, dimensions.y); // center
                sf::Vector2f s2 = cam.project(pBonded, dimensions.x, dimensions.y);  // end

                if (s1.x <= -9999 || s2.x <= -9999)
                    continue;
                if ((pCentral - pBonded).length() > 5.f)
                    continue;
                if ((pCentral - cam.eye()).length() > 200.f)
                    continue;

                sf::Vector2f dir = s2 - s1;
                float len = dir.length();
                if (len < 2.0f)
                    continue;
                dir /= len;

                sf::Vector2f perp{-dir.y, dir.x};

                uint8_t lines = static_cast<uint8_t>(bond.type);
                float shrink = 2.f;   // pixels
                float spacing = 10.f; // pixels between lines

                std::vector<sf::Vertex> verts;
                verts.reserve(lines * 2);

                for (int32_t i = 0; i < lines; ++i)
                {
                    float offset = spacing * (i - (lines - 1) * 0.5f);
                    sf::Vector2f off = perp * offset / (pCentral - cam.eye()).length();

                    sf::Vector2f start = s1 + dir * shrink + off;
                    sf::Vector2f end = s2 - dir * shrink + off;

                    verts.emplace_back(start, sf::Color::White);
                    verts.emplace_back(end, sf::Color::White);
                }

                window.getWindow().draw(verts.data(), verts.size(), sf::PrimitiveType::Lines);
            }
        }

        void universe::drawCylinder(core::window_t &window,
                                    size_t i, size_t j,
                                    int32_t segments,
                                    float radius,
                                    sf::Color color)
        {
            const sf::Vector3f a = pos(i);
            const sf::Vector3f b = pos(j);

            sf::Vector3f dir = b - a;
            float length = dir.length();
            if (length < 0.01f)
                return;
            dir /= length;

            sf::Vector2f s1 = project(window, a);
            sf::Vector2f s2 = project(window, b);
            if (s1.x < -9999 || s2.x < -9999)
                return;

            sf::Vector3f up = (std::abs(dir.z) < 0.9f)
                                  ? sf::Vector3f(0, 0, 1)
                                  : sf::Vector3f(0, 1, 0);
            sf::Vector3f right = dir.cross(up).normalized();
            up = right.cross(dir).normalized();

            sf::VertexArray cylinder(sf::PrimitiveType::TriangleStrip, (segments + 1) * 2);

            for (int32_t k = 0; k <= segments; ++k)
            {
                float angle = 2.0f * M_PI * k / segments;
                float c = std::cos(angle);
                float s = std::sin(angle);

                sf::Vector3f offset = (right * c + up * s) * radius;

                sf::Vector3f p1 = a + offset;
                sf::Vector3f p2 = b + offset;

                sf::Vector2f screen1 = project(window, p1);
                sf::Vector2f screen2 = project(window, p2);

                if (screen1.x < -1000.f || screen2.x < -1000.f)
                    continue;

                float depth1 = (p1 - cam.eye()).length();
                float depth2 = (p2 - cam.eye()).length();
                float avg_depth = (depth1 + depth2) * 0.5f;

                uint8_t alpha = static_cast<uint8_t>(255 * std::exp(-avg_depth * 0.003f));
                color.a = alpha;

                cylinder[k * 2 + 0] = sf::Vertex(screen1, color);
                cylinder[k * 2 + 1] = sf::Vertex(screen2, color);
            }

            window.getWindow().draw(cylinder);
        }

        void universe::drawReactiveBonds(core::window_t &window)
        {
            for (auto &bond : reactive_bonds)
            {
                sf::Vector2f dimensions = window.getWindow().getView().getSize();

                const sf::Vector3f &pCentral = pos(bond.i); // pB
                const sf::Vector3f &pBonded = pos(bond.j);  // pA

                sf::Vector2f s1 = cam.project(pCentral, dimensions.x, dimensions.y); // center
                sf::Vector2f s2 = cam.project(pBonded, dimensions.x, dimensions.y);  // end

                if (s1.x <= -9999 || s2.x <= -9999)
                    continue;
                if ((pCentral - pBonded).length() > 8.f)
                    continue;
                if ((pCentral - cam.eye()).length() > 200.f)
                    continue;

                sf::Vector2f dir = (s2 - s1).normalized();
                float len = dir.length();
                dir /= len;

                sf::Vector2f perp{-dir.y, dir.x};

                if (bond.type != BondType::PARTIAL)
                {
                    uint8_t lines = static_cast<uint8_t>(bond.type);
                    float shrink = 2.f;   // pixels
                    float spacing = 10.f; // pixels between lines

                    std::vector<sf::Vertex> verts;
                    verts.reserve(lines * 2);

                    for (int32_t i = 0; i < lines; ++i)
                    {
                        float offset = spacing * (i - (lines - 1) * 0.5f);
                        sf::Vector2f off = perp * offset / (pCentral - cam.eye()).length();

                        sf::Vector2f start = s1 + dir * shrink + off;
                        sf::Vector2f end = s2 - dir * shrink + off;

                        verts.emplace_back(start, sf::Color::White);
                        verts.emplace_back(end, sf::Color::White);
                    }

                    window.getWindow().draw(verts.data(), verts.size(), sf::PrimitiveType::Lines);
                }
            }
        }

        void universe::drawHydrogenBond(core::window_t &window, size_t H)
        {
            if (timeStep == 0 || cells.size() == 0)
                return;

            if (data.q[H] > -0.2f && data.q[H] < 0.2f)
                return;

            size_t D = SIZE_MAX;
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
            constexpr float MIN_COS_ANGLE = 0.9f;
            constexpr float STRONG_ENERGY = -300.0f;
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

            for (size_t d = 0; d < 14; ++d)
            {
                int32_t iix = ix + offsets[d][0];
                int32_t iiy = iy + offsets[d][1];
                int32_t iiz = iz + offsets[d][2];

                int32_t cell_id = iix + cx * (iiy + cy * iiz);

                if (cell_id < -1 || cell_id >= (cx * cy * cz))
                    continue;

                for (size_t A : cells[cell_id])
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

            size_t vert_idx = 0;

            for (size_t ix = 0; ix < cx; ++ix)
                for (size_t iy = 0; iy < cy; ++iy)
                    for (size_t iz = 0; iz < cz; ++iz)
                    {
                        size_t cell_id = ix + cx * (iy + cy * iz);
                        bool has_atoms = !cells[cell_id].empty();

                        sf::Color color = has_atoms ? hotColor : cellColor;

                        float x0 = ix * box.x;
                        float y0 = iy * box.y;
                        float z0 = iz * box.z;
                        float x1 = x0 + box.x;
                        float y1 = y0 + box.y;
                        float z1 = z0 + box.z;

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

                            sf::Vector3f camToA = c[edges[e][0]] - cam.eye();
                            sf::Vector3f camToB = c[edges[e][1]] - cam.eye();
                            if (camToA.normalized().dot(camToB.normalized()) < 0.2f)
                                continue;

                            grid[vert_idx++] = sf::Vertex(a, color);
                            grid[vert_idx++] = sf::Vertex(b, color);
                        }
                    }

            window.getWindow().draw(grid);
        }

        void universe::log(size_t step)
        {
        }

        size_t universe::createAtom(sf::Vector3f p, sf::Vector3f v, uint8_t ZIndex, uint8_t numNeutrons, uint8_t numElectron, int32_t chirality)
        {
            atom newAtom{};
            newAtom.ZIndex = ZIndex;

            emplace_vel(v);
            emplace_pos(p);

            std::pair<float, float> constants = constants::getAtomConstants(ZIndex);

            newAtom.sigma = constants.first;
            newAtom.epsilon = constants.second;
            newAtom.radius = constants.first / 1.3f;
            newAtom.electrons = numElectron;
            newAtom.NCount = numNeutrons;
            newAtom.mass = ZIndex * MASS_PROTON + numNeutrons * MASS_NEUTRON + numElectron * MASS_ELECTRON;
            newAtom.chirality = chirality;
            newAtom.bondCount = 0;

            atoms.emplace_back(std::move(newAtom));
            data.q.emplace_back(ZIndex - numElectron);
            data.fx.resize(atoms.size());
            data.fy.resize(atoms.size());
            data.fz.resize(atoms.size());
            data.bond_orders.resize(atoms.size());

            return atoms.size() - 1;
        }

        void universe::createBond(size_t idx1, size_t idx2, BondType type)
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

        size_t universe::createSubset(const def_subset &nSub, const size_t baseAtom, const size_t baseSubset)
        {
            subset nSubset{};
            nSubset.mainAtomIdx = nSub.mainAtomIdx + baseAtom;
            nSubset.bondedSubsetIdx = nSub.bondedSubset == SIZE_MAX ? SIZE_MAX : nSub.bondedSubset + baseSubset;
            nSubset.bondingSubsetIdx = nSub.bondingSubset == SIZE_MAX ? SIZE_MAX : nSub.bondingSubset + baseSubset;

            std::vector<size_t> connected{nSub.connectedIdx};
            std::vector<size_t> hydrogens{nSub.hydrogensIdx};
            std::vector<uint8_t> neighbourZs;
            neighbourZs.reserve(nSub.connectedIdx.size() + nSub.hydrogensIdx.size());

            for (size_t i = 0; i < connected.size(); ++i)
            {
                const size_t bondedAtom = connected[i];

                connected[i] += baseAtom;
                neighbourZs.emplace_back(atoms[bondedAtom].ZIndex);
            }

            for (size_t h = 0; h < hydrogens.size(); ++h)
            {
                const size_t bondedAtom = hydrogens[h];

                hydrogens[h] += baseAtom;
                neighbourZs.emplace_back(1);
            }

            nSubset.connectedIdx = std::move(connected);
            nSubset.hydrogenIdx = std::move(hydrogens);
            subsets.emplace_back(std::move(nSubset));

            return subsets.size() - 1; // Index
        }

        void universe::createMolecule(molecule_structure structure, sf::Vector3f pos, sf::Vector3f vel)
        {
            size_t baseAtomIndex = atoms.size();

            molecule nMolecule{};

            for (size_t i = 0; i < structure.atoms.size(); ++i)
            {
                const def_atom &a = structure.atoms[i];
                if (pos.z != 0)
                    structure.positions[i].z += 0.01f * i;
                createAtom(structure.positions[i] + pos, vel, a.ZIndex, a.NIndex, a.ZIndex - a.charge, a.chirality);
                data.q[i] += structure.atoms[i].charge;

                nMolecule.atomIdx.emplace_back(baseAtomIndex + i);
            }

            if (!react)
            {
                size_t baseBondIndex = bonds.size();
                for (size_t b = 0; b < structure.bonds.size(); ++b)
                {
                    const def_bond &db = structure.bonds[b];
                    size_t central = baseAtomIndex + db.centralAtomIdx;
                    size_t bonded = baseAtomIndex + db.bondingAtomIdx;
                    createBond(bonded, central, db.type);
                    nMolecule.bondIdx.emplace_back(baseBondIndex + b);
                }
            }

            size_t baseSubset = subsets.size();
            for (size_t s = 0; s < structure.subsets.size(); ++s)
            {
                createSubset(structure.subsets[s], baseAtomIndex, baseSubset);
                nMolecule.subsetIdx.emplace_back(baseSubset + s);
            }

            size_t baseAngle = angles.size();
            for (size_t a = 0; a < structure.angles.size(); ++a)
            {
                angle angle = structure.angles[a];
                angle.A += baseAtomIndex;
                angle.B += baseAtomIndex;
                angle.C += baseAtomIndex;

                angles.emplace_back(angle);
                nMolecule.angles.emplace_back(baseAngle + a);
            }

            size_t baseDihedral = dihedral_angles.size();
            for (size_t a = 0; a < structure.dihedral_angles.size(); ++a)
            {
                dihedral_angle angle = structure.dihedral_angles[a];
                angle.A += baseAtomIndex;
                angle.B += baseAtomIndex;
                angle.C += baseAtomIndex;
                angle.D += baseAtomIndex;

                dihedral_angles.emplace_back(angle);
                nMolecule.dihedrals.emplace_back(baseDihedral + a);
            }

            nMolecule.name = moleculeName(nMolecule.subsetIdx);

            if (!react)
            {
                balanceMolecularCharges(subsets[baseSubset]);
                rebuildBondTopology();
            }

            if (structure.atoms.size() == 3 && structure.atoms[0].ZIndex == 8 && structure.atoms[1].ZIndex == 1 && structure.atoms[2].ZIndex == 1)
                nMolecule.water = true;

            molecules.emplace_back(std::move(nMolecule));
        }

        void universe::balanceMolecularCharges(subset &mol)
        {
            auto it = std::find_if(subsets.begin(), subsets.end(), [&](const subset &s)
                                   { return &s == &mol; });
            if (it == subsets.end())
                return;

            size_t startIdx = std::distance(subsets.begin(), it);

            std::vector<size_t> atomsToBalance;
            std::unordered_set<size_t> visited;
            std::queue<size_t> q;
            q.push(startIdx);
            visited.insert(startIdx);

            while (!q.empty())
            {
                size_t currentIdx = q.front();
                q.pop();
                const auto &currentSubset = subsets[currentIdx];

                atomsToBalance.push_back(currentSubset.mainAtomIdx);
                atomsToBalance.insert(atomsToBalance.end(),
                                      currentSubset.connectedIdx.begin(), currentSubset.connectedIdx.end());
                atomsToBalance.insert(atomsToBalance.end(),
                                      currentSubset.hydrogenIdx.begin(), currentSubset.hydrogenIdx.end());

                size_t nextIdx = currentSubset.bondedSubsetIdx;
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

            std::vector<std::vector<size_t>> neighborLists(atoms.size());
            for (const auto &bond : bonds)
            {
                size_t a = bond.bondedAtom;
                size_t b = bond.centralAtom;
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
            for (size_t i = 0; i < atomsToBalance.size(); ++i)
            {
                size_t idx = atomsToBalance[i];
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
                for (size_t i = 0; i < atomsToBalance.size(); ++i)
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
                for (size_t idx : atomsToBalance)
                {
                    data.q[idx] += adjustment;
                }
            }
        }

        void universe::boundCheck(size_t i)
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

        float universe::ljPot(size_t i, size_t j)
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

        sf::Vector3f universe::ljForce(size_t i, size_t j)
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

        sf::Vector3f universe::coulombForce(size_t i, size_t j, sf::Vector3f &dr_vec)
        {
            float dr = dr_vec.length();

            if (dr < EPSILON || dr > COULOMB_CUTOFF)
                return {0.f, 0.f, 0.f};

            float qq = data.q[i] * data.q[j];
            if (qq == 0.f)
                return {0.f, 0.f, 0.f};

            float forceMag = COULOMB_K * data.q[i] * data.q[j] / dr;
            return -forceMag * dr_vec / dr;
        }

        void universe::calcBondForces()
        {
            for (size_t i = 0; i < bonds.size(); ++i)
            {
                bond &bond = bonds[i];

                size_t idx1 = bond.bondedAtom;
                size_t idx2 = bond.centralAtom;
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
            // size_t count = 0;

            for (const angle &ang : angles)
            {
                size_t i = ang.A; // left atom
                size_t j = ang.B; // central atom
                size_t k = ang.C; // right atom

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
            // size_t count = 0;

            for (size_t ix = 0; ix < cx; ++ix)
                for (size_t iy = 0; iy < cy; ++iy)
                    for (size_t iz = 0; iz < cz; ++iz)
                    {
                        size_t cell_id = ix + cx * (iy + cy * iz);

                        const std::vector<size_t> &cell = cells[cell_id];

                        for (int32_t dz = -1; dz <= 1; ++dz)
                            for (int32_t dy = -1; dy <= 1; ++dy)
                                for (int32_t dx = -1; dx <= 1; ++dx)
                                {
                                    int32_t n_ix = (ix + dx + cx) % cx;
                                    int32_t n_iy = (iy + dy + cy) % cy;
                                    int32_t n_iz = (iz + dz + cz) % cz;
                                    size_t neighbour_id = getCellID(n_ix, n_iy, n_iz);

                                    const std::vector<size_t> &neighbour_cell = cells[neighbour_id];

                                    for (size_t ii = 0; ii < cell.size(); ++ii)
                                    {
                                        const size_t &i = cell[ii];
                                        for (size_t jj = 0; jj < neighbour_cell.size(); ++jj)
                                        {
                                            const size_t &j = neighbour_cell[jj];

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
            for (size_t d = 0; d < dihedral_angles.size(); ++d)
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

                float torque_magnitude = -d_angle.K * d_angle.periodicity *
                                         std::sin(d_angle.periodicity * phi - d_angle.K);

                sf::Vector3f axis = (pc - pb).normalized();

                sf::Vector3f rA = (pa - (pb.length() > 0.f ? pb.normalized() : sf::Vector3f{0.f, 0.f, 0.f}));
                sf::Vector3f torqueA = axis.cross(rA).normalized() * torque_magnitude;

                sf::Vector3f rD = (pd - pc).normalized();
                sf::Vector3f torqueD = axis.cross(rD).normalized() * (-torque_magnitude);

                add_force(d_angle.A, torqueA * 0.5f);
                add_force(d_angle.B, torqueA * 0.5f - torqueD * 0.5f);
                add_force(d_angle.C, torqueD * 0.5f);
                add_force(d_angle.D, -torqueD * 0.5f);
            }
        }

        void universe::calcElectrostaticForces()
        {
            for (size_t ix = 0; ix < cx; ++ix)
                for (size_t iy = 0; iy < cy; ++iy)
                    for (size_t iz = 0; iz < cz; ++iz)
                    {
                        size_t cell_id = getCellID(ix, iy, iz);

                        const std::vector<size_t> &cell = cells[cell_id];

                        for (size_t ii = 0; ii < cell.size(); ++ii)
                        {
                            size_t i = cell[ii];
                            if (data.q[i] == 0.f)
                                continue;
                            for (size_t jj = ii + 1; jj < cell.size(); ++jj)
                            {
                                const size_t &j = cell[jj];

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

                        for (size_t d = 1; d < 14; ++d)
                        {
                            int32_t n_ix = ix + offsets[d][0], n_iy = iy + offsets[d][1], n_iz = iz + offsets[d][2];
                            size_t neighbour_id = getCellID(n_ix, n_iy, n_iz);

                            const std::vector<size_t> &neighbour_cell = cells[neighbour_id];

                            for (size_t ii = 0; ii < cell.size(); ++ii)
                            {
                                const size_t &i = cell[ii];
                                if (data.q[i] == 0.f)
                                    continue;
                                for (size_t jj = 0; jj < neighbour_cell.size(); ++jj)
                                {
                                    const size_t &j = neighbour_cell[jj];

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

        std::vector<sf::Vector3f> universe::processCellUnbonded(size_t ix, size_t iy, size_t iz)
        {
            thread_local size_t local_count = 0;
            thread_local float local_virial = 0.f;

            size_t cellID = getCellID(ix, iy, iz);
            const auto &cell = cells[cellID];
            thread_local std::vector<sf::Vector3f> local_forces(atoms.size(), {0, 0, 0});

            if (local_forces.size() != atoms.size()) 
                local_forces.assign(atoms.size(), {0.0f, 0.0f, 0.0f});

            std::fill(local_forces.begin(), local_forces.end(), sf::Vector3f{0,0,0});

            for (size_t ii = 0; ii < cell.size(); ++ii)
            {
                const size_t i = cell[ii];

                for (size_t jj = ii + 1; jj < cell.size(); ++jj)
                {
                    const size_t j = cell[jj];
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
                        size_t neighbor_id = getCellID(n_ix, n_iy, n_iz);

                        if (neighbor_id == cellID)
                            continue;

                        const auto &neighbor_cell = cells[neighbor_id];

                        for (size_t ii = 0; ii < cell.size(); ++ii)
                        {
                            size_t i = cell[ii];

                            for (size_t jj = 0; jj < neighbor_cell.size(); ++jj)
                            {
                                size_t j = neighbor_cell[jj];

                                if (j <= i)
                                    continue;

                                sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                                float r = dr.length();

                                if (r > CELL_CUTOFF)
                                    continue;

                                if (areBonded(i, j))
                                    continue;

                                ++local_count;

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
            total_count.fetch_add(local_count, std::memory_order_relaxed);

            return local_forces;
        }

        void universe::calcBondedForcesParallel()
        {
            const int32_t n_threads = std::max(1, static_cast<int32_t>(std::thread::hardware_concurrency()));
            std::vector<std::future<std::vector<sf::Vector3f>>> futures;

            auto make_task = [this](auto&& func) {
                return [this, func](size_t start, size_t end) -> std::vector<sf::Vector3f> {
                    thread_local std::vector<sf::Vector3f> local_forces;
                    if (local_forces.size() != atoms.size()) {
                        local_forces.assign(atoms.size(), {0.0f, 0.0f, 0.0f});
                    }
                    std::fill(local_forces.begin(), local_forces.end(), sf::Vector3f{0,0,0});

                    func(start, end, local_forces);

                    return local_forces; 
                };
            };

            auto bond_func = [this](size_t start, size_t end, std::vector<sf::Vector3f>& lf) 
            {
                for (size_t i = start; i < end; ++i) {
                    const bond& b = bonds[i];
                    size_t a = b.bondedAtom;
                    size_t c = b.centralAtom;

                    sf::Vector3f dr = minImageVec(pos(c) - pos(a));
                    float len = dr.length();
                    if (len <= EPSILON) continue;

                    float dl = len - b.equilibriumLength;
                    sf::Vector3f force = (b.k * dl / len) * dr;

                    lf[a] += force;
                    lf[c] -= force;
                }
            };

            for (int32_t t = 0; t < n_threads; ++t) 
            {
                size_t start = bonds.size() * t / n_threads;
                size_t end   = bonds.size() * (t + 1) / n_threads;
                futures.emplace_back(std::async(std::launch::async, make_task(bond_func), start, end));
            }

            auto angle_func = [this](size_t start, size_t end, std::vector<sf::Vector3f>& lf) 
            {
                for (size_t a = start; a < end; ++a) {
                    const angle& ang = angles[a];
                    size_t i = ang.A, j = ang.B, k = ang.C;

                    sf::Vector3f r_ji = minImageVec(pos(i) - pos(j));
                    sf::Vector3f r_jk = minImageVec(pos(k) - pos(j));
                    float len_ji = r_ji.length();
                    float len_jk = r_jk.length();
                    if (len_ji < EPSILON || len_jk < EPSILON) continue;

                    sf::Vector3f u_ji = r_ji / len_ji;
                    sf::Vector3f u_jk = r_jk / len_jk;
                    float cos_theta = std::clamp(u_ji.dot(u_jk), -1.0f, 1.0f);
                    float sin_theta = std::sqrt(std::max(1.0f - cos_theta * cos_theta, 0.0f));
                    if (sin_theta < 1e-6f) sin_theta = 1e-6f;

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
                size_t start = angles.size() * t / n_threads;
                size_t end   = angles.size() * (t + 1) / n_threads;
                futures.emplace_back(std::async(std::launch::async, make_task(angle_func), start, end));
            }

            auto dihedral_func = [this](size_t start, size_t end, std::vector<sf::Vector3f>& lf) 
            {
                for (size_t d = start; d < end; ++d) {
                    const dihedral_angle& da = dihedral_angles[d];
                    float phi = calculateDihedral(pos(da.A), pos(da.B), pos(da.C), pos(da.D));

                    float target = da.rad;
                    if (target == 0.0f && da.periodicity == 1) {
                        int32_t chi = atoms[da.B].chirality ? atoms[da.B].chirality : atoms[da.C].chirality;
                        if (chi == 1)  target = (phi < M_PI) ? 1.047f : 5.236f;
                        if (chi == 2)  target = (phi < M_PI) ? 5.236f : 1.047f;
                    }

                    float diff = phi - target;
                    while (diff > M_PI) diff -= 2.0f * M_PI;
                    while (diff < -M_PI) diff += 2.0f * M_PI;

                    float torque = -da.K * da.periodicity * std::sin(da.periodicity * phi);

                    sf::Vector3f axis = (pos(da.C) - pos(da.B)).normalized();
                    sf::Vector3f rA = pos(da.A) - pos(da.B);
                    sf::Vector3f rD = pos(da.D) - pos(da.C);

                    sf::Vector3f tA = axis.cross(rA).normalized() * torque;
                    sf::Vector3f tD = axis.cross(rD).normalized() * (-torque);

                    lf[da.A] += tA * 0.5f;
                    lf[da.B] += tA * 0.5f - tD * 0.5f;
                    lf[da.C] += tD * 0.5f;
                    lf[da.D] -= tD * 0.5f;
                }
            };

            for (int32_t t = 0; t < n_threads; ++t) 
            {
                size_t start = dihedral_angles.size() * t / n_threads;
                size_t end   = dihedral_angles.size() * (t + 1) / n_threads;
                futures.emplace_back(std::async(std::launch::async, make_task(dihedral_func), start, end));
            }

            for (auto& fut : futures) {
                auto local_f = fut.get();
                for (size_t i = 0; i < atoms.size(); ++i) 
                    add_force(i, local_f[i]);
            }
        }

        void universe::calcUnbondedForcesParallel()
        {
            const int32_t n_threads = std::thread::hardware_concurrency();

            std::vector<std::future<std::vector<sf::Vector3f>>> futures;

            auto worker = [this](size_t start_flat, size_t end_flat) -> std::vector<sf::Vector3f>
            {
                std::vector<sf::Vector3f> thread_forces(atoms.size(), {0,0,0});

                for (size_t flat = start_flat; flat < end_flat; ++flat) 
                {
                    size_t iz = flat % cz;
                    size_t iy = (flat / cz) % cy;
                    size_t ix = flat / (cz * cy);

                    auto cell_forces = processCellUnbonded(ix, iy, iz);
                    for (size_t i = 0; i < atoms.size(); ++i) 
                        thread_forces[i] += cell_forces[i];
                }
                return thread_forces;
            };

            for (int32_t t = 0; t < n_threads; ++t) 
            {
                size_t start = (cells.size() * t) / n_threads;
                size_t end   = (cells.size() * (t + 1)) / n_threads;

                futures.push_back(std::async(std::launch::async, worker, start, end));
            }

            for (auto &fut : futures)
            {
                std::vector<sf::Vector3f> local_f = fut.get();
                for (size_t i = 0; i < atoms.size(); ++i)
                    add_force(i, local_f[i]);
            }

            std::cout << "Non Bonded Calcs: " << total_count.load() << std::endl;
        }

        void universe::calcUnbondedForces()
        {
            for (size_t ix = 0; ix < cx; ++ix)
                for (size_t iy = 0; iy < cy; ++iy)
                    for (size_t iz = 0; iz < cz; ++iz)
                    {
                        size_t cell_id = getCellID(ix, iy, iz);
                        const auto &cell = cells[cell_id];

                        for (size_t ii = 0; ii < cell.size(); ++ii)
                        {
                            size_t i = cell[ii];

                            for (size_t jj = ii + 1; jj < cell.size(); ++jj)
                            {
                                size_t j = cell[jj];

                                sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                                float r2 = dr.lengthSquared();
                                if (r2 > CELL_CUTOFF)
                                    continue;

                                if (areBonded(i, j))
                                    continue;

                                ++total_count;

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
                                    size_t neighbor_id = getCellID(n_ix, n_iy, n_iz);

                                    if (neighbor_id == cell_id)
                                        continue;

                                    const auto &neighbor_cell = cells[neighbor_id];

                                    for (size_t ii = 0; ii < cell.size(); ++ii)
                                    {
                                        size_t i = cell[ii];

                                        for (size_t jj = 0; jj < neighbor_cell.size(); ++jj)
                                        {
                                            size_t j = neighbor_cell[jj];

                                            if (j <= i)
                                                continue;

                                            sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                                            float r = dr.length();

                                            if (r > CELL_CUTOFF)
                                                continue;
                                            if (areBonded(i, j))
                                                continue;

                                            ++total_count;

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

        std::vector<float> total_bo{};
        std::vector<float> total_bo_uncorrected{};

        float taper(float r, float cutoff)
        {
            if (r >= cutoff)
                return 0.0f;
            float x = r / cutoff;
            return 1.f - powf(x, 7);
        }

        float universe::calculateUncorrectedBondOrder(size_t i, size_t j)
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

        float universe::calculateBondEnergy(size_t i, size_t j, float bo_sigma, float bo_pi, float bo_pp)
        {
            auto &parami = constants::getParams(atoms[i].ZIndex);
            auto &paramj = constants::getParams(atoms[j].ZIndex);

            return -(sqrt(parami.De_sigma * paramj.De_sigma) * bo_sigma + sqrt(parami.De_pi * paramj.De_pi) * bo_pi + sqrt(parami.De_pp * paramj.De_pp) * bo_pp);
        }

        float taper7(float x)
        {
            if (x >= 1.f)
                return 0.f;
            float x3 = x * x * x;
            float x4 = x3 * x;
            float x5 = x4 * x;
            return 1.f - 10.f * x3 + 35.f * x4 - 50.f * x5 + 35.f * x3 * x3 * x3 - 10.f * x3 * x3 * x3 * x3;
        }

        void universe::processReactivePair(size_t i, size_t j, float cutoff, float vis_thresh)
        {
            sf::Vector3f dr_vec = minImageVec(pos(j) - pos(i));
            float r = dr_vec.length();
            if (r >= cutoff || r < EPSILON)
                return;

            uint8_t Zi = atoms[i].ZIndex;
            uint8_t Zj = atoms[j].ZIndex;

            auto &parami = constants::getParams(Zi);
            auto &paramj = constants::getParams(Zj);
            auto pair_param = constants::getPairReaxParams(Zi, Zj);
            float gamma_ij = std::sqrtf(parami.gamma * paramj.gamma);
            float r_vdw_ij = std::sqrtf(parami.r_vdw * paramj.r_vdw);
            float D_ij = std::sqrtf(parami.D_vdw * paramj.D_vdw);
            float gamma3_ij = powf(1.0f / gamma_ij, 3.f);
            float S = std::cbrtf(r * r * r + gamma3_ij);
            float S_vdw = std::cbrtf(1.f + powf(gamma_ij * r, 3.f));
            constexpr float alpha = 12.f;
            constexpr float coulomb_constant = COULOMB_K * JOULE_TO_CAL;

            // Average parameters for pair
            float p_bo1_ij = 0.5f * (parami.p_bo1 + paramj.p_bo1);
            float p_bo2_ij = 0.5f * (parami.p_bo2 + paramj.p_bo2);
            float p_bo3_ij = 0.5f * (parami.p_bo3 + paramj.p_bo3);
            float p_bo4_ij = 0.5f * (parami.p_bo4 + paramj.p_bo4);
            float p_bo5_ij = 0.5f * (parami.p_bo5 + paramj.p_bo5);
            float p_bo6_ij = 0.5f * (parami.p_bo6 + paramj.p_bo6);
            float p_boc1_ij = 0.5f * (parami.p_boc1 + paramj.p_boc1);
            float p_boc2_ij = 0.5f * (parami.p_boc2 + paramj.p_boc2);
            float p_boc3_ij = 0.5f * (parami.p_boc3 + paramj.p_boc3);
            float p_boc4_ij = 0.5f * (parami.p_boc4 + paramj.p_boc4);
            float p_boc5_ij = 0.5f * (parami.p_boc5 + paramj.p_boc5);
            float p_be1_ij = 0.5f * (parami.p_be1 + paramj.p_be1);
            float p_be2_ij = 0.5f * (parami.p_be2 + paramj.p_be2);
            float r0_sigma_ij = 0.5f * (parami.r0_sigma + paramj.r0_sigma);
            float r0_pi_ij = 0.5f * (parami.r0_pi + paramj.r0_pi);
            float r0_pp_ij = 0.5f * (parami.r0_pp + paramj.r0_pp);
            float De_sig = std::sqrtf(parami.De_sigma * paramj.De_sigma);
            float De_pi = std::sqrtf(parami.De_pi * paramj.De_pi);
            float De_pp = std::sqrtf(parami.De_pp * paramj.De_pp);

            float BO_sigma = std::exp(p_bo1_ij * std::pow(r / r0_sigma_ij, p_bo2_ij));
            float BO_pi = 0.0f;
            float BO_pp = 0.0f;

            if (parami.r0_pi > 0.f && paramj.r0_pi > 0.f)
                BO_pi = std::exp(p_bo3_ij * std::pow(r / r0_pi_ij, p_bo4_ij));
            if (parami.r0_pp > 0.f && paramj.r0_pp > 0.f)
                BO_pp = std::exp(p_bo5_ij * std::pow(r / r0_pp_ij, p_bo6_ij));

            float BO_uncorr = BO_sigma + BO_pi + BO_pp;

            float val_i = constants::getUsualBonds(Zi);
            float val_j = constants::getUsualBonds(Zj);
            float lonePairs_i = constants::lonePairsBO(Zi, total_bo_uncorrected[i]);
            float lonePairs_j = constants::lonePairsBO(Zj, total_bo_uncorrected[j]);

            float Delta_i = total_bo_uncorrected[i] - val_i + lonePairs_i;
            float Delta_j = total_bo_uncorrected[j] - val_j + lonePairs_j;

            float Delta_boc_i = Delta_i;
            float Delta_boc_j = Delta_j;

            float f1 = 1.0f;
            if (Delta_i > 0.0f || Delta_j > 0.0f)
            {
                float f1_i = (Delta_i > 0.0f) ? 1.0f / (1.0f + std::exp(p_boc1_ij * Delta_i)) : 1.0f;
                float f1_j = (Delta_j > 0.0f) ? 1.0f / (1.0f + std::exp(p_boc1_ij * Delta_j)) : 1.0f;
                f1 = 0.5f * (f1_i + f1_j);
            }

            float f2 = std::exp(-p_boc1_ij * Delta_i) + std::exp(-p_boc1_ij * Delta_j);
            float f3 = -1.0f / p_boc2_ij * std::log(0.5f * (std::exp(-p_boc2_ij * Delta_i) + std::exp(-p_boc2_ij * Delta_j)));

            float f4 = 1.0f;
            float f5 = 1.0f;

            if (Delta_i > -1.0f && BO_pi > 0.f)
                f4 = 1.0f / (1.0f + std::exp(-p_boc5_ij * (p_boc4_ij + Delta_i)));
            if (Delta_j > -1.0f && BO_pi > 0.f)
                f5 = 1.0f / (1.0f + std::exp(-p_boc5_ij * (p_boc4_ij + Delta_j)));

            // Corrected bond orders
            float BO_sigma_corr = BO_sigma * f1;
            float BO_pi_corr = BO_pi * f1 * f1 * f4 * f5;
            float BO_pp_corr = BO_pp * f1 * f1 * f4 * f5;
            float BO_corr = BO_sigma_corr + BO_pi_corr + BO_pp_corr;
            total_bo[i] += BO_corr;
            total_bo[j] += BO_corr;

            float dE_bond_dr = 0.0f;
            if (BO_corr > 1e-4f)
            {
                // Uncorrected derivatives
                float dBO_sigma_dr = BO_sigma * (p_bo1_ij * p_bo2_ij) * std::pow(r / r0_sigma_ij, p_bo2_ij - 1.0f) / (r0_sigma_ij + 1e-9f);
                float dBO_pi_dr = (BO_pi > 1e-12f) ? BO_pi * (p_bo3_ij * p_bo4_ij) * std::pow(r / r0_pi_ij, p_bo4_ij - 1.0f) / (r0_pi_ij + 1e-9f) : 0.0f;
                float dBO_pp_dr = (BO_pp > 1e-12f) ? BO_pp * (p_bo5_ij * p_bo6_ij) * std::pow(r / r0_pp_ij, p_bo6_ij - 1.0f) / (r0_pp_ij + 1e-9f) : 0.0f;

                float dBO_sigma_corr_dr = dBO_sigma_dr * f1;
                float dBO_pi_corr_dr = dBO_pi_dr * f1 * f1 * f4 * f5;
                float dBO_pp_corr_dr = dBO_pp_dr * f1 * f1 * f4 * f5;

                float pow_bo_sigma = std::pow(BO_sigma_corr, p_be2_ij);
                float exp_term = std::exp(p_be1_ij * (1.0f - pow_bo_sigma));
                float dE_sigma_dBO = -De_sig * exp_term - De_sig * BO_sigma_corr * exp_term * p_be1_ij * (-p_be2_ij) * std::pow(BO_sigma_corr, p_be2_ij - 1.0f);
                dE_bond_dr = dE_sigma_dBO * dBO_sigma_corr_dr - De_pi * dBO_pi_corr_dr - De_pp * dBO_pp_corr_dr;
            }

            if (BO_corr > 0.1f)
            {
                BondType type = BondType::SINGLE;
                if (BO_corr > 1.4f)
                    type = BondType::DOUBLE;
                if (BO_corr > 2.1f)
                    type = BondType::TRIPLE;
                reactive_bonds.emplace_back(i, j, BO_corr, type);
            }

            // Van der Waals
            float rho = r / (S_vdw * r_vdw_ij); // Use S_vdw or S
            float exp12 = std::expf(alpha * (1.f - rho));
            float exp6 = std::expf(alpha / 2.f * (1.f - rho));
            float E_vdW = D_ij * (exp12 - 2.f * exp6);

            float E_core = 0.0f, dE_core_dr = 0.0f;
            constexpr float R_CORE = 0.8f;
            if (r < R_CORE)
            {
                float ratio = r / R_CORE;
                float f = 1.0f - ratio;
                E_core = 4000.0f * f * f;
                dE_core_dr = -4000.0f * f / R_CORE;
            }

            // Full drho/dr
            float dS_vdw_dr = powf(gamma_ij, 3.f) * r * r / (S_vdw * S_vdw);
            float drho_dr = 1.0f / (S_vdw * r_vdw_ij) - rho / S_vdw * dS_vdw_dr;
            float dE_vdW_dr = D_ij * alpha * drho_dr * (exp6 - exp12) + dE_core_dr;

            // Taper
            float x = r / REACTION_CUTOFF;
            float T = taper7(x);
            float dT_dr = (r < REACTION_CUTOFF) ? -7.0f * powf(1.0f - x, 6) / REACTION_CUTOFF : 0.0f;

            // Coulomb
            float E_coul = coulomb_constant * data.q[i] * data.q[j] / S;
            float dE_coul_dr = coulomb_constant * data.q[i] * data.q[j] * (-r * r / powf(S, 4.f));
            float E_nonbond = T * (E_vdW + E_coul);
            float dE_nonbond_dr = dT_dr * (E_vdW + E_coul) + T * (dE_vdW_dr + dE_coul_dr);

            sf::Vector3f unit = dr_vec / r;
            sf::Vector3f F_bond = dE_bond_dr * unit;
            sf::Vector3f F_nonbond = dE_nonbond_dr * unit;
            sf::Vector3f F_total = F_bond + F_nonbond;

            add_force(i, F_total);
            add_force(j, -F_total); // Newton's 3rd law
        }

        void universe::handleReactiveForces()
        {
            reactive_bonds.clear();

            total_bo.resize(atoms.size());
            total_bo_uncorrected.resize(atoms.size());
            std::fill(total_bo.begin(), total_bo.end(), 0.0f);
            std::fill(total_bo_uncorrected.begin(), total_bo_uncorrected.end(), 0.0f);

            if (timeStep % REACTION_UPDATE_Q == 0)
            {
                data.q.resize(atoms.size());
                std::fill(data.q.begin(), data.q.end(), 0.0f);
            }

            for (size_t ix = 0; ix < cx; ++ix)
                for (size_t iy = 0; iy < cy; ++iy)
                    for (size_t iz = 0; iz < cz; ++iz)
                    {
                        size_t current_id = getCellID(ix, iy, iz);
                        auto &currentCell = cells[current_id];

                        for (size_t ii = 0; ii < currentCell.size(); ++ii)
                        {
                            size_t i = currentCell[ii];
                            for (size_t jj = ii + 1; jj < currentCell.size(); ++jj)
                            {
                                size_t j = currentCell[jj];

                                if (j <= i)
                                    continue;

                                sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                                float r = dr.length();

                                if (r >= REACTION_CUTOFF || r < 0.1f)
                                    continue;

                                float uncorr = calculateUncorrectedBondOrder(i, j);

                                total_bo_uncorrected[i] += uncorr * 0.5f;
                                total_bo_uncorrected[j] += uncorr * 0.5f;
                            }
                        }

                        for (int32_t dz = -1; dz <= 1; ++dz)
                            for (int32_t dy = -1; dy <= 1; ++dy)
                                for (int32_t dx = -1; dx <= 1; ++dx)
                                {
                                    if (dz == 0 && dy == 0 && dx == 0)
                                        continue;

                                    int32_t n_ix = ix + dx, n_iy = iy + dy, n_iz = iz + dz;
                                    size_t neighbour_id = getCellID(n_ix, n_iy, n_iz);

                                    if (neighbour_id == current_id)
                                        continue;

                                    auto &neighbourCell = cells[neighbour_id];

                                    for (size_t ii = 0; ii < currentCell.size(); ++ii)
                                    {
                                        size_t i = currentCell[ii];
                                        for (size_t jj = 0; jj < neighbourCell.size(); ++jj)
                                        {
                                            size_t j = neighbourCell[jj];

                                            if (j <= i)
                                                continue;

                                            sf::Vector3f dr = minImageVec(pos(j) - pos(i));
                                            float r = dr.length();

                                            if (r >= REACTION_CUTOFF || r < 0.1f)
                                                continue;

                                            float uncorr = calculateUncorrectedBondOrder(i, j);

                                            total_bo_uncorrected[i] += uncorr * 0.5f;
                                            total_bo_uncorrected[j] += uncorr * 0.5f;
                                        }
                                    }
                                }
                    }

            for (size_t ix = 0; ix < cx; ++ix)
                for (size_t iy = 0; iy < cy; ++iy)
                    for (size_t iz = 0; iz < cz; ++iz)
                    {
                        size_t current_id = getCellID(ix, iy, iz);
                        auto &currentCell = cells[current_id];

                        for (int32_t dz = -1; dz <= 1; ++dz)
                            for (int32_t dy = -1; dy <= 1; ++dy)
                                for (int32_t dx = -1; dx <= 1; ++dx)
                                {
                                    int32_t x = ix + dx, y = iy + dy, z = iz + dz;
                                    size_t neighbour_id = getCellID(x, y, z);
                                    auto &neighbourCell = cells[neighbour_id];

                                    for (size_t ii = 0; ii < currentCell.size(); ++ii)
                                    {
                                        size_t i = currentCell[ii];
                                        for (size_t jj = 0; jj < neighbourCell.size(); ++jj)
                                        {
                                            size_t j = neighbourCell[jj];

                                            if (j <= i)
                                                continue;

                                            processReactivePair(i, j, REACTION_CUTOFF);
                                        }
                                    }
                                }
                    }

            float sumQ = std::accumulate(data.q.begin(), data.q.end(), 0.f);
            if (sumQ > std::numeric_limits<float>::epsilon())
            {
                float correction = -sumQ / static_cast<float>(data.q.size());
                for (float &q : data.q)
                {
                    q += correction;
                }
            }
        }

        void universe::update(float targetTemperature, float targetPressure)
        {
            data.old_fx.resize(atoms.size());
            data.old_fy.resize(atoms.size());
            data.old_fz.resize(atoms.size());
            std::swap(data.fx, data.old_fx);
            std::swap(data.fy, data.old_fy);
            std::swap(data.fz, data.old_fz);

            data.fx.assign(atoms.size(), 0.0f);
            data.fy.assign(atoms.size(), 0.0f);
            data.fz.assign(atoms.size(), 0.0f);

            total_virial = 0.0f;

            if (!react)
            {
                calcBondedForcesParallel();
                calcUnbondedForcesParallel();
            }
            else
                handleReactiveForces();

            setPressure(targetPressure);

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector3f a_old = old_force(i) / atoms[i].mass;
                sf::Vector3f a_new = force(i) / atoms[i].mass;

                if (gravity)
                    add_vel(i, sf::Vector3f(0.f, 0.f, -mag_gravity * DT));

                add_pos(i, vel(i) * DT + 0.5f * (a_old + a_new) * DT * DT);
                boundCheck(i);
            }

            std::fill(data.fx.begin(), data.fx.end(), 0.f);
            std::fill(data.fy.begin(), data.fy.end(), 0.f);
            std::fill(data.fz.begin(), data.fz.end(), 0.f);

            if (!react)
            {
                calcBondedForcesParallel();
                calcUnbondedForcesParallel();
            }
            else
                handleReactiveForces();

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector3f accel = (old_force(i) + force(i)) * (0.5f * DT / atoms[i].mass);
                add_vel(i, accel);
            }

            buildCells();
            setTemperature(targetTemperature);

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

            mu = std::clamp(mu, 0.8f, 1.2f);

            float scale = std::cbrt(mu);

            box.z *= scale;

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                data.z[i] *= scale;
            }
        }

        void universe::setTemperature(float kelvin)
        {
            if (timeStep % THERMOSTAT_INTERVAL != 0)
                return;

            float avg_KE = calculateKineticEnergy() / atoms.size();
            temp = (2.0f / 3.0f) * (avg_KE * KB);
            float lambda = sqrtf(kelvin / (temp + 1e-5f));

            if (timeStep > 100)
            {
                constexpr float tau = 0.6f;
                lambda = 1.0f + (lambda - 1.0f) / tau;
            }

            for (size_t i = 0; i < data.vx.size(); ++i)
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
            for (size_t i = 0; i < atoms.size(); ++i)
                kinetic_energy += 0.5f * atoms[i].mass * vel(i).lengthSquared();

            return kinetic_energy;
        }

        void universe::buildCells()
        {
            const float cutoff = react ? REACTION_CUTOFF : CELL_CUTOFF;
            if ((cutoff - box.x) > EPSILON || timeStep == 0)
            {
                cx = static_cast<size_t>(std::ceil(cutoff / box.x));
                cy = static_cast<size_t>(std::ceil(cutoff / box.y));
                cz = static_cast<size_t>(std::ceil(cutoff / box.z));
            }

            size_t ncells = cx * cy * cz;
            if (cells.size() != ncells)
                cells.assign(ncells, {});
            else
                for (auto &c : cells)
                    c.clear();

            const float inv = 1.0f / box.x;

            for (size_t i = 0; i < atoms.size(); ++i)
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

                size_t id = static_cast<size_t>(ix) + cx * (static_cast<size_t>(iy) + cy * static_cast<size_t>(iz));
                cells[id].emplace_back(i);
            }
        }

        float universe::calculateAtomTemperature(size_t i)
        {
            float ke = 0.5f * atoms[i].mass * vel(i).lengthSquared();
            return (2.0f / 3.0f) * ke * KB;
        }

        // Camera
        sf::Vector2f universe::project(core::window_t &window, const sf::Vector3f &p) const
        {
            auto size = window.getWindow().getView().getSize();
            return cam.project(p, size.x, size.y);
        }

        void universe::handleCamera(bool leftDown, bool rightDown, const sf::Vector2i &mousePos, float wheelDelta, const std::vector<sf::Keyboard::Key> &keys)
        {
            if (leftDown)
            {
                sf::Vector2i delta = mousePos - lastMouse;
                cam.azimuth -= delta.x * 0.3f;
                cam.elevation = std::clamp(cam.elevation - delta.y * 0.3f,
                                           -89.f, 89.f);
            }

            sf::Vector3f forward = cam.distance == 0.f ? sf::Vector3f(1.f, 0.f, 0.f) : (cam.target - cam.eye()).normalized();
            sf::Vector3f right = forward.cross({0, 0, 1}).normalized();
            sf::Vector3f up = right.cross(forward);

            if (rightDown)
            {
                sf::Vector2i delta = mousePos - lastMouse;
                cam.target += sf::Vector3f{right.x * delta.x, right.y, right.z} * cam.distance * 0.002f +
                              sf::Vector3f{up.x, up.y * delta.y, up.z} * cam.distance * 0.002f;
            }

            if (wheelDelta != 0.0f)
            {
                float base_speed = 1.5f;
                float speed = base_speed * (cam.distance / 50.0f);

                if (wheelDelta > 0)
                {
                    cam.distance *= std::pow(0.9f, speed * wheelDelta);
                }
                else
                {
                    cam.distance *= std::pow(1.1f, speed * -wheelDelta);
                }
                cam.distance = std::clamp(cam.distance, 10.0f, 500.0f);
            }

            for (auto k : keys)
            {
                if (k == sf::Keyboard::Key::W)
                    cam.target += forward * 0.5f;
                if (k == sf::Keyboard::Key::S)
                    cam.target -= forward * 0.5f;
                if (k == sf::Keyboard::Key::A)
                    cam.target -= right * 0.5f;
                if (k == sf::Keyboard::Key::D)
                    cam.target += right * 0.5f;
                if (k == sf::Keyboard::Key::Q)
                    cam.target -= up * 0.5f; // Down
                if (k == sf::Keyboard::Key::E)
                    cam.target += up * 0.5f; // Up
            }

            lastMouse = mousePos;
        }

        // TO DO
        std::string universe::moleculeName(const std::vector<size_t> &subsetIdx)
        {
            std::map<uint8_t, size_t> ZIndices;

            for (int32_t i = 0; i < subsetIdx.size(); ++i)
            {
                auto &s = subsets[subsetIdx[i]];
                size_t centralAtom = s.mainAtomIdx;
                ++ZIndices[atoms[centralAtom].ZIndex];
                ZIndices[1] += s.hydrogenIdx.size();
            }

            std::unordered_map<size_t, std::string> prefixes{
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
            constants::reaxParams[1] = {
                .p_boc1 = 6.500f,
                .p_boc2 = 6.500f,
                .p_boc3 = 0.000f,
                .p_boc4 = 0.000f,
                .p_boc5 = 0.000f,

                .p_bo1 = -0.040f,
                .p_bo2 = 5.000f,
                .p_bo3 = -0.040f,
                .p_bo4 = 5.000f,
                .p_bo5 = 0.000f,
                .p_bo6 = 0.000f,

                .p_be1 = 0.000f,
                .p_be2 = 1.000f,

                .r0_sigma = 0.960f,

                .r0_pi = 0.000f,
                .r0_pp = 0.000f,

                .De_sigma = 103.50f,
                .De_pi = 0.0f,
                .De_pp = 0.0f,

                .gamma = 1.100f,
                .r_vdw = 1.850f,
                .D_vdw = 0.015f,
            };

            constants::reaxParams[6] = {
                .p_boc1 = 1.540f,
                .p_boc2 = 6.500f,
                .p_bo1 = -0.200f,
                .p_bo2 = 5.000f,
                .p_bo3 = -0.100f,
                .p_bo4 = 8.000f,
                .p_bo5 = -0.200f,
                .p_bo6 = 6.000f,
                .r0_sigma = 1.515f,
                .r0_pi = 1.350f,
                .r0_pp = 1.200f,
                .De_sigma = 95.0f,
                .De_pi = 115.0f,
                .De_pp = 10.0f,
                .gamma = 1.60f,
                .r_vdw = 1.992f,
                .D_vdw = 0.006f,
            };

            constants::reaxParams[7] = {
                .p_boc1 = 2.0f,
                .p_boc2 = 6.5f,
                .p_boc3 = 1.2f,
                .p_boc4 = 1.5f,
                .p_boc5 = 7.0f,
                .p_bo1 = -0.18f,
                .p_bo2 = 6.0f,
                .p_bo3 = -0.18f,
                .p_bo4 = 7.0f,
                .p_bo5 = -0.18f,
                .p_bo6 = 7.0f,
                .p_be1 = 0.1f,
                .p_be2 = 1.0f,
                .r0_sigma = 1.45f,
                .r0_pi = 1.30f,
                .r0_pp = 1.10f,
                .De_sigma = 45.0f,
                .De_pi = 105.0f,
                .De_pp = 140.0f,
                .gamma = 1.85f,
                .r_vdw = 1.95f,
                .D_vdw = 0.012f,
            };

            // O (Z = 8)
            constants::reaxParams[8] = {
                .p_boc1 = 2.0f,
                .p_boc2 = 7.0f,
                .p_boc3 = 1.0f,
                .p_boc4 = 2.0f,
                .p_boc5 = 7.0f,
                .p_bo1 = -0.15f,
                .p_bo2 = 6.5f,
                .p_bo3 = -0.15f,
                .p_bo4 = 7.0f,
                .p_bo5 = 0.0f,
                .p_bo6 = 0.0f,
                .p_be1 = 0.2f,
                .p_be2 = 1.0f,
                .r0_sigma = 1.42f,
                .r0_pi = 1.20f,
                .r0_pp = 0.0f,
                .De_sigma = 40.0f,
                .De_pi = 110.0f,
                .De_pp = 0.0f,
                .gamma = 2.40f,
                .r_vdw = 1.82f,
                .D_vdw = 0.018f,
            };

            // F (Z = 9)
            constants::reaxParams[9] = {
                .p_boc1 = 3.0f,
                .p_boc2 = 7.5f,
                .p_boc3 = 0.0f,
                .p_boc4 = 0.0f,
                .p_boc5 = 0.0f,
                .p_bo1 = -0.15f,
                .p_bo2 = 7.0f,
                .p_bo3 = 0.0f,
                .p_bo4 = 0.0f,
                .p_bo5 = 0.0f,
                .p_bo6 = 0.0f,
                .p_be1 = 0.0f,
                .p_be2 = 1.0f,
                .r0_sigma = 1.35f,
                .r0_pi = 0.0f,
                .r0_pp = 0.0f,
                .De_sigma = 110.0f,
                .De_pi = 0.0f,
                .De_pp = 0.0f,
                .gamma = 2.90f,
                .r_vdw = 1.750f,
                .D_vdw = 0.080f,
            };

            // Na (Z = 11)
            constants::reaxParams[11] = {
                .p_boc1 = 1.0f,
                .p_boc2 = 8.0f,
                .p_boc3 = 0.0f,
                .p_boc4 = 0.0f,
                .p_boc5 = 0.0f,
                .p_bo1 = -0.08f,
                .p_bo2 = 7.5f,
                .p_bo3 = 0.0f,
                .p_bo4 = 0.0f,
                .p_bo5 = 0.0f,
                .p_bo6 = 0.0f,
                .p_be1 = 0.0f,
                .p_be2 = 1.0f,
                .r0_sigma = 2.30f,
                .r0_pi = 0.0f,
                .r0_pp = 0.0f,
                .De_sigma = 20.0f,
                .De_pi = 0.0f,
                .De_pp = 0.0f,
                .gamma = 1.20f,
                .r_vdw = 2.300f,
                .D_vdw = 0.020f,
            };

            // 12 = Magnesium
            constants::reaxParams[12] = {
                .p_boc1 = 1.2f, .p_boc2 = 8.0f, .p_boc3 = 0.0f, .p_boc4 = 0.0f, .p_boc5 = 0.0f, .p_bo1 = -0.1f, .p_bo2 = 7.0f, .p_bo3 = 0.0f, .p_bo4 = 0.0f, .p_bo5 = 0.0f, .p_bo6 = 0.0f, .r0_sigma = 2.20f, .r0_pi = 0.0f, .r0_pp = 0.0f, .De_sigma = 30.0f, .De_pi = 0.0f, .De_pp = 0.0f, .gamma = 1.300f, .r_vdw = 2.200f, .D_vdw = 0.030f};

            // 13 = Aluminium
            constants::reaxParams[13] = {
                .p_boc1 = 1.5f, .p_boc2 = 7.5f, .p_boc3 = 7.5f, .p_boc4 = 8.0f, .p_boc5 = 7.5f, .p_bo1 = -0.12f, .p_bo2 = 6.5f, .p_bo3 = 0.0f, .p_bo4 = 0.0f, .p_bo5 = 0.0f, .p_bo6 = 0.0f, .r0_sigma = 2.60f, .r0_pi = 0.0f, .r0_pp = 0.0f, .De_sigma = 70.0f, .De_pi = 0.0f, .De_pp = 0.0f, .gamma = 1.400f, .r_vdw = 2.600f, .D_vdw = 0.050f};

            // 14 = Silicon
            constants::reaxParams[14] = {
                .p_boc1 = 1.8f,
                .p_boc2 = 6.8f,
                .p_boc3 = 6.8f,
                .p_boc4 = 7.5f,
                .p_boc5 = 6.8f,
                .p_bo1 = -0.15f,
                .p_bo2 = 6.0f,
                .p_bo3 = -0.15f,
                .p_bo4 = 6.5f,
                .p_bo5 = 0.0f,
                .p_bo6 = 0.0f,
                .p_be1 = 0.0f,
                .p_be2 = 1.0f,
                .r0_sigma = 2.35f,
                .r0_pi = 2.10f,
                .r0_pp = 0.0f,
                .De_sigma = 75.0f,
                .De_pi = 50.0f,
                .De_pp = 0.0f,
                .gamma = 1.60f,
                .r_vdw = 2.350f,
                .D_vdw = 0.075f,
            };

            // 15 = Phosphorus
            constants::reaxParams[15] = {
                .p_boc1 = 2.0f, .p_boc2 = 7.0f, .p_boc3 = 7.0f, .p_boc4 = 8.0f, .p_boc5 = 7.0f, .p_bo1 = -0.16f, .p_bo2 = 6.5f, .p_bo3 = -0.16f, .p_bo4 = 7.0f, .p_bo5 = 0.0f, .p_bo6 = 0.0f, .p_be1 = 0.15f, .p_be2 = 1.0f, .r0_sigma = 2.20f, .r0_pi = 1.90f, .r0_pp = 0.0f, .De_sigma = 60.0f, .De_pi = 80.0f, .De_pp = 0.0f, .gamma = 1.800f, .r_vdw = 2.200f, .D_vdw = 0.060f};

            // 16 = Sulfur
            constants::reaxParams[16] = {
                .p_boc1 = 1.8f, .p_boc2 = 7.0f, .p_boc3 = 7.0f, .p_boc4 = 8.0f, .p_boc5 = 7.0f, .p_bo1 = -0.15f, .p_bo2 = 6.5f, .p_bo3 = -0.15f, .p_bo4 = 7.0f, .p_bo5 = 0.0f, .p_bo6 = 0.0f, .r0_sigma = 2.05f, .r0_pi = 1.80f, .r0_pp = 0.0f, .De_sigma = 65.0f, .De_pi = 90.0f, .De_pp = 0.0f, .gamma = 2.000f, .r_vdw = 2.050f, .D_vdw = 0.090f};

            // 17 = Chlorine
            constants::reaxParams[17] = {
                .p_boc1 = 2.5f, .p_boc2 = 7.5f, .p_boc3 = 0.0f, .p_boc4 = 0.0f, .p_boc5 = 0.0f, .p_bo1 = -0.14f, .p_bo2 = 7.0f, .p_bo3 = 0.0f, .p_bo4 = 0.0f, .p_bo5 = 0.0f, .p_bo6 = 0.0f, .r0_sigma = 1.80f, .r0_pi = 0.0f, .r0_pp = 0.0f, .De_sigma = 60.0f, .De_pi = 0.0f, .De_pp = 0.0f, .gamma = 2.300f, .r_vdw = 1.980f, .D_vdw = 0.100f};

            // 26 = Iron
            constants::reaxParams[26] = {
                .p_boc1 = 1.0f,
                .p_boc2 = 8.0f,
                .p_boc3 = 8.0f,
                .p_boc4 = 9.0f,
                .p_boc5 = 8.0f,
                .p_bo1 = -0.10f,
                .p_bo2 = 7.0f,
                .p_bo3 = 0.0f,
                .p_bo4 = 0.0f,
                .p_bo5 = 0.0f,
                .p_bo6 = 0.0f,
                .p_be1 = 0.0f,
                .p_be2 = 1.0f,
                .r0_sigma = 2.50f,
                .r0_pi = 0.0f,
                .r0_pp = 0.0f,
                .De_sigma = 100.0f,
                .De_pi = 0.0f,
                .De_pp = 0.0f,
                .gamma = 1.40f,
                .r_vdw = 2.500f,
                .D_vdw = 0.200f,
            };
        }
        // loading and saving scenes

        void universe::saveScene(const std::filesystem::path path)
        {
            nlohmann::json scene{};
        }

        void universe::loadScene(const std::filesystem::path path)
        {
            nlohmann::json nScene{};

            std::ifstream file(path);
            if (!file.is_open() || path.extension() != ".json")
            {
                std::cerr << "[Simulation] Cannot open: " << path << '\n';
                return;
            }

            file >> nScene;
        }
    } // namespace fun
} // namespace sim
