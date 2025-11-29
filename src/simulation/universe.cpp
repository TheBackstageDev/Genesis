#include "universe.hpp"

#include <iostream>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <queue>

#include "constants.hpp"

namespace sim
{
    namespace fun
    {

        universe::universe(float universeSize, float cell_size, bool react)
            : boxSize(universeSize), cell_size(cell_size), react(react)
        {
        }

        universe::universe(const std::filesystem::path path)
        {
        }

        void universe::drawBox(core::window_t &window)
        {
            const std::array<sf::Vector3f, 8> corners =
                {{{0.f, 0.f, 0.f},
                  {boxSize, 0.f, 0.f},
                  {boxSize, boxSize, 0.f},
                  {0.f, boxSize, 0.f},
                  {0.f, 0.f, boxSize},
                  {boxSize, 0.f, boxSize},
                  {boxSize, boxSize, boxSize},
                  {0.f, boxSize, boxSize}}};

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

            if (!react)
                drawBonds(window);
            else
                drawReactiveBonds(window);

            for (size_t i = 0; i < atoms.size(); ++i)
            {
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

        void universe::drawChargeField(core::window_t &window)
        {
            if (!(cx > 0 && cy > 0 && cz > 0)) return;
            
            sf::RenderStates states;
            states.blendMode = sf::BlendAlpha;

            constexpr int32_t GRID = 50;
            float step = boxSize / GRID;

            const float isovalue_pos = 1.f;
            const float isovalue_neg = -1.f;
            const float thickness = 0.35f;

            struct Blob 
            {
                sf::Vector2f screenPos;
                float        depth;
                float        size;
                sf::Color    color;
            };

            std::vector<Blob> blobs;
            blobs.reserve(GRID*GRID*GRID);

            for (int iz = 0; iz < GRID; ++iz)
            for (int iy = 0; iy < GRID; ++iy)
            for (int ix = 0; ix < GRID; ++ix)
            {
                sf::Vector3f worldPos(
                    (ix + 0.5f) * step,
                    (iy + 0.5f) * step,
                    (iz + 0.5f) * step);

                sf::Vector2f screen = project(window, worldPos);
                if (screen.x < -500) continue;

                float density = 0.0f;

                for (size_t d = 0; d < 14; ++d)
                {
                    auto& cell = cells[getCellID(ix + offsets[d][0], iy + offsets[d][1], iz + offsets[d][2])];
                    for (size_t i = 0; i < cell.size(); ++i)
                    {
                        sf::Vector3f dr = minImageVec(worldPos - pos(cell[i]));
                        float r2 = dr.lengthSquared();
                        density += data.q[cell[i]] * std::exp(-r2 / (1.2f * 1.2f));
                    }
                }

                if (std::abs(density) < 0.1f) continue;

                float dist_to_pos = std::abs(density - isovalue_pos);
                float dist_to_neg = std::abs(density - isovalue_neg);
                float dist = std::min(dist_to_pos, dist_to_neg);

                if (dist > thickness) continue;

                sf::Color col;
                if (density <= -1.f)       col = sf::Color(220,  20,  60, 180);   // deep red
                else if (density <= -0.5f)  col = sf::Color(255,  80,  80, 180);   // red
                else if (density <= -0.2f) col = sf::Color(255, 147, 147, 150);   // mild red / light pink
                else if (density <  0.2f)  col = sf::Color(100, 240, 100, 200);   // green
                else if (density <  0.5f)  col = sf::Color(135, 206, 255, 160);   // light sky blue
                else if (density <  1.f)   col = sf::Color(70, 130, 255, 180);   // medium blue
                else                        col = sf::Color(30,  70, 200, 180);   // deep royal blue

                uint8_t alpha = static_cast<uint8_t>(267 * (1.0f - dist / thickness));
                col.a = alpha;

                float distance = (worldPos - cam.eye()).length();
                float size = 1.8f * step * GRID / distance;

                blobs.emplace_back(screen, distance, size, col);
            }

            std::sort(blobs.begin(), blobs.end(), [=](const Blob& a, const Blob& b)
            {
                return a.depth > b.depth;
            });

            sf::VertexArray billboards(sf::PrimitiveType::Triangles);
            billboards.resize(blobs.size() * 6);

            size_t idx = 0;
            for (const auto& b : blobs)
            {
                sf::Vector2f p = b.screenPos;
                float s = b.size;
                sf::Color c = b.color;

                billboards[idx + 0] = {p + sf::Vector2f(-s,-s), c};
                billboards[idx + 1] = {p + sf::Vector2f( s,-s), c};
                billboards[idx + 2] = {p + sf::Vector2f( s, s), c};

                billboards[idx + 3] = {p + sf::Vector2f(-s,-s), c};
                billboards[idx + 4] = {p + sf::Vector2f( s, s), c};
                billboards[idx + 5] = {p + sf::Vector2f(-s, s), c};

                idx += 6;
            }
            window.getWindow().draw(billboards, states);
        }

        void universe::drawBonds(core::window_t &window)
        {
            sf::Vector2f dimensions = window.getWindow().getView().getSize();

            for (size_t b = 0; b < bonds.size(); ++b)
            {
                const bond &bond = bonds[b];
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
            for (auto &pair : reactive_bonds)
            {
                reactive_bond &bond = pair.second;

                sf::Vector2f dimensions = window.getWindow().getView().getSize();

                const sf::Vector3f &pCentral = pos(bond.i); // pB
                const sf::Vector3f &pBonded = pos(bond.j);  // pA

                sf::Vector2f s1 = cam.project(pCentral, dimensions.x, dimensions.y); // center
                sf::Vector2f s2 = cam.project(pBonded, dimensions.x, dimensions.y);  // end

                if (s1.x <= -9999 || s2.x <= -9999)
                    continue;
                if ((pCentral - pBonded).length() > 5.f)
                    continue;
                if ((pCentral - cam.eye()).length() > 200.f)
                    continue;

                sf::Vector2f dir = (s2 - s1).normalized();
                float len = dir.length();
                if (len < 2.0f)
                    continue;
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
            constexpr float MIN_COS_ANGLE = 0.8f;
            constexpr float STRONG_ENERGY = -100.0f;
            constexpr float WEAK_ENERGY = -5.0f;

            sf::Vector3f p = pos(H);
            p.x -= boxSize * std::floor(p.x * (1.0f / boxSize));
            p.y -= boxSize * std::floor(p.y * (1.0f / boxSize));
            p.z -= boxSize * std::floor(p.z * (1.0f / boxSize));

            int32_t ix = static_cast<int32_t>(p.x * 1 / cell_size);
            int32_t iy = static_cast<int32_t>(p.y * 1 / cell_size);
            int32_t iz = static_cast<int32_t>(p.z * 1 / cell_size);

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
            if (cx == 0 || cell_size <= 0.f)
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

                        float x0 = ix * cell_size;
                        float y0 = iy * cell_size;
                        float z0 = iz * cell_size;
                        float x1 = x0 + cell_size;
                        float y1 = y0 + cell_size;
                        float z1 = z0 + cell_size;

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

            int8_t bondCount = static_cast<int8_t>(type);

            atoms[idx1].bondCount += bondCount;
            atoms[idx2].bondCount += bondCount;

            float EN1 = constants::getElectronegativity(atoms[idx1].ZIndex);
            float EN2 = constants::getElectronegativity(atoms[idx2].ZIndex);
            float deltaEN = std::abs(EN1 - EN2);

            if (deltaEN > 0.3f) // Significant electronegativity difference
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
            data.x[i] = std::fmod(data.x[i] + boxSize, boxSize);
            data.y[i] = std::fmod(data.y[i] + boxSize, boxSize);
            data.z[i] = std::fmod(data.z[i] + boxSize, boxSize);
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

            const float sigma = (sigma_i + sigma_j) / 2.0f;

            if (dr < sigma * CUTOFF && dr > EPSILON)
            {
                if (areBonded(i, j))
                    return sf::Vector3f{0.f, 0.f, 0.f};

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

        sf::Vector3f universe::ljGrad(size_t i)
        {
            sf::Vector3f gradient({0.f, 0.f, 0.f});

            for (size_t j = 0; j < atoms.size(); ++j)
            {
                if (j == i)
                    continue;

                gradient += ljForce(i, j);
            }

            return gradient;
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

            if (areBonded(i, j))
                return sf::Vector3f{0.f, 0.f, 0.f};

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

                /* if (delta_r > bond.equilibriumLength * 2.f)
                {
                    breakBond(idx1, idx2);
                    continue;
                } */

                float bond_k = constants::getBondHarmonicConstantFromEnergy(atoms[idx1].ZIndex, atoms[idx2].ZIndex, bond.type);
                float force_magnitude = bond_k * delta_r;

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
                float k_angle = 2.0f * ang.K;

                sf::Vector3f dtheta_dri = (cos_theta * u_ji - u_jk) / (len_ji * sin_theta);
                sf::Vector3f dtheta_drk = (cos_theta * u_jk - u_ji) / (len_jk * sin_theta);

                sf::Vector3f F_i = -k_angle * delta_theta * dtheta_dri;
                sf::Vector3f F_k = -k_angle * delta_theta * dtheta_drk;
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

                        for (size_t ii = 0; ii < cell.size(); ++ii)
                        {
                            size_t i = cell[ii];
                            for (size_t jj = ii + 1; jj < cell.size(); ++jj)
                            {
                                size_t j = cell[jj];

                                if (j <= i)
                                    continue;

                                sf::Vector3f f = ljForce(i, j);
                                add_force(i, f);
                                add_force(j, -f);
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
                                for (size_t jj = 0; jj < neighbour_cell.size(); ++jj)
                                {
                                    const size_t &j = neighbour_cell[jj];

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

            sf::Vector3f n1 = v1.cross(v2).normalized();
            sf::Vector3f n2 = v2.cross(v3).normalized();

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

                sf::Vector3f rA = (pa - pb.normalized());
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

        void universe::calcUnbondedForces()
        {
            calcLjForces();
            calcElectrostaticForces();
        }

        void universe::update(float targetTemperature)
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

            if (!react)
            {
                calcBondedForces();
                calcUnbondedForces();
            }
            else
                handleReactiveForces();

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector3f a_old = old_force(i) / atoms[i].mass;
                sf::Vector3f a_new = force(i) / atoms[i].mass;

                add_pos(i, vel(i) * DT + 0.5f * (a_old + a_new) * DT * DT);
                boundCheck(i);
            }

            std::fill(data.fx.begin(), data.fx.end(), 0.f);
            std::fill(data.fy.begin(), data.fy.end(), 0.f);
            std::fill(data.fz.begin(), data.fz.end(), 0.f);

            if (!react)
            {
                calcBondedForces();
                calcUnbondedForces();
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

        void universe::setTemperature(float kelvin)
        {
            if (timeStep % THERMOSTAT_INTERVAL == 0)
            {
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
            if ((CELL_CUTOFF - cell_size) > EPSILON || timeStep == 0)
            {
                cx = static_cast<size_t>(std::ceil(boxSize / cell_size));
                cy = cx;
                cz = cx;
            }

            size_t ncells = cx * cy * cz;
            if (cells.size() != ncells)
                cells.assign(ncells, {});
            else
                for (auto &c : cells)
                    c.clear();

            const float inv = 1.0f / cell_size;

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                sf::Vector3f p = pos(i);
                p.x -= boxSize * std::floor(p.x / boxSize);
                p.y -= boxSize * std::floor(p.y / boxSize);
                p.z -= boxSize * std::floor(p.z / boxSize);

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

        // Reactions

        float universe::calculateAtomTemperature(size_t i)
        {
            float ke = 0.5f * atoms[i].mass * vel(i).lengthSquared();
            return (2.0f / 3.0f) * ke * KB;
        }

        float universe::calculateBondEnergy(size_t i, size_t j, const BondType type)
        {
            float r = minImageVec(pos(i) - pos(j)).length();
            float r0 = constants::getBondLength(atoms[i].ZIndex, atoms[j].ZIndex, type);
            float k = constants::getBondHarmonicConstantFromEnergy(atoms[i].ZIndex, atoms[j].ZIndex, type);
            float dr = r - r0;
            return 0.5f * k * dr * dr;
        }

        float universe::calculateNonBondedEnergy(size_t i, size_t j)
        {
            float r_ij = (pos(i) - pos(j)).length();
            if (r_ij < EPSILON)
                return 1e6f;

            float coulombPotential = COULOMB_K * (data.q[i] * data.q[j]) / r_ij;
            float lennardjonesPotential = ljPot(i, j);

            return lennardjonesPotential + coulombPotential;
        }

        void universe::calcReactiveAngleForces()
        {
            for (size_t ix = 0; ix < cx; ++ix)
                for (size_t iy = 0; iy < cy; ++iy)
                    for (size_t iz = 0; iz < cz; ++cz)
                    {
                    }
        }

        void universe::assignCharges()
        {
            std::vector<float> new_q(atoms.size(), 0.0f);

            for (const auto &[key, rb] : reactive_bonds)
            {
                if (rb.strength < 0.3f)
                    continue;

                size_t i = rb.i, j = rb.j;
                float bo = rb.strength;

                float chi_i = constants::getElectronegativity(atoms[i].ZIndex);
                float chi_j = constants::getElectronegativity(atoms[j].ZIndex);
                float delta_chi = chi_j - chi_i;

                float ionic_char = 1.0f - std::exp(-0.25f * delta_chi * delta_chi);
                float dq = 0.9f * ionic_char * bo;

                if (chi_j > chi_i)
                {
                    new_q[i] += dq;
                    new_q[j] -= dq;
                }
                else
                {
                    new_q[i] -= dq;
                    new_q[j] += dq;
                }
            }

            for (size_t i = 0; i < atoms.size(); ++i)
            {
                data.q[i] = 0.7f * data.q[i] + 0.3f * new_q[i];
            }
        }

        void universe::calculateBondOrder()
        {
            reactive_bonds.clear();

            size_t N = atoms.size();

            std::vector<float> uncorrected_bo(atoms.size(), 0.0f);
            for (size_t ix = 0; ix < cx; ++ix)
                for (size_t iy = 0; iy < cy; ++iy)
                    for (size_t iz = 0; iz < cz; ++iz)
                    {
                        size_t currentCell = getCellID(ix, iy, iz);

                        for (int32_t d = 0; d < 14; ++d)
                        {
                            int32_t x = ix + offsets[d][0], y = iy + offsets[d][1], z = iz + offsets[d][2];
                            size_t neighbourCell = getCellID(x, y, z);

                            for (size_t ii = 0; ii < cells[currentCell].size(); ++ii)
                            {
                                size_t i = cells[currentCell][ii];
                                for (size_t jj = 0; jj < cells[neighbourCell].size(); ++jj)
                                {
                                    size_t j = cells[neighbourCell][jj];
                                    if (currentCell == neighbourCell && j <= i)
                                        continue;

                                    sf::Vector3f dr_vec = minImageVec(pos(j) - pos(i));
                                    float r = dr_vec.length();
                                    if (r < 0.3f || r > REACTION_CUT_BO)
                                        continue;

                                    uint8_t Z1 = atoms[i].ZIndex;
                                    uint8_t Z2 = atoms[j].ZIndex;

                                    float r0 = constants::covalent_radius[Z1] + constants::covalent_radius[Z2];
                                    if (r0 < 0.1f)
                                        continue;

                                    // Note: Change the exponents later for a table
                                    float bo_raw = 0.0f;
                                    float ratio = r / r0;

                                    float bo_sigma = std::exp(-1.5f * std::pow(ratio - 1.0f, 2.0f));
                                    float bo_pi = 0.f;

                                    if ((Z1 == 6 && Z2 == 6) || (Z1 == 6 && Z2 == 7) || (Z1 == 7 && Z2 == 7))
                                    {
                                        bo_pi = 0.7f * std::exp(-3.0f * std::pow(ratio - 0.9f, 2.0f));
                                    }

                                    float bo_total_uncorrected = bo_sigma + bo_pi;
                                    if (bo_total_uncorrected < 0.01f)
                                        continue;

                                    uint64_t key = (uint64_t(i) << 32) | j;
                                    reactive_bonds[key] = {i, j, bo_sigma, bo_pi, 0.f, 0.0f, BondType::SINGLE};

                                    uncorrected_bo[i] += bo_total_uncorrected;
                                    uncorrected_bo[j] += bo_total_uncorrected;
                                }
                            }
                        }
                    }

            std::vector<float> Delta(atoms.size(), 0.0f);
            for (size_t i = 0; i < atoms.size(); ++i)
            {
                float ideal = constants::getUsualBonds(atoms[i].ZIndex);
                Delta[i] = uncorrected_bo[i] - ideal;
            }

            // OverCoordination
            for (auto &[key, rb] : reactive_bonds)
            {
                size_t i = rb.i, j = rb.j;

                float f_c_i = (Delta[i] <= 0) ? 1.0f : (Delta[i] >= 2.0f) ? 0.0f
                                                                          : 0.5f * (1.0f + std::cos(M_PI * Delta[i]));

                float f_c_j = (Delta[j] <= 0) ? 1.0f : (Delta[j] >= 2.0f) ? 0.0f
                                                                          : 0.5f * (1.0f + std::cos(M_PI * Delta[j]));

                float total_bo_raw = rb.bo_sigma + rb.bo_pi + rb.bo_pipi;
                rb.strength = total_bo_raw * f_c_i * f_c_j;

                if (rb.strength > 1.7f)
                    rb.type = BondType::TRIPLE;
                else if (rb.strength > 1.2f)
                    rb.type = BondType::DOUBLE;
                else if (rb.strength > 0.7f)
                    rb.type = BondType::SINGLE;
                // else                         rb.type = BondType::PARTIAL;
            }
        }

        void universe::handleReactiveForces()
        {
            calculateBondOrder();
            if (timeStep % 8 == 0)
                assignCharges();

            for (const auto &[key, rb] : reactive_bonds)
            {
                if (rb.strength < 0.1f)
                    continue;

                size_t i = rb.i, j = rb.j;
                sf::Vector3f dr_vec = minImageVec(pos(j) - pos(i));
                float r = dr_vec.length();
                if (r < EPSILON)
                    continue;

                uint8_t Z1 = atoms[i].ZIndex;
                uint8_t Z2 = atoms[j].ZIndex;

                float r0_single = constants::getBondLength(Z1, Z2, rb.type);
                float De_single = constants::getBondEnergy(Z1, Z2, rb.type);
                float beta = 1.5f;

                float p = 2.0f;
                float r0 = r0_single * (1.0f + (1.0f - std::pow(rb.strength, p)) * 0.12f);
                float De = De_single * rb.strength;

                float exp_term = std::exp(-beta * (r - r0));
                float force_mag = -2.0f * De * beta * (1.5f - exp_term) * exp_term;

                sf::Vector3f dir = dr_vec / r;
                sf::Vector3f force = force_mag * dir;

                add_force(i, force);
                add_force(j, -force);
            }

            calcLjForces();
            calcElectrostaticForces();
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
                size_t centralAtom = subsets[subsetIdx[i]].mainAtomIdx;
                ++ZIndices[atoms[centralAtom].ZIndex];
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
                    name += prefixes[ZIndices[i]] + constants::getAtomName(i) + " ";
                }
            }
            else // organic nomeclature
            {
            }

            return name;
        }

        // loading and saving scenes

        void universe::saveScene(const std::filesystem::path path)
        {
            nlohmann::json scene{};
        }

        void universe::loadScene(const std::filesystem::path path)
        {
        }
    } // namespace fun
} // namespace sim
