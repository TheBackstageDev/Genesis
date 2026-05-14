#pragma once

#include <fstream>
#include <filesystem>
#include <vector>

namespace sim
{
    struct reaxxAtomParams 
    {
        double radius, valence, chi, J;
        double vdwA, vdwB;
    };

    struct reaxxBondParams 
    {
        double r0, De;
        double pbo1, pbo2, pbo3, pbo4;
        double pbe1, pbe2;
    };

    std::pair<std::vector<reaxxAtomParams>, std::vector<reaxxBondParams>> getReaxxParams(std::filesystem::path path);
} // namespace sim::io
