#include "reaxx.hpp"

namespace sim
{
    std::pair<std::vector<reaxxAtomParams>, std::vector<reaxxBondParams>> getReaxxParams(std::filesystem::path path)
    {
        std::ifstream file(path);
        if (path.extension() != ".ff" || !file) throw std::runtime_error("Cannot open ffield file");

        std::string line;
        int natoms, nbonds;

        std::getline(file, line);
        std::istringstream(line) >> natoms;
        std::vector<reaxxAtomParams> atoms(natoms);

        for (int i = 0; i < natoms; i++) {
            std::getline(file, line);
            std::istringstream iss(line);
            iss >> atoms[i].radius >> atoms[i].valence
                >> atoms[i].chi >> atoms[i].J
                >> atoms[i].vdwA >> atoms[i].vdwB;
        }

        std::getline(file, line);
        std::istringstream(line) >> nbonds;
        std::vector<reaxxBondParams> bonds(nbonds);

        for (int i = 0; i < nbonds; i++) {
            std::getline(file, line);
            std::istringstream iss(line);
            iss >> bonds[i].r0 >> bonds[i].De
                >> bonds[i].pbo1 >> bonds[i].pbo2
                >> bonds[i].pbo3 >> bonds[i].pbo4
                >> bonds[i].pbe1 >> bonds[i].pbe2;
        }

        return {atoms, bonds};
    }
} // namespace sim
