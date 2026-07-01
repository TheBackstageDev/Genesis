#include "parameter_table.hpp"

#include <fstream>
#include <iostream>
#include "simulation/constants.hpp"

namespace sim
{
    void parameter_table::save(const std::filesystem::path& file)
    {
    
    }

    void parameter_table::loadLJ(const std::filesystem::path& file)
    {
        std::ifstream in(file);

        if (!in)
        {
            std::cerr << "parameters file unable to be loaded! " << file << "\n";  
            return;
        }

        std::string line;
        std::getline(in, line);

        if (line != "Z,element,sigma,epsilon")
        {
            in.close();
            std::cerr << "parameters file doesn't have Z, element, sigma, epsilon. Unable to be loaded! " << file << "\n";
            return;
        }

        m_universalljparams.clear();
        m_universalljparams.resize(119); // 118 atoms
        
        while (std::getline(in, line)) 
        {
            std::stringstream ss(line);
            std::string token;

            int32_t Z;
            std::string element;
            float sigma, epsilon;

            std::getline(ss, token, ','); Z = std::stoi(token);
            std::getline(ss, element, ',');
            std::getline(ss, token, ','); sigma = std::stof(token);
            std::getline(ss, token, ','); epsilon = std::stof(token);

            m_universalljparams[Z] = LJParams{sigma, epsilon};
        }

        m_ljparams = m_universalljparams;
    }

    void parameter_table::loadMorse(const std::filesystem::path& file)
    {
        std::ifstream in(file);
        if (!in)
        {
            std::cerr << "Morse parameters file unable to be loaded! " << file << "\n";
            return;
        }

        std::string line;
        std::getline(in, line);

        if (line != "bond,D,alpha,r0")
        {
            in.close();
            std::cerr << "Morse file doesn't have bond,D,alpha,r0 header. Unable to be loaded! " << file << "\n";
            return;
        }

        m_universalmorseparams.clear();

        while (std::getline(in, line))
        {
            std::stringstream ss(line);
            std::string bondStr, token;
            float D, alpha, r0;

            std::getline(ss, bondStr, ',');
            std::getline(ss, token, ','); D = std::stof(token);
            std::getline(ss, token, ','); alpha = std::stof(token);
            std::getline(ss, token, ','); r0 = std::stof(token);

            char order = '-';
            if (bondStr.find('=') != std::string::npos) order = '=';
            else if (bondStr.find('#') != std::string::npos) order = '#';

            std::string atom1, atom2;
            size_t pos = bondStr.find_first_of("-=#");
            atom1 = bondStr.substr(0, pos);
            atom2 = bondStr.substr(pos+1, bondStr.size() - pos - 1);

            uint32_t Zi = static_cast<uint32_t>(constants::symbolToZ(atom1));
            uint32_t Zj = static_cast<uint32_t>(constants::symbolToZ(atom2));

            uint64_t hash = bondHash(Zi, Zj, order);
            m_universalmorseparams[hash] = MorseParams{D, r0, alpha, order};
        }

        m_morseparams = m_universalmorseparams;
    }
} // namespace sim
