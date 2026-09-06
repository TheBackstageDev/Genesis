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

    void parameter_table::loadTersoff(const std::filesystem::path& file)
    {
        std::ifstream in(file);
        if (!in)
        {
            std::cerr << "Tersoff parameters file unable to be loaded! " << file << "\n";
            return;
        }

        std::string line;
        std::getline(in, line);

        if (line != "element1,element2,element3,m,gamma,lambda3,c,d,cos_theta0,n,beta,lambda2,B,R,D,lambda1,A")
        {
            in.close();
            std::cerr << "Tersoff file doesn't have correct header. Unable to be loaded! " << file << "\n";
            return;
        }

        m_universaltersoffParams.clear();
        m_universaltersoffParamsFlat.clear();

        while (std::getline(in, line)) 
        {
            auto pos = line.find('#');
            if (pos != std::string::npos) line = line.substr(0, pos);

            if (line.empty()) continue;

            std::stringstream ss(line);
            std::string token;
            TersoffParams tp;

            std::getline(ss, token, ','); // element1
            std::getline(ss, token, ','); // element2
            std::getline(ss, token, ','); // element3

            std::getline(ss, token, ','); tp.m       = std::stof(token);
            std::getline(ss, token, ','); tp.gamma   = std::stof(token);
            std::getline(ss, token, ','); tp.lambda3 = std::stof(token);
            std::getline(ss, token, ','); tp.c       = std::stof(token);
            std::getline(ss, token, ','); tp.d       = std::stof(token);
            std::getline(ss, token, ','); tp.h = std::stof(token);
            std::getline(ss, token, ','); tp.n       = std::stof(token);
            std::getline(ss, token, ','); tp.beta    = std::stof(token);
            std::getline(ss, token, ','); tp.lambda2 = std::stof(token);
            std::getline(ss, token, ','); tp.B       = std::stof(token);
            std::getline(ss, token, ','); tp.R       = std::stof(token);
            std::getline(ss, token, ','); tp.D       = std::stof(token);
            std::getline(ss, token, ','); tp.lambda1 = std::stof(token);
            std::getline(ss, token, ','); tp.A       = std::stof(token);

            m_universaltersoffParams.emplace_back(tp);

            m_universaltersoffParamsFlat.push_back(tp.m);
            m_universaltersoffParamsFlat.push_back(tp.gamma);
            m_universaltersoffParamsFlat.push_back(tp.lambda3);
            m_universaltersoffParamsFlat.push_back(tp.c);
            m_universaltersoffParamsFlat.push_back(tp.d);
            m_universaltersoffParamsFlat.push_back(tp.h);
            m_universaltersoffParamsFlat.push_back(tp.n);
            m_universaltersoffParamsFlat.push_back(tp.beta);
            m_universaltersoffParamsFlat.push_back(tp.lambda2);
            m_universaltersoffParamsFlat.push_back(tp.B);
            m_universaltersoffParamsFlat.push_back(tp.R);
            m_universaltersoffParamsFlat.push_back(tp.D);
            m_universaltersoffParamsFlat.push_back(tp.lambda1);
            m_universaltersoffParamsFlat.push_back(tp.A);
        }
    }
} // namespace sim
