#pragma once

#include "common.hpp"
#include "lennard_jones.hpp"
#include "morse.hpp"

#include <filesystem>

namespace sim
{
    struct parameter_table_create_info
    {
        std::filesystem::path universalLj{};
        std::filesystem::path universalMorse{};
    };

    class parameter_table 
    {
    public:
        parameter_table(const parameter_table_create_info info)
        {
            loadLJ(info.universalLj);
            loadMorse(info.universalMorse);
        }

        void save(const std::filesystem::path& path);
        void loadLJ(const std::filesystem::path& file);
        void loadMorse(const std::filesystem::path& file);

        const LJParams& lj(uint32_t atomType) const 
        { 
            if (m_ljparams.empty() && m_customParamsLJChosen) return LJParams(2.0, 0.5);
            if (m_universalljparams.empty() && !m_customParamsLJChosen) return LJParams(2.0, 0.5);
            
            return m_customParamsLJChosen ? m_ljparams.at(atomType) : m_universalljparams.at(atomType); 
        }
        
        const MorseParams& morse(uint32_t Zi, uint32_t Zj, char bondType) const 
        { 
            uint64_t hash = bondHash(Zi, Zj, bondType);
            return m_morseparams.at(hash); 
        }

        void editLJ(uint32_t atomType, const LJParams& newParams)
        {
            m_ljparams[atomType] = newParams;
        }

        void editMorse(uint32_t Zi, uint32_t Zj, const MorseParams& newParams)
        {
            uint64_t hash = bondHash(Zi, Zj, newParams.order);
            m_morseparams[hash] = newParams;
        }

        bool customLJ() { return m_customParamsLJChosen; }
        bool customMorse() { return m_customParamsMorseChosen; }

    private:
        uint64_t bondHash(uint32_t a, uint32_t b, char bondType) const
        {
            if (a > b) std::swap(a, b);

            uint64_t orderCode = 0;
            switch (bondType) {
                case '-': orderCode = 1; break;
                case '=': orderCode = 2; break;
                case '#': orderCode = 3; break;
                default:  orderCode = 0; break;
            }

            return (static_cast<uint64_t>(a) << 32) |
                (static_cast<uint64_t>(b) << 16) |
                orderCode;
        }


        bool m_customParamsLJChosen = false;
        bool m_customParamsMorseChosen = false;

        std::vector<LJParams> m_universalljparams;
        std::unordered_map<uint64_t, MorseParams> m_universalmorseparams;

        // Custom
        std::vector<LJParams> m_ljparams;
        std::unordered_map<uint64_t, MorseParams> m_morseparams;
        //std::vector<std::vector<LJParams>> m_ljpair;
    };
} // namespace sim
