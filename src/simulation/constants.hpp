#pragma once

#include <vector>
#include <map>
#include <string>

namespace sim
{
    namespace fun { enum class BondType { NONE, SINGLE, DOUBLE, TRIPLE, QUADRUPLE }; };
}; // namespace sim


namespace constants
{
#define M_PI 3.1415926535
#define RADIAN M_PI / 180
#define DEGREE 180 / M_PI

#define MASS_PROTON 1.0007     // Daltons
#define MASS_NEUTRON 1.0008    // Daltons
#define MASS_ELECTRON 1 / 1337 // Daltons

#define EPSILON 0.1f
#define DT 0.001f // ps
#define MULT_FACTOR 1.f
#define ANGSTROM 1e20f
#define PICOSECOND 1e24f

#define VERLET_SKIN 1.f
#define CUTOFF 2.5f
#define COULOMB_CUTOFF 8.f * MULT_FACTOR
#define REACTION_CUTOFF 3.0f * MULT_FACTOR

#define AVOGADRO 6.02214076e26f                                   // conversion from Daltons to Kg
#define BOLTZMAN_CONSTANT 1.380649e-23f                           // Boltzman Constant m^2 kg s^-2 K^-1
#define KB (BOLTZMAN_CONSTANT * AVOGADRO * ANGSTROM) / PICOSECOND // A^2 D ps^-2 K^-1

#define REBUILD_THRESHOLD 2.f * MULT_FACTOR
#define THERMOSTAT_INTERVAL 2
#define REACTION_INTERVAL 50

#define COULOMB_K 1389.3546f // kJ·mol⁻¹· Å ·e⁻²
#define BOND_K 50000.f        // Harmonic force constant
#define ANGLE_K 5000.f       // J/mol/rad² for angular potential
#define BOND_LENGTH_FACTOR 1.1f

    inline sf::Color getElementColor(uint8_t Z)
    {
        switch (Z)
        {
            case  1: return sf::Color(233, 233, 233);    // H  - Hydrogen     (White)
            case  2: return sf::Color(217, 255, 255);    // He - Helium       (Light Cyan)
            case  3: return sf::Color(204, 128, 255);    // Li - Lithium      (Violet)
            case  4: return sf::Color(178, 255, 255);    // Be - Beryllium    (Light Green)
            case  5: return sf::Color(128, 128, 128);    // B  - Boron        (Dark Gray)
            case  6: return sf::Color(144, 144, 144);    // C  - Carbon       (Gray)
            case  7: return sf::Color( 48,  80, 248);    // N  - Nitrogen     (Blue)
            case  8: return sf::Color(255,  13,  13);    // O  - Oxygen       (Red)
            case  9: return sf::Color(  0, 255, 255);    // F  - Fluorine     (Cyan)
            case 10: return sf::Color(179, 227, 245);    // Ne - Neon         (Pale Blue)
            case 11: return sf::Color( 42,  82, 190);    // Na - Sodium       (Blue)
            case 12: return sf::Color(  0, 130,  70);    // Mg - Magnesium    (Dark Green)
            case 13: return sf::Color(194, 194, 194);    // Al - Aluminum     (Silver)
            case 14: return sf::Color( 99,  99,  99);    // Si - Silicon      (Gray)
            case 15: return sf::Color(255, 165,   0);    // P  - Phosphorus   (Orange)
            case 16: return sf::Color(255, 255,   0);    // S  - Sulfur       (Yellow)
            case 17: return sf::Color( 50, 205,  50);    // Cl - Chlorine     (Green)
            case 18: return sf::Color(128, 128, 128);    // Ar - Argon        (Dark Gray)
            case 19: return sf::Color(160,  82,  45);    // K  - Potassium    (Brown)
            case 20: return sf::Color(135, 206, 235);    // Ca - Calcium      (Sky Blue)
            case 21: return sf::Color(138, 138, 138);    // Sc - Scandium     (Gray)
            case 22: return sf::Color(  0, 100,   0);    // Ti - Titanium     (Dark Green)
            case 23: return sf::Color(148,   0, 211);    // V  - Vanadium     (Purple)
            case 24: return sf::Color( 70, 130, 180);    // Cr - Chromium     (Steel Blue)
            case 25: return sf::Color( 25,  25, 112);    // Mn - Manganese    (Midnight Blue)
            case 26: return sf::Color(194,  78,  84);    // Fe - Iron         (Rust Red)
            default: return sf::Color(200, 200, 200);    // Unknown
        }
    }

    inline uint8_t symbolToZ(const std::string &sym)
    {
        static const std::map<std::string, uint8_t> elements = {
            {"H", 1},
            {"He", 2},
            {"Li", 3},
            {"Be", 4},
            {"B", 5},
            {"C", 6},
            {"N", 7},
            {"O", 8},
            {"F", 9},
            {"Ne", 10},
            {"Na", 11},
            {"Mg", 12},
            {"Al", 13},
            {"Si", 14},
            {"P", 15},
            {"S", 16},
            {"Cl", 17},
            {"Ar", 18},
            {"K", 19},
            {"Ca", 20},
            {"Sc", 21},
            {"Ti", 22},
            {"V", 23},
            {"Cr", 24},
            {"Mn", 25},
            {"Fe", 26},
        };
        auto it = elements.find(sym);
        return it != elements.end() ? it->second : 0; 
    }

    inline uint8_t getValenceElectrons(uint8_t ZIndex)
    {
        switch (ZIndex)
        {
        case 1:
            return 1; // H
        case 2:
            return 2; // He
        case 3:
            return 1; // Li
        case 4:
            return 2; // Be
        case 5:
            return 3; // B
        case 6:
            return 4; // C
        case 7:
            return 5; // N
        case 8:
            return 6; // O
        case 9:
            return 7; // F
        case 10:
            return 8; // Ne
        case 11:
            return 1; // Na
        case 12:
            return 2; // Mg
        case 13:
            return 3; // Al
        case 14:
            return 4; // Si
        case 15:
            return 5; // P
        case 16:
            return 6; // S
        case 17:
            return 7; // Cl
        case 18:
            return 8; // Ar
        case 19:
            return 1; // K
        case 20:
            return 2; // Ca
        case 21:
            return 3; // Sc
        case 22:
            return 4; // Ti
        case 23:
            return 5; // V
        case 24:
            return 6; // Cr
        case 25:
            return 7; // Mn
        case 26:
            return 8; // Fe
        default:
            return 0;
        }
    }

    inline uint8_t getUsualBonds(uint8_t ZIndex)
    {
        uint8_t valenceShell = getValenceElectrons(ZIndex);

        // special cases
        if (ZIndex == 1)
            return 1;
        if (valenceShell == 8 || ZIndex == 2)
            return 0; // noble gases

        // octet rule
        uint8_t bonds = 8 - valenceShell;

        // exceptions
        switch (ZIndex)
        {
        case 4:
            return 2; // Be
        case 5:
            return 3; // B
        case 11:
            return 1; // Na
        case 15:
            return 3; // P
        case 16:
            return 2; // S
        case 21:
            return 3; // Sc
        case 22:
            return 4; // Ti
        case 23:
            return 5; // V
        case 24:
            return 3; // Cr
        case 25:
            return 2; // Mn
        case 26:
            return 2; // Fe
        }

        return bonds;
    }

    inline float getAngles(uint8_t centralZIndex, const std::vector<uint8_t> &neighborZs, const std::vector<sim::fun::BondType> &types)
    {
        size_t bond_count = neighborZs.size();
        if (bond_count < 2)
            return 0.0f; // No angle if fewer than 2 neighbors

        uint32_t totalBondOrder = 0;
        for (const auto& t : types)
        {
            totalBondOrder += static_cast<uint32_t>(t);
        }

        const uint32_t valence = getValenceElectrons(centralZIndex);

        const uint32_t electronsInBonds = totalBondOrder * 2;
        int32_t remaining = valence - electronsInBonds;
        if (remaining < 0) remaining = 0;

        const uint32_t totalPairs = bond_count + remaining;

        float ideal_angle = 0.0f;
        switch (totalPairs)
        {
        case 2:
            ideal_angle = M_PI;
            break; // Linear
        case 3:
            ideal_angle = 120.0f * RADIAN;
            break; // Trigonal planar
        case 4:
            ideal_angle = 109.5f * RADIAN; // Tetrahedral base
            if (remaining == 1 && bond_count == 3)
                ideal_angle = 107.0f * RADIAN;
            if (remaining == 2 && bond_count == 2)
                ideal_angle = 104.5f * RADIAN;
            break;
        case 5:
            ideal_angle = 120.0f * RADIAN;
            break; // Trigonal bipyramidal (equatorial)
        case 6:
            ideal_angle = 90.0f * RADIAN;
            break; // Octahedral
        }

        return ideal_angle;
    }

    inline std::pair<float, float> getAtomConstants(uint32_t ZIndex)
    {
        std::pair<float, float> constants; // {sigma (Å), epsilon (kJ/mol)}
        switch (ZIndex)
        {
        case 1:
            constants = {2.50f, 0.1255f};
            break; // H
        case 2:
            constants = {2.58f, 0.0870f};
            break; // He
        case 3:
            constants = {1.82f, 0.1100f};
            break; // Li
        case 4:
            constants = {2.75f, 0.2510f};
            break; // Be
        case 5:
            constants = {3.64f, 0.3347f};
            break; // B
        case 6:
            constants = {3.55f, 0.2929f};
            break; // C
        case 7:
            constants = {3.25f, 0.7113f};
            break; // N
        case 8:
            constants = {3.03f, 0.7113f};
            break; // O
        case 9:
            constants = {2.95f, 0.6276f};
            break; // F
        case 10:
            constants = {2.78f, 0.1632f};
            break; // Ne
        case 11:
            constants = {2.43f, 0.4184f};
            break; // Na
        case 12:
            constants = {3.43f, 0.4184f};
            break; // Mg
        case 13:
            constants = {4.01f, 0.5021f};
            break; // Al
        case 14:
            constants = {3.83f, 0.8368f};
            break; // Si
        case 15:
            constants = {3.74f, 0.8368f};
            break; // P
        case 16:
            constants = {3.56f, 1.0460f};
            break; // S
        case 17:
            constants = {3.47f, 1.2552f};
            break; // Cl
        case 18:
            constants = {3.30f, 0.4184f};
            break; // Ar
        case 19:
            constants = {2.75f, 0.4184f};
            break; // K
        case 20:
            constants = {3.50f, 0.5021f};
            break; // Ca
        case 21:
            constants = {3.80f, 0.7113f};
            break; // Sc
        case 22:
            constants = {3.75f, 0.8368f};
            break; // Ti
        case 23:
            constants = {3.70f, 0.9205f};
            break; // V
        case 24:
            constants = {3.65f, 1.0460f};
            break; // Cr
        case 25:
            constants = {3.60f, 1.1715f};
            break; // Mn
        case 26:
            constants = {3.55f, 1.2970f};
            break; // Fe
        default:
            constants = {2.50f, 0.1255f};
            break; // H fallback
        }

        constants.first *= MULT_FACTOR;
        constants.second *= MULT_FACTOR;
        return constants;
    }

    inline const std::map<uint8_t, float> electronegativity =
        {
            {1, 2.20f},  // H
            {2, 0.00f},  // He
            {3, 0.98f},  // Li
            {4, 1.57f},  // Be
            {5, 2.04f},  // B
            {6, 2.55f},  // C
            {7, 3.04f},  // N
            {8, 3.44f},  // O
            {9, 3.98f},  // F
            {10, 0.00f}, // Ne
            {11, 0.93f}, // Na
            {12, 1.31f}, // Mg
            {13, 1.61f}, // Al
            {14, 1.90f}, // Si
            {15, 2.19f}, // P
            {16, 2.58f}, // S
            {17, 3.16f}, // Cl
            {18, 0.00f}, // Ar
            {19, 0.82f}, // K
            {20, 1.00f}, // Ca
            {21, 1.36f}, // Sc
            {22, 1.54f}, // Ti
            {23, 1.63f}, // V
            {24, 1.66f}, // Cr
            {25, 1.55f}, // Mn
            {26, 1.83f}  // Fe
    };

    inline std::string getAtomLetter(uint32_t ZIndex)
    {
        switch (ZIndex)
        {
        case 1:
            return "H"; // Hydrogen
        case 2:
            return "He"; // Helium
        case 3:
            return "Li"; // Lithium
        case 4:
            return "Be"; // Beryllium
        case 5:
            return "B"; // Boron
        case 6:
            return "C"; // Carbon
        case 7:
            return "N"; // Nitrogen
        case 8:
            return "O"; // Oxygen
        case 9:
            return "F"; // Fluorine
        case 10:
            return "Ne"; // Neon
        case 11:
            return "Na"; // Sodium
        case 12:
            return "Mg"; // Magnesium
        case 13:
            return "Al"; // Aluminum
        case 14:
            return "Si"; // Silicon
        case 15:
            return "P"; // Phosphorus
        case 16:
            return "S"; // Sulfur
        case 17:
            return "Cl"; // Chlorine
        case 18:
            return "Ar"; // Argon
        case 19:
            return "K"; // Potassium
        case 20:
            return "Ca"; // Calcium
        case 21:
            return "Sc"; // Scandium
        case 22:
            return "Ti"; // Titanium
        case 23:
            return "V"; // Vanadium
        case 24:
            return "Cr"; // Chromium
        case 25:
            return "Mn"; // Manganese
        case 26:
            return "Fe"; // Iron
        default:
            return "H";
        }
    }

    inline float getBondLength(uint8_t ZIndex1, uint8_t ZIndex2, sim::fun::BondType type)
    {
        // Fallback: average covalent radii (scaled)
        float base = (getAtomConstants(ZIndex1).first / MULT_FACTOR +
                    getAtomConstants(ZIndex2).first / MULT_FACTOR) / 2.3f * BOND_LENGTH_FACTOR;
        float baseBase = base;

        if (ZIndex1 > ZIndex2) std::swap(ZIndex1, ZIndex2);
        else std::swap(ZIndex2, ZIndex1);

        switch (type)
        {
        case sim::fun::BondType::SINGLE:
            // --- Carbon bonds ---
            if (ZIndex1 == 6 && ZIndex2 == 6)  base = 1.54f; // C–C (sp³)
            if (ZIndex1 == 6 && ZIndex2 == 1)  base = 1.09f; // C–H
            if (ZIndex1 == 6 && ZIndex2 == 7)  base = 1.47f; // C–N
            if (ZIndex1 == 6 && ZIndex2 == 8)  base = 1.43f; // C–O
            if (ZIndex1 == 6 && ZIndex2 == 9)  base = 1.35f; // C–F
            if (ZIndex1 == 6 && ZIndex2 == 15) base = 1.85f; // C–P
            if (ZIndex1 == 6 && ZIndex2 == 16) base = 1.82f; // C–S
            if (ZIndex1 == 6 && ZIndex2 == 17) base = 1.77f; // C–Cl
            if (ZIndex1 == 6 && ZIndex2 == 35) base = 1.94f; // C–Br
            if (ZIndex1 == 6 && ZIndex2 == 53) base = 2.14f; // C–I

            // --- Nitrogen bonds ---
            if (ZIndex1 == 7 && ZIndex2 == 7)  base = 1.45f; // N–N
            if (ZIndex1 == 7 && ZIndex2 == 8)  base = 1.46f; // N–O
            if (ZIndex1 == 7 && ZIndex2 == 1)  base = 1.01f; // N–H

            // --- Oxygen bonds ---
            if (ZIndex1 == 8 && ZIndex2 == 8)  base = 1.48f; // O–O
            if (ZIndex1 == 8 && ZIndex2 == 1)  base = 0.96f; // O–H
            if (ZIndex1 == 8 && ZIndex2 == 16) base = 1.57f; // O–S

            // --- Halogens ---
            if (ZIndex1 == 9 && ZIndex2 == 1)   base = 0.92f; // F–H
            if (ZIndex1 == 17 && ZIndex2 == 1)  base = 1.27f; // Cl–H
            if (ZIndex1 == 35 && ZIndex2 == 1)  base = 1.41f; // Br–H
            if (ZIndex1 == 53 && ZIndex2 == 1)  base = 1.61f; // I–H

            // --- Phosphorus ---
            if (ZIndex1 == 15 && ZIndex2 == 15) base = 2.22f; // P–P
            if (ZIndex1 == 15 && ZIndex2 == 8)  base = 1.60f; // P–O
            if (ZIndex1 == 15 && ZIndex2 == 1)  base = 1.42f; // P–H

            // --- Sulfur ---
            if (ZIndex1 == 16 && ZIndex2 == 16) base = 2.05f; // S–S
            if (ZIndex1 == 16 && ZIndex2 == 1)  base = 1.34f; // S–H

            // --- Metals (coordination) ---
            if (ZIndex1 == 26 && ZIndex2 == 6)  base = 1.92f; // Fe–C
            if (ZIndex1 == 26 && ZIndex2 == 8)  base = 1.62f; // Fe–O
            if (ZIndex1 == 28 && ZIndex2 == 6)  base = 1.88f; // Ni–C
            if (ZIndex1 == 29 && ZIndex2 == 6)  base = 1.90f; // Cu–C
            if (ZIndex1 == 30 && ZIndex2 == 6)  base = 1.95f; // Zn–C

            // --- Aromatic (special case) ---
            if ((ZIndex1 == 6 && ZIndex2 == 6) && 
                (getAtomConstants(ZIndex1).second || getAtomConstants(ZIndex2).second))
                base = 1.39f; // C–C aromatic (benzene)

            break;

        case sim::fun::BondType::DOUBLE:
            if (ZIndex1 == 6 && ZIndex2 == 6)  base = 1.34f; // C=C
            if (ZIndex1 == 6 && ZIndex2 == 7)  base = 1.27f; // C=N
            if (ZIndex1 == 6 && ZIndex2 == 8)  base = 1.22f; // C=O
            if (ZIndex1 == 6 && ZIndex2 == 16) base = 1.61f; // C=S
            if (ZIndex1 == 7 && ZIndex2 == 7)  base = 1.24f; // N=N
            if (ZIndex1 == 7 && ZIndex2 == 8)  base = 1.28f; // N=O
            if (ZIndex1 == 8 && ZIndex2 == 8)  base = 1.21f; // O=O
            if (ZIndex1 == 8 && ZIndex2 == 16) base = 1.43f; // S=O
            if (ZIndex1 == 15 && ZIndex2 == 8) base = 1.44f; // P=O
            break;

        case sim::fun::BondType::TRIPLE:
            if (ZIndex1 == 6 && ZIndex2 == 6)  base = 1.20f; // C≡C
            if (ZIndex1 == 6 && ZIndex2 == 7)  base = 1.16f; // C≡N
            if (ZIndex1 == 7 && ZIndex2 == 7)  base = 1.10f; // N≡N
            break;

        case sim::fun::BondType::QUADRUPLE:
            // Rare (e.g. Cr≡Cr), fallback
            base *= 0.55f;
            break;

        default:
            break;
        }

        if (!std::isnan(base) && base != baseBase)
            base *= MULT_FACTOR;

        return base;
    }
}; // namespace constants