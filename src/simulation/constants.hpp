#pragma once

#include <vector>
#include <map>
#include <string>
#include <SFML/Graphics.hpp>

namespace sim
{
    namespace fun { enum class BondType { NONE, SINGLE, DOUBLE, TRIPLE, QUADRUPLE, PARTIAL }; }; // Partial is for bonds that have partial pi bonds
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

#define JOULE_TO_CAL 1/4.164f

#define VERLET_SKIN 1.f
#define CUTOFF 2.5f
#define COULOMB_CUTOFF 12.f * MULT_FACTOR
#define REACTION_CUTOFF 2.f * MULT_FACTOR

#define CELL_CUTOFF COULOMB_CUTOFF + VERLET_SKIN

#define AVOGADRO 6.02214076e26f                                   // conversion from Daltons to Kg
#define BOLTZMAN_CONSTANT 1.380649e-23f                           // Boltzman Constant m^2 kg s^-2 K^-1
#define KB (BOLTZMAN_CONSTANT * AVOGADRO * ANGSTROM) / PICOSECOND // A^2 D ps^-2 K^-1

#define REBUILD_THRESHOLD 2.5f * MULT_FACTOR
#define THERMOSTAT_INTERVAL 2

#define COULOMB_K 1389.3546f // kJ·mol⁻¹· Å ·e⁻²
#define BOND_K 340000.f        // Harmonic force constant
#define ANGLE_K 12000.f       // J/mol/rad² for angular potential
#define BOND_LENGTH_FACTOR 1.f

#define REACTION_INTERVAL 3 // time steps
#define REACTION_STRETCH_FACTOR 1.5f
#define REACTION_FORMING_FACTOR 1.3f
#define REACTION_CUT_BO 2.0f
#define REACTION_VELOCITY_THRESHOLD 0.05f 

#define COUNT_ATOMS 118

    inline constexpr std::array<float, 119> covalent_radius = 
    {
        0.00f,  // Z=0 (dummy)
        0.31f,  // H
        0.28f,  // He
        1.28f,  // Li
        0.96f,  // Be
        0.84f,  // B
        0.76f,  // C
        0.71f,  // N
        0.66f,  // O
        0.57f,  // F
        0.58f,  // Ne
        1.66f,  // Na
        1.41f,  // Mg
        1.21f,  // Al
        1.11f,  // Si
        1.07f,  // P
        1.05f,  // S
        1.02f,  // Cl
        1.00f,  // Ar
        2.03f,  // K
        1.76f,  // Ca
        1.70f,  // Sc
        1.60f,  // Ti
        1.50f,  // V
        1.42f,  // Cr
        1.46f,  // Mn
        1.48f,  // Fe
        1.40f,  // Co
        1.37f,  // Ni
        1.39f,  // Cu
        1.25f,  // Zn
        1.46f,  // Ga
        1.21f,  // Ge
        1.07f,  // As
        1.05f,  // Se
        1.02f,  // Br
        1.00f,  // Kr
        2.20f,  // Rb
        1.95f,  // Sr
        1.90f,  // Y
        1.75f,  // Zr
        1.65f,  // Nb
        1.55f,  // Mo
        1.45f,  // Tc
        1.40f,  // Ru
        1.39f,  // Rh
        1.35f,  // Pd
        1.30f,  // Ag
        1.25f,  // Cd
        1.35f,  // In
        1.22f,  // Sn
        1.20f,  // Sb
        1.19f,  // Te
        1.20f,  // I
        1.20f,  // Xe
        2.44f,  // Cs
        2.15f,  // Ba
        2.05f,  // La
        2.00f,  // Hf (Lu=1.75, but grouped)
        1.90f,  // Ta
        1.85f,  // W
        1.80f,  // Re
        1.75f,  // Os
        1.75f,  // Ir
        1.75f,  // Pt
        1.55f,  // Au
        1.45f,  // Hg
        1.60f,  // Tl
        1.50f,  // Pb
        1.45f,  // Bi
        1.45f,  // Po
        1.40f,  // At
        1.40f,  // Rn
        2.60f,  // Fr
        2.25f,  // Ra
        0.00f,  // Ac

        1.80f,  // Ce 1.80, Pr 1.80, Nd 1.80, Pm 1.80, Sm 1.80, Eu 1.80, Gd 1.80,
        1.80f,  // Tb 1.75, Dy 1.75, Ho 1.75, Er 1.75, Tm 1.75, Yb 1.75, Lu 1.75

        1.70f,  // Th 1.75, Pa 1.70, U 1.65, Np 1.65, Pu 1.65, Am 1.65, Cm 1.65,
        1.65f,  // Bk 1.65, Cf 1.65, Es 1.65, Fm 1.65, Md 1.65, No 1.65

        1.60f,  // Rf 1.60, Db 1.60, Sg 1.60, Bh 1.60, Hs 1.60, Mt 1.50,
        1.50f,  // Ds 1.40, Rg 1.35, Cn 1.30, Nh 1.30, Fl 1.30, Mc 1.30, Lv 1.30,
        1.30f,  // Ts 1.30, Og 1.30
    };

    inline sf::Color getElementColor(uint8_t Z)
    {
        switch (Z)
        {
            case  1: return sf::Color(210, 210, 210);    // H  - Hydrogen     (White)
            case  2: return sf::Color(217, 255, 255);    // He - Helium       (Light Cyan)
            case  3: return sf::Color(204, 128, 255);    // Li - Lithium      (Violet)
            case  4: return sf::Color(178, 255, 255);    // Be - Beryllium    (Light Green)
            case  5: return sf::Color(128, 128, 128);    // B  - Boron        (Dark Gray)
            case  6: return sf::Color(144, 144, 144);    // C  - Carbon       (Gray)
            case  7: return sf::Color( 48,  80, 248);    // N  - Nitrogen     (Blue)
            case  8: return sf::Color(255,  13,  13);    // O  - Oxygen       (Red)
            case  9: return sf::Color(  0, 233, 233);    // F  - Fluorine     (Cyan)
            case 10: return sf::Color(179, 227, 245);    // Ne - Neon         (Pale Blue)
            case 11: return sf::Color( 42,  82, 190);    // Na - Sodium       (Blue)
            case 12: return sf::Color(  0, 130,  70);    // Mg - Magnesium    (Dark Green)
            case 13: return sf::Color(194, 194, 194);    // Al - Aluminum     (Silver)
            case 14: return sf::Color( 99,  99,  99);    // Si - Silicon      (Gray)
            case 15: return sf::Color(255, 165,   0);    // P  - Phosphorus   (Orange)
            case 16: return sf::Color(200, 200,   1);    // S  - Sulfur       (Yellow)
            case 17: return sf::Color( 50, 205,  50);    // Cl - Chlorine     (Green)
            case 18: return sf::Color(128, 128, 128);    // Ar - Argon        (Dark Gray)
            case 19: return sf::Color(160,  82,  45);    // K  - Potassium    (Brown)
            case 20: return sf::Color(135, 206, 235);    // Ca - Calcium      (Sky Blue)
            case 21: return sf::Color(108, 108, 100);    // Sc - Scandium     (Gray)
            case 22: return sf::Color(  0, 100,   0);    // Ti - Titanium     (Dark Green)
            case 23: return sf::Color(148,   0, 211);    // V  - Vanadium     (Purple)
            case 24: return sf::Color( 70, 130, 180);    // Cr - Chromium     (Steel Blue)
            case 25: return sf::Color( 25,  25, 112);    // Mn - Manganese    (Midnight Blue)
            case 26: return sf::Color(194,  78,  84);    // Fe - Iron         (Rust Red)
            case  36: return sf::Color(50, 200, 50); // Kr
            case  54: return sf::Color(  0, 120, 120); // Xe 
            case 79: return sf::Color(222, 222, 1); // Au
            case  86: return sf::Color( 80,  80,  80); // Rn (dark gray)
            case 92: return sf::Color(100, 200, 0); // U
            default:
            {
                if (Z >= 21 && Z <= 30)  return sf::Color(200, 150, 150); // transition metals reddish
                if (Z >= 31 && Z <= 48)  return sf::Color(150, 200, 150); // post-transition greenish
                if (Z >= 57 && Z <= 71)  return sf::Color(255, 200, 100); // lanthanides orange
                if (Z >= 89 && Z <= 103) return sf::Color(255, 150, 200); // actinides pink
                return sf::Color(200, 200, 200); // unknown
            }
        }
    }

    inline uint8_t symbolToZ(const std::string& sym)
    {
        static const std::map<std::string, uint8_t> table = {
            {"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},{"O",8},{"F",9},{"Ne",10},
            {"Na",11},{"Mg",12},{"Al",13},{"Si",14},{"P",15},{"S",16},{"Cl",17},{"Ar",18},{"K",19},{"Ca",20},
            {"Sc",21},{"Ti",22},{"V",23},{"Cr",24},{"Mn",25},{"Fe",26},{"Co",27},{"Ni",28},{"Cu",29},{"Zn",30},
            {"Ga",31},{"Ge",32},{"As",33},{"Se",34},{"Br",35},{"Kr",36},
            {"Rb",37},{"Sr",38},{"Y",39},{"Zr",40},{"Nb",41},{"Mo",42},{"Tc",43},{"Ru",44},{"Rh",45},{"Pd",46},{"Ag",47},{"Cd",48},
            {"In",49},{"Sn",50},{"Sb",51},{"Te",52},{"I",53},{"Xe",54},
            {"Cs",55},{"Ba",56},{"La",57},{"Ce",58},{"Pr",59},{"Nd",60},{"Pm",61},{"Sm",62},{"Eu",63},{"Gd",64},
            {"Tb",65},{"Dy",66},{"Ho",67},{"Er",68},{"Tm",69},{"Yb",70},{"Lu",71},
            {"Hf",72},{"Ta",73},{"W",74},{"Re",75},{"Os",76},{"Ir",77},{"Pt",78},{"Au",79},{"Hg",80},
            {"Tl",81},{"Pb",82},{"Bi",83},{"Po",84},{"At",85},{"Rn",86},
            {"Fr",87},{"Ra",88},{"Ac",89},{"Th",90},{"Pa",91},{"U",92},{"Np",93},{"Pu",94},{"Am",95},{"Cm",96},
            {"Bk",97},{"Cf",98},{"Es",99},{"Fm",100},{"Md",101},{"No",102},{"Lr",103},
            {"Rf",104},{"Db",105},{"Sg",106},{"Bh",107},{"Hs",108},{"Mt",109},{"Ds",110},{"Rg",111},{"Cn",112},
            {"Nh",113},{"Fl",114},{"Mc",115},{"Lv",116},{"Ts",117},{"Og",118}
        };
        auto it = table.find(sym);
        return (it != table.end()) ? it->second : 0;
    }

    inline uint8_t getValenceElectrons(uint8_t ZIndex)
    {
        switch (ZIndex)
        {
            // Period 1
            case   1: return 1; // H  : 1s¹
            case   2: return 2; // He : 1s²

            // Period 2
            case   3: return 1; // Li : [He] 2s¹
            case   4: return 2; // Be : [He] 2s²
            case   5: return 3; // B  : [He] 2s²2p¹
            case   6: return 4; // C  : [He] 2s²2p²
            case   7: return 5; // N  : [He] 2s²2p³
            case   8: return 6; // O  : [He] 2s²2p⁴
            case   9: return 7; // F  : [He] 2s²2p⁵
            case  10: return 8; // Ne : [He] 2s²2p⁶

            // Period 3
            case  11: return 1; // Na : [Ne] 3s¹
            case  12: return 2; // Mg : [Ne] 3s²
            case  13: return 3; // Al : [Ne] 3s²3p¹
            case  14: return 4; // Si : [Ne] 3s²3p²
            case  15: return 5; // P  : [Ne] 3s²3p³
            case  16: return 6; // S  : [Ne] 3s²3p⁴
            case  17: return 7; // Cl : [Ne] 3s²3p⁵
            case  18: return 8; // Ar : [Ne] 3s²3p⁶

            // Period 4
            case  19: return 1; // K  : [Ar] 4s¹
            case  20: return 2; // Ca : [Ar] 4s²
            case  21: return 3; // Sc : [Ar] 4s²3d¹ → valence = 3
            case  22: return 4; // Ti : [Ar] 4s²3d² → valence = 4
            case  23: return 5; // V  : [Ar] 4s²3d³ → valence = 5
            case  24: return 6; // Cr : [Ar] 4s¹3d⁵ → valence = 6
            case  25: return 7; // Mn : [Ar] 4s²3d⁵ → valence = 7
            case  26: return 8; // Fe : [Ar] 4s²3d⁶ → valence = 8
            case  27: return 9; // Co : [Ar] 4s²3d⁷ → valence = 9
            case  28: return 10;// Ni : [Ar] 4s²3d⁸ → valence = 10
            case  29: return 11;// Cu : [Ar] 4s¹3d¹⁰ → valence = 11
            case  30: return 12;// Zn : [Ar] 4s²3d¹⁰ → valence = 12
            case  31: return 3; // Ga : [Ar] 4s²3d¹⁰4p¹ → valence = 3
            case  32: return 4; // Ge : ...4p² → 4
            case  33: return 5; // As : ...4p³ → 5
            case  34: return 6; // Se : ...4p⁴ → 6
            case  35: return 7; // Br : ...4p⁵ → 7
            case  36: return 8; // Kr : ...4p⁶ → 8

            // Period 5
            case  37: return 1; // Rb : 5s¹
            case  38: return 2; // Sr : 5s²
            case  39: return 3; // Y  : 5s²4d¹
            case  40: return 4; // Zr : 5s²4d²
            case  41: return 5; // Nb : 5s¹4d⁴
            case  42: return 6; // Mo : 5s¹4d⁵
            case  43: return 7; // Tc : 5s²4d⁵
            case  44: return 8; // Ru : 5s¹4d⁷
            case  45: return 9; // Rh : 5s¹4d⁸
            case  46: return 10;// Pd : 4d¹⁰ (5s⁰)
            case  47: return 11;// Ag : 5s¹4d¹⁰
            case  48: return 12;// Cd : 5s²4d¹⁰
            case  49: return 3; // In : 5s²5p¹
            case  50: return 4; // Sn : 5s²5p²
            case  51: return 5; // Sb : 5s²5p³
            case  52: return 6; // Te : 5s²5p⁴
            case  53: return 7; // I  : 5s²5p⁵
            case  54: return 8; // Xe : 5s²5p⁶

            // Period 6 (includes f-block)
            case  55: return 1; // Cs : 6s¹
            case  56: return 2; // Ba : 6s²
            case  57: return 3; // La : 6s²5d¹
            case  72: return 4; // Hf : 6s²5d²
            case  73: return 5; // Ta : 6s²5d³
            case  74: return 6; // W  : 6s²5d⁴
            case  75: return 7; // Re : 6s²5d⁵
            case  76: return 8; // Os : 6s²5d⁶
            case  77: return 9; // Ir : 6s²5d⁷
            case  78: return 10;// Pt : 6s¹5d⁹
            case  79: return 11;// Au : 6s¹5d¹⁰
            case  80: return 12;// Hg : 6s²5d¹⁰
            case  81: return 3; // Tl : 6s²6p¹
            case  82: return 4; // Pb : 6s²6p²
            case  83: return 5; // Bi : 6s²6p³
            case  84: return 6; // Po : 6s²6p⁴
            case  85: return 7; // At : 6s²6p⁵
            case  86: return 8; // Rn : 6s²6p⁶

            // Period 7
            case  87: return 1; // Fr : 7s¹
            case  88: return 2; // Ra : 7s²
            case  89: return 3; // Ac : 7s²6d¹
            case 104: return 4; // Rf : ~7s²6d²
            case 105: return 5; // Db
            case 106: return 6; // Sg
            case 107: return 7; // Bh
            case 108: return 8; // Hs
            case 109: return 9; // Mt
            case 110: return 10;// Ds
            case 111: return 11;// Rg
            case 112: return 12;// Cn
            case 113: return 3; // Nh : ~7s²7p¹
            case 114: return 4; // Fl : ~7s²7p²
            case 115: return 5; // Mc
            case 116: return 6; // Lv
            case 117: return 7; // Ts
            case 118: return 8; // Og

            default:
                if ((ZIndex >= 58 && ZIndex <= 71) || (ZIndex >= 90 && ZIndex <= 103))
                    return 3;
                return 0;
        }
    }

    inline uint8_t getUsualBonds(uint8_t ZIndex)
    {
        if (ZIndex == 1)  return 1; // H
        if (ZIndex == 2 || ZIndex == 10 || ZIndex == 18 || 
            ZIndex == 36 || ZIndex == 54 || ZIndex == 86 || ZIndex == 118) 
            return 0; // Noble gases

        uint8_t ve = getValenceElectrons(ZIndex);
        if (ve == 0) return 0;

        uint8_t bonds = 8 - ve;

        switch (ZIndex)
        {
            case   1: return 1;  // H
            case   4: return 2;  // Be (BeF₂, BeCl₂)
            case   5: return 3;  // B  (BF₃, BH₃)
            case   6: return 4;  // C  (always 4)
            case   7: return 3;  // N  (NH₃, rarely 4)
            case   8: return 2;  // O  (H₂O, ethers)
            case   9: return 1;  // F  (always 1)
            case  13: return 3;  // Al (AlCl₃)
            case  14: return 4;  // Si (SiO₂, silanes)
            case  15: return 3;  // P  (PCl₃, PH₃), can do 5 but usual is 3
            case  16: return 2;  // S  (H₂S, R₂S), can do 4 or 6 but usual is 2
            case  17: return 1;  // Cl (HCl, RCl), can do 3,5,7 but usual is 1
            case  35: return 1;  // Br
            case  53: return 1;  // I

            case  21: return 3;  // Sc
            case  22: return 4;  // Ti (TiO₂, TiCl₄)
            case  23: return 5;  // V  (VO₄³⁻, VCl₅)
            case  24: return 3;  // Cr (CrO₄²⁻, Cr³⁺ common)
            case  25: return 2;  // Mn (MnO₄⁻ is 7, but Mn²⁺ common → 2 bonds in simple view)
            case  26: return 2;  // Fe (Fe²⁺/³⁺, many complexes)
            case  27: return 2;  // Co
            case  28: return 2;  // Ni (Ni²⁺, Ni(0) in carbonyls = 4, but 2 common)
            case  29: return 1;  // Cu (Cu⁺/Cu²⁺)
            case  30: return 2;  // Zn (always 2+)

            case  50: return 4;  // Sn (SnCl₄)
            case  82: return 4;  // Pb (PbO₂, Pb(IV))

            default:
                if (bonds > 4) bonds = 4;
                break;
        }

        // Final clamp
        return (bonds > 8) ? 4 : bonds;
    }

    inline float getAngles(uint8_t centralZIndex, const std::vector<uint8_t> &neighborZs, const std::vector<sim::fun::BondType> &types)
    {
        const size_t bond_count = types.size();
        if (bond_count < 2) return 0.0f;

        const uint32_t domains_from_bonds = bond_count;

        uint32_t electrons_in_bonds = 0;
        for (const auto& t : types) 
            electrons_in_bonds += static_cast<uint32_t>(t);

        const uint32_t valence = getValenceElectrons(centralZIndex);
        const int32_t lone_pair_electrons = valence - electrons_in_bonds;
        if (lone_pair_electrons < 0) 
            return 0.0f;
        
        const uint32_t lone_pairs = lone_pair_electrons / 2;
        const uint32_t total_domains = domains_from_bonds + lone_pairs;

        float ideal_angle = 0.0f;

        switch (total_domains) {
            case 2:
                ideal_angle = 180.0f * RADIAN; // Linear: AX₂
                break;

            case 3:
                ideal_angle = 120.0f * RADIAN; // Trigonal planar: AX₃
                break;

            case 4:
                if (bond_count == 4) {
                    ideal_angle = 109.5f * RADIAN; // AX₄
                } else if (bond_count == 3) {
                    ideal_angle = 107.0f * RADIAN; // AX₃E
                } else if (bond_count == 2) {
                    ideal_angle = 104.5f * RADIAN; // AX₂E₂
                }
                break;

            case 5:
                ideal_angle = 120.0f * RADIAN; // AX₅ equatorial
                break;

            case 6:
                ideal_angle = 90.0f * RADIAN; // AX₆ octahedral
                break;

            default:
                ideal_angle = 109.5f * RADIAN; // Fallback
                break;
        }

        return ideal_angle;
    }

    inline std::pair<float, float> getAtomConstants(uint32_t ZIndex)
    {
        std::pair<float, float> constants; // {sigma (Å), epsilon (kJ/mol)}
        switch (ZIndex)
        {
            case   1: constants = {2.886f, 0.184f}; break; // H  Hydrogen
            case   2: constants = {2.362f, 0.0844f}; break; // He Helium
            case   3: constants = {3.345f, 0.335f}; break; // Li Lithium
            case   4: constants = {3.051f, 0.418f}; break; // Be Beryllium
            case   5: constants = {3.660f, 0.251f}; break; // B  Boron
            case   6: constants = {3.500f, 0.276f}; break; // C  Carbon
            case   7: constants = {3.300f, 0.460f}; break; // N  Nitrogen
            case   8: constants = {3.100f, 0.585f}; break; // O  Oxygen
            case   9: constants = {2.950f, 0.502f}; break; // F  Fluorine
            case  10: constants = {2.800f, 0.205f}; break; // Ne Neon
            case  11: constants = {3.830f, 0.502f}; break; // Na Sodium
            case  12: constants = {3.580f, 0.544f}; break; // Mg Magnesium
            case  13: constants = {4.010f, 0.335f}; break; // Al Aluminium
            case  14: constants = {3.900f, 0.544f}; break; // Si Silicon
            case  15: constants = {3.800f, 0.711f}; break; // P  Phosphorus
            case  16: constants = {3.700f, 0.879f}; break; // S  Sulfur
            case  17: constants = {3.520f, 1.046f}; break; // Cl Chlorine
            case  18: constants = {3.410f, 0.996f}; break; // Ar Argon
            case  19: constants = {4.230f, 0.335f}; break; // K  Potassium
            case  20: constants = {3.830f, 0.418f}; break; // Ca Calcium
            case  21: constants = {3.400f, 0.084f}; break; // Sc Scandium
            case  22: constants = {3.350f, 0.084f}; break; // Ti Titanium
            case  23: constants = {3.300f, 0.084f}; break; // V  Vanadium
            case 24: constants = {3.250f, 0.084f}; break; // Cr Chromium
            case  25: constants = {3.200f, 0.084f}; break; // Mn Manganese
            case  26: constants = {2.935f, 0.105f}; break; // Fe Iron
            case  27: constants = {2.880f, 0.088f}; break; // Co Cobalt
            case  28: constants = {2.750f, 0.071f}; break; // Ni Nickel
            case  29: constants = {2.650f, 0.100f}; break; // Cu Copper
            case  30: constants = {3.120f, 0.167f}; break; // Zn Zinc
            case  31: constants = {4.200f, 0.251f}; break; // Ga Gallium
            case  32: constants = {4.050f, 0.335f}; break; // Ge Germanium
            case  33: constants = {3.950f, 0.418f}; break; // As Arsenic
            case  34: constants = {3.850f, 0.502f}; break; // Se Selenium
            case  35: constants = {3.700f, 0.879f}; break; // Br Bromine
            case  36: constants = {3.710f, 1.420f}; break; // Kr Krypton
            case  37: constants = {4.600f, 0.335f}; break; // Rb Rubidium
            case  38: constants = {4.300f, 0.418f}; break; // Sr Strontium
            case  39: constants = {3.900f, 0.084f}; break; // Y  Yttrium
            case  40: constants = {3.800f, 0.084f}; break; // Zr Zirconium
            case  41: constants = {3.750f, 0.084f}; break; // Nb Niobium
            case  42: constants = {3.700f, 0.084f}; break; // Mo Molybdenum
            case  43: constants = {3.650f, 0.084f}; break; // Tc Technetium
            case  44: constants = {3.600f, 0.084f}; break; // Ru Ruthenium
            case  45: constants = {3.550f, 0.084f}; break; // Rh Rhodium
            case  46: constants = {3.500f, 0.084f}; break; // Pd Palladium
            case  47: constants = {3.450f, 0.084f}; break; // Ag Silver
            case  48: constants = {3.400f, 0.167f}; break; // Cd Cadmium
            case  49: constants = {4.400f, 0.251f}; break; // In Indium
            case  50: constants = {4.300f, 0.335f}; break; // Sn Tin
            case  51: constants = {4.200f, 0.418f}; break; // Sb Antimony
            case  52: constants = {4.100f, 0.502f}; break; // Te Tellurium
            case  53: constants = {4.000f, 0.879f}; break; // I  Iodine
            case  54: constants = {4.050f, 1.920f}; break; // Xe Xenon
            case  55: constants = {4.900f, 0.335f}; break; // Cs Caesium
            case  56: constants = {4.600f, 0.418f}; break; // Ba Barium
            case  57: constants = {4.000f, 0.100f}; break; // La Lanthanum
            case  58: constants = {3.950f, 0.100f}; break; // Ce Cerium
            case  59: constants = {3.930f, 0.100f}; break; // Pr Praseodymium
            case  60: constants = {3.900f, 0.100f}; break; // Nd Neodymium
            case  61: constants = {3.880f, 0.100f}; break; // Pm Promethium
            case  62: constants = {3.860f, 0.100f}; break; // Sm Samarium
            case  63: constants = {3.840f, 0.100f}; break; // Eu Europium
            case  64: constants = {3.820f, 0.100f}; break; // Gd Gadolinium
            case  65: constants = {3.800f, 0.100f}; break; // Tb Terbium
            case  66: constants = {3.780f, 0.100f}; break; // Dy Dysprosium
            case  67: constants = {3.760f, 0.100f}; break; // Ho Holmium
            case  68: constants = {3.740f, 0.100f}; break; // Er Erbium
            case  69: constants = {3.720f, 0.100f}; break; // Tm Thulium
            case  70: constants = {3.700f, 0.100f}; break; // Yb Ytterbium
            case  71: constants = {3.680f, 0.100f}; break; // Lu Lutetium
            case  72: constants = {3.650f, 0.084f}; break; // Hf Hafnium
            case  73: constants = {3.620f, 0.084f}; break; // Ta Tantalum
            case  74: constants = {3.590f, 0.084f}; break; // W  Tungsten
            case  75: constants = {3.570f, 0.084f}; break; // Re Rhenium
            case  76: constants = {3.550f, 0.084f}; break; // Os Osmium
            case  77: constants = {3.530f, 0.084f}; break; // Ir Iridium
            case  78: constants = {3.510f, 0.084f}; break; // Pt Platinum
            case  79: constants = {3.490f, 0.084f}; break; // Au Gold
            case  80: constants = {3.470f, 0.167f}; break; // Hg Mercury
            case  81: constants = {4.450f, 0.251f}; break; // Tl Thallium
            case  82: constants = {4.400f, 0.335f}; break; // Pb Lead
            case  83: constants = {4.350f, 0.418f}; break; // Bi Bismuth
            case  84: constants = {4.300f, 0.502f}; break; // Po Polonium (estimated)
            case  85: constants = {4.250f, 0.628f}; break; // At Astatine (estimated)
            case  86: constants = {4.200f, 1.255f}; break; // Rn Radon (estimated)
            case  87: constants = {5.000f, 0.335f}; break; // Fr Francium (estimated)
            case  88: constants = {4.700f, 0.418f}; break; // Ra Radium
            case  89: constants = {4.100f, 0.100f}; break; // Ac Actinium
            case  90: constants = {4.050f, 0.100f}; break; // Th Thorium
            case  91: constants = {4.000f, 0.100f}; break; // Pa Protactinium
            case  92: constants = {3.950f, 0.100f}; break; // U  Uranium
            case  93: constants = {3.900f, 0.100f}; break; // Np Neptunium
            case 94: constants = {3.850f, 0.100f}; break; // Pu Plutonium
            case  95: constants = {3.850f, 0.100f}; break; // Am Americium
            case  96: constants = {3.850f, 0.100f}; break; // Cm Curium
            case  97: constants = {3.850f, 0.100f}; break; // Bk Berkelium
            case  98: constants = {3.850f, 0.100f}; break; // Cf Californium
            case  99: constants = {3.500f, 1.255f}; break; // Es (Einsteinium)
            case 100: constants = {3.450f, 1.380f}; break; // Fm (Fermium)
            case 101: constants = {3.400f, 1.464f}; break; // Md (Mendelevium)
            case 102: constants = {3.350f, 1.548f}; break; // No (Nobelium)
            case 103: constants = {3.300f, 1.130f}; break; // Lr (Lawrencium)
            case 104: constants = {3.700f, 0.084f}; break; // Rf Rutherfordium
            case 105: constants = {3.680f, 0.084f}; break; // Db Dubnium
            case 106: constants = {3.660f, 0.084f}; break; // Sg Seaborgium
            case 107: constants = {3.640f, 0.084f}; break; // Bh Bohrium
            case 108: constants = {3.620f, 0.084f}; break; // Hs Hassium
            case 109: constants = {3.600f, 0.084f}; break; // Mt Meitnerium
            case 110: constants = {3.580f, 0.084f}; break; // Ds Darmstadtium
            case 111: constants = {3.560f, 0.084f}; break; // Rg Roentgenium
            case 112: constants = {3.540f, 0.167f}; break; // Cn Copernicium
            case 113: constants = {4.500f, 0.251f}; break; // Nh Nihonium
            case 114: constants = {4.450f, 0.335f}; break; // Fl Flerovium
            case 115: constants = {4.400f, 0.418f}; break; // Mc Moscovium
            case 116: constants = {4.350f, 0.502f}; break; // Lv Livermorium
            case 117: constants = {4.300f, 0.628f}; break; // Ts Tennessine
            case 118: constants = {4.250f, 0.836f}; break; // Og Oganesson
        default:
            constants = {2.50f, 0.1255f};
            break; // H fallback
        }

        constants.first *= MULT_FACTOR;
        constants.second *= MULT_FACTOR;

        return constants;
    }

    inline float getElectronegativity(uint8_t Z)
    {
        static const std::map<uint8_t, float> en = {
            {1,2.20f},{3,0.98f},{4,1.57f},{5,2.04f},{6,2.55f},{7,3.04f},{8,3.44f},{9,3.98f},
            {11,0.93f},{12,1.31f},{13,1.61f},{14,1.90f},{15,2.19f},{16,2.58f},{17,3.16f},
            {19,0.82f},{20,1.00f},{21,1.36f},{22,1.54f},{23,1.63f},{24,1.66f},{25,1.55f},{26,1.83f},{27,1.88f},{28,1.91f},{29,1.90f},{30,1.65f},
            {31,1.81f},{32,2.01f},{33,2.18f},{34,2.55f},{35,2.96f},{36,3.00f},
            {37,0.82f},{38,0.95f},{39,1.22f},{40,1.33f},{41,1.60f},{42,2.16f},{43,1.90f},{44,2.20f},{45,2.28f},{46,2.20f},{47,1.93f},{48,1.69f},
            {49,1.78f},{50,1.96f},{51,2.05f},{52,2.10f},{53,2.66f},{54,2.60f},
            {55,0.79f},{56,0.89f},{57,1.10f},{72,1.30f},{73,1.50f},{74,2.36f},{75,1.90f},{76,2.20f},{77,2.20f},{78,2.28f},{79,2.54f},{80,2.00f},
            {81,1.62f},{82,2.33f},{83,2.02f},{84,2.00f},{85,2.20f}
        };
        auto it = en.find(Z);
        return (it != en.end()) ? it->second : 0.0f;
    }

    inline std::string getAtomName(uint8_t Z)
    {
        static const std::vector<std::string> names = {
            "", "Hydrogen","Helium","Lithium","Beryllium","Boron","Carbon","Nitrogen","Oxygen","Fluorine","Neon",
            "Sodium","Magnesium","Aluminium","Silicon","Phosphorus","Sulfur","Chlorine","Argon","Potassium","Calcium",
            "Scandium","Titanium","Vanadium","Chromium","Manganese","Iron","Cobalt","Nickel","Copper","Zinc",
            "Gallium","Germanium","Arsenic","Selenium","Bromine","Krypton",
            "Rubidium","Strontium","Yttrium","Zirconium","Niobium","Molybdenum","Technetium","Ruthenium","Rhodium","Palladium","Silver","Cadmium",
            "Indium","Tin","Antimony","Tellurium","Iodine","Xenon",
            "Caesium","Barium","Lanthanum","Cerium","Praseodymium","Neodymium","Promethium","Samarium","Europium","Gadolinium",
            "Terbium","Dysprosium","Holmium","Erbium","Thulium","Ytterbium","Lutetium",
            "Hafnium","Tantalum","Tungsten","Rhenium","Osmium","Iridium","Platinum","Gold","Mercury",
            "Thallium","Lead","Bismuth","Polonium","Astatine","Radon",
            "Francium","Radium","Actinium","Thorium","Protactinium","Uranium","Neptunium","Plutonium","Americium","Curium",
            "Berkelium","Californium","Einsteinium","Fermium","Mendelevium","Nobelium","Lawrencium",
            "Rutherfordium","Dubnium","Seaborgium","Bohrium","Hassium","Meitnerium","Darmstadtium","Roentgenium","Copernicium",
            "Nihonium","Flerovium","Moscovium","Livermorium","Tennessine","Oganesson"
        };
        return (Z < names.size()) ? names[Z] : "Unknown";
    }

    inline std::string getAtomLetter(uint8_t Z)
    {
        static const char* letters[] = {
            nullptr,"H","He","Li","Be","B","C","N","O","F","Ne",
            "Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
            "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
            "Ga","Ge","As","Se","Br","Kr",
            "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
            "In","Sn","Sb","Te","I","Xe",
            "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd",
            "Tb","Dy","Ho","Er","Tm","Yb","Lu",
            "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
            "Tl","Pb","Bi","Po","At","Rn",
            "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm",
            "Bk","Cf","Es","Fm","Md","No","Lr",
            "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn",
            "Nh","Fl","Mc","Lv","Ts","Og"
        };
        return (Z < 119) ? letters[Z] : "?";
    }

    inline float getBondLength(uint8_t ZIndex1, uint8_t ZIndex2, sim::fun::BondType type)
    {
        // Fallback: average covalent radii (scaled)
        float base = (getAtomConstants(ZIndex1).first / MULT_FACTOR +
                    getAtomConstants(ZIndex2).first / MULT_FACTOR) / 5.f * BOND_LENGTH_FACTOR;
        float baseBase = base;

        switch (type)
        {
        case sim::fun::BondType::SINGLE:
            // --- Carbon bonds ---
            if ((ZIndex1 == 6 && ZIndex2 == 6) || (ZIndex1 == 6 && ZIndex2 == 6))  base = 1.54f; // C–C
            else if ((ZIndex1 == 6 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 6))  base = 1.09f; // C–H
            else if ((ZIndex1 == 6 && ZIndex2 == 7) || (ZIndex1 == 7 && ZIndex2 == 6))  base = 1.47f; // C–N
            else if ((ZIndex1 == 6 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 6))  base = 1.43f; // C–O
            else if ((ZIndex1 == 6 && ZIndex2 == 9) || (ZIndex1 == 9 && ZIndex2 == 6))  base = 1.35f; // C–F
            else if ((ZIndex1 == 6 && ZIndex2 == 15) || (ZIndex1 == 15 && ZIndex2 == 6)) base = 1.85f; // C–P
            else if ((ZIndex1 == 6 && ZIndex2 == 16) || (ZIndex1 == 16 && ZIndex2 == 6)) base = 1.82f; // C–S
            else if ((ZIndex1 == 6 && ZIndex2 == 17) || (ZIndex1 == 17 && ZIndex2 == 6)) base = 1.77f; // C–Cl
            else if ((ZIndex1 == 6 && ZIndex2 == 35) || (ZIndex1 == 35 && ZIndex2 == 6)) base = 1.94f; // C–Br
            else if ((ZIndex1 == 6 && ZIndex2 == 53) || (ZIndex1 == 53 && ZIndex2 == 6)) base = 2.14f; // C–I

            // --- Nitrogen bonds ---
            else if ((ZIndex1 == 7 && ZIndex2 == 7)) base = 1.45f; // N–N
            else if ((ZIndex1 == 7 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 7)) base = 1.46f; // N–O
            else if ((ZIndex1 == 7 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 7)) base = 1.01f; // N–H

            // --- Oxygen bonds ---
            else if ((ZIndex1 == 8 && ZIndex2 == 8)) base = 1.48f; // O–O
            else if ((ZIndex1 == 8 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 8)) base = 0.96f; // O–H
            else if ((ZIndex1 == 8 && ZIndex2 == 16) || (ZIndex1 == 16 && ZIndex2 == 8)) base = 1.57f; // O–S

            // --- Halogens ---
            else if ((ZIndex1 == 9 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 9))   base = 0.92f; // F–H
            else if ((ZIndex1 == 17 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 17)) base = 1.27f; // Cl–H
            else if ((ZIndex1 == 17 && ZIndex2 == 5) || (ZIndex1 == 5 && ZIndex2 == 17)) base = 1.81f; // Cl-B
            else if ((ZIndex1 == 35 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 35)) base = 1.41f; // Br–H
            else if ((ZIndex1 == 53 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 53)) base = 1.61f; // I–H

            // --- Phosphorus ---
            else if ((ZIndex1 == 15 && ZIndex2 == 15)) base = 2.22f; // P–P
            else if ((ZIndex1 == 15 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 15)) base = 1.60f; // P–O
            else if ((ZIndex1 == 15 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 15)) base = 1.42f; // P–H
            else if ((ZIndex1 == 15 && ZIndex2 == 9) || (ZIndex1 == 9 && ZIndex2 == 15)) base = 1.561f; // P–F

            // --- Sulfur ---
            else if ((ZIndex1 == 16 && ZIndex2 == 16)) base = 2.05f; // S–S
            else if ((ZIndex1 == 16 && ZIndex2 == 1) || (ZIndex1 == 1 && ZIndex2 == 16)) base = 1.34f; // S–H

            // --- Metals ---
            else if ((ZIndex1 == 26 && ZIndex2 == 6) || (ZIndex1 == 6 && ZIndex2 == 26)) base = 1.92f; // Fe–C
            else if ((ZIndex1 == 26 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 26)) base = 1.62f; // Fe–O
            else if ((ZIndex1 == 28 && ZIndex2 == 6) || (ZIndex1 == 6 && ZIndex2 == 28)) base = 1.88f; // Ni–C
            else if ((ZIndex1 == 29 && ZIndex2 == 6) || (ZIndex1 == 6 && ZIndex2 == 29)) base = 1.90f; // Cu–C
            else if ((ZIndex1 == 30 && ZIndex2 == 6) || (ZIndex1 == 6 && ZIndex2 == 30)) base = 1.95f; // Zn–C

            // --- Aromatic C–C ---
            else if ((ZIndex1 == 6 && ZIndex2 == 6) &&
                    (getAtomConstants(ZIndex1).second || getAtomConstants(ZIndex2).second))
                base = 1.39f;

            break;

        case sim::fun::BondType::DOUBLE:
            if ((ZIndex1 == 6 && ZIndex2 == 6)) base = 1.34f; // C=C
            else if ((ZIndex1 == 6 && ZIndex2 == 7) || (ZIndex1 == 7 && ZIndex2 == 6)) base = 1.27f; // C=N
            else if ((ZIndex1 == 6 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 6)) base = 1.22f; // C=O
            else if ((ZIndex1 == 6 && ZIndex2 == 16) || (ZIndex1 == 16 && ZIndex2 == 6)) base = 1.61f; // C=S
            else if ((ZIndex1 == 7 && ZIndex2 == 7)) base = 1.24f; // N=N
            else if ((ZIndex1 == 7 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 7)) base = 1.28f; // N=O
            else if ((ZIndex1 == 8 && ZIndex2 == 8)) base = 1.21f; // O=O
            else if ((ZIndex1 == 8 && ZIndex2 == 16) || (ZIndex1 == 16 && ZIndex2 == 8)) base = 1.43f; // S=O
            else if ((ZIndex1 == 15 && ZIndex2 == 8) || (ZIndex1 == 8 && ZIndex2 == 15)) base = 1.44f; // P=O
            break;

        case sim::fun::BondType::TRIPLE:
            if ((ZIndex1 == 6 && ZIndex2 == 6)) base = 1.20f; // C≡C
            else if ((ZIndex1 == 6 && ZIndex2 == 7) || (ZIndex1 == 7 && ZIndex2 == 6)) base = 1.16f; // C≡N
            else if ((ZIndex1 == 7 && ZIndex2 == 7)) base = 1.10f; // N≡N
            break;

        case sim::fun::BondType::QUADRUPLE:
            base *= 0.55f;
            break;

        default:
            break;
        }

        if (base != baseBase)
            base *= MULT_FACTOR;

        return base;
    }

    inline float getBondEnergy(uint8_t Z1, uint8_t Z2, sim::fun::BondType type)
    {
        if (Z1 > Z2) std::swap(Z1, Z2);

        static const std::map<std::pair<uint8_t, uint8_t>, float> singleBondEnergy = 
        {
            // C–X
            {{6,6},  348.0f},   // C–C  (alkane)
            {{6,1},  413.0f},   // C–H
            {{6,7},  305.0f},   // C–N  (amine)
            {{6,8},  358.0f},   // C–O  (alcohol)
            {{6,9},  485.0f},   // C–F
            {{6,15}, 272.0f},   // C–P
            {{6,16}, 272.0f},   // C–S  (thioether)
            {{6,17}, 327.0f},   // C–Cl
            {{6,35}, 285.0f},   // C–Br
            {{6,53}, 234.0f},   // C–I

            // N–X
            {{7,7},  163.0f},   // N–N  (hydrazine)
            {{7,1},  391.0f},   // N–H
            {{7,8},  201.0f},   // N–O  (hydroxylamine)

            // O–X
            {{8,8},  146.0f},   // O–O  (peroxide)
            {{8,1},  463.0f},   // O–H
            {{8,16}, 265.0f},   // O–S  (sulfoxide)

            // Halogens
            {{1,9},   565.0f},  // H–F
            {{1,17},  431.0f},  // H–Cl
            {{1,35},  366.0f},  // H–Br
            {{1,53},  298.0f},  // H–I
            {{9,9},   159.0f},  // F–F
            {{17,17}, 243.0f},  // Cl–Cl
            {{35,35}, 193.0f},  // Br–Br
            {{53,53}, 151.0f},  // I–I

            // P, S, Si
            {{15,15}, 201.0f},  // P–P
            {{15,1},  322.0f},  // P–H
            {{15,8},  360.0f},  // P–O
            {{16,16}, 226.0f},  // S–S  (disulfide)
            {{14,14}, 226.0f},  // Si–Si
            {{14,1},  318.0f},  // Si–H
            {{14,8},  452.0f},  // Si–O

            // Metal–carbon (common in catalysis)
            {{26,6},  180.0f},  // Fe–C  (typical organometallic)
            {{27,6},  200.0f},  // Co–C
            {{28,6},  220.0f},  // Ni–C
            {{29,6},  220.0f},  // Cu–C
            {{44,6},  280.0f},  // Ru–C
            {{45,6},  310.0f},  // Rh–C
            {{46,6},  340.0f},  // Pd–C
            {{77,6},  380.0f},  // Ir–C
            {{78,6},  400.0f},  // Pt–C
            {{79,6},  250.0f},  // Au–C

            // Other common
            {{5,1},   389.0f},  // B–H
            {{5,6},   435.0f},  // B–C
            {{13,8},  502.0f},  // Al–O
        };

        float base = 0.0f;

        auto key = std::make_pair(Z1, Z2);
        auto it = singleBondEnergy.find(key);
        if (it != singleBondEnergy.end()) 
            base = it->second;

        switch (type) 
        {
            case sim::fun::BondType::SINGLE:
                if (base > 0.0f) return base;

            case sim::fun::BondType::DOUBLE:
                if (Z1 == 6 && Z2 == 6) base = 614.0f;     // C=C
                if (Z1 == 6 && Z2 == 8) base = 799.0f;     // C=O (carbonyl)
                if (Z1 == 6 && Z2 == 7) base = 615.0f;     // C=N
                if (Z1 == 7 && Z2 == 8) base = 607.0f;     // N=O
                if (Z1 == 8 && Z2 == 8) base = 498.0f;     // O=O (molecular oxygen)
                if (Z1 == 6 && Z2 == 16) base = 536.0f;    // C=S
                if (Z1 == 15 && Z2 == 8) base = 544.0f;    // P=O

                if (base > 0.0f) return base * 1.65f;      // average double/single ratio
                break;

            case sim::fun::BondType::TRIPLE:
                if (Z1 == 6 && Z2 == 6) base = 837.0f;     // C≡C
                if (Z1 == 6 && Z2 == 7) base = 887.0f;     // C≡N
                if (Z1 == 7 && Z2 == 7) base = 941.0f;     // N≡N
                if (base > 0.0f) return base * 2.40f;
                break;

            case sim::fun::BondType::QUADRUPLE:
                if (base > 0.0f) return base * 3.2f;
                break;

            default:
                break;
        }

        if (base > 0.f)
            return base; // J/Mol

        if (base == 0.0f) 
        {
            float DA = 0.0f, DB = 0.0f;
            switch (Z1) 
            {
                case 1: DA = 436.0f; break;   // H–H
                case 6: DA = 348.0f; break;   // C–C
                case 7: DA = 163.0f; break;   // N–N
                case 8: DA = 146.0f; break;   // O–O
                case 9: DA = 159.0f; break;   // F–F
                case 14: DA = 226.0f; break;  // Si–Si
                case 15: DA = 201.0f; break;  // P–P
                case 16: DA = 226.0f; break;  // S–S
                case 17: DA = 243.0f; break;  // Cl–Cl
                case 35: DA = 193.0f; break;  // Br–Br
                case 53: DA = 151.0f; break;  // I–I
            }
            switch (Z2) 
            {
                case 1: DB = 436.0f; break;
                case 6: DB = 348.0f; break;
                case 7: DB = 163.0f; break;
                case 8: DB = 146.0f; break;
                case 9: DB = 159.0f; break;
                case 14: DB = 226.0f; break;
                case 15: DB = 201.0f; break;
                case 16: DB = 226.0f; break;
                case 17: DB = 243.0f; break;
                case 35: DB = 193.0f; break;
                case 53: DB = 151.0f; break;
            }

            float chiA = getElectronegativity(Z1);
            float chiB = getElectronegativity(Z2);
            float ionic_term = (chiA > 0.0f && chiB > 0.0f) ? 96.5f * (chiA - chiB) * (chiA - chiB) : 0.0f;

            if (DA > 0.0f && DB > 0.0f) 
            {
                base = 0.5f * (DA + DB) + ionic_term;
            } 
            else 
            {
                float r1 = covalent_radius[Z1];
                float r2 = covalent_radius[Z2];
                if (r1 > 0.0f && r2 > 0.0f) 
                {
                    float avg_r = (r1 + r2) * 0.5f;
                    base = 400.0f * std::exp(-avg_r / 1.0f);
                } 
                else 
                {
                    base = 250.0f; 
                }
            }
        }

        switch (type) 
        {
            case sim::fun::BondType::DOUBLE: base *= 1.65f; break;
            case sim::fun::BondType::TRIPLE: base *= 2.40f; break;
            case sim::fun::BondType::QUADRUPLE: base *= 3.20f; break;
            default: break;
        }

        return base;  // J/mol
    }

    inline float getBondHarmonicConstantFromEnergy(uint8_t Z1, uint8_t Z2, sim::fun::BondType type)
    {
        float D = getBondEnergy(Z1, Z2, type);        // kJ/mol
        float r0 = getBondLength(Z1, Z2, type);       // Å

        float D_aJ = D * 1000.0f / AVOGADRO / 4.184f;
        float k = (D / 348.0f) * BOND_K;

        return k * static_cast<float>(type);
    }
}; // namespace constants