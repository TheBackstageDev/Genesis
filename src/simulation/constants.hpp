#pragma once

#include <vector>
#include <map>
#include <string>
#include <SFML/Graphics.hpp>
#include <cmath>

namespace sim
{
    namespace fun { enum class BondType : uint32_t { NONE, SINGLE, DOUBLE, TRIPLE, QUADRUPLE }; };
}; // namespace sim


namespace constants
{
#define M_PI 3.14159265f
#define RADIAN M_PI / 180
#define DEGREE 180 / M_PI

#define MASS_PROTON 1.0007     // Daltons
#define MASS_NEUTRON 1.0008    // Daltons
#define MASS_ELECTRON 1 / 1337 // Daltons

#define EPSILON 0.001f
#define FEMTOSECOND 0.001f // ps
#define MULT_FACTOR 1.f
#define ANGSTROM 1e20f
#define PICOSECOND 1e24f

#define JOULE_TO_CAL 1/4.164f
#define PRESSURE_CONVERSION 16387.9f // converts (kcal/mol)/Å³ → bar

#define VERLET_SKIN 2.f
#define CUTOFF 4.f
#define COULOMB_CUTOFF 12.f * MULT_FACTOR

#define CELL_CUTOFF 12 + VERLET_SKIN

#define AVOGADRO 6.02214076e26f                                   // conversion from Daltons to Kg
#define BOLTZMAN_CONSTANT 1.380649e-23f                           // Boltzman Constant m^2 kg s^-2 K^-1
#define KB (BOLTZMAN_CONSTANT * AVOGADRO * ANGSTROM) / PICOSECOND // A^2 D ps^-2 K^-1

#define REBUILD_THRESHOLD 2.5f * MULT_FACTOR
#define THERMOSTAT_INTERVAL 5
#define BAROSTAT_INTERVAL 5

#define COULOMB_K 1389.354576f // kJ·mol⁻¹· Å ·e⁻²
#define BOND_K 340000.f        // Harmonic force constant
#define ANGLE_K 12000.f       // J/mol/rad² for angular potential
#define BOND_LENGTH_FACTOR 1.f

#define REACTION_INTERVAL 3 // time steps
#define REACTION_CUTOFF 10.f * MULT_FACTOR
#define REACTION_UPDATE_Q 2
#define REACTION_CUT_BO 0.3f

#define COUNT_ATOMS 118

    inline constexpr std::array<float, 119> covalent_radius = {
        0.00f,   // Z=0 (dummy)

        0.32f,   // 1  H   (Pyykkö 32 pm)
        0.46f,   // 2  He  (estimated for rare cases)

        1.33f,   // 3  Li
        0.98f,   // 4  Be
        0.85f,   // 5  B
        0.75f,   // 6  C
        0.71f,   // 7  N
        0.63f,   // 8  O
        0.64f,   // 9  F
        0.58f,   // 10 Ne

        1.55f,   // 11 Na
        1.39f,   // 12 Mg
        1.26f,   // 13 Al
        1.16f,   // 14 Si
        1.11f,   // 15 P
        1.03f,   // 16 S
        0.99f,   // 17 Cl
        0.96f,   // 18 Ar

        1.96f,   // 19 K
        1.71f,   // 20 Ca
        1.48f,   // 21 Sc
        1.36f,   // 22 Ti
        1.34f,   // 23 V
        1.22f,   // 24 Cr
        1.19f,   // 25 Mn
        1.16f,   // 26 Fe
        1.11f,   // 27 Co
        1.10f,   // 28 Ni
        1.12f,   // 29 Cu
        1.18f,   // 30 Zn
        1.24f,   // 31 Ga
        1.21f,   // 32 Ge
        1.21f,   // 33 As
        1.16f,   // 34 Se
        1.14f,   // 35 Br
        1.17f,   // 36 Kr

        2.10f,   // 37 Rb
        1.85f,   // 38 Sr
        1.63f,   // 39 Y
        1.54f,   // 40 Zr
        1.47f,   // 41 Nb
        1.38f,   // 42 Mo
        1.28f,   // 43 Tc
        1.25f,   // 44 Ru
        1.25f,   // 45 Rh
        1.20f,   // 46 Pd
        1.28f,   // 47 Ag
        1.36f,   // 48 Cd
        1.42f,   // 49 In
        1.40f,   // 50 Sn
        1.40f,   // 51 Sb
        1.36f,   // 52 Te
        1.33f,   // 53 I
        1.31f,   // 54 Xe

        2.32f,   // 55 Cs
        1.96f,   // 56 Ba
        1.80f,   // 57 La
        1.63f,   // 58 Ce
        1.76f,   // 59 Pr
        1.74f,   // 60 Nd
        1.73f,   // 61 Pm
        1.72f,   // 62 Sm
        1.68f,   // 63 Eu
        1.69f,   // 64 Gd
        1.68f,   // 65 Tb
        1.67f,   // 66 Dy
        1.66f,   // 67 Ho
        1.65f,   // 68 Er
        1.64f,   // 69 Tm
        1.70f,   // 70 Yb
        1.62f,   // 71 Lu

        1.52f,   // 72 Hf
        1.46f,   // 73 Ta
        1.37f,   // 74 W
        1.31f,   // 75 Re
        1.29f,   // 76 Os
        1.22f,   // 77 Ir
        1.23f,   // 78 Pt
        1.24f,   // 79 Au
        1.33f,   // 80 Hg
        1.44f,   // 81 Tl
        1.51f,   // 82 Pb
        1.45f,   // 83 Bi
        1.47f,   // 84 Po
        1.42f,   // 85 At
        1.42f,   // 86 Rn

        2.23f,   // 87 Fr
        2.01f,   // 88 Ra

        1.86f,   // 89 Ac
        1.75f,   // 90 Th
        1.69f,   // 91 Pa
        1.70f,   // 92 U
        1.71f,   // 93 Np
        1.72f,   // 94 Pu
        1.66f,   // 95 Am
        1.66f,   // 96 Cm
        1.68f,   // 97 Bk
        1.68f,   // 98 Cf
        1.65f,   // 99 Es
        1.67f,   // 100 Fm
        1.73f,   // 101 Md
        1.76f,   // 102 No
        1.61f,   // 103 Lr

        1.57f,   // 104 Rf
        1.49f,   // 105 Db
        1.43f,   // 106 Sg
        1.41f,   // 107 Bh
        1.34f,   // 108 Hs
        1.29f,   // 109 Mt
        1.28f,   // 110 Ds
        1.21f,   // 111 Rg
        1.22f,   // 112 Cn
        1.36f,   // 113 Nh
        1.62f,   // 114 Fl
        1.75f,   // 115 Mc
        1.65f,   // 116 Lv
        1.57f,   // 117 Ts
        1.57f    // 118 Og
    };

    constexpr uint8_t NEUTRON_COUNTS[119] = 
    {
        0,   // 0: dummy
        0,   // 1  H  → ¹H (protium)
        2,   // 2  He → ⁴He
        4,   // 3  Li → ⁷Li (most abundant)
        5,   // 4  Be → ⁹Be
        6,   // 5  B  → ¹¹B (slightly more abundant)
        6,   // 6  C  → ¹²C (by definition for atomic mass unit)
        7,   // 7  N  → ¹⁴N
        8,   // 8  O  → ¹⁶O
        10,  // 9  F  → ¹⁹F
        10,  // 10 Ne → ²⁰Ne
        12,  // 11 Na → ²³Na
        14,  // 12 Mg → ²⁴Mg
        14,  // 13 Al → ²⁷Al
        14,  // 14 Si → ²⁸Si
        16,  // 15 P  → ³¹P
        16,  // 16 S  → ³²S
        18,  // 17 Cl → ³⁵Cl (most abundant)
        20,  // 18 Ar → ⁴⁰Ar
        20,  // 19 K  → ³⁹K
        20,  // 20 Ca → ⁴⁰Ca
        24,  // 21 Sc → ⁴⁵Sc
        26,  // 22 Ti → ⁴⁸Ti
        28,  // 23 V  → ⁵¹V
        28,  // 24 Cr → ⁵²Cr
        30,  // 25 Mn → ⁵⁵Mn
        30,  // 26 Fe → ⁵⁶Fe
        32,  // 27 Co → ⁵⁹Co
        30,  // 28 Ni → ⁵⁸Ni (most abundant)
        34,  // 29 Cu → ⁶³Cu
        34,  // 30 Zn → ⁶⁴Zn
        39,  // 31 Ga → ⁶⁹Ga (slightly more)
        36,  // 32 Ge → ⁷⁴Ge (most abundant)
        42,  // 33 As → ⁷⁵As
        42,  // 34 Se → ⁸⁰Se (most abundant)
        44,  // 35 Br → ⁷⁹Br (very close to ⁸¹Br)
        48,  // 36 Kr → ⁸⁴Kr
        48,  // 37 Rb → ⁸⁵Rb
        50,  // 38 Sr → ⁸⁸Sr
        50,  // 39 Y  → ⁸⁹Y
        50,  // 40 Zr → ⁹⁰Zr
        52,  // 41 Nb → ⁹³Nb
        54,  // 42 Mo → ⁹⁸Mo
        55,  // 43 Tc → ⁹⁸Tc (longest-lived; no stable)
        58,  // 44 Ru → ¹⁰²Ru (most abundant)
        58,  // 45 Rh → ¹⁰³Rh
        64,  // 46 Pd → ¹⁰⁶Pd (most abundant)
        64,  // 47 Ag → ¹⁰⁷Ag (close to ¹⁰⁹Ag)
        66,  // 48 Cd → ¹¹²Cd
        66,  // 49 In → ¹¹⁵In
        70,  // 50 Sn → ¹²⁰Sn (most abundant)
        74,  // 51 Sb → ¹²¹Sb (close to ¹²³Sb)
        76,  // 52 Te → ¹³⁰Te
        77,  // 53 I  → ¹²⁷I
        78,  // 54 Xe → ¹³²Xe (most abundant)
        81,  // 55 Cs → ¹³³Cs
        82,  // 56 Ba → ¹³⁸Ba
        82,  // 57 La → ¹³⁹La
        82,  // 58 Ce → ¹⁴⁰Ce
        84,  // 59 Pr → ¹⁴¹Pr
        90,  // 60 Nd → ¹⁴²Nd (most abundant)
        94,  // 61 Pm → ¹⁴⁵Pm (longest-lived; no stable)
        90,  // 62 Sm → ¹⁵²Sm
        94,  // 63 Eu → ¹⁵³Eu (close to ¹⁵¹Eu)
        98,  // 64 Gd → ¹⁵⁸Gd (most abundant)
        98,  // 65 Tb → ¹⁵⁹Tb
        102, // 66 Dy → ¹⁶⁴Dy
        104, // 67 Ho → ¹⁶⁵Ho
        104, // 68 Er → ¹⁶⁶Er
        106, // 69 Tm → ¹⁶⁹Tm
        104, // 70 Yb → ¹⁷⁴Yb
        106, // 71 Lu → ¹⁷⁵Lu (¹⁷⁶Lu is long-lived but low abundance)
        110, // 72 Hf → ¹⁸⁰Hf
        110, // 73 Ta → ¹⁸¹Ta
        110, // 74 W  → ¹⁸⁴W
        114, // 75 Re → ¹⁸⁷Re
        114, // 76 Os → ¹⁹²Os (most abundant)
        116, // 77 Ir → ¹⁹³Ir
        118, // 78 Pt → ¹⁹⁵Pt (most abundant)
        118, // 79 Au → ¹⁹⁷Au
        120, // 80 Hg → ²⁰²Hg (most abundant)
        122, // 81 Tl → ²⁰⁵Tl
        122, // 82 Pb → ²⁰⁸Pb (most abundant / very stable)
        124, // 83 Bi → ²⁰⁹Bi (very long-lived; effectively stable)
        126, // 84 Po → ²¹⁰Po (longest-lived)
        136, // 85 At → ²¹⁰At or ²¹¹At (~short)
        136, // 86 Rn → ²²²Rn
        138, // 87 Fr → ²²³Fr
        138, // 88 Ra → ²²⁶Ra
        138, // 89 Ac → ²²⁷Ac
        140, // 90 Th → ²³²Th
        140, // 91 Pa → ²³¹Pa
        146, // 92 U  → ²³⁸U
        144, // 93 Np → ²³⁷Np
        150, // 94 Pu → ²⁴⁴Pu (longest-lived)
        148, // 95 Am → ²⁴³Am
        150, // 96 Cm → ²⁴⁷Cm or ²⁴⁸Cm
        150, // 97 Bk → ²⁴⁷Bk
        150, // 98 Cf → ²⁵¹Cf
        152, // 99 Es → ²⁵²Es
        157, // 100 Fm → ²⁵⁷Fm
        157, // 101 Md → ²⁵⁸Md
        157, // 102 No → ²⁵⁹No
        157, // 103 Lr → ²⁶⁰Lr or ²⁶²Lr
        160, // 104 Rf → ²⁶¹Rf or ²⁶⁷Rf (very short-lived)
        160, // 105 Db → ~²⁶⁸Db
        160, // 106 Sg → ~²⁶⁹Sg
        160, // 107 Bh → ~²⁷⁰Bh
        160, // 108 Hs → ~²⁶⁹Hs or ²⁷⁰Hs
        160, // 109 Mt → ~²⁷⁸Mt
        160, // 110 Ds → ~²⁸¹Ds
        160, // 111 Rg → ~²⁸⁰Rg or ²⁸²Rg
        160, // 112 Cn → ~²⁸⁵Cn
        160, // 113 Nh → ~²⁸⁶Nh
        160, // 114 Fl → ~²⁸⁹Fl
        160, // 115 Mc → ~²⁸⁹Mc or ²⁹⁰Mc
        160, // 116 Lv → ~²⁹³Lv
        160, // 117 Ts → ~²⁹⁴Ts
        160   // 118 Og → ~²⁹⁴Og
    };

    inline constexpr float VDW_RADII[119] = 
    {
        0.00f,   // Z=0 (dummy)

        1.10f,   //  1  H 
        1.40f,   //  2  He

        1.81f,   //  3  Li
        1.53f,   //  4  Be
        1.92f,   //  5  B
        1.70f,   //  6  C
        1.55f,   //  7  N
        1.52f,   //  8  O
        1.47f,   //  9  F
        1.54f,   // 10 Ne

        2.27f,   // 11 Na
        1.73f,   // 12 Mg
        1.84f,   // 13 Al
        2.10f,   // 14 Si
        1.80f,   // 15 P
        1.80f,   // 16 S
        1.75f,   // 17 Cl
        1.88f,   // 18 Ar

        2.75f,   // 19 K
        2.31f,   // 20 Ca
        2.11f,   // 21 Sc
        2.00f,   // 22 Ti
        2.00f,   // 23 V
        2.00f,   // 24 Cr
        2.00f,   // 25 Mn
        2.00f,   // 26 Fe
        2.00f,   // 27 Co
        1.63f,   // 28 Ni 
        1.40f,   // 29 Cu 
        1.39f,   // 30 Zn
        1.87f,   // 31 Ga
        2.11f,   // 32 Ge
        1.85f,   // 33 As
        1.90f,   // 34 Se
        1.83f,   // 35 Br
        2.02f,   // 36 Kr

        3.03f,   // 37 Rb
        2.49f,   // 38 Sr
        2.00f,   // 39 Y 
        2.00f,   // 40 Zr
        2.00f,   // 41 Nb
        2.00f,   // 42 Mo
        2.00f,   // 43 Tc
        2.00f,   // 44 Ru
        2.00f,   // 45 Rh
        1.63f,   // 46 Pd 
        1.72f,   // 47 Ag
        1.58f,   // 48 Cd
        1.93f,   // 49 In
        2.17f,   // 50 Sn
        2.06f,   // 51 Sb
        2.06f,   // 52 Te
        1.98f,   // 53 I
        2.16f,   // 54 Xe

        3.43f,   // 55 Cs
        2.68f,   // 56 Ba

        2.35f,   // 57 La
        2.30f,   // 58 Ce
        2.30f,   // 59 Pr
        2.30f,   // 60 Nd
        2.30f,   // 61 Pm
        2.30f,   // 62 Sm
        2.30f,   // 63 Eu
        2.30f,   // 64 Gd
        2.30f,   // 65 Tb
        2.30f,   // 66 Dy
        2.30f,   // 67 Ho
        2.30f,   // 68 Er
        2.30f,   // 69 Tm
        2.30f,   // 70 Yb
        2.30f,   // 71 Lu

        2.23f,   // 72 Hf   
        2.00f,   // 73 Ta
        2.00f,   // 74 W
        2.00f,   // 75 Re
        2.00f,   // 76 Os
        2.00f,   // 77 Ir
        1.75f,   // 78 Pt   
        1.66f,   // 79 Au
        1.55f,   // 80 Hg   
        1.96f,   // 81 Tl
        2.02f,   // 82 Pb
        2.07f,   // 83 Bi
        1.97f,   // 84 Po
        2.02f,   // 85 At
        2.20f,   // 86 Rn

        3.48f,   // 87 Fr
        2.83f,   // 88 Ra

        2.40f,   // 89 Ac
        2.40f,   // 90 Th
        2.30f,   // 91 Pa
        2.30f,   // 92 U
        2.30f,   // 93 Np
        2.30f,   // 94 Pu
        2.30f,   // 95 Am
        2.30f,   // 96 Cm
        2.30f,   // 97 Bk
        2.30f,   // 98 Cf
        2.30f,   // 99 Es
        2.30f,   //100 Fm
        2.30f,   //101 Md
        2.30f,   //102 No
        2.30f,   //103 Lr

        2.20f,   //104 Rf
        2.10f,   //105 Db
        2.00f,   //106 Sg
        2.00f,   //107 Bh
        2.00f,   //108 Hs
        2.00f,   //109 Mt
        2.00f,   //110 Ds
        2.00f,   //111 Rg
        2.00f,   //112 Cn
        2.00f,   //113 Nh
        2.00f,   //114 Fl
        2.00f,   //115 Mc
        2.00f,   //116 Lv
        2.00f,   //117 Ts
        2.00f    //118 Og
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
            case  6: return sf::Color(138, 138, 138);    // C  - Carbon       (Gray)
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
            case 18: return sf::Color(115, 3, 192);    // Ar - Argon        (Dark Purple)
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

    struct ReaxParams 
    {
        float p_boc1, p_boc2, p_boc3, p_boc4, p_boc5;  // For corrections
        float p_bo1, p_bo2, p_bo3, p_bo4, p_bo5, p_bo6;  // Exponents for sigma/pi/pp
        float p_be1 = 0.f, p_be2 = 10.f;
        float r0_sigma, r0_pi, r0_pp;  // Equilibrium distances
        float De_sigma, De_pi, De_pp;  // Dissociation energies

        float gamma;      // shielding (Å⁻¹)
        float r_vdw;      // van der Waals radius (Å)
        float D_vdw;     // van der Waals depth (kcal/mol)
    };

    struct PairReaxParams 
    {
        float De_sigma, De_pi, De_pp;
    };

    static std::map<uint8_t, ReaxParams> reaxParams;

    static ReaxParams& getParams(uint8_t Z)
    {
        if (!reaxParams.count(Z))
            return reaxParams[6];

        return reaxParams[Z];
    }

    static PairReaxParams getPairReaxParams(uint8_t Zi, uint8_t Zj) 
    {
        const auto& pi = constants::getParams(Zi);
        const auto& pj = constants::getParams(Zj);
        PairReaxParams p;
        p.De_sigma = std::sqrt(pi.De_sigma * pj.De_sigma);
        p.De_pi    = std::sqrt(pi.De_pi    * pj.De_pi);
        p.De_pp    = std::sqrt(pi.De_pp    * pj.De_pp);
        return p;
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

    inline float lonePairsBO(uint8_t ZIndex, float BO)
    {
        return (getValenceElectrons(ZIndex) - BO) / 2;
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

                if (centralZIndex == 15)
        {
            ideal_angle = 0.0f;
        }

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

    inline float getBondOrderExponent(uint8_t Zi, uint8_t Zj, sim::fun::BondType type = sim::fun::BondType::SINGLE)
    {
        if (Zi > Zj) std::swap(Zi, Zj);

        float exponent = 4.f;

        if (type == sim::fun::BondType::SINGLE)
        {
            if (Zi == 1 && Zj == 1)   exponent = 4.0f;   // H-H
            if (Zi == 1 && Zj == 6)   exponent = 3.8f;   // C-H
            if (Zi == 1 && Zj == 7)   exponent = 4.1f;   // N-H
            if (Zi == 1 && Zj == 8)   exponent = 4.5f;   // O-H
            if (Zi == 1 && Zj == 9)   exponent = 5.0f;   // F-H
            if (Zi == 1 && Zj == 17)  exponent = 4.2f;   // Cl-H
            if (Zi == 1 && Zj == 35)  exponent = 4.0f;   // Br-H
    
            if (Zi == 6 && Zj == 6)   exponent = 3.8f;   // C-C
            if (Zi == 6 && Zj == 7)   exponent = 4.0f;   // C-N
            if (Zi == 6 && Zj == 8)   exponent = 4.2f;   // C-O
            if (Zi == 6 && Zj == 9)   exponent = 4.8f;   // C-F
            if (Zi == 6 && Zj == 17)  exponent = 4.0f;   // C-Cl
    
            if (Zi == 7 && Zj == 7)   exponent = 4.2f;   // N-N
            if (Zi == 7 && Zj == 8)   exponent = 4.4f;   // N-O
            if (Zi == 7 && Zj == 9)   exponent = 5.0f;   // N-F
    
            if (Zi == 8 && Zj == 8)   exponent = 4.5f;   // O-O
            if (Zi == 8 && Zj == 17)  exponent = 4.3f;   // O-Cl
    
            if (Zi == 9 && Zj == 9)   exponent = 5.2f;   // F-F
            if (Zi == 17 && Zj == 17) exponent = 4.0f;   // Cl-Cl
            if (Zi == 35 && Zj == 35) exponent = 3.8f;   // Br-Br
    
            // === METALS / OTHER ===
            if (Zi == 13 && Zj == 13) exponent = 3.5f;   // Al-Al
            if (Zi == 14 && Zj == 14) exponent = 3.6f;   // Si-Si
            if (Zi == 15 && Zj == 15) exponent = 3.5f;   // P-P
            if (Zi == 16 && Zj == 16) exponent = 3.6f;   // S-S
        }

        if (type == sim::fun::BondType::DOUBLE)
        {
            if (Zi == 6 && Zj == 6)   exponent = 4.8f;   // C=C
            if (Zi == 6 && Zj == 8)   exponent = 5.2f;   // C=O
            if (Zi == 7 && Zj == 8)   exponent = 5.0f;   // N=O
            if (Zi == 8 && Zj == 8)   exponent = 5.5f;   // O=O
            exponent = 5.0f;
        }

        if (type ==  sim::fun:: BondType::TRIPLE)
        {
            if (Zi == 6 && Zj == 6)   exponent = 5.5f;   // C≡C
            if (Zi == 7 && Zj == 7)   exponent = 5.8f;   // N≡N
            exponent = 5.5f;
        }

        return exponent;
    }

    inline std::pair<float, float> getAtomConstants(uint32_t ZIndex)
    {
        std::pair<float, float> constants; // {sigma (Å), epsilon (kJ/mol)}
        switch (ZIndex)
        {
            case   1: constants = {2.532f, 0.184f}; break; // H  Hydrogen
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

    inline float getBondLength(uint8_t Z1, uint8_t Z2, sim::fun::BondType type)
    {
        float base = (getAtomConstants(Z1).first + getAtomConstants(Z2).first) 
                    / (4.0f * MULT_FACTOR) * BOND_LENGTH_FACTOR;
        float original = base;

        #define BONDED(a, b) ((Z1 == (a) && Z2 == (b)) || (Z1 == (b) && Z2 == (a)))

        switch (type)
        {
            case sim::fun::BondType::SINGLE:
            
            if      (BONDED(6,1))  base = 1.09f;   // C–H
            else if (BONDED(6,3))  base = 1.82f;   // C–Li 
            else if (BONDED(6,4))  base = 1.70f;   // C–Be 
            else if (BONDED(6,5))  base = 1.56f;   // C–B 
            else if (BONDED(6,6))  base = 1.54f;   // C–C (sp³)
            else if (BONDED(6,7))  base = 1.47f;   // C–N
            else if (BONDED(6,8))  base = 1.43f;   // C–O
            else if (BONDED(6,9))  base = 1.35f;   // C–F
            else if (BONDED(6,11)) base = 2.16f;   // C–Na
            else if (BONDED(6,12)) base = 1.97f;   // C–Mg
            else if (BONDED(6,13)) base = 1.97f;   // C–Al 
            else if (BONDED(6,14)) base = 1.87f;   // C–Si
            else if (BONDED(6,15)) base = 1.85f;   // C–P
            else if (BONDED(6,16)) base = 1.82f;   // C–S
            else if (BONDED(6,17)) base = 1.77f;   // C–Cl
            else if (BONDED(6,33)) base = 1.98f;   // C–As 
            else if (BONDED(6,34)) base = 1.98f;   // C–Se 
            else if (BONDED(6,35)) base = 1.94f;   // C–Br
            else if (BONDED(6,50)) base = 2.15f;   // C–Sn 
            else if (BONDED(6,51)) base = 2.16f;   // C–Sb 
            else if (BONDED(6,52)) base = 2.21f;   // C–Te 
            else if (BONDED(6,53)) base = 2.14f;   // C–I
            else if (BONDED(6,82)) base = 2.16f;   // C–Pb 

            // Nitrogen
            else if (BONDED(7,1))  base = 1.01f;   // N–H
            else if (BONDED(7,5))  base = 1.36f;   // N–B 
            else if (BONDED(7,7))  base = 1.45f;   // N–N
            else if (BONDED(7,8))  base = 1.43f;   // N–O
            else if (BONDED(7,9))  base = 1.39f;   // N–F 
            else if (BONDED(7,13)) base = 2.05f;   // N–Al 
            else if (BONDED(7,14)) base = 1.74f;   // N–Si 
            else if (BONDED(7,15)) base = 1.77f;   // N–P 
            else if (BONDED(7,16)) base = 1.68f;   // N–S 
            else if (BONDED(7,17)) base = 1.87f;   // N–Cl 

            // Oxygen
            else if (BONDED(8,1))  base = 0.96f;   // O–H
            else if (BONDED(8,5))  base = 1.36f;   // O–B 
            else if (BONDED(8,8))  base = 1.48f;   // O–O (peroxide)
            else if (BONDED(8,9))  base = 1.42f;   // O–F 
            else if (BONDED(8,13)) base = 1.75f;   // O–Al 
            else if (BONDED(8,14)) base = 1.63f;   // O–Si
            else if (BONDED(8,15)) base = 1.61f;   // O–P
            else if (BONDED(8,16)) base = 1.57f;   // O–S

            // Halogens
            else if (BONDED(9,9))   base = 1.42f;   // F–F
            else if (BONDED(9,17))  base = 1.63f;   // F–Cl
            else if (BONDED(17,17)) base = 1.99f;   // Cl–Cl
            else if (BONDED(35,35)) base = 2.28f;   // Br–Br
            else if (BONDED(53,53)) base = 2.66f;   // I–I

            // Phosphorus
            else if (BONDED(15,1))  base = 1.42f;   // P–H
            else if (BONDED(15,5))  base = 1.96f;   // P–B 
            else if (BONDED(15,6))  base = 1.85f;   // P–C
            else if (BONDED(15,7))  base = 1.77f;   // P–N 
            else if (BONDED(15,8))  base = 1.61f;   // P–O (P(III))
            else if (BONDED(15,9))  base = 1.56f;   // P–F
            else if (BONDED(15,13)) base = 2.27f;   // P–Al 
            else if (BONDED(15,14)) base = 2.27f;   // P–Si 
            else if (BONDED(15,15)) base = 2.20f;   // P–P
            else if (BONDED(15,16)) base = 2.10f;   // P–S
            else if (BONDED(15,17)) base = 2.01f;   // P–Cl
            else if (BONDED(15,33)) base = 2.30f;   // P–As 
            else if (BONDED(15,34)) base = 2.31f;   // P–Se 
            else if (BONDED(15,35)) base = 2.18f;   // P–Br
            else if (BONDED(15,53)) base = 2.40f;   // P–I

            // Sulfur
            else if (BONDED(16,16)) base = 2.05f;   // S–S
            else if (BONDED(16,5))  base = 1.82f;   // S–B
            else if (BONDED(16,14)) base = 2.15f;   // S–Si
            else if (BONDED(16,1))  base = 1.34f;   // S–H
            else if (BONDED(16,17)) base = 1.99f;   // S–Cl

            // Silicon
            else if (BONDED(14,14)) base = 2.35f;   // Si–Si
            else if (BONDED(14,1))  base = 1.48f;   // Si–H
            else if (BONDED(14,17)) base = 2.05f;   // Si–Cl

            // Boron
            else if (BONDED(5,5))   base = 1.76f;   // B–B
            else if (BONDED(5,1))   base = 1.19f;   // B–H
            else if (BONDED(5,9))   base = 1.30f;   // B–F
            else if (BONDED(5,17))  base = 1.81f;   // B–Cl

            // Metals
            else if (BONDED(3,3))   base = 2.67f;   // Li–Li
            else if (BONDED(11,11)) base = 3.08f;   // Na–Na
            else if (BONDED(12,12)) base = 2.82f;   // Mg–Mg
            else if (BONDED(13,13)) base = 2.86f;   // Al–Al
            else if (BONDED(26,6))  base = 1.92f;   // Fe–C
            else if (BONDED(26,8))  base = 1.62f;   // Fe–O
            else if (BONDED(28,6))  base = 1.88f;   // Ni–C
            else if (BONDED(29,6))  base = 1.90f;   // Cu–C
            else if (BONDED(46,6))  base = 2.05f;   // Pd–C
            else if (BONDED(78,6))  base = 2.00f;   // Pt–C
            else if (BONDED(47,1)) base = 1.62f;   // Ag–H
            else if (BONDED(48,1)) base = 1.66f;   // Cd–H
            else if (BONDED(79,1)) base = 1.52f;   // Au–H
            else if (BONDED(80,1)) base = 1.74f;   // Hg–H

            else if (BONDED(6,6))   base = 1.39f;

            // Hydrides

            if (BONDED(3,1))  base = 1.60f;   // Li–H
            else if (BONDED(4,1))  base = 1.34f;   // Be–H (in BeH2 polymeric ~1.34 Å)
            else if (BONDED(11,1)) base = 1.96f;   // Na–H
            else if (BONDED(12,1)) base = 1.89f;   // Mg–H
            else if (BONDED(19,1)) base = 2.32f;   // K–H
            else if (BONDED(20,1)) base = 2.04f;   // Ca–H
            else if (BONDED(21,1)) base = 1.92f;   // Sc–H 
            else if (BONDED(22,1)) base = 1.70f;   // Ti–H
            else if (BONDED(23,1)) base = 1.68f;   // V–H
            else if (BONDED(24,1)) base = 1.66f;   // Cr–H
            else if (BONDED(25,1)) base = 1.62f;   // Mn–H
            else if (BONDED(26,1)) base = 1.60f;   // Fe–H
            else if (BONDED(27,1)) base = 1.58f;   // Co–H
            else if (BONDED(28,1)) base = 1.56f;   // Ni–H
            else if (BONDED(29,1)) base = 1.55f;   // Cu–H
            else if (BONDED(30,1)) base = 1.60f;   // Zn–H
            else if (BONDED(37,1)) base = 2.45f;   // Rb–H
            else if (BONDED(38,1)) base = 2.20f;   // Sr–H
            else if (BONDED(39,1)) base = 2.05f;   // Y–H
            else if (BONDED(40,1)) base = 1.85f;   // Zr–H
            else if (BONDED(41,1)) base = 1.82f;   // Nb–H
            else if (BONDED(42,1)) base = 1.80f;   // Mo–H
            else if (BONDED(43,1)) base = 1.78f;   // Tc–H 
            else if (BONDED(44,1)) base = 1.75f;   // Ru–H
            else if (BONDED(45,1)) base = 1.72f;   // Rh–H
            else if (BONDED(46,1)) base = 1.70f;   // Pd–H
            else if (BONDED(47,1)) base = 1.62f;   // Ag–H
            else if (BONDED(48,1)) base = 1.66f;   // Cd–H
            else if (BONDED(55,1)) base = 2.60f;   // Cs–H
            else if (BONDED(56,1)) base = 2.35f;   // Ba–H
            else if (BONDED(57,1)) base = 2.06f;   // La–H
            else if (BONDED(58,1)) base = 2.04f;   // Ce–H
            else if (BONDED(59,1)) base = 2.02f;   // Pr–H
            else if (BONDED(60,1)) base = 2.01f;   // Nd–H

            if (BONDED(61,1)) base = 2.00f;   // Pm–H 
            else if (BONDED(62,1)) base = 1.99f;   // Sm–H
            else if (BONDED(63,1)) base = 1.98f;   // Eu–H
            else if (BONDED(64,1)) base = 1.97f;   // Gd–H
            else if (BONDED(65,1)) base = 1.96f;   // Tb–H
            else if (BONDED(66,1)) base = 1.95f;   // Dy–H
            else if (BONDED(67,1)) base = 1.94f;   // Ho–H
            else if (BONDED(68,1)) base = 1.93f;   // Er–H
            else if (BONDED(69,1)) base = 1.92f;   // Tm–H
            else if (BONDED(70,1)) base = 1.91f;   // Yb–H
            else if (BONDED(71,1)) base = 1.90f;   // Lu–H
            else if (BONDED(72,1)) base = 1.88f;   // Hf–H
            else if (BONDED(73,1)) base = 1.86f;   // Ta–H
            else if (BONDED(74,1)) base = 1.84f;   // W–H
            else if (BONDED(75,1)) base = 1.82f;   // Re–H
            else if (BONDED(76,1)) base = 1.80f;   // Os–H
            else if (BONDED(77,1)) base = 1.78f;   // Ir–H
            else if (BONDED(78,1)) base = 1.76f;   // Pt–H
            else if (BONDED(79,1)) base = 1.52f;   // Au–H
            else if (BONDED(80,1)) base = 1.74f;   // Hg–H
            else if (BONDED(81,1)) base = 1.70f;   // Tl–H
            else if (BONDED(82,1)) base = 1.90f;   // Pb–H
            else if (BONDED(83,1)) base = 1.85f;   // Bi–H
            else if (BONDED(84,1)) base = 1.80f;   // Po–H 
            else if (BONDED(85,1)) base = 1.75f;   // At–H 
            else if (BONDED(87,1)) base = 2.80f;   // Fr–H 
            else if (BONDED(88,1)) base = 2.50f;   // Ra–H 
            else if (BONDED(89,1)) base = 2.20f;   // Ac–H
            else if (BONDED(90,1)) base = 2.00f;   // Th–H
            else if (BONDED(92,1)) base = 1.95f;   // U–H

            // Lanthanides
            else if (BONDED(57,1)) base = 2.06f;   // La–H
            else if (BONDED(64,1)) base = 1.97f;   // Gd–H

            break;

        case sim::fun::BondType::DOUBLE:
            if      (BONDED(6,5))  base = 1.42f;   // C=B 
            else if (BONDED(6,6))  base = 1.34f;   // C=C
            else if (BONDED(6,7))  base = 1.31f;   // C=N
            else if (BONDED(6,8))  base = 1.22f;   // C=O
            else if (BONDED(6,15)) base = 1.67f;   // C=P
            else if (BONDED(6,16)) base = 1.61f;   // C=S
            else if (BONDED(7,7))  base = 1.24f;   // N=N
            else if (BONDED(7,8))  base = 1.28f;   // N=O
            else if (BONDED(7,15)) base = 1.54f;   // N=P
            else if (BONDED(8,8))  base = 1.21f;   // O=O
            else if (BONDED(8,15)) base = 1.48f;   // O=P
            else if (BONDED(8,16)) base = 1.43f;   // O=S
            else if (BONDED(15,15)) base = 1.89f;  // P=P
            else if (BONDED(15,16)) base = 1.95f;  // P=S
            break;

        case sim::fun::BondType::TRIPLE:
            if      (BONDED(6,6))  base = 1.20f;   // C≡C
            else if (BONDED(6,7))  base = 1.16f;   // C≡N
            else if (BONDED(6,15)) base = 1.50f;   // C≡P
            else if (BONDED(7,7))  base = 1.10f;   // N≡N
            else if (BONDED(7,15)) base = 1.49f;   // N≡P
            else if (BONDED(6,8))  base = 1.13f;   // C≡O (CO)
            break;

        case sim::fun::BondType::QUADRUPLE:
            base *= 0.55f;
            break;

        default:
            break;
        }

        #undef BONDED

        if (base != original)
            base *= MULT_FACTOR;

        return base;
    }

    inline float getBondEnergy(uint8_t Z1, uint8_t Z2, sim::fun::BondType type)
    {
        #define BONDED(a, b) ((Z1 == (a) && Z2 == (b)) || (Z1 == (b) && Z2 == (a)))

        static const std::map<std::pair<uint8_t, uint8_t>, float> singleBondEnergy = {
            // Carbon
            {{6,1},  413.0f},   // C–H
            {{6,6},  348.0f},   // C–C
            {{6,7},  305.0f},   // C–N
            {{6,8},  358.0f},   // C–O
            {{6,9},  485.0f},   // C–F
            {{6,15}, 272.0f},   // C–P
            {{6,16}, 272.0f},   // C–S
            {{6,17}, 327.0f},   // C–Cl

            // Nitrogen
            {{7,1},  391.0f},   // N–H
            {{7,7},  163.0f},   // N–N
            {{7,8},  201.0f},   // N–O

            // Oxygen
            {{8,1},  463.0f},   // O–H
            {{8,8},  146.0f},   // O–O

            // Halogens
            {{1,9},  565.0f},   // H–F
            {{1,17}, 431.0f},   // H–Cl
            {{9,9},  159.0f},   // F–F
            {{17,17},243.0f},   // Cl–Cl

            // Phosphorus
            {{15,1},  322.0f},   // P–H
            {{15,6},  264.0f},   // P–C
            {{15,7},  285.0f},   // P–N
            {{15,8},  360.0f},   // P–O
            {{15,15}, 201.0f},   // P–P
            {{15,9},  490.0f},   // P–F
            {{15,17}, 325.0f},   // P–Cl
            {{15,16}, 230.0f},   // P–S

            // Others
            {{14,8},  452.0f},   // Si–O
            {{16,16}, 226.0f},   // S–S
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
            break;

        case sim::fun::BondType::DOUBLE:

            if      (BONDED(6,6))  base = 614.0f;   // C=C
            else if (BONDED(6,7))  base = 607.0f;   // C=N
            else if (BONDED(6,8))  base = 799.0f;   // C=O
            else if (BONDED(7,8))  base = 607.0f;   // N=O
            else if (BONDED(8,8))  base = 498.0f;   // O=O
            else if (BONDED(15,8)) base = 545.0f;   // P=O
            else if (BONDED(6,16)) base = 536.0f;   // C=S
            else if (base > 0.0f)  base *= 1.65f;
            break;

        case sim::fun::BondType::TRIPLE:
            if      (BONDED(6,6))  base = 837.0f;   // C≡C
            else if (BONDED(6,7))  base = 887.0f;   // C≡N
            else if (BONDED(7,7))  base = 941.0f;   // N≡N
            else if (base > 0.0f)  base *= 2.40f;
            break;

        case sim::fun::BondType::QUADRUPLE:
            if (base > 0.0f) base *= 3.2f;
            break;

        default:
            break;
        }

        if (base <= 0.0f)
        {
            float DA = 0.0f, DB = 0.0f;

            if (BONDED(1,1))  DA = DB = 436.0f;
            else if (BONDED(6,6)) DA = DB = 348.0f;
            else if (BONDED(7,7)) DA = DB = 163.0f;
            else if (BONDED(8,8)) DA = DB = 146.0f;
            else if (BONDED(9,9)) DA = DB = 159.0f;
            else if (BONDED(15,15)) DA = DB = 201.0f;
            else if (BONDED(16,16)) DA = DB = 226.0f;
            else if (BONDED(17,17)) DA = DB = 243.0f;

            float chiA = getElectronegativity(Z1);
            float chiB = getElectronegativity(Z2);
            float ionic_term = (chiA > 0.0f && chiB > 0.0f) ? 96.5f * (chiA - chiB) * (chiA - chiB) : 0.0f;

            if (DA > 0.0f && DB > 0.0f)
                base = 0.5f * (DA + DB) + ionic_term;
            else
            {
                float r1 = covalent_radius[Z1];
                float r2 = covalent_radius[Z2];
                if (r1 > 0.0f && r2 > 0.0f)
                    base = 400.0f * std::exp(- (r1 + r2) / 2.0f);
                else
                    base = 250.0f;
            }
        }

        #undef BONDED

        return base;  // in kJ/mol
    }

    inline float getAngleHarmonicConstant(uint8_t ZA, uint8_t ZB, uint8_t ZC)
    {
        float K = 300.0f;

    #define ANY(a, b) ((ZA == (a) && ZC == (b)) || (ZA == (b) && ZC == (a)))

        if (ZB == 6)          // ───── CARBON ─────
        {
            if      (ANY(1, 1)) K = 310.0f;   // H-C-H        sp³
            else if (ANY(1, 6)) K = 360.0f;   // H-C-C        sp³
            else if (ANY(6, 6)) K = 520.0f;   // C-C-C        alkane
            else if (ANY(1, 7)) K = 350.0f;   // H-C-N
            else if (ANY(1, 8)) K = 380.0f;   // H-C-O
            else if (ANY(1, 9)) K = 400.0f;   // H-C-F
            else if (ANY(1,17)) K = 420.0f;   // H-C-Cl
            else if (ANY(6, 8)) K = 670.0f;   // C-C=O        carbonyl
            else if (ANY(8, 8)) K = 670.0f;   // O=C-O        carboxyl, ester
            else if (ANY(6, 7)) K = 550.0f;   // C-C-N
            else if (ANY(7, 8)) K = 650.0f;   // N-C=O        amide (not peptide)
            else if (ANY(6,16)) K = 580.0f;   // C-C-S
            else if (ANY(7, 7)) K = 600.0f;   // N-C-N        guanidinium-like
        }
        else if (ZB == 7)     // ───── NITROGEN ─────
        {
            if      (ANY(1, 1)) K = 370.0f;   // H-N-H        ammonia, amines
            else if (ANY(1, 6)) K = 380.0f;   // H-N-C
            else if (ANY(6, 6)) K = 460.0f;   // C-N-C        tertiary amine
            else if (ANY(6, 8)) K = 850.0f;   // C-N-C=O      PEPTIDE BOND
            else if (ANY(1, 8)) K = 400.0f;   // H-N-O
            else if (ANY(6, 7)) K = 700.0f;   // C-N-C        imine, guanidine
            else if (ANY(7, 7)) K = 750.0f;   // N-N-N        (azides, rare)
            else if (ANY(6, 7)) K = 680.0f;   // C-N=C        imidazole, His
            else if (ANY(7, 6)) K = 720.0f;   // N=C-N        imidazole ring, arginine
        }
        else if (ZB == 8)     // ───── OXYGEN ─────
        {
            if      (ANY(1, 1)) K = 460.0f;   // H-O-H        water
            else if (ANY(1, 6)) K = 460.0f;   // H-O-C        alcohols
            else if (ANY(6, 6)) K = 545.0f;   // C-O-C        ethers
            else if (ANY(1,15)) K = 450.0f;   // H-O-P        phosphoric acid
            else if (ANY(8,16)) K = 700.0f;   // O-S-O        sulfate (extra stiff)
        }
        else if (ZB == 15)  // Phosphorus central atom
        {
            K = 520.0f;

            if (ANY(1,1))     K = 460.0f;   // H–P–H    (phosphine PH₃)
            if (ANY(1,6))     K = 480.0f;   // H–P–C
            if (ANY(6,6))     K = 520.0f;   // C–P–C    (typical alkylphosphines)
            if (ANY(6,7))     K = 550.0f;   // C–P–N
            if (ANY(1,17))    K = 520.0f;   // H–P–Cl

            if (ANY(8,8))     K = 720.0f;   // O–P–O    (phosphate esters, phosphoric acid)
            if (ANY(8,6))     K = 680.0f;   // O–P–C    (phosphonates)
            if (ANY(8,7))     K = 700.0f;   // O–P–N    (phosphoramidates)
            if (ANY(8,9))     K = 780.0f;   // O–P–F    (phosphoryl fluorides)
            if (ANY(8,17))    K = 720.0f;   // O–P–Cl   (POCl₃, phosphorochloridates)

            if (ANY(8,8) && (ZA == 8 || ZC == 8))
                K = 820.0f;
        }
        else if (ZB == 16)    // ───── SULFUR ─────
        {
            if      (ANY(1, 1)) K = 420.0f;   // H-S-H        thiols
            else if (ANY(6, 6)) K = 545.0f;   // C-S-C        disulfide, Met
            else if (ANY(8, 8)) K = 670.0f;   // O=S=O        sulfate
            else if (ANY(6, 8)) K = 620.0f;   // C-S=O
            else if (ANY(16,16)) K = 800.0f;  // S-S-S        (elemental sulfur, rare)
        }
        else if (ZB == 14)    // Silicon
        {
            if      (ANY(6, 6)) K = 450.0f;
            else if (ANY(8, 8)) K = 550.0f;
            else if (ANY(1, 8)) K = 480.0f;
        }
        else if (ZB == 5)     // Boron
        {
            if (ANY(8, 8)) K = 600.0f;   // O-B-O
            else K = 550.0f;
        }

        else if (ZB == 6 && ANY(7, 7)) K = 700.0f;  // N=C(N)
        else if (ZB == 7 && ANY(6, 6)) K = 720.0f;  // N=C-N
        else if (ZB == 7 && ANY(6, 8)) K = 900.0f;  // C=N-C=O

        // Halogens
        else if (ZB == 9  || ZB == 17 || ZB == 35 || ZB == 53)
            K = 600.0f;

        // Metals
        else if (ZB == 12 || ZB == 20 || ZB == 26 || ZB == 29 || ZB == 30) // Mg, Ca, Fe, Cu, Zn
            K = 350.0f;

    #undef ANY

        return K * 1000.f;
    }

    inline float getBondHarmonicConstantFromEnergy(uint8_t Z1, uint8_t Z2, sim::fun::BondType type)
    {
        float D = getBondEnergy(Z1, Z2, type);        // kJ/mol
        float r0 = getBondLength(Z1, Z2, type);       // Å

        float D_aJ = D * 1000.0f / AVOGADRO / 4.184f;
        float k = (D / 348.0f) * BOND_K;

        return k * static_cast<float>(type);
    }

    inline float getBondHarmonicConstantFromEnergy(float D, float r0, sim::fun::BondType type)
    {
        float D_aJ = D * 1000.0f / AVOGADRO / 4.184f;
        float k = (D / 348.0f) * BOND_K;

        return k * static_cast<float>(type);
    }
}; // namespace constants
