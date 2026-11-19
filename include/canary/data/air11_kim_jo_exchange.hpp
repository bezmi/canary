#ifndef CANARY_DATA_AIR11_KIM_JO_EXCHANGE_HPP
#define CANARY_DATA_AIR11_KIM_JO_EXCHANGE_HPP

// Generated from database/data/air11_kim_jo_exchange.json by
// tools/species_database/exchange_to_header.py; do not edit.
// Reference: V-T: Landau-Teller with Millikan-White constants and limiting
// cross-sections from J. G. Kim and S. M. Jo, Modification of chemical-kinetic
// parameters for 11-air species in re-entry flows, Int. J. Heat Mass Transf.
// 169, 2021, Table 1; E-T: elastic electron-heavy exchange with the
// cross-section fits of Gnoffo, Gupta and Shinn (NASA TP-2867, 1989) and
// Coulomb cross-sections; V-E: Landau-Teller rate with the relaxation time of
// J.-H. Lee, AIAA Paper 85-1035, 1985, Eq. (19) with the rate coefficients of
// its Table 2, fitted by scripts/electron_vibration_relaxation/fit.py; C-V:
// impartial dissociation (Treanor and Marrone, Phys. Fluids 5, 1962; Knab,
// Fruhauf and Messerschmid, J. Thermophys. Heat Transf. 9(2), 1995), with the
// dissociation energies 2 Hf(N), 2 Hf(O) and Hf(N) + Hf(O) - Hf(NO) from the
// heats of formation of the NASA CEA thermo.inp and the vibrational
// temperatures omega_e h c / k_B with omega_e of the NIST Chemistry WebBook
// (Huber and Herzberg); C-E: electron-impact ionisation (Gnoffo, Gupta and
// Shinn 1989) for the two electron-impact ionisation reactions of the
// mechanism, with the ionisation energies of the species file

#include <array>
#include <canary/data/air11.hpp>
#include <canary/data/air11_kim_jo.hpp>
#include <canary/energy_exchange.hpp>
#include <canary/species.hpp>
#include <canary/units.hpp>
#include <concepts>
#include <cstddef>
#include <string_view>

namespace canary
{

// The 30 vibration-translation pairs and 21 chemistry-vibration
// couplings of air11_kim_jo_exchange.json in file order, in SI units per mole.
inline constexpr EnergyExchangeInput<30, 21> air11_kim_jo_exchange {
  .vibration_translation =
    {
      {
        {
          .species = "N2",
          .partner = "N2",
          .millikan_white =
            {
              .a = 221.53,
              .b = 0.029,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 1.7e-23,
              .reference_temperature = 0.0,
              .exponent = 0.55,
            },
        },
        {
          .species = "N2",
          .partner = "O2",
          .millikan_white =
            {
              .a = 228.76,
              .b = 0.0295,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 4.5e-23,
              .reference_temperature = 0.0,
              .exponent = 0.5,
            },
        },
        {
          .species = "N2",
          .partner = "NO",
          .millikan_white =
            {
              .a = 225.3,
              .b = 0.0293,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 4.5e-23,
              .reference_temperature = 0.0,
              .exponent = 0.5,
            },
        },
        {
          .species = "N2",
          .partner = "N",
          .millikan_white =
            {
              .a = 184.44,
              .b = 0.0389,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 1e-22,
              .reference_temperature = 0.0,
              .exponent = 0.52,
            },
        },
        {
          .species = "N2",
          .partner = "O",
          .millikan_white =
            {
              .a = 72.4,
              .b = 0.015,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 4.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.3,
            },
        },
        {
          .species = "N2",
          .partner = "N2+",
          .millikan_white =
            {
              .a = 221.53,
              .b = 0.029,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 1.7e-23,
              .reference_temperature = 0.0,
              .exponent = 0.55,
            },
        },
        {
          .species = "N2",
          .partner = "O2+",
          .millikan_white =
            {
              .a = 228.76,
              .b = 0.0295,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 4.5e-23,
              .reference_temperature = 0.0,
              .exponent = 0.5,
            },
        },
        {
          .species = "N2",
          .partner = "NO+",
          .millikan_white =
            {
              .a = 225.3,
              .b = 0.0293,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 4.5e-23,
              .reference_temperature = 0.0,
              .exponent = 0.5,
            },
        },
        {
          .species = "N2",
          .partner = "N+",
          .millikan_white =
            {
              .a = 180.0,
              .b = 0.0262,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 4.5e-23,
              .reference_temperature = 0.0,
              .exponent = 0.5,
            },
        },
        {
          .species = "N2",
          .partner = "O+",
          .millikan_white =
            {
              .a = 188.89,
              .b = 0.0268,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 4.5e-23,
              .reference_temperature = 0.0,
              .exponent = 0.5,
            },
        },
        {
          .species = "O2",
          .partner = "N2",
          .millikan_white =
            {
              .a = 131.32,
              .b = 0.0295,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 7.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.28,
            },
        },
        {
          .species = "O2",
          .partner = "O2",
          .millikan_white =
            {
              .a = 130.15,
              .b = 0.0303,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 5e-24,
              .reference_temperature = 0.0,
              .exponent = 0.7,
            },
        },
        {
          .species = "O2",
          .partner = "NO",
          .millikan_white =
            {
              .a = 133.71,
              .b = 0.0298,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 7.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.28,
            },
        },
        {
          .species = "O2",
          .partner = "N",
          .millikan_white =
            {
              .a = 72.4,
              .b = 0.015,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 7.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.28,
            },
        },
        {
          .species = "O2",
          .partner = "O",
          .millikan_white =
            {
              .a = -27.32,
              .b = 0.1465,
            },
          .correction =
            {
              .model = RelaxationCorrection::none,
              .cross_section = 0.0,
              .reference_temperature = 0.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "O2",
          .partner = "N2+",
          .millikan_white =
            {
              .a = 131.32,
              .b = 0.0295,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 7.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.28,
            },
        },
        {
          .species = "O2",
          .partner = "O2+",
          .millikan_white =
            {
              .a = 130.15,
              .b = 0.0303,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 5e-24,
              .reference_temperature = 0.0,
              .exponent = 0.7,
            },
        },
        {
          .species = "O2",
          .partner = "NO+",
          .millikan_white =
            {
              .a = 133.71,
              .b = 0.0298,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 7.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.28,
            },
        },
        {
          .species = "O2",
          .partner = "N+",
          .millikan_white =
            {
              .a = 106.06,
              .b = 0.0265,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 7.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.28,
            },
        },
        {
          .species = "O2",
          .partner = "O+",
          .millikan_white =
            {
              .a = 110.97,
              .b = 0.0271,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 7.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.28,
            },
        },
        {
          .species = "NO",
          .partner = "N2",
          .millikan_white =
            {
              .a = 49.5,
              .b = 0.042,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1985,
              .cross_section = 5.8e-20,
              .reference_temperature = 0.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "NO",
          .partner = "O2",
          .millikan_white =
            {
              .a = 49.5,
              .b = 0.042,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1985,
              .cross_section = 5.8e-20,
              .reference_temperature = 0.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "NO",
          .partner = "NO",
          .millikan_white =
            {
              .a = 49.5,
              .b = 0.042,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1985,
              .cross_section = 5.8e-20,
              .reference_temperature = 0.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "NO",
          .partner = "N",
          .millikan_white =
            {
              .a = 49.5,
              .b = 0.042,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1985,
              .cross_section = 5.8e-20,
              .reference_temperature = 0.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "NO",
          .partner = "O",
          .millikan_white =
            {
              .a = 49.5,
              .b = 0.042,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1985,
              .cross_section = 5.8e-20,
              .reference_temperature = 0.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "NO",
          .partner = "N2+",
          .millikan_white =
            {
              .a = 175.67,
              .b = 0.0293,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 8.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.35,
            },
        },
        {
          .species = "NO",
          .partner = "O2+",
          .millikan_white =
            {
              .a = 181.6,
              .b = 0.0298,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 8.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.35,
            },
        },
        {
          .species = "NO",
          .partner = "NO+",
          .millikan_white =
            {
              .a = 178.76,
              .b = 0.0295,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 8.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.35,
            },
        },
        {
          .species = "NO",
          .partner = "N+",
          .millikan_white =
            {
              .a = 142.62,
              .b = 0.0264,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 8.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.35,
            },
        },
        {
          .species = "NO",
          .partner = "O+",
          .millikan_white =
            {
              .a = 149.08,
              .b = 0.027,
            },
          .correction =
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 8.5e-22,
              .reference_temperature = 0.0,
              .exponent = 0.35,
            },
        },
      },
    },
  .chemistry_vibration =
    {
      {
        {
          .reaction = 0,
          .species = "N2",
          .dissociation_energy = 945360.0,
          .vibrational_temperature = 3393.4559799744534,
        },
        {
          .reaction = 1,
          .species = "N2",
          .dissociation_energy = 945360.0,
          .vibrational_temperature = 3393.4559799744534,
        },
        {
          .reaction = 2,
          .species = "N2",
          .dissociation_energy = 945360.0,
          .vibrational_temperature = 3393.4559799744534,
        },
        {
          .reaction = 3,
          .species = "N2",
          .dissociation_energy = 945360.0,
          .vibrational_temperature = 3393.4559799744534,
        },
        {
          .reaction = 4,
          .species = "N2",
          .dissociation_energy = 945360.0,
          .vibrational_temperature = 3393.4559799744534,
        },
        {
          .reaction = 5,
          .species = "N2",
          .dissociation_energy = 945360.0,
          .vibrational_temperature = 3393.4559799744534,
        },
        {
          .reaction = 6,
          .species = "N2",
          .dissociation_energy = 945360.0,
          .vibrational_temperature = 3393.4559799744534,
        },
        {
          .reaction = 7,
          .species = "N2",
          .dissociation_energy = 945360.0,
          .vibrational_temperature = 3393.4559799744534,
        },
        {
          .reaction = 8,
          .species = "N2",
          .dissociation_energy = 945360.0,
          .vibrational_temperature = 3393.4559799744534,
        },
        {
          .reaction = 9,
          .species = "N2",
          .dissociation_energy = 945360.0,
          .vibrational_temperature = 3393.4559799744534,
        },
        {
          .reaction = 15,
          .species = "O2",
          .dissociation_energy = 498350.006,
          .vibrational_temperature = 2273.5451503935737,
        },
        {
          .reaction = 16,
          .species = "O2",
          .dissociation_energy = 498350.006,
          .vibrational_temperature = 2273.5451503935737,
        },
        {
          .reaction = 17,
          .species = "O2",
          .dissociation_energy = 498350.006,
          .vibrational_temperature = 2273.5451503935737,
        },
        {
          .reaction = 18,
          .species = "O2",
          .dissociation_energy = 498350.006,
          .vibrational_temperature = 2273.5451503935737,
        },
        {
          .reaction = 19,
          .species = "O2",
          .dissociation_energy = 498350.006,
          .vibrational_temperature = 2273.5451503935737,
        },
        {
          .reaction = 10,
          .species = "O2",
          .dissociation_energy = 498350.006,
          .vibrational_temperature = 2273.5451503935737,
        },
        {
          .reaction = 11,
          .species = "O2",
          .dissociation_energy = 498350.006,
          .vibrational_temperature = 2273.5451503935737,
        },
        {
          .reaction = 12,
          .species = "O2",
          .dissociation_energy = 498350.006,
          .vibrational_temperature = 2273.5451503935737,
        },
        {
          .reaction = 13,
          .species = "O2",
          .dissociation_energy = 498350.006,
          .vibrational_temperature = 2273.5451503935737,
        },
        {
          .reaction = 14,
          .species = "O2",
          .dissociation_energy = 498350.006,
          .vibrational_temperature = 2273.5451503935737,
        },
        {
          .reaction = 20,
          .species = "NO",
          .dissociation_energy = 630583.693,
          .vibrational_temperature = 2739.48872584259,
        },
      },
    },
};

// The electron exchanges of air11_kim_jo_exchange.json: 10 elastic
// electron-heavy pairs, 1 vibration-electron exchanges and 2
// electron-impact ionisations (ionisation energies in J/mol from the
// species file).
inline constexpr ElectronExchangeInput<10, 1, 2> air11_kim_jo_exchange_electrons {
  .electron_translation =
    {
      {
        {
          .partner = "N2",
          .cross_section = ElectronCrossSection::gnoffo_neutral,
          .a = 7.5e-20,
          .b = 5.5e-24,
          .c = -1e-28,
          .maximum_temperature = 30000.0,
        },
        {
          .partner = "O2",
          .cross_section = ElectronCrossSection::gnoffo_neutral,
          .a = 2e-20,
          .b = 6e-24,
          .c = 0.0,
          .maximum_temperature = 30000.0,
        },
        {
          .partner = "N",
          .cross_section = ElectronCrossSection::gnoffo_neutral,
          .a = 5e-20,
          .b = 0.0,
          .c = 0.0,
          .maximum_temperature = 30000.0,
        },
        {
          .partner = "O",
          .cross_section = ElectronCrossSection::gnoffo_neutral,
          .a = 1.2e-20,
          .b = 1.7e-24,
          .c = -2e-29,
          .maximum_temperature = 30000.0,
        },
        {
          .partner = "NO",
          .cross_section = ElectronCrossSection::gnoffo_neutral,
          .a = 1e-19,
          .b = 0.0,
          .c = 0.0,
          .maximum_temperature = 30000.0,
        },
        {
          .partner = "N+",
          .cross_section = ElectronCrossSection::coulomb,
          .a = 0.0,
          .b = 0.0,
          .c = 0.0,
          .maximum_temperature = 0.0,
        },
        {
          .partner = "O+",
          .cross_section = ElectronCrossSection::coulomb,
          .a = 0.0,
          .b = 0.0,
          .c = 0.0,
          .maximum_temperature = 0.0,
        },
        {
          .partner = "N2+",
          .cross_section = ElectronCrossSection::coulomb,
          .a = 0.0,
          .b = 0.0,
          .c = 0.0,
          .maximum_temperature = 0.0,
        },
        {
          .partner = "O2+",
          .cross_section = ElectronCrossSection::coulomb,
          .a = 0.0,
          .b = 0.0,
          .c = 0.0,
          .maximum_temperature = 0.0,
        },
        {
          .partner = "NO+",
          .cross_section = ElectronCrossSection::coulomb,
          .a = 0.0,
          .b = 0.0,
          .c = 0.0,
          .maximum_temperature = 0.0,
        },
      },
    },
  .vibration_electron =
    {
      {
        {
          .species = "N2",
          .relaxation_time =
            {
              .reference_temperature = 10000.0,
              .constant = -24.91199,
              .exponent = 4.475475,
              .inverse_coefficients =
                {
                  {
                    3.192608,
                    -0.08318516,
                    0.002015418,
                  },
                },
            },
        },
      },
    },
  .electron_impact_ionisation =
    {
      {
        {
          .reaction = 38,
          .species = "O",
          .ionisation_energy = 1313942.5595485026,
        },
        {
          .reaction = 39,
          .species = "N",
          .ionisation_energy = 1402330.360173364,
        },
      },
    },
};

// Energy exchange of air11_kim_jo_exchange over the species and modes of Air11ThreeTemperature.
template <std::floating_point P, UnitSystem U>
using Air11KimJoExchange = MultiTemperatureExchange<P,
                                                    U,
                                                    air11_species.size(),
                                                    Air11ThreeTemperature<P, U>::mode_count,
                                                    30,
                                                    21,
                                                    10,
                                                    1,
                                                    2>;

namespace detail
{
inline constexpr auto air11_kim_jo_exchange_species_names = [] {
    std::array<std::string_view, air11_species.size()> names {};
    for (const SpeciesIndex i : species_indices(0, names.size())) {
        names[i.value()] = air11_species[i.value()].name;
    }
    return names;
}();
}  // namespace detail

// Air11KimJoExchange model of air11_kim_jo_exchange for the species of air11_species, in
// their order, with their species data, the reactions of mechanism and the
// modes of gas, in the unit system of the scales. See
// make_multi_temperature_exchange.
template <std::floating_point P, UnitSystem U>
Air11KimJoExchange<P, U>
  make_air11_kim_jo_exchange(const SpeciesData<P, U, air11_species.size()>& species,
                             const Air11KimJo<P, U>& mechanism,
                             const Air11ThreeTemperature<P, U>& gas,
                             const UnitScales<P, U>& scales)
{
    return make_multi_temperature_exchange(air11_kim_jo_exchange,
                                           air11_kim_jo_exchange_electrons,
                                           detail::air11_kim_jo_exchange_species_names,
                                           species,
                                           mechanism,
                                           gas,
                                           scales);
}

}  // namespace canary

#endif  // CANARY_DATA_AIR11_KIM_JO_EXCHANGE_HPP
