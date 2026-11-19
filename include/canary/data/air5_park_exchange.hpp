#ifndef CANARY_DATA_AIR5_PARK_EXCHANGE_HPP
#define CANARY_DATA_AIR5_PARK_EXCHANGE_HPP

// Generated from database/data/air5_park_exchange.json by
// tools/species_database/exchange_to_header.py; do not edit.
// Reference: V-T: Landau-Teller with Millikan-White a, b from C. Park, J.
// Thermophys. Heat Transf. 7(3), 1993, Table 1, and Park's high-temperature
// correction (ibid. p. 387); C-V: impartial dissociation (Treanor and Marrone,
// Phys. Fluids 5, 1962; Knab, Fruhauf and Messerschmid, J. Thermophys. Heat
// Transf. 9(2), 1995). From GDTk (gdtk-uq/gdtk bf559c3f)
// examples/lmr/3D/cylinder-larc/eem-park-5sp-2T.inp via prep-kinetics with
// gm-air-5sp-2T and rm-park-5sp-2T.

#include <array>
#include <canary/data/air5.hpp>
#include <canary/data/air5_park.hpp>
#include <canary/energy_exchange.hpp>
#include <canary/species.hpp>
#include <canary/units.hpp>
#include <concepts>
#include <cstddef>
#include <string_view>

namespace canary
{

// The 15 vibration-translation pairs and 3 chemistry-vibration
// couplings of air5_park_exchange.json in file order, in SI units per mole.
inline constexpr EnergyExchangeInput<15, 3> air5_park_exchange {
  .vibration_translation =
    {
      {
        {
          .species = "N2",
          .partner = "N",
          .millikan_white =
            {
              .a = 180.0,
              .b = 0.0262,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
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
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "N2",
          .partner = "N2",
          .millikan_white =
            {
              .a = 221.0,
              .b = 0.029,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "N2",
          .partner = "O2",
          .millikan_white =
            {
              .a = 229.0,
              .b = 0.0295,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "N2",
          .partner = "NO",
          .millikan_white =
            {
              .a = 225.0,
              .b = 0.0293,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
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
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "O2",
          .partner = "O",
          .millikan_white =
            {
              .a = 47.7,
              .b = 0.059,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "O2",
          .partner = "N2",
          .millikan_white =
            {
              .a = 134.0,
              .b = 0.0295,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "O2",
          .partner = "O2",
          .millikan_white =
            {
              .a = 138.0,
              .b = 0.03,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
            },
        },
        {
          .species = "O2",
          .partner = "NO",
          .millikan_white =
            {
              .a = 136.0,
              .b = 0.0298,
            },
          .correction =
            {
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
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
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
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
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
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
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
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
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
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
              .model = RelaxationCorrection::park_1993,
              .cross_section = 3e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
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
          .vibrational_temperature = 3393.44,
        },
        {
          .reaction = 1,
          .species = "O2",
          .dissociation_energy = 498350.006,
          .vibrational_temperature = 2273.53,
        },
        {
          .reaction = 2,
          .species = "NO",
          .dissociation_energy = 630583.693,
          .vibrational_temperature = 2739.7,
        },
      },
    },
};

// Energy exchange of air5_park_exchange over the species of air5_species.
template <std::floating_point P, UnitSystem U>
using Air5ParkExchange = EnergyExchange<P, U, air5_species.size(), 15, 3>;

namespace detail
{
inline constexpr auto air5_park_exchange_species_names = [] {
    std::array<std::string_view, air5_species.size()> names {};
    for (const SpeciesIndex i : species_indices(0, names.size())) {
        names[i.value()] = air5_species[i.value()].name;
    }
    return names;
}();
}  // namespace detail

// Air5ParkExchange model of air5_park_exchange for the species of air5_species, in
// their order, with their species data, and the reactions of mechanism, in the
// unit system of the scales. See make_energy_exchange.
template <std::floating_point P, UnitSystem U>
Air5ParkExchange<P, U>
  make_air5_park_exchange(const SpeciesData<P, U, air5_species.size()>& species,
                          const Air5Park<P, U>& mechanism,
                          const UnitScales<P, U>& scales)
{
    return make_energy_exchange(air5_park_exchange,
                                detail::air5_park_exchange_species_names,
                                species,
                                mechanism,
                                scales);
}

}  // namespace canary

#endif  // CANARY_DATA_AIR5_PARK_EXCHANGE_HPP
