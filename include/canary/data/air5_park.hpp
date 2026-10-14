#ifndef CANARY_DATA_AIR5_PARK_HPP
#define CANARY_DATA_AIR5_PARK_HPP

// Generated from database/data/air5_park.json by
// tools/species_database/reactions_to_header.py; do not edit.
// Reference: C. Park, Review of chemical kinetic problems of future NASA
// missions, I: Earth entries, J. Thermophys. Heat Transf. 7(3), 1993, reduced
// to five species; from GDTk (gdtk-uq/gdtk bf559c3f)
// examples/lmr/3D/cylinder-larc/rm-park-5sp-2T.inp

#include <array>
#include <concepts>
#include <cstddef>
#include <canary/data/air5.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/species.hpp>
#include <canary/units.hpp>
#include <string_view>

namespace canary
{

// The 5 reactions of air5_park.json in file order, in SI units per mole,
// with the temperature limits of rate evaluation in K.
inline constexpr MechanismInput<5, 5> air5_park {
  .temperature_lower = 300.0,
  .temperature_upper = 50000.0,
  .reactions =
    {
      {
        {
          .equation = "N2 + M <=> N + N + M",
          .reactants =
            {
              .count = 1,
              .entries =
                {
                  {
                    {
                      .name = "N2",
                      .value = 1.0,
                    },
                  },
                },
            },
          .products =
            {
              .count = 1,
              .entries =
                {
                  {
                    {
                      .name = "N",
                      .value = 2.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 7000000000000000.0,
              .n = -1.6,
              .activation_temperature = 113200.0,
              .translational_exponent = 0.5,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .third_body = true,
          .efficiencies =
            {
              .count = 5,
              .entries =
                {
                  {
                    {
                      .name = "N2",
                      .value = 1.0,
                    },
                    {
                      .name = "O2",
                      .value = 1.0,
                    },
                    {
                      .name = "NO",
                      .value = 1.0,
                    },
                    {
                      .name = "N",
                      .value = 4.2857,
                    },
                    {
                      .name = "O",
                      .value = 4.2857,
                    },
                  },
                },
            },
        },
        {
          .equation = "O2 + M <=> O + O + M",
          .reactants =
            {
              .count = 1,
              .entries =
                {
                  {
                    {
                      .name = "O2",
                      .value = 1.0,
                    },
                  },
                },
            },
          .products =
            {
              .count = 1,
              .entries =
                {
                  {
                    {
                      .name = "O",
                      .value = 2.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 2000000000000000.0,
              .n = -1.5,
              .activation_temperature = 59500.0,
              .translational_exponent = 0.5,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .third_body = true,
          .efficiencies =
            {
              .count = 5,
              .entries =
                {
                  {
                    {
                      .name = "N2",
                      .value = 1.0,
                    },
                    {
                      .name = "O2",
                      .value = 1.0,
                    },
                    {
                      .name = "NO",
                      .value = 1.0,
                    },
                    {
                      .name = "N",
                      .value = 5.0,
                    },
                    {
                      .name = "O",
                      .value = 5.0,
                    },
                  },
                },
            },
        },
        {
          .equation = "NO + M <=> N + O + M",
          .reactants =
            {
              .count = 1,
              .entries =
                {
                  {
                    {
                      .name = "NO",
                      .value = 1.0,
                    },
                  },
                },
            },
          .products =
            {
              .count = 2,
              .entries =
                {
                  {
                    {
                      .name = "N",
                      .value = 1.0,
                    },
                    {
                      .name = "O",
                      .value = 1.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 5000000000.0,
              .n = 0.0,
              .activation_temperature = 75500.0,
              .translational_exponent = 0.5,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .third_body = true,
          .efficiencies =
            {
              .count = 5,
              .entries =
                {
                  {
                    {
                      .name = "N2",
                      .value = 1.0,
                    },
                    {
                      .name = "O2",
                      .value = 1.0,
                    },
                    {
                      .name = "NO",
                      .value = 1.0,
                    },
                    {
                      .name = "N",
                      .value = 20.0,
                    },
                    {
                      .name = "O",
                      .value = 20.0,
                    },
                  },
                },
            },
        },
        {
          .equation = "NO + O <=> O2 + N",
          .reactants =
            {
              .count = 2,
              .entries =
                {
                  {
                    {
                      .name = "NO",
                      .value = 1.0,
                    },
                    {
                      .name = "O",
                      .value = 1.0,
                    },
                  },
                },
            },
          .products =
            {
              .count = 2,
              .entries =
                {
                  {
                    {
                      .name = "O2",
                      .value = 1.0,
                    },
                    {
                      .name = "N",
                      .value = 1.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 8400000.0,
              .n = 0.0,
              .activation_temperature = 19450.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .third_body = false,
          .efficiencies =
            {
              .count = 0,
              .entries =
                {
                  {},
                },
            },
        },
        {
          .equation = "N2 + O <=> NO + N",
          .reactants =
            {
              .count = 2,
              .entries =
                {
                  {
                    {
                      .name = "N2",
                      .value = 1.0,
                    },
                    {
                      .name = "O",
                      .value = 1.0,
                    },
                  },
                },
            },
          .products =
            {
              .count = 2,
              .entries =
                {
                  {
                    {
                      .name = "NO",
                      .value = 1.0,
                    },
                    {
                      .name = "N",
                      .value = 1.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 640000000000.0,
              .n = -1.0,
              .activation_temperature = 38400.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .third_body = false,
          .efficiencies =
            {
              .count = 0,
              .entries =
                {
                  {},
                },
            },
        },
      },
    },
};

// Mechanism of air5_park over the species of air5_species.
template <std::floating_point P, UnitSystem U>
using Air5Park = ReactionMechanism<P, U, air5_species.size(), 5>;

namespace detail
{
inline constexpr auto air5_park_species_names = [] {
    std::array<std::string_view, air5_species.size()> names {};
    for (const SpeciesIndex i : species_indices(0, names.size())) {
        names[i.value()] = air5_species[i.value()].name;
    }
    return names;
}();
}  // namespace detail

// Air5Park model of air5_park for the species of air5_species, in
// their order, with their species data, in the unit system of the scales,
// with standard pressure p0. See make_reaction_mechanism.
template <std::floating_point P, UnitSystem U>
Air5Park<P, U> make_air5_park(const SpeciesData<P, U, air5_species.size()>& species,
                              const UnitScales<P, U>& scales,
                              Pressure<P, U> standard_pressure)
{
    return make_reaction_mechanism(air5_park,
                                   detail::air5_park_species_names,
                                   species,
                                   scales,
                                   standard_pressure);
}

// As above, with the standard pressure 1 bar (1e5 Pa).
template <std::floating_point P, UnitSystem U>
Air5Park<P, U> make_air5_park(const SpeciesData<P, U, air5_species.size()>& species,
                              const UnitScales<P, U>& scales)
{ return make_reaction_mechanism(air5_park, detail::air5_park_species_names, species, scales); }

// As above for the thermo model gas with several energy modes: rates'
// other temperatures resolve to its modes. See make_reaction_mechanism.
template <std::floating_point P, UnitSystem U, typename G>
Air5Park<P, U> make_air5_park(const SpeciesData<P, U, air5_species.size()>& species,
                              const UnitScales<P, U>& scales,
                              Pressure<P, U> standard_pressure,
                              const G& gas)
{
    return make_reaction_mechanism(air5_park,
                                   detail::air5_park_species_names,
                                   species,
                                   scales,
                                   standard_pressure,
                                   gas);
}

}  // namespace canary

#endif  // CANARY_DATA_AIR5_PARK_HPP
