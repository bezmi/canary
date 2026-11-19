#ifndef CANARY_DATA_AIR5_KIM_JO_HPP
#define CANARY_DATA_AIR5_KIM_JO_HPP

// Generated from database/data/air5_kim_jo.json by
// tools/species_database/reactions_to_header.py; do not edit.
// Reference: Kim and Jo, Int. J. Heat Mass Transf. 169 (2021) (dissociation)
// and Park (1993) (NO exchange), as given in GDTk
// examples/lmr/2D/shock-tunnel-nozzles/five-species-air/rr-kim-air5-1T.inp
// (gdtk-uq/gdtk commit bf559c3f), converted to SI

#include <array>
#include <canary/data/air5.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/species.hpp>
#include <canary/units.hpp>
#include <concepts>
#include <cstddef>
#include <string_view>

namespace canary
{

// The 13 reactions of air5_kim_jo.json in file order, in SI units per mole,
// with the temperature limits of rate evaluation in K.
inline constexpr MechanismInput<13, 5> air5_kim_jo {
  .temperature_lower = 300.0,
  .temperature_upper = 30000.0,
  .reactions =
    {
      {
        {
          .equation = "N2 + N2 <=> N + N + N2",
          .reactants =
            {
              .count = 1,
              .entries =
                {
                  {
                    {
                      .name = "N2",
                      .value = 2.0,
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
                      .value = 2.0,
                    },
                    {
                      .name = "N2",
                      .value = 1.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 121600000000000.0,
              .n = -1.214,
              .activation_temperature = 113200.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .equilibrium_temperature = TranslationTemperatureInput {},
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
          .equation = "N2 + O2 <=> N + N + O2",
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
                      .name = "O2",
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
                      .value = 2.0,
                    },
                    {
                      .name = "O2",
                      .value = 1.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 7000000000000000.0,
              .n = -1.6,
              .activation_temperature = 113200.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .equilibrium_temperature = TranslationTemperatureInput {},
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
          .equation = "N2 + NO <=> N + N + NO",
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
                      .value = 2.0,
                    },
                    {
                      .name = "NO",
                      .value = 1.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 7000000000000000.0,
              .n = -1.6,
              .activation_temperature = 113200.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .equilibrium_temperature = TranslationTemperatureInput {},
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
          .equation = "N2 + N <=> N + N + N",
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
                      .name = "N",
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
                      .value = 3.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 359100000000000.0,
              .n = -1.226,
              .activation_temperature = 113200.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .equilibrium_temperature = TranslationTemperatureInput {},
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
          .equation = "N2 + O <=> N + N + O",
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
                      .name = "N",
                      .value = 2.0,
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
              .A = 3e+16,
              .n = -1.6,
              .activation_temperature = 113200.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .equilibrium_temperature = TranslationTemperatureInput {},
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
          .equation = "O2 + N2 <=> O + O + N2",
          .reactants =
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
                      .name = "N2",
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
                      .name = "O",
                      .value = 2.0,
                    },
                    {
                      .name = "N2",
                      .value = 1.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 3354000000.0,
              .n = -0.2726,
              .activation_temperature = 59500.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .equilibrium_temperature = TranslationTemperatureInput {},
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
          .equation = "O2 + O2 <=> O + O + O2",
          .reactants =
            {
              .count = 1,
              .entries =
                {
                  {
                    {
                      .name = "O2",
                      .value = 2.0,
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
                      .name = "O",
                      .value = 2.0,
                    },
                    {
                      .name = "O2",
                      .value = 1.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 1.117e+19,
              .n = -2.585,
              .activation_temperature = 59500.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .equilibrium_temperature = TranslationTemperatureInput {},
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
          .equation = "O2 + NO <=> O + O + NO",
          .reactants =
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
                      .name = "O",
                      .value = 2.0,
                    },
                    {
                      .name = "NO",
                      .value = 1.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 3354000000.0,
              .n = -0.2726,
              .activation_temperature = 59500.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .equilibrium_temperature = TranslationTemperatureInput {},
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
          .equation = "O2 + N <=> O + O + N",
          .reactants =
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
          .products =
            {
              .count = 2,
              .entries =
                {
                  {
                    {
                      .name = "O",
                      .value = 2.0,
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
              .A = 1e+16,
              .n = -1.5,
              .activation_temperature = 59500.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .equilibrium_temperature = TranslationTemperatureInput {},
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
          .equation = "O2 + O <=> O + O + O",
          .reactants =
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
                      .name = "O",
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
                      .value = 3.0,
                    },
                  },
                },
            },
          .forward =
            {
              .A = 3000000000000000.0,
              .n = -1.5,
              .activation_temperature = 59500.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .equilibrium_temperature = TranslationTemperatureInput {},
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
              .A = 1450000000.0,
              .n = 0.0,
              .activation_temperature = 75200.0,
            },
          .backward_model = BackwardRateModel::equilibrium,
          .backward =
            {
              .A = 0.0,
              .n = 0.0,
              .activation_temperature = 0.0,
            },
          .equilibrium_temperature = TranslationTemperatureInput {},
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
                      .value = 0.6648275862068966,
                    },
                    {
                      .name = "N",
                      .value = 0.6648275862068966,
                    },
                    {
                      .name = "O",
                      .value = 0.6648275862068966,
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
          .equilibrium_temperature = TranslationTemperatureInput {},
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
          .equilibrium_temperature = TranslationTemperatureInput {},
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

// Mechanism of air5_kim_jo over the species of air5_species.
template <std::floating_point P, UnitSystem U>
using Air5KimJo = ReactionMechanism<P, U, air5_species.size(), 13>;

namespace detail
{
inline constexpr auto air5_kim_jo_species_names = [] {
    std::array<std::string_view, air5_species.size()> names {};
    for (const SpeciesIndex i : species_indices(0, names.size())) {
        names[i.value()] = air5_species[i.value()].name;
    }
    return names;
}();
}  // namespace detail

// Air5KimJo model of air5_kim_jo for the species of air5_species, in
// their order, in the unit system of the scales,
// with standard pressure p0. See make_reaction_mechanism.
template <std::floating_point P, UnitSystem U>
Air5KimJo<P, U> make_air5_kim_jo(const UnitScales<P, U>& scales, Pressure<P, U> standard_pressure)
{
    return make_reaction_mechanism(air5_kim_jo,
                                   reaction_species(air5_species),
                                   scales,
                                   standard_pressure);
}

// As above, with the standard pressure 1 bar (1e5 Pa).
template <std::floating_point P, UnitSystem U>
Air5KimJo<P, U> make_air5_kim_jo(const UnitScales<P, U>& scales)
{ return make_reaction_mechanism(air5_kim_jo, reaction_species(air5_species), scales); }

// As above for the thermo model gas with several energy modes: rates'
// other temperatures resolve to its modes. See make_reaction_mechanism.
template <std::floating_point P, UnitSystem U, typename G>
Air5KimJo<P, U>
  make_air5_kim_jo(const UnitScales<P, U>& scales, Pressure<P, U> standard_pressure, const G& gas)
{
    return make_reaction_mechanism(air5_kim_jo,
                                   reaction_species(air5_species),
                                   scales,
                                   standard_pressure,
                                   gas);
}

}  // namespace canary

#endif  // CANARY_DATA_AIR5_KIM_JO_HPP
