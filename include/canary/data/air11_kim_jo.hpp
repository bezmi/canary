// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_DATA_AIR11_KIM_JO_HPP
#define CANARY_DATA_AIR11_KIM_JO_HPP

// Generated from database/data/air11_kim_jo.json by
// tools/species_database/reactions_to_header.py; do not edit.
// Reference: Dissociation by heavy particles: J. G. Kim and S. M. Jo,
// Modification of chemical-kinetic parameters for 11-air species in re-entry
// flows, Int. J. Heat Mass Transf. 169, 2021, Table 3; N2 + e-, exchange,
// associative ionisation, charge exchange and electron-impact ionisation: C.
// Park, Review of chemical-kinetic problems of future NASA missions, I: Earth
// entries, J. Thermophys. Heat Transf. 7(3), 1993, Table 2; rate coefficients
// converted from cm^3/(mol s) to SI; third-body efficiencies of NO + M are
// ratios of the tabulated rate coefficients. Controlling temperatures: sqrt(T
// T_v) for dissociation by heavy particles; sqrt(T T_e) for N2 + e- and T_e
// for electron-impact ionisation, with T_e the temperature of the model's
// electron mode. The cited papers have one temperature for vibration and
// electrons, for which these are their forms

#include <array>
#include <canary/data/air11.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/species.hpp>
#include <canary/units.hpp>
#include <concepts>
#include <cstddef>
#include <string_view>

namespace canary
{

// The 40 reactions of air11_kim_jo.json in file order, in SI units per mole,
// with the temperature limits of rate evaluation in K.
inline constexpr MechanismInput<40, 11>
  air11_kim_jo {
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"N2"},
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"N2"},
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"N2"},
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"N2"},
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"N2"},
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
            .equation = "N2 + N2+ <=> N + N + N2+",
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
                        .name = "N2+",
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
                        .name = "N2+",
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"N2"},
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
            .equation = "N2 + O2+ <=> N + N + O2+",
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
                        .name = "O2+",
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
                        .name = "O2+",
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"N2"},
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
            .equation = "N2 + NO+ <=> N + N + NO+",
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
                        .name = "NO+",
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
                        .name = "NO+",
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"N2"},
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
            .equation = "N2 + N+ <=> N + N + N+",
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
                        .name = "N+",
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
                        .name = "N+",
                        .value = 1.0,
                      },
                    },
                  },
              },
            .forward =
              {
                .A = 359100000000000.0,
                .n = -1.226,
                .activation_temperature = 113200.0,
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"N2"},
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
            .equation = "N2 + O+ <=> N + N + O+",
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
                        .name = "O+",
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
                        .name = "O+",
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"N2"},
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"O2"},
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"O2"},
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"O2"},
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"O2"},
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"O2"},
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
            .equation = "O2 + N2+ <=> O + O + N2+",
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
                        .name = "N2+",
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
                        .name = "N2+",
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"O2"},
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
            .equation = "O2 + O2+ <=> O + O + O2+",
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
                        .name = "O2+",
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
                        .name = "O2+",
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"O2"},
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
            .equation = "O2 + NO+ <=> O + O + NO+",
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
                        .name = "NO+",
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
                        .name = "NO+",
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"O2"},
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
            .equation = "O2 + N+ <=> O + O + N+",
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
                        .name = "N+",
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
                        .name = "N+",
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"O2"},
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
            .equation = "O2 + O+ <=> O + O + O+",
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
                        .name = "O+",
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
                        .name = "O+",
                        .value = 1.0,
                      },
                    },
                  },
              },
            .forward =
              {
                .A = 3000000000000000.0,
                .n = -1.5,
                .activation_temperature = 59500.0,
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"O2"},
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
                .translational_exponent = 0.5,
                .other_temperature = VibrationTemperatureInput {"NO"},
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
                .count = 11,
                .entries =
                  {
                    {
                      {
                        .name = "N2",
                        .value = 1.0,
                      },
                      {
                        .name = "N2+",
                        .value = 1.0,
                      },
                      {
                        .name = "O2",
                        .value = 1.0,
                      },
                      {
                        .name = "O2+",
                        .value = 1.0,
                      },
                      {
                        .name = "NO",
                        .value = 0.6648275862068966,
                      },
                      {
                        .name = "NO+",
                        .value = 0.6648275862068966,
                      },
                      {
                        .name = "N",
                        .value = 0.6648275862068966,
                      },
                      {
                        .name = "N+",
                        .value = 0.6648275862068966,
                      },
                      {
                        .name = "O",
                        .value = 0.6648275862068966,
                      },
                      {
                        .name = "O+",
                        .value = 0.6648275862068966,
                      },
                      {
                        .name = "e-",
                        .value = 0.0,
                      },
                    },
                  },
              },
          },
          {
            .equation = "N2 + e- <=> N + N + e-",
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
                        .name = "e-",
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
                        .name = "e-",
                        .value = 1.0,
                      },
                    },
                  },
              },
            .forward =
              {
                .A = 1.2e+19,
                .n = -1.6,
                .activation_temperature = 113200.0,
                .translational_exponent = 0.5,
                .other_temperature = FreeElectronTemperatureInput {},
              },
            .backward_model = BackwardRateModel::equilibrium,
            .backward =
              {
                .A = 0.0,
                .n = 0.0,
                .activation_temperature = 0.0,
              },
            .equilibrium_temperature = FreeElectronTemperatureInput {},
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
          {
            .equation = "N + O <=> NO+ + e-",
            .reactants =
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
            .products =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "NO+",
                        .value = 1.0,
                      },
                      {
                        .name = "e-",
                        .value = 1.0,
                      },
                    },
                  },
              },
            .forward =
              {
                .A = 880.0,
                .n = 1.0,
                .activation_temperature = 31900.0,
              },
            .backward_model = BackwardRateModel::equilibrium,
            .backward =
              {
                .A = 0.0,
                .n = 0.0,
                .activation_temperature = 0.0,
              },
            .equilibrium_temperature = FreeElectronTemperatureInput {},
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
            .equation = "O + O <=> O2+ + e-",
            .reactants =
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
            .products =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "O2+",
                        .value = 1.0,
                      },
                      {
                        .name = "e-",
                        .value = 1.0,
                      },
                    },
                  },
              },
            .forward =
              {
                .A = 0.00071,
                .n = 2.7,
                .activation_temperature = 80600.0,
              },
            .backward_model = BackwardRateModel::equilibrium,
            .backward =
              {
                .A = 0.0,
                .n = 0.0,
                .activation_temperature = 0.0,
              },
            .equilibrium_temperature = FreeElectronTemperatureInput {},
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
            .equation = "N + N <=> N2+ + e-",
            .reactants =
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
            .products =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "N2+",
                        .value = 1.0,
                      },
                      {
                        .name = "e-",
                        .value = 1.0,
                      },
                    },
                  },
              },
            .forward =
              {
                .A = 44.0,
                .n = 1.5,
                .activation_temperature = 67500.0,
              },
            .backward_model = BackwardRateModel::equilibrium,
            .backward =
              {
                .A = 0.0,
                .n = 0.0,
                .activation_temperature = 0.0,
              },
            .equilibrium_temperature = FreeElectronTemperatureInput {},
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
            .equation = "NO+ + O <=> N+ + O2",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "NO+",
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
                        .name = "N+",
                        .value = 1.0,
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
                .A = 1000000.0,
                .n = 0.5,
                .activation_temperature = 77200.0,
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
            .equation = "N+ + N2 <=> N2+ + N",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "N+",
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
                        .name = "N2+",
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
                .A = 1000000.0,
                .n = 0.5,
                .activation_temperature = 12200.0,
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
            .equation = "O2+ + N <=> N+ + O2",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "O2+",
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
                        .name = "N+",
                        .value = 1.0,
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
                .A = 87000000.0,
                .n = 0.14,
                .activation_temperature = 28600.0,
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
            .equation = "O+ + NO <=> N+ + O2",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "O+",
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
                        .name = "N+",
                        .value = 1.0,
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
                .A = 0.14,
                .n = 1.9,
                .activation_temperature = 26600.0,
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
            .equation = "O2+ + N2 <=> N2+ + O2",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "O2+",
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
                        .name = "N2+",
                        .value = 1.0,
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
                .A = 9900000.0,
                .n = 0.0,
                .activation_temperature = 40700.0,
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
            .equation = "O2+ + O <=> O+ + O2",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "O2+",
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
                        .name = "O+",
                        .value = 1.0,
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
                .A = 4000000.0,
                .n = -0.09,
                .activation_temperature = 18000.0,
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
            .equation = "NO+ + N <=> O+ + N2",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "NO+",
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
                        .name = "O+",
                        .value = 1.0,
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
                .A = 34000000.0,
                .n = -1.08,
                .activation_temperature = 12800.0,
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
            .equation = "NO+ + O2 <=> O2+ + NO",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "NO+",
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
                        .name = "O2+",
                        .value = 1.0,
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
                .A = 24000000.0,
                .n = 0.41,
                .activation_temperature = 32600.0,
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
            .equation = "NO+ + O <=> O2+ + N",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "NO+",
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
                        .name = "O2+",
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
                .A = 7200000.0,
                .n = 0.29,
                .activation_temperature = 48600.0,
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
            .equation = "O+ + N2 <=> N2+ + O",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "O+",
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
                        .name = "N2+",
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
                .A = 910000.0,
                .n = 0.36,
                .activation_temperature = 22800.0,
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
            .equation = "NO+ + N <=> N2+ + O",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "NO+",
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
                        .name = "N2+",
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
                .A = 72000000.0,
                .n = 0.0,
                .activation_temperature = 35500.0,
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
            .equation = "O + e- <=> O+ + e- + e-",
            .reactants =
              {
                .count = 2,
                .entries =
                  {
                    {
                      {
                        .name = "O",
                        .value = 1.0,
                      },
                      {
                        .name = "e-",
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
                        .name = "O+",
                        .value = 1.0,
                      },
                      {
                        .name = "e-",
                        .value = 2.0,
                      },
                    },
                  },
              },
            .forward =
              {
                .A = 3.9e+27,
                .n = -3.78,
                .activation_temperature = 158500.0,
                .translational_exponent = 0.0,
                .other_temperature = FreeElectronTemperatureInput {},
              },
            .backward_model = BackwardRateModel::equilibrium,
            .backward =
              {
                .A = 0.0,
                .n = 0.0,
                .activation_temperature = 0.0,
              },
            .equilibrium_temperature = FreeElectronTemperatureInput {},
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
            .equation = "N + e- <=> N+ + e- + e-",
            .reactants =
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
                        .name = "e-",
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
                        .name = "N+",
                        .value = 1.0,
                      },
                      {
                        .name = "e-",
                        .value = 2.0,
                      },
                    },
                  },
              },
            .forward =
              {
                .A = 2.5e+28,
                .n = -3.82,
                .activation_temperature = 168600.0,
                .translational_exponent = 0.0,
                .other_temperature = FreeElectronTemperatureInput {},
              },
            .backward_model = BackwardRateModel::equilibrium,
            .backward =
              {
                .A = 0.0,
                .n = 0.0,
                .activation_temperature = 0.0,
              },
            .equilibrium_temperature = FreeElectronTemperatureInput {},
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

// Mechanism of air11_kim_jo over the species of air11_species.
template <std::floating_point P, UnitSystem U>
using Air11KimJo = ReactionMechanism<P, U, air11_species.size(), 40>;

namespace detail
{
inline constexpr auto air11_kim_jo_species_names = [] {
    std::array<std::string_view, air11_species.size()> names {};
    for (const SpeciesIndex i : species_indices(0, names.size())) {
        names[i.value()] = air11_species[i.value()].name;
    }
    return names;
}();
}  // namespace detail

// Air11KimJo model of air11_kim_jo for the species of air11_species, in
// their order, in the unit system of the scales,
// with standard pressure p0. See make_reaction_mechanism.
template <std::floating_point P, UnitSystem U>
Air11KimJo<P, U> make_air11_kim_jo(const UnitScales<P, U>& scales, Pressure<P, U> standard_pressure)
{
    return make_reaction_mechanism(air11_kim_jo,
                                   reaction_species(air11_species),
                                   scales,
                                   standard_pressure);
}

// As above, with the standard pressure 1 bar (1e5 Pa).
template <std::floating_point P, UnitSystem U>
Air11KimJo<P, U> make_air11_kim_jo(const UnitScales<P, U>& scales)
{ return make_reaction_mechanism(air11_kim_jo, reaction_species(air11_species), scales); }

// As above for the thermo model gas with several energy modes: rates'
// other temperatures resolve to its modes. See make_reaction_mechanism.
template <std::floating_point P, UnitSystem U, typename G>
Air11KimJo<P, U>
  make_air11_kim_jo(const UnitScales<P, U>& scales, Pressure<P, U> standard_pressure, const G& gas)
{
    return make_reaction_mechanism(air11_kim_jo,
                                   reaction_species(air11_species),
                                   scales,
                                   standard_pressure,
                                   gas);
}

}  // namespace canary

#endif  // CANARY_DATA_AIR11_KIM_JO_HPP
