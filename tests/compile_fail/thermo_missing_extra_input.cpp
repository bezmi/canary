// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// A model whose closure requires an extra input cannot close without it.
#include "test_models.hpp"

#include <canary/thermo.hpp>

namespace
{
constexpr canary::test::GuessRequiredGas gas {
  {canary::SpecificGasConstant<double, canary::si> {287.0},
   canary::SpecificHeat<double, canary::si> {717.5}}};
constexpr canary::Density<double, canary::si> rho {1.2};
constexpr canary::SpecificEnergy<double, canary::si> e {2.15e5};
}  // namespace

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::thermo_state(gas, rho, e);
#endif
auto good = canary::thermo_state(
  gas, rho, e, canary::test::TemperatureGuess {canary::Temperature<double, canary::si> {300.0}});
