// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// An SI quantity cannot be passed to a model in scaled units.
#include <canary/models/calorically_perfect.hpp>
#include <canary/thermo.hpp>

namespace
{
const auto gas = canary::make_calorically_perfect(
  canary::SpecificGasConstant<double, canary::scaled> {1.0}, canary::HeatCapacityRatio {1.4});
}  // namespace

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::thermo_state(
  gas, canary::Density<double, canary::si> {1.0}, canary::SpecificEnergy<double, canary::si> {1.0});
#endif
auto good = canary::thermo_state(gas,
                                 canary::Density<double, canary::scaled> {1.0},
                                 canary::SpecificEnergy<double, canary::scaled> {1.0});
