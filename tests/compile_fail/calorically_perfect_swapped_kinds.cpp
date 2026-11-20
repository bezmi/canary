// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// make_calorically_perfect takes its parameters in a fixed order of kinds.
#include <canary/models/calorically_perfect.hpp>

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::make_calorically_perfect(
  canary::HeatCapacityRatio {1.4}, canary::SpecificGasConstant<double, canary::si> {287.0});
#endif
auto good = canary::make_calorically_perfect(
  canary::SpecificGasConstant<double, canary::si> {287.0}, canary::HeatCapacityRatio {1.4});
