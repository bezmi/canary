// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// Two densities do not determine a state.
#include "test_models.hpp"

#include <canary/thermo.hpp>

namespace
{
constexpr canary::test::IdealGas gas {{canary::SpecificGasConstant<double, canary::si> {287.0},
                                       canary::SpecificHeat<double, canary::si> {717.5}}};
constexpr canary::Density<double, canary::si> rho {1.2};
}  // namespace

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::thermo_state(gas, rho, rho);
#endif
auto good = canary::thermo_state(gas, rho, canary::Temperature<double, canary::si> {300.0});
