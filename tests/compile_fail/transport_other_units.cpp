// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// A viscosity law in SI units does not accept a scaled temperature.
#include <canary/transport/sutherland.hpp>

namespace
{
constexpr canary::Sutherland<double, canary::si> law {
  .mu_0 = canary::DynamicViscosity<double, canary::si> {1.716e-5},
  .t_0 = canary::Temperature<double, canary::si> {273.15},
  .s = canary::Temperature<double, canary::si> {110.4},
};
}  // namespace

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::viscosity(law, canary::Temperature<double, canary::scaled> {1.0});
#endif
auto good = canary::viscosity(law, canary::Temperature<double, canary::si> {300.0});
