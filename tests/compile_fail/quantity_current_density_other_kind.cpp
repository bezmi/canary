// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// The current density needs a momentum density, not a mass density.
#include <canary/quantity.hpp>

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::current_density(canary::SpecificCharge<double, canary::si> {1.0},
                                   canary::Density<double, canary::si> {1.0});
#endif
auto good = canary::current_density(canary::SpecificCharge<double, canary::si> {1.0},
                                    canary::MomentumDensity<double, canary::si> {1.0});
