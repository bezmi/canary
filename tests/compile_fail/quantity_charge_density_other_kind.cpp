// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// The charge density needs the charge per unit mass, not the charge of one
// particle.
#include <canary/quantity.hpp>

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::charge_density(canary::Density<double, canary::si> {1.0},
                                  canary::Charge<double, canary::si> {1.0});
#endif
auto good = canary::charge_density(canary::Density<double, canary::si> {1.0},
                                   canary::SpecificCharge<double, canary::si> {1.0});
