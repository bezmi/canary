// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// The mixture transport models take the mole fractions from a
// NormalisedComposition; the mixture value of the calorically perfect
// mixture is not one.
#include <canary/models/calorically_perfect.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/transport/mixture.hpp>

#ifdef CANARY_COMPILE_FAIL
static_assert(
  canary::
    HasNormalisedComposition<canary::CaloricallyPerfectMixture<double, canary::si, 5>, double, 5>);
#endif
static_assert(
  canary::
    HasNormalisedComposition<canary::ThermallyPerfectMixture<double, canary::si, 5, 4>, double, 5>);
