// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// The Gupta-Yos model needs the isobaric specific heat of every species; the
// calorically perfect mixture does not give it.
#include <canary/models/calorically_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/transport/mixture.hpp>

#ifdef CANARY_COMPILE_FAIL
static_assert(
  canary::HasSpeciesIsobaricSpecificHeat<canary::CaloricallyPerfectMixture<double, canary::si, 5>,
                                         double,
                                         5>);
#endif
static_assert(
  canary::HasSpeciesIsobaricSpecificHeat<canary::TwoTemperatureMixture<double, canary::si, 5, 4>,
                                         double,
                                         5>);
