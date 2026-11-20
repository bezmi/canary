// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// The single-temperature model has no vibrational-electronic energy per
// species.
#include <canary/energy_exchange.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>

#ifdef CANARY_COMPILE_FAIL
static_assert(
  canary::HasVibrationalElectronicEnergy<canary::ThermallyPerfectMixture<double, canary::si, 5, 4>,
                                         double,
                                         5>);
#endif
static_assert(
  canary::HasVibrationalElectronicEnergy<canary::TwoTemperatureMixture<double, canary::si, 5, 4>,
                                         double,
                                         5>);
