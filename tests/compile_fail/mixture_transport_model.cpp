// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// Gupta-Yos transport is a transport model of the two-temperature mixture,
// with one conductivity per mode, but not of the calorically perfect
// mixture.
#include <canary/data/air5.hpp>
#include <canary/models/calorically_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/transport/transport.hpp>

#ifdef CANARY_COMPILE_FAIL
static_assert(canary::TransportModel<canary::Air5GuptaYosTransport<double, canary::si>,
                                     canary::CaloricallyPerfectMixture<double, canary::si, 5>,
                                     double>);
#endif
static_assert(canary::TransportModel<canary::Air5GuptaYosTransport<double, canary::si>,
                                     canary::Air5TwoTemperature<double, canary::si>,
                                     double>);
