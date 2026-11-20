// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// vibration_translation_source needs a model with the
// vibrational-electronic energy of every species; the single-temperature
// model has none.
#include <array>
#include <canary/data/air5.hpp>
#include <canary/data/air5_park.hpp>
#include <canary/data/air5_park_exchange.hpp>
#include <canary/energy_exchange.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/units.hpp>
#include <cstddef>
#include <mdspan/mdspan.hpp>

namespace
{
using Extents = Kokkos::dextents<std::size_t, 1>;
const canary::UnitScales<double, canary::si> scales {};
const auto built = canary::make_air5_two_temperature<double>(scales);
const auto single = canary::make_air5<double>(scales).model;
const auto mechanism = canary::make_air5_park(scales);
const auto exchange = canary::make_air5_park_exchange(built.species, mechanism, scales);
const std::array<double, 5> y {0.767, 0.233, 0.0, 0.0, 0.0};
const auto mix = canary::mixture(
  built.model,
  canary::as_kind<canary::kind::mass_fraction>(Kokkos::mdspan<const double, Extents>(y.data(), 5)));
const canary::Density<double, canary::si> rho {0.1};
const canary::Temperature<double, canary::si> temp {5000.0};
const canary::Temperature<double, canary::si> temp_ve {1000.0};
}  // namespace

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::vibration_translation_source(exchange, single, mix, rho, temp, temp_ve);
#endif
auto good = canary::vibration_translation_source(exchange, built.model, mix, rho, temp, temp_ve);
