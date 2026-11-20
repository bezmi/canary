// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// vibration_translation_source of a multi-temperature exchange needs a
// model with the energy of every species in each mode; the
// single-temperature model has none.
#include <array>
#include <cstddef>
#include <canary/data/air11.hpp>
#include <canary/data/air11_kim_jo.hpp>
#include <canary/data/air11_kim_jo_exchange.hpp>
#include <canary/energy_exchange.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/units.hpp>
#include <mdspan/mdspan.hpp>

namespace
{
using Extents = Kokkos::dextents<std::size_t, 1>;
const canary::UnitScales<double, canary::si> scales {};
const auto built = canary::make_air11_three_temperature<double>(scales);
const auto single = canary::make_air11<double>(scales).model;
const auto mechanism =
  canary::make_air11_kim_jo(scales, canary::Pressure<double, canary::si> {101325.0}, built.model);
const auto exchange =
  canary::make_air11_kim_jo_exchange(built.species, mechanism, built.model, scales);
const std::array<double, 11> y {0.767, 0.233, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const auto mix = canary::mixture(built.model,
                                 canary::as_kind<canary::kind::mass_fraction>(
                                   Kokkos::mdspan<const double, Extents>(y.data(), 11)));
using Mode = canary::ModeTemperature<double, canary::si>;
using ModeEnergy = canary::SpecificEnergy<double, canary::si>;
const canary::MultiTemperatureState<double, canary::si, 2> state {
  .rho = canary::Density<double, canary::si> {0.1},
  .e = canary::SpecificEnergy<double, canary::si> {0.0},
  .p = canary::Pressure<double, canary::si> {0.0},
  .temp = canary::Temperature<double, canary::si> {5000.0},
  .mode_e = {ModeEnergy {0.0}, ModeEnergy {0.0}},
  .mode_temp = {Mode {1000.0}, Mode {2000.0}},
};
}  // namespace

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::vibration_translation_source(exchange, single, mix, state);
#endif
auto good = canary::vibration_translation_source(exchange, built.model, mix, state);
