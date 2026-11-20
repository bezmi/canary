// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// A mode energy density view in SI units cannot be passed to a
// multi-temperature model in scaled units.
#include <array>
#include <canary/data/air5.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/thermo.hpp>
#include <canary/units.hpp>
#include <cstddef>
#include <mdspan/mdspan.hpp>

namespace
{
using Extents = Kokkos::dextents<std::size_t, 1>;
const auto gas =
  canary::make_air5_four_temperature<double>(canary::UnitScales<double, canary::scaled> {{
                                               .mass = 4.65e-26,
                                               .density = 1.0e-3,
                                               .temperature = 1000.0,
                                               .length = 1.0,
                                             }})
    .model;
const std::array<double, 5> y {0.767, 0.233, 0.0, 0.0, 0.0};
const auto mix = canary::mixture(
  gas,
  canary::as_kind<canary::kind::mass_fraction>(Kokkos::mdspan<const double, Extents>(y.data(), 5)));
const std::array<double, 3> rho_e_m {1.0, 1.0, 0.0};
}  // namespace

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::thermo_state(gas,
                                mix,
                                canary::Density<double, canary::scaled> {1.0},
                                canary::EnergyDensity<double, canary::scaled> {1.0},
                                canary::as_kind<canary::kind::mode_energy_density, canary::si>(
                                  Kokkos::mdspan<const double, Extents>(rho_e_m.data(), 3)));
#endif
auto good = canary::thermo_state(gas,
                                 mix,
                                 canary::Density<double, canary::scaled> {1.0},
                                 canary::EnergyDensity<double, canary::scaled> {1.0},
                                 canary::as_kind<canary::kind::mode_energy_density, canary::scaled>(
                                   Kokkos::mdspan<const double, Extents>(rho_e_m.data(), 3)));
