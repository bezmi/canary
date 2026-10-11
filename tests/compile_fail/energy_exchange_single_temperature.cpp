// vibration_translation_source needs a model with the
// vibrational-electronic energy of every species; the single-temperature
// model has none.
#include <array>
#include <cstddef>
#include <gaslib/data/air5.hpp>
#include <gaslib/data/air5_park.hpp>
#include <gaslib/data/air5_park_exchange.hpp>
#include <gaslib/energy_exchange.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/species.hpp>
#include <gaslib/units.hpp>
#include <mdspan/mdspan.hpp>

namespace
{
using Extents = Kokkos::dextents<std::size_t, 1>;
const gaslib::UnitScales<double, gaslib::si> scales {};
const auto built = gaslib::make_air5_two_temperature<double>(scales);
const auto single = gaslib::make_air5<double>(scales).model;
const auto mechanism = gaslib::make_air5_park(built.species, scales);
const auto exchange = gaslib::make_air5_park_exchange(built.species, mechanism, scales);
const std::array<double, 5> y {0.767, 0.233, 0.0, 0.0, 0.0};
const auto mix = gaslib::mixture(
  built.model,
  gaslib::as_kind<gaslib::kind::mass_fraction>(Kokkos::mdspan<const double, Extents>(y.data(), 5)));
const gaslib::Density<double, gaslib::si> rho {0.1};
const gaslib::Temperature<double, gaslib::si> temp {5000.0};
const gaslib::Temperature<double, gaslib::si> temp_ve {1000.0};
}  // namespace

#ifdef GASLIB_COMPILE_FAIL
auto bad = gaslib::vibration_translation_source(exchange, single, mix, rho, temp, temp_ve);
#endif
auto good = gaslib::vibration_translation_source(exchange, built.model, mix, rho, temp, temp_ve);
