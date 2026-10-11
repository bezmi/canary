// vibration_translation_source of a multi-temperature exchange needs a
// model with the energy of every species in each mode; the
// single-temperature model has none.
#include <array>
#include <cstddef>
#include <gaslib/data/air11.hpp>
#include <gaslib/data/air11_kim_jo.hpp>
#include <gaslib/data/air11_kim_jo_exchange.hpp>
#include <gaslib/energy_exchange.hpp>
#include <gaslib/models/multi_temperature.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/species.hpp>
#include <gaslib/state.hpp>
#include <gaslib/units.hpp>
#include <mdspan/mdspan.hpp>

namespace
{
using Extents = Kokkos::dextents<std::size_t, 1>;
const gaslib::UnitScales<double, gaslib::si> scales {};
const auto built = gaslib::make_air11_three_temperature<double>(scales);
const auto single = gaslib::make_air11<double>(scales).model;
const auto mechanism = gaslib::make_air11_kim_jo(
  built.species, scales, gaslib::Pressure<double, gaslib::si> {101325.0}, built.model);
const auto exchange =
  gaslib::make_air11_kim_jo_exchange(built.species, mechanism, built.model, scales);
const std::array<double, 11> y {0.767, 0.233, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const auto mix = gaslib::mixture(built.model,
                                 gaslib::as_kind<gaslib::kind::mass_fraction>(
                                   Kokkos::mdspan<const double, Extents>(y.data(), 11)));
using Mode = gaslib::ModeTemperature<double, gaslib::si>;
using ModeEnergy = gaslib::SpecificEnergy<double, gaslib::si>;
const gaslib::MultiTemperatureState<double, gaslib::si, 2> state {
  .rho = gaslib::Density<double, gaslib::si> {0.1},
  .e = gaslib::SpecificEnergy<double, gaslib::si> {0.0},
  .p = gaslib::Pressure<double, gaslib::si> {0.0},
  .temp = gaslib::Temperature<double, gaslib::si> {5000.0},
  .mode_e = {ModeEnergy {0.0}, ModeEnergy {0.0}},
  .mode_temp = {Mode {1000.0}, Mode {2000.0}},
};
}  // namespace

#ifdef GASLIB_COMPILE_FAIL
auto bad = gaslib::vibration_translation_source(exchange, single, mix, state);
#endif
auto good = gaslib::vibration_translation_source(exchange, built.model, mix, state);
