// An SI quantity cannot be passed to a thermally perfect model in scaled
// units.
#include <array>
#include <cstddef>
#include <gaslib/data/air5.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/species.hpp>
#include <gaslib/thermo.hpp>
#include <gaslib/units.hpp>
#include <mdspan/mdspan.hpp>

namespace
{
const auto gas = gaslib::make_air5<double>(gaslib::UnitScales<double, gaslib::scaled> {
                                             .mass = 4.65e-26,
                                             .density = 1.0e-3,
                                             .temperature = 1000.0,
                                             .length = 1.0,
                                           })
                   .model;
const std::array<double, 5> y {0.767, 0.233, 0.0, 0.0, 0.0};
const auto mix =
  gaslib::mixture(gas,
                  gaslib::as_kind<gaslib::kind::mass_fraction>(
                    Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data(), 5)));
}  // namespace

#ifdef GASLIB_COMPILE_FAIL
auto bad = gaslib::thermo_state(gas,
                                mix,
                                gaslib::Density<double, gaslib::si> {1.0},
                                gaslib::SpecificEnergy<double, gaslib::si> {1.0});
#endif
auto good = gaslib::thermo_state(gas,
                                 mix,
                                 gaslib::Density<double, gaslib::scaled> {1.0},
                                 gaslib::SpecificEnergy<double, gaslib::scaled> {1.0});
