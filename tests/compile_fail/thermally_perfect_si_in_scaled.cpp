// An SI quantity cannot be passed to a thermally perfect model in scaled
// units.
#include <array>
#include <cstddef>
#include <canary/data/air5.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/thermo.hpp>
#include <canary/units.hpp>
#include <mdspan/mdspan.hpp>

namespace
{
const auto gas = canary::make_air5<double>(canary::UnitScales<double, canary::scaled> {
                                             .mass = 4.65e-26,
                                             .density = 1.0e-3,
                                             .temperature = 1000.0,
                                             .length = 1.0,
                                           })
                   .model;
const std::array<double, 5> y {0.767, 0.233, 0.0, 0.0, 0.0};
const auto mix =
  canary::mixture(gas,
                  canary::as_kind<canary::kind::mass_fraction>(
                    Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data(), 5)));
}  // namespace

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::thermo_state(gas,
                                mix,
                                canary::Density<double, canary::si> {1.0},
                                canary::SpecificEnergy<double, canary::si> {1.0});
#endif
auto good = canary::thermo_state(gas,
                                 mix,
                                 canary::Density<double, canary::scaled> {1.0},
                                 canary::SpecificEnergy<double, canary::scaled> {1.0});
