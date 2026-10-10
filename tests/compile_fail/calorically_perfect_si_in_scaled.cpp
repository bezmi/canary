// An SI quantity cannot be passed to a model in scaled units.
#include <gaslib/models/calorically_perfect.hpp>
#include <gaslib/thermo.hpp>

namespace
{
const auto gas = gaslib::make_calorically_perfect(
  gaslib::SpecificGasConstant<double, gaslib::scaled> {1.0}, gaslib::HeatCapacityRatio {1.4});
}  // namespace

#ifdef GASLIB_COMPILE_FAIL
auto bad = gaslib::thermo_state(
  gas, gaslib::Density<double, gaslib::si> {1.0}, gaslib::SpecificEnergy<double, gaslib::si> {1.0});
#endif
auto good = gaslib::thermo_state(gas,
                                 gaslib::Density<double, gaslib::scaled> {1.0},
                                 gaslib::SpecificEnergy<double, gaslib::scaled> {1.0});
