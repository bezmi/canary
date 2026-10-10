// A model whose closure requires an extra input cannot close without it.
#include "test_models.hpp"

#include <gaslib/thermo.hpp>

namespace
{
constexpr gaslib::test::GuessRequiredGas gas {
  {gaslib::SpecificGasConstant<double, gaslib::si> {287.0},
   gaslib::SpecificHeat<double, gaslib::si> {717.5}}};
constexpr gaslib::Density<double, gaslib::si> rho {1.2};
constexpr gaslib::SpecificEnergy<double, gaslib::si> e {2.15e5};
}  // namespace

#ifdef GASLIB_COMPILE_FAIL
auto bad = gaslib::thermo_state(gas, rho, e);
#endif
auto good = gaslib::thermo_state(
  gas, rho, e, gaslib::test::TemperatureGuess {gaslib::Temperature<double, gaslib::si> {300.0}});
