// Two densities do not determine a state.
#include "test_models.hpp"

#include <gaslib/thermo.hpp>

namespace
{
constexpr gaslib::test::IdealGas gas {{gaslib::SpecificGasConstant<double, gaslib::si> {287.0},
                                       gaslib::SpecificHeat<double, gaslib::si> {717.5}}};
constexpr gaslib::Density<double, gaslib::si> rho {1.2};
}  // namespace

#ifdef GASLIB_COMPILE_FAIL
auto bad = gaslib::thermo_state(gas, rho, rho);
#endif
auto good = gaslib::thermo_state(gas, rho, gaslib::Temperature<double, gaslib::si> {300.0});
