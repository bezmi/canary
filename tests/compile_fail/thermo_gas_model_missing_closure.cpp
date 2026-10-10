// A model without the (p, T) closure is not a GasModel.
#include "test_models.hpp"

#include <gaslib/thermo.hpp>

namespace
{
template <gaslib::GasModel<double> G> int use(const G& /*gas*/) { return 0; }
constexpr gaslib::test::IdealParameters air {
  gaslib::SpecificGasConstant<double, gaslib::si> {287.0},
  gaslib::SpecificHeat<double, gaslib::si> {717.5}};
}  // namespace

#ifdef GASLIB_COMPILE_FAIL
int bad = use(gaslib::test::NoPressureTemperatureGas {air});
#endif
int good = use(gaslib::test::IdealGas {air});
