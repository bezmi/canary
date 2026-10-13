// A model without the (p, T) closure is not a GasModel.
#include "test_models.hpp"

#include <canary/thermo.hpp>

namespace
{
template <canary::GasModel<double> G> int use(const G& /*gas*/) { return 0; }
constexpr canary::test::IdealParameters air {
  canary::SpecificGasConstant<double, canary::si> {287.0},
  canary::SpecificHeat<double, canary::si> {717.5}};
}  // namespace

#ifdef CANARY_COMPILE_FAIL
int bad = use(canary::test::NoPressureTemperatureGas {air});
#endif
int good = use(canary::test::IdealGas {air});
