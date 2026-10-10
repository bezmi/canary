// make_calorically_perfect takes its parameters in a fixed order of kinds.
#include <gaslib/models/calorically_perfect.hpp>

#ifdef GASLIB_COMPILE_FAIL
auto bad = gaslib::make_calorically_perfect(
  gaslib::HeatCapacityRatio {1.4}, gaslib::SpecificGasConstant<double, gaslib::si> {287.0});
#endif
auto good = gaslib::make_calorically_perfect(
  gaslib::SpecificGasConstant<double, gaslib::si> {287.0}, gaslib::HeatCapacityRatio {1.4});
