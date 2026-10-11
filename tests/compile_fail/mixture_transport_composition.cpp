// The mixture transport models take the mole fractions from a
// NormalisedComposition; the mixture value of the calorically perfect
// mixture is not one.
#include <gaslib/models/calorically_perfect.hpp>
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/transport/mixture.hpp>

#ifdef GASLIB_COMPILE_FAIL
static_assert(
  gaslib::
    HasNormalisedComposition<gaslib::CaloricallyPerfectMixture<double, gaslib::si, 5>, double, 5>);
#endif
static_assert(
  gaslib::
    HasNormalisedComposition<gaslib::ThermallyPerfectMixture<double, gaslib::si, 5, 4>, double, 5>);
