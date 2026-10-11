// The Gupta-Yos model needs the isobaric specific heat of every species; the
// calorically perfect mixture does not give it.
#include <gaslib/models/calorically_perfect.hpp>
#include <gaslib/models/two_temperature.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/transport/mixture.hpp>

#ifdef GASLIB_COMPILE_FAIL
static_assert(
  gaslib::HasSpeciesIsobaricSpecificHeat<gaslib::CaloricallyPerfectMixture<double, gaslib::si, 5>,
                                         double,
                                         5>);
#endif
static_assert(
  gaslib::HasSpeciesIsobaricSpecificHeat<gaslib::TwoTemperatureMixture<double, gaslib::si, 5, 4>,
                                         double,
                                         5>);
