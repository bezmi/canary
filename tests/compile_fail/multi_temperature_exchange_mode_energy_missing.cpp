// The single-temperature model has no energy modes and no energy per mode
// of each species.
#include <gaslib/energy_exchange.hpp>
#include <gaslib/models/multi_temperature.hpp>
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/quantity.hpp>

#ifdef GASLIB_COMPILE_FAIL
static_assert(
  gaslib::HasModeEnergies<gaslib::ThermallyPerfectMixture<double, gaslib::si, 11, 4>, double, 11>);
#endif
static_assert(
  gaslib::HasModeEnergies<gaslib::MultiTemperatureMixture<double, gaslib::si, 11, 2>, double, 11>);
