// The single-temperature model has no energy modes and no energy per mode
// of each species.
#include <canary/energy_exchange.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/quantity.hpp>

#ifdef CANARY_COMPILE_FAIL
static_assert(
  canary::HasModeEnergies<canary::ThermallyPerfectMixture<double, canary::si, 11, 4>, double, 11>);
#endif
static_assert(
  canary::HasModeEnergies<canary::MultiTemperatureMixture<double, canary::si, 11, 2>, double, 11>);
