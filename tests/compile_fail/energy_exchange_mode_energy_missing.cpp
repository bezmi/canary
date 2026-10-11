// The single-temperature model has no vibrational-electronic energy per
// species.
#include <gaslib/energy_exchange.hpp>
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/models/two_temperature.hpp>
#include <gaslib/quantity.hpp>

#ifdef GASLIB_COMPILE_FAIL
static_assert(
  gaslib::HasVibrationalElectronicEnergy<gaslib::ThermallyPerfectMixture<double, gaslib::si, 5, 4>,
                                         double,
                                         5>);
#endif
static_assert(
  gaslib::HasVibrationalElectronicEnergy<gaslib::TwoTemperatureMixture<double, gaslib::si, 5, 4>,
                                         double,
                                         5>);
