// The Gupta-Yos mode conductivity needs a model whose single mode holds the
// vibrational and electronic energy of every species. The two-temperature
// mixture declares that mode; a multi-temperature mixture with one mode does
// not.
#include <canary/models/multi_temperature.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/transport/mixture.hpp>
#include <canary/units.hpp>

#ifdef CANARY_COMPILE_FAIL
static_assert(
  canary::HasVibrationalElectronicMode<canary::MultiTemperatureMixture<double, canary::si, 11, 1>>);
#endif
static_assert(
  canary::HasVibrationalElectronicMode<canary::TwoTemperatureMixture<double, canary::si, 5, 4>>);
