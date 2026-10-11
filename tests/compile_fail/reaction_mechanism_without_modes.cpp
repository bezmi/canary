// make_reaction_mechanism with a gas model resolves the other temperatures
// of the rates to the model's energy modes; the single-temperature model has
// none.
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/models/two_temperature.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/reactions.hpp>

#ifdef GASLIB_COMPILE_FAIL
static_assert(gaslib::HasEnergyModes<gaslib::ThermallyPerfectMixture<double, gaslib::si, 5, 4>>);
#endif
static_assert(gaslib::HasEnergyModes<gaslib::TwoTemperatureMixture<double, gaslib::si, 5, 4>>);
