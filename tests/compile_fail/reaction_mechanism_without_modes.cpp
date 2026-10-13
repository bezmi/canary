// make_reaction_mechanism with a gas model resolves the other temperatures
// of the rates to the model's energy modes; the single-temperature model has
// none.
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>

#ifdef CANARY_COMPILE_FAIL
static_assert(canary::HasEnergyModes<canary::ThermallyPerfectMixture<double, canary::si, 5, 4>>);
#endif
static_assert(canary::HasEnergyModes<canary::TwoTemperatureMixture<double, canary::si, 5, 4>>);
