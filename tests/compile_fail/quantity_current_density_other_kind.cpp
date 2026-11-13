// The current density needs a momentum density, not a mass density.
#include <canary/quantity.hpp>

#ifdef CANARY_COMPILE_FAIL
auto bad = canary::current_density(canary::SpecificCharge<double, canary::si> {1.0},
                                   canary::Density<double, canary::si> {1.0});
#endif
auto good = canary::current_density(canary::SpecificCharge<double, canary::si> {1.0},
                                    canary::MomentumDensity<double, canary::si> {1.0});
