// Quantities of different kinds cannot be added.
#include <canary/quantity.hpp>

#ifdef CANARY_COMPILE_FAIL
auto bad =
  canary::Pressure<double, canary::si> {1.0} + canary::Temperature<double, canary::si> {1.0};
#endif
auto good = canary::Pressure<double, canary::si> {1.0} + canary::Pressure<double, canary::si> {1.0};
