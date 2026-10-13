// Quantities in different unit systems cannot be added.
#include <canary/quantity.hpp>

#ifdef CANARY_COMPILE_FAIL
auto bad =
  canary::Density<double, canary::si> {1.0} + canary::Density<double, canary::scaled> {1.0};
#endif
auto good =
  canary::Density<double, canary::scaled> {1.0} + canary::Density<double, canary::scaled> {1.0};
