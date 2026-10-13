// A product of different kinds exists only as a relation function.
#include <canary/quantity.hpp>

#ifdef CANARY_COMPILE_FAIL
auto bad =
  canary::Density<double, canary::si> {1.0} * canary::Temperature<double, canary::si> {1.0};
#endif
auto good = canary::Density<double, canary::si> {1.0} * 2.0;
