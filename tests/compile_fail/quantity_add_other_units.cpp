// Quantities in different unit systems cannot be added.
#include <gaslib/quantity.hpp>

#ifdef GASLIB_COMPILE_FAIL
auto bad =
  gaslib::Density<double, gaslib::si> {1.0} + gaslib::Density<double, gaslib::scaled> {1.0};
#endif
auto good =
  gaslib::Density<double, gaslib::scaled> {1.0} + gaslib::Density<double, gaslib::scaled> {1.0};
