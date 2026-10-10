// Quantities of different kinds cannot be added.
#include <gaslib/quantity.hpp>

#ifdef GASLIB_COMPILE_FAIL
auto bad =
  gaslib::Pressure<double, gaslib::si> {1.0} + gaslib::Temperature<double, gaslib::si> {1.0};
#endif
auto good = gaslib::Pressure<double, gaslib::si> {1.0} + gaslib::Pressure<double, gaslib::si> {1.0};
