// A product of different kinds exists only as a relation function.
#include <gaslib/quantity.hpp>

#ifdef GASLIB_COMPILE_FAIL
auto bad =
  gaslib::Density<double, gaslib::si> {1.0} * gaslib::Temperature<double, gaslib::si> {1.0};
#endif
auto good = gaslib::Density<double, gaslib::si> {1.0} * 2.0;
