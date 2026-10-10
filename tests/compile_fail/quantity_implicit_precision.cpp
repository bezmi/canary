// A double quantity does not convert implicitly to a float quantity.
#include <gaslib/quantity.hpp>

#ifdef GASLIB_COMPILE_FAIL
gaslib::Pressure<float, gaslib::si> bad = gaslib::Pressure<double, gaslib::si> {1.0};
#endif
gaslib::Pressure<float, gaslib::si> good {gaslib::Pressure<double, gaslib::si> {1.0}};
