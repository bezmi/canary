// Kinds with an amount of substance exist only in SI.
#include <gaslib/quantity.hpp>

#ifdef GASLIB_COMPILE_FAIL
gaslib::Quantity<gaslib::kind::molar_mass, double, gaslib::scaled> bad {0.029};
#endif
gaslib::MolarMass<double> good {0.029};
