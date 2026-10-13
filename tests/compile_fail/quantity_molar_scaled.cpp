// Kinds with an amount of substance exist only in SI.
#include <canary/quantity.hpp>

#ifdef CANARY_COMPILE_FAIL
canary::Quantity<canary::kind::molar_mass, double, canary::scaled> bad {0.029};
#endif
canary::MolarMass<double> good {0.029};
