// A double quantity does not convert implicitly to a float quantity.
#include <canary/quantity.hpp>

#ifdef CANARY_COMPILE_FAIL
canary::Pressure<float, canary::si> bad = canary::Pressure<double, canary::si> {1.0};
#endif
canary::Pressure<float, canary::si> good {canary::Pressure<double, canary::si> {1.0}};
