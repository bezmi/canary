// A type with arithmetic but no ScalarTraits specialisation is not a state
// scalar.
#include "scalar_types.hpp"

#ifdef GASLIB_COMPILE_FAIL
bool bad = gaslib::math::isfinite(gaslib::test::NoTraits {1.0});
#endif
bool good = gaslib::math::isfinite(1.0);
