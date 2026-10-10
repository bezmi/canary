// math::exp needs exp declared next to a non-floating scalar type; this type
// provides only sqrt.
#include "scalar_types.hpp"

#ifdef GASLIB_COMPILE_FAIL
gaslib::test::OnlySqrt bad = gaslib::math::exp(gaslib::test::OnlySqrt {1.0});
#endif
gaslib::test::OnlySqrt good = gaslib::math::sqrt(gaslib::test::OnlySqrt {1.0});
