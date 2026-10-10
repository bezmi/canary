// A pointer is not a state scalar.
#include <gaslib/scalar.hpp>

template <gaslib::Scalar T> T twice(T x) { return x + x; }

#ifdef GASLIB_COMPILE_FAIL
int* bad = twice(static_cast<int*>(nullptr));
#endif
double good = twice(1.0);
