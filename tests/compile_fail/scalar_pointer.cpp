// A pointer is not a state scalar.
#include <canary/scalar.hpp>

template <canary::Scalar T> T twice(T x) { return x + x; }

#ifdef CANARY_COMPILE_FAIL
int* bad = twice(static_cast<int*>(nullptr));
#endif
double good = twice(1.0);
