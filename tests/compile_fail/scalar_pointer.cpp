// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// A pointer is not a state scalar.
#include <canary/scalar.hpp>

template <canary::Scalar T> T twice(T x) { return x + x; }

#ifdef CANARY_COMPILE_FAIL
int* bad = twice(static_cast<int*>(nullptr));
#endif
double good = twice(1.0);
