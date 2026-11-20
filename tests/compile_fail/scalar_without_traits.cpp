// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// A type with arithmetic but no ScalarTraits specialisation is not a state
// scalar.
#include "scalar_types.hpp"

#ifdef CANARY_COMPILE_FAIL
bool bad = canary::math::isfinite(canary::test::NoTraits {1.0});
#endif
bool good = canary::math::isfinite(1.0);
