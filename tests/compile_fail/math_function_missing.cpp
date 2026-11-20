// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// math::exp needs exp declared next to a non-floating scalar type; this type
// provides only sqrt.
#include "scalar_types.hpp"

#ifdef CANARY_COMPILE_FAIL
canary::test::OnlySqrt bad = canary::math::exp(canary::test::OnlySqrt {1.0});
#endif
canary::test::OnlySqrt good = canary::math::sqrt(canary::test::OnlySqrt {1.0});
