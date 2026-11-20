// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// Kinds with an amount of substance exist only in SI.
#include <canary/quantity.hpp>

#ifdef CANARY_COMPILE_FAIL
canary::Quantity<canary::kind::molar_mass, double, canary::scaled> bad {0.029};
#endif
canary::MolarMass<double> good {0.029};
