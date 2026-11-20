// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// Conversion keeps the kind: an energy density does not become a pressure.
#include <canary/units.hpp>

const canary::UnitScales<double, canary::scaled> scales {{1.0, 1.0, 1.0, 1.0}};

#ifdef CANARY_COMPILE_FAIL
canary::Pressure<double, canary::scaled> bad =
  canary::convert<canary::scaled>(canary::EnergyDensity<double, canary::si> {1.0}, scales);
#endif
canary::Pressure<double, canary::scaled> good =
  canary::convert<canary::scaled>(canary::Pressure<double, canary::si> {1.0}, scales);
