// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// The one-species closure from (rho, rho_e) through canary, and the same
// arithmetic and validity check on plain doubles, as separate functions for comparing their
// generated code.

#include <cmath>
#include <canary/models/calorically_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/state.hpp>

// The functions have external linkage so that the compiler emits them.
// NOLINTBEGIN(misc-use-internal-linkage)
struct RawClosure {
    double rho;
    double e;
    double p;
    double temp;
    unsigned char status;  // 0 valid, 3 invalid, as canary::Status
};

struct RawGas {
    double gas_constant;
    double gamma;
};

canary::Closure<double, canary::si>
  canary_close(const canary::CaloricallyPerfect<double, canary::si>& gas,
               canary::Density<double, canary::si> rho,
               canary::EnergyDensity<double, canary::si> rho_e)
{ return close(gas, rho, rho_e); }

RawClosure raw_close(const RawGas& gas, double rho, double rho_e)
{
    const double p = (gas.gamma - 1.0) * rho_e;
    const double temp = p / (rho * gas.gas_constant);
    const auto positive = [](double x) { return std::isfinite(x) && x > 0.0; };
    const bool valid = positive(rho) && positive(rho_e) && positive(temp);
    return {
      .rho = rho,
      .e = rho_e / rho,
      .p = p,
      .temp = temp,
      .status = static_cast<unsigned char>(valid ? 0 : 3),
    };
}
// NOLINTEND(misc-use-internal-linkage)
