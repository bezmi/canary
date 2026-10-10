// The one-species closure from (rho, rho_e) through gaslib, and the same
// arithmetic and validity check on plain doubles, as separate functions for comparing their
// generated code.

#include <cmath>
#include <gaslib/models/calorically_perfect.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/state.hpp>

// The functions have external linkage so that the compiler emits them.
// NOLINTBEGIN(misc-use-internal-linkage)
struct RawClosure {
    double rho;
    double e;
    double p;
    double temp;
    unsigned char status;  // 0 valid, 3 invalid, as gaslib::Status
};

struct RawGas {
    double gas_constant;
    double gamma;
};

gaslib::Closure<double, gaslib::si>
  gaslib_close(const gaslib::CaloricallyPerfect<double, gaslib::si>& gas,
               gaslib::Density<double, gaslib::si> rho,
               gaslib::EnergyDensity<double, gaslib::si> rho_e)
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
