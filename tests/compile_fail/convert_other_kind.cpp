// Conversion keeps the kind: an energy density does not become a pressure.
#include <gaslib/units.hpp>

constexpr gaslib::UnitScales<double, gaslib::scaled> scales {1.0, 1.0, 1.0, 1.0};

#ifdef GASLIB_COMPILE_FAIL
gaslib::Pressure<double, gaslib::scaled> bad =
  gaslib::convert<gaslib::scaled>(gaslib::EnergyDensity<double, gaslib::si> {1.0}, scales);
#endif
gaslib::Pressure<double, gaslib::scaled> good =
  gaslib::convert<gaslib::scaled>(gaslib::Pressure<double, gaslib::si> {1.0}, scales);
