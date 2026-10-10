// A viscosity law in SI units does not accept a scaled temperature.
#include <gaslib/transport/sutherland.hpp>

namespace
{
constexpr gaslib::Sutherland<double, gaslib::si> law {
  .mu_0 = gaslib::DynamicViscosity<double, gaslib::si> {1.716e-5},
  .t_0 = gaslib::Temperature<double, gaslib::si> {273.15},
  .s = gaslib::Temperature<double, gaslib::si> {110.4},
};
}  // namespace

#ifdef GASLIB_COMPILE_FAIL
auto bad = gaslib::viscosity(law, gaslib::Temperature<double, gaslib::scaled> {1.0});
#endif
auto good = gaslib::viscosity(law, gaslib::Temperature<double, gaslib::si> {300.0});
