#include "dual.hpp"
#include "test_suites.hpp"

#include <boost/ut.hpp>
#include <concepts>
#include <gaslib/quantity.hpp>
#include <gaslib/state.hpp>
#include <gaslib/status.hpp>
#include <type_traits>

namespace
{

using gaslib::scaled;
using gaslib::si;
using D = gaslib::test::Dual<double>;

// A state without a temperature member.
struct NoTemperature {
    using value_type = double;
    using units = si;
    gaslib::Density<double, si> rho;
    gaslib::SpecificEnergy<double, si> e;
    gaslib::Pressure<double, si> p;
};

}  // namespace

static_assert(gaslib::ThermoStateOf<gaslib::ThermoState<double, si>, double, si>);
static_assert(gaslib::ThermoStateOf<gaslib::ThermoState<D, scaled>, D, scaled>);
static_assert(!gaslib::ThermoStateOf<gaslib::ThermoState<double, si>, double, scaled>);
static_assert(!gaslib::ThermoStateOf<gaslib::ThermoState<double, si>, float, si>);
static_assert(!gaslib::ThermoStateOf<NoTemperature, double, si>);
static_assert(std::is_trivially_copyable_v<gaslib::ThermoState<double, si>>);
static_assert(std::is_trivially_copyable_v<gaslib::Closure<double, si>>);
static_assert(std::is_trivially_copyable_v<gaslib::Properties<float, scaled>>);
static_assert(sizeof(gaslib::ThermoState<double, si>) == 4 * sizeof(double));
static_assert(
  std::same_as<gaslib::Closure<double, si>, gaslib::ClosureOf<gaslib::ThermoState<double, si>>>);

void test_state()
{
    using namespace boost::ut;

    "a closure holds the state members and the status"_test = [] {
        const gaslib::Closure<double, si> c {
          .state =
            {
              .rho = gaslib::Density<double, si> {1.2},
              .e = gaslib::SpecificEnergy<double, si> {2.0e5},
              .p = gaslib::Pressure<double, si> {1.0e5},
              .temp = gaslib::Temperature<double, si> {290.0},
            },
          .status = gaslib::Status::Valid,
        };
        expect(c.state.rho.value() == 1.2_d);
        expect(c.state.temp.value() == 290.0_d);
        expect(c.status == gaslib::Status::Valid);
    };
}
