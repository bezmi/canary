#include "dual.hpp"
#include "test_suites.hpp"

#include <boost/ut.hpp>
#include <concepts>
#include <canary/quantity.hpp>
#include <canary/state.hpp>
#include <canary/status.hpp>
#include <type_traits>

namespace
{

using canary::scaled;
using canary::si;
using D = canary::test::Dual<double>;

// A state without a temperature member.
struct NoTemperature {
    using value_type = double;
    using units = si;
    canary::Density<double, si> rho;
    canary::SpecificEnergy<double, si> e;
    canary::Pressure<double, si> p;
};

}  // namespace

static_assert(canary::ThermoStateOf<canary::ThermoState<double, si>, double, si>);
static_assert(canary::ThermoStateOf<canary::ThermoState<D, scaled>, D, scaled>);
static_assert(!canary::ThermoStateOf<canary::ThermoState<double, si>, double, scaled>);
static_assert(!canary::ThermoStateOf<canary::ThermoState<double, si>, float, si>);
static_assert(!canary::ThermoStateOf<NoTemperature, double, si>);
static_assert(std::is_trivially_copyable_v<canary::ThermoState<double, si>>);
static_assert(std::is_trivially_copyable_v<canary::Closure<double, si>>);
static_assert(std::is_trivially_copyable_v<canary::Properties<float, scaled>>);
static_assert(canary::ThermoStateOf<canary::MultiTemperatureState<double, si, 1>, double, si>);
static_assert(canary::ThermoStateOf<canary::MultiTemperatureState<D, scaled, 2>, D, scaled>);
static_assert(!canary::ThermoStateOf<canary::MultiTemperatureState<double, si, 1>, double, scaled>);
static_assert(std::is_trivially_copyable_v<canary::MultiTemperatureState<double, si, 1>>);
static_assert(sizeof(canary::MultiTemperatureState<double, si, 2>) == 8 * sizeof(double));
static_assert(sizeof(canary::ThermoState<double, si>) == 4 * sizeof(double));
static_assert(
  std::same_as<canary::Closure<double, si>, canary::ClosureOf<canary::ThermoState<double, si>>>);

void test_state()
{
    using namespace boost::ut;

    "a closure holds the state members and the status"_test = [] {
        const canary::Closure<double, si> c {
          .state =
            {
              .rho = canary::Density<double, si> {1.2},
              .e = canary::SpecificEnergy<double, si> {2.0e5},
              .p = canary::Pressure<double, si> {1.0e5},
              .temp = canary::Temperature<double, si> {290.0},
            },
          .status = canary::Status::Valid,
        };
        expect(c.state.rho.value() == 1.2_d);
        expect(c.state.temp.value() == 290.0_d);
        expect(c.status == canary::Status::Valid);
    };
}
