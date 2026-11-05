#include "dual.hpp"
#include "test_models.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <canary/detail/closure.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/status.hpp>
#include <canary/thermo.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <mdspan/mdspan.hpp>

namespace
{

using canary::si;
using canary::Status;
using D = canary::test::Dual<double>;
using Rho = canary::Density<double, si>;
using E = canary::SpecificEnergy<double, si>;
using RhoE = canary::EnergyDensity<double, si>;
using P = canary::Pressure<double, si>;
using Temp = canary::Temperature<double, si>;

// Air-like parameters: R = 287, cv = 717.5.
constexpr canary::test::IdealParameters air {
  .gas_constant = canary::SpecificGasConstant<double, si> {287.0},
  .cv = canary::SpecificHeat<double, si> {717.5},
};

bool identical(const canary::Closure<double, si>& a, const canary::Closure<double, si>& b)
{
    return a.state.rho == b.state.rho && a.state.e == b.state.e && a.state.p == b.state.p &&
           a.state.temp == b.state.temp && a.status == b.status &&
           a.diagnostic.reason == b.diagnostic.reason && a.diagnostic.mode == b.diagnostic.mode;
}

using ModeView = canary::ModeEnergyDensityView<double, si>;
using ModeTempView = canary::ModeTemperatureView<double, si>;

}  // namespace

static_assert(canary::StateVariable<Rho> && canary::StateVariable<P> &&
              canary::StateVariable<Temp>);
static_assert(canary::StateVariable<E> && canary::StateVariable<RhoE> &&
              canary::StateVariable<ModeView>);
static_assert(canary::StateVariable<ModeTempView>);
static_assert(canary::detail::state_rank<ModeTempView>() == 4);
static_assert(canary::detail::state_rank<ModeTempView>() == canary::detail::state_rank<ModeView>());
// A mode temperature view goes after both state variables.
static_assert(canary::detail::argument_order<ModeTempView, Temp, Rho>() ==
              std::array<std::size_t, 3> {2, 1, 0});
static_assert(canary::detail::argument_order<P, ModeTempView, Temp>() ==
              std::array<std::size_t, 3> {0, 2, 1});
static_assert(canary::ModeRangeOf<ModeTempView, canary::kind::mode_temperature, double, si>);
static_assert(!canary::ModeRangeOf<ModeTempView, canary::kind::mode_temperature, float, si>);
static_assert(
  !canary::ModeRangeOf<ModeTempView, canary::kind::mode_temperature, double, canary::scaled>);
static_assert(!canary::ModeRangeOf<ModeView, canary::kind::mode_temperature, double, si>);
using TypedModes = canary::IndexedView<canary::ModeIndex, ModeTempView>;
using WrongModes = canary::IndexedView<canary::SpeciesIndex, ModeTempView>;
static_assert(canary::ModeRangeOf<TypedModes, canary::kind::mode_temperature, double, si>);
static_assert(!canary::ModeRangeOf<WrongModes, canary::kind::mode_temperature, double, si>);
static_assert(
  !canary::ModeRangeOf<TypedModes, canary::kind::mode_temperature, double, canary::scaled>);
static_assert(canary::detail::state_rank<TypedModes>() == 4);
static_assert(canary::detail::argument_order<TypedModes, Temp, Rho>() ==
              std::array<std::size_t, 3> {2, 1, 0});
static_assert(!canary::StateVariable<canary::MassFraction<double>>);
static_assert(!canary::StateVariable<canary::test::TemperatureGuess>);

static_assert(canary::GasModel<canary::test::IdealGas, double>);
static_assert(canary::GasModel<canary::test::IdealGas, D>);
static_assert(canary::GasModel<canary::test::IdealGas, float>);
static_assert(canary::GasModel<canary::test::GuessGas, double>);
static_assert(!canary::GasModel<canary::test::NoPressureTemperatureGas, double>);
static_assert(!canary::GasModel<canary::test::GuessRequiredGas, double>);
static_assert(!canary::GasModel<canary::test::IdealGas, int>);
static_assert(canary::GasModel<canary::test::TwoSpeciesGas, double>);
static_assert(canary::MixtureGasModel<canary::test::TwoSpeciesGas, double>);
static_assert(!canary::MixtureGasModel<canary::test::IdealGas, double>);
static_assert(!canary::MultiTemperatureGasModel<canary::test::IdealGas, double>);
static_assert(!canary::MultiTemperatureGasModel<canary::test::TwoSpeciesGas, double>);

static_assert(canary::ClosesFrom<canary::test::IdealGas, double, Rho, E>);
static_assert(canary::ClosesFrom<canary::test::IdealGas, double, E, Rho>);
static_assert(canary::ClosesFrom<canary::test::IdealGas, double, RhoE, Rho>);
static_assert(!canary::ClosesFrom<canary::test::IdealGas, double, Rho, Rho>);
static_assert(!canary::ClosesFrom<canary::test::IdealGas, double, E, RhoE>);
static_assert(!canary::ClosesFrom<canary::test::NoPressureTemperatureGas, double, Temp, P>);

static_assert(canary::AcceptsInput<canary::test::GuessGas, canary::test::TemperatureGuess>);
static_assert(canary::AcceptsInput<canary::test::GuessRequiredGas, canary::test::TemperatureGuess>);
static_assert(!canary::AcceptsInput<canary::test::IdealGas, canary::test::TemperatureGuess>);
static_assert(!canary::AcceptsInput<canary::test::GuessGas, double>);

// combine keeps the more severe check in either order and the first check on a
// tie, with its reason and mode.
constexpr canary::detail::ClosureCheck valid_check {};
constexpr canary::detail::ClosureCheck invalid_density {
  .status = Status::Invalid,
  .diagnostic = {.reason = canary::ClosureReason::NonfiniteDensity, .mode = {}},
};
constexpr canary::detail::ClosureCheck invalid_mode_energy {
  .status = Status::Invalid,
  .diagnostic = {.reason = canary::ClosureReason::AbsentModeEnergy,
                 .mode = canary::OptionalModeIndex {canary::ModeIndex {1}}},
};
constexpr canary::detail::ClosureCheck not_converged {
  .status = Status::NotConverged,
  .diagnostic = {.reason = canary::ClosureReason::IterationLimit, .mode = {}},
};
constexpr bool same(canary::detail::ClosureCheck a, canary::detail::ClosureCheck b)
{
    return a.status == b.status && a.diagnostic.reason == b.diagnostic.reason &&
           a.diagnostic.mode == b.diagnostic.mode;
}
static_assert(same(canary::detail::combine(valid_check, valid_check), valid_check));
static_assert(same(canary::detail::combine(valid_check, invalid_mode_energy), invalid_mode_energy));
static_assert(same(canary::detail::combine(invalid_mode_energy, valid_check), invalid_mode_energy));
static_assert(same(canary::detail::combine(not_converged, invalid_density), invalid_density));
static_assert(same(canary::detail::combine(invalid_density, not_converged), invalid_density));
static_assert(same(canary::detail::combine(invalid_density, invalid_mode_energy), invalid_density));
static_assert(same(canary::detail::combine(invalid_mode_energy, invalid_density),
                   invalid_mode_energy));

void test_thermo()
{
    using namespace boost::ut;
    using canary::thermo_state;

    "finite_check is Valid for a finite value part and Invalid with the reason otherwise"_test =
      [] {
          using canary::detail::finite_check;
          const auto reason = canary::ClosureReason::NonfiniteEnergy;
          const double nan = std::numeric_limits<double>::quiet_NaN();
          const double inf = std::numeric_limits<double>::infinity();
          expect(same(finite_check(1.0, reason), valid_check));
          expect(same(finite_check(D {1.0, nan}, reason), valid_check));
          for (const double bad : {nan, inf, -inf}) {
              const canary::detail::ClosureCheck expected {
                .status = Status::Invalid,
                .diagnostic = {.reason = reason, .mode = {}},
              };
              expect(same(finite_check(bad, reason), expected));
              expect(same(finite_check(D {bad, 1.0}, reason), expected));
          }
      };

    const canary::test::IdealGas gas {air};
    const Rho rho {1.2};
    const E e {2.15e5};
    const RhoE rho_e {1.2 * 2.15e5};
    const P p {101325.0};
    const Temp temp {300.0};

    "both argument orders give identical results for every pair"_test = [&] {
        expect(identical(thermo_state(gas, rho, e), thermo_state(gas, e, rho)));
        expect(identical(thermo_state(gas, rho, rho_e), thermo_state(gas, rho_e, rho)));
        expect(identical(thermo_state(gas, rho, p), thermo_state(gas, p, rho)));
        expect(identical(thermo_state(gas, rho, temp), thermo_state(gas, temp, rho)));
        expect(identical(thermo_state(gas, p, temp), thermo_state(gas, temp, p)));
    };

    "the generic energy-density closure equals the specific-energy closure at rho_e / rho"_test =
      [&] {
          expect(identical(thermo_state(gas, rho, rho_e),
                           thermo_state(gas, rho, canary::specific_energy(rho_e, rho))));
      };

    "a model's own energy-density closure is selected over the generic one"_test = [&] {
        const canary::test::OwnEnergyDensityGas own {air};
        expect(thermo_state(own, rho, rho_e).status == Status::CompositionViolation);
        expect(thermo_state(own, rho_e, rho).status == Status::CompositionViolation);
        expect(thermo_state(own, rho, e).status == Status::Valid);
    };

    "an extra input reaches the model's overload that takes it, after the state variables"_test =
      [&] {
          const canary::test::GuessGas guess_gas {air};
          const canary::test::TemperatureGuess guess {Temp {290.0}};
          expect(thermo_state(guess_gas, rho, e, guess).status == Status::NotConverged);
          expect(thermo_state(guess_gas, e, rho, guess).status == Status::NotConverged);
          expect(thermo_state(guess_gas, rho, e).status == Status::Valid);
          // The generic energy-density closure passes the extra input on.
          expect(thermo_state(guess_gas, rho_e, rho, guess).status == Status::NotConverged);
      };

    "closures are consistent: (rho, e) -> (p, T) -> back"_test = [&] {
        const auto c = thermo_state(gas, rho, e);
        const auto from_rho_p = thermo_state(gas, rho, c.state.p);
        const auto from_rho_t = thermo_state(gas, rho, c.state.temp);
        const auto from_p_t = thermo_state(gas, c.state.p, c.state.temp);
        const double tol = 4.0 * 0x1p-52;
        expect(std::abs(from_rho_p.state.e.value() - e.value()) <= tol * e.value());
        expect(std::abs(from_rho_t.state.e.value() - e.value()) <= tol * e.value());
        expect(std::abs(from_p_t.state.rho.value() - rho.value()) <= tol * rho.value());
        expect(std::abs(from_p_t.state.e.value() - e.value()) <= tol * e.value());
    };

    "dual-number states pass through the front end with their derivative"_test = [&] {
        // d p / d rho at constant e is R T.
        const auto c = thermo_state(gas,
                                    canary::SpecificEnergy<D, si> {D(2.15e5)},
                                    canary::Density<D, si> {D(1.2, 1.0)});
        expect(c.state.p.value().dual() == 287.0 * c.state.temp.value().real());
    };

    "a float state with double parameters closes and agrees with the double result"_test = [&] {
        const auto c = thermo_state(gas,
                                    canary::Density<float, si> {1.2F},
                                    canary::SpecificEnergy<float, si> {2.15e5F});
        static_assert(std::same_as<decltype(c), const canary::Closure<float, si>>);
        const auto reference = thermo_state(gas, Rho {1.2F}, E {2.15e5F});
        const double tol = 2.0 * 0x1p-23;
        expect(std::abs(c.state.p.value() - reference.state.p.value()) <=
               tol * reference.state.p.value());
        expect(std::abs(c.state.temp.value() - reference.state.temp.value()) <=
               tol * reference.state.temp.value());
    };

    "a mixture model takes its mixture value before the state variables"_test = [&] {
        const canary::test::TwoSpeciesGas two {.species = {air, air}};
        const std::array y {0.25, 0.75};
        const auto mix = canary::test::mixture(
          two,
          canary::as_kind<canary::kind::mass_fraction>(
            Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data(), 2)));
        const auto c = thermo_state(two, mix, e, rho);
        expect(identical(c, thermo_state(two, mix, rho, e)));
        expect(identical(thermo_state(two, mix, rho_e, rho), thermo_state(two, mix, rho, e)));
        expect(std::abs(c.state.p.value() - thermo_state(gas, rho, e).state.p.value()) <=
               4.0 * 0x1p-52 * c.state.p.value());
    };
}
