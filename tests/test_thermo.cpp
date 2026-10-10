#include "dual.hpp"
#include "test_models.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <gaslib/quantity.hpp>
#include <gaslib/species.hpp>
#include <gaslib/state.hpp>
#include <gaslib/status.hpp>
#include <gaslib/thermo.hpp>
#include <mdspan/mdspan.hpp>

namespace
{

using gaslib::si;
using gaslib::Status;
using D = gaslib::test::Dual<double>;
using Rho = gaslib::Density<double, si>;
using E = gaslib::SpecificEnergy<double, si>;
using RhoE = gaslib::EnergyDensity<double, si>;
using P = gaslib::Pressure<double, si>;
using Temp = gaslib::Temperature<double, si>;

// Air-like parameters: R = 287, cv = 717.5.
constexpr gaslib::test::IdealParameters air {
  .gas_constant = gaslib::SpecificGasConstant<double, si> {287.0},
  .cv = gaslib::SpecificHeat<double, si> {717.5},
};

bool identical(const gaslib::Closure<double, si>& a, const gaslib::Closure<double, si>& b)
{
    return a.state.rho == b.state.rho && a.state.e == b.state.e && a.state.p == b.state.p &&
           a.state.temp == b.state.temp && a.status == b.status;
}

using ModeView = gaslib::ModeEnergyDensityView<double, si>;
using ModeTempView = gaslib::ModeTemperatureView<double, si>;

}  // namespace

static_assert(gaslib::StateVariable<Rho> && gaslib::StateVariable<P> &&
              gaslib::StateVariable<Temp>);
static_assert(gaslib::StateVariable<E> && gaslib::StateVariable<RhoE> &&
              gaslib::StateVariable<ModeView>);
static_assert(gaslib::StateVariable<ModeTempView>);
static_assert(gaslib::detail::state_rank<ModeTempView>() == 4);
static_assert(gaslib::detail::state_rank<ModeTempView>() == gaslib::detail::state_rank<ModeView>());
// A mode temperature view goes after both state variables.
static_assert(gaslib::detail::argument_order<ModeTempView, Temp, Rho>() ==
              std::array<std::size_t, 3> {2, 1, 0});
static_assert(gaslib::detail::argument_order<P, ModeTempView, Temp>() ==
              std::array<std::size_t, 3> {0, 2, 1});
static_assert(gaslib::ModeRangeOf<ModeTempView, gaslib::kind::mode_temperature, double, si>);
static_assert(!gaslib::ModeRangeOf<ModeTempView, gaslib::kind::mode_temperature, float, si>);
static_assert(
  !gaslib::ModeRangeOf<ModeTempView, gaslib::kind::mode_temperature, double, gaslib::scaled>);
static_assert(!gaslib::ModeRangeOf<ModeView, gaslib::kind::mode_temperature, double, si>);
static_assert(!gaslib::StateVariable<gaslib::MassFraction<double>>);
static_assert(!gaslib::StateVariable<gaslib::test::TemperatureGuess>);

static_assert(gaslib::GasModel<gaslib::test::IdealGas, double>);
static_assert(gaslib::GasModel<gaslib::test::IdealGas, D>);
static_assert(gaslib::GasModel<gaslib::test::IdealGas, float>);
static_assert(gaslib::GasModel<gaslib::test::GuessGas, double>);
static_assert(!gaslib::GasModel<gaslib::test::NoPressureTemperatureGas, double>);
static_assert(!gaslib::GasModel<gaslib::test::GuessRequiredGas, double>);
static_assert(!gaslib::GasModel<gaslib::test::IdealGas, int>);
static_assert(gaslib::GasModel<gaslib::test::TwoSpeciesGas, double>);
static_assert(gaslib::MixtureGasModel<gaslib::test::TwoSpeciesGas, double>);
static_assert(!gaslib::MixtureGasModel<gaslib::test::IdealGas, double>);
static_assert(!gaslib::MultiTemperatureGasModel<gaslib::test::IdealGas, double>);
static_assert(!gaslib::MultiTemperatureGasModel<gaslib::test::TwoSpeciesGas, double>);

static_assert(gaslib::ClosesFrom<gaslib::test::IdealGas, double, Rho, E>);
static_assert(gaslib::ClosesFrom<gaslib::test::IdealGas, double, E, Rho>);
static_assert(gaslib::ClosesFrom<gaslib::test::IdealGas, double, RhoE, Rho>);
static_assert(!gaslib::ClosesFrom<gaslib::test::IdealGas, double, Rho, Rho>);
static_assert(!gaslib::ClosesFrom<gaslib::test::IdealGas, double, E, RhoE>);
static_assert(!gaslib::ClosesFrom<gaslib::test::NoPressureTemperatureGas, double, Temp, P>);

static_assert(gaslib::AcceptsInput<gaslib::test::GuessGas, gaslib::test::TemperatureGuess>);
static_assert(gaslib::AcceptsInput<gaslib::test::GuessRequiredGas, gaslib::test::TemperatureGuess>);
static_assert(!gaslib::AcceptsInput<gaslib::test::IdealGas, gaslib::test::TemperatureGuess>);
static_assert(!gaslib::AcceptsInput<gaslib::test::GuessGas, double>);

void test_thermo()
{
    using namespace boost::ut;
    using gaslib::thermo_state;

    const gaslib::test::IdealGas gas {air};
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
                           thermo_state(gas, rho, gaslib::specific_energy(rho_e, rho))));
      };

    "a model's own energy-density closure is selected over the generic one"_test = [&] {
        const gaslib::test::OwnEnergyDensityGas own {air};
        expect(thermo_state(own, rho, rho_e).status == Status::CompositionViolation);
        expect(thermo_state(own, rho_e, rho).status == Status::CompositionViolation);
        expect(thermo_state(own, rho, e).status == Status::Valid);
    };

    "an extra input reaches the model's overload that takes it, after the state variables"_test =
      [&] {
          const gaslib::test::GuessGas guess_gas {air};
          const gaslib::test::TemperatureGuess guess {Temp {290.0}};
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
                                    gaslib::SpecificEnergy<D, si> {D(2.15e5)},
                                    gaslib::Density<D, si> {D(1.2, 1.0)});
        expect(c.state.p.value().dual() == 287.0 * c.state.temp.value().real());
    };

    "a float state with double parameters closes and agrees with the double result"_test = [&] {
        const auto c = thermo_state(gas,
                                    gaslib::Density<float, si> {1.2F},
                                    gaslib::SpecificEnergy<float, si> {2.15e5F});
        static_assert(std::same_as<decltype(c), const gaslib::Closure<float, si>>);
        const auto reference = thermo_state(gas, Rho {1.2F}, E {2.15e5F});
        const double tol = 2.0 * 0x1p-23;
        expect(std::abs(c.state.p.value() - reference.state.p.value()) <=
               tol * reference.state.p.value());
        expect(std::abs(c.state.temp.value() - reference.state.temp.value()) <=
               tol * reference.state.temp.value());
    };

    "a mixture model takes its mixture value before the state variables"_test = [&] {
        const gaslib::test::TwoSpeciesGas two {.species = {air, air}};
        const std::array y {0.25, 0.75};
        const auto mix = gaslib::test::mixture(
          two,
          gaslib::as_kind<gaslib::kind::mass_fraction>(
            Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data(), 2)));
        const auto c = thermo_state(two, mix, e, rho);
        expect(identical(c, thermo_state(two, mix, rho, e)));
        expect(identical(thermo_state(two, mix, rho_e, rho), thermo_state(two, mix, rho, e)));
        expect(std::abs(c.state.p.value() - thermo_state(gas, rho, e).state.p.value()) <=
               4.0 * 0x1p-52 * c.state.p.value());
    };
}
