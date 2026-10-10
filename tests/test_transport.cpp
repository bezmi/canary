#include "dual.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <gaslib/models/calorically_perfect.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/species.hpp>
#include <gaslib/state.hpp>
#include <gaslib/thermo.hpp>
#include <gaslib/transport/constant.hpp>
#include <gaslib/transport/power_law.hpp>
#include <gaslib/transport/prandtl.hpp>
#include <gaslib/transport/sutherland.hpp>
#include <gaslib/transport/transport.hpp>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <type_traits>

namespace
{

using gaslib::scaled;
using gaslib::si;
using D = gaslib::test::Dual<double>;
using Mu = gaslib::DynamicViscosity<double, si>;
using Temp = gaslib::Temperature<double, si>;

constexpr double eps = std::numeric_limits<double>::epsilon();

bool within_ulp(double a, double b, double ulps)
{ return std::abs(a - b) <= ulps * eps * std::abs(b); }

// Air: mu_0 = 1.716e-5 Pa s at T_0 = 273.15 K, S = 110.4 K.
constexpr gaslib::Sutherland<double, si> sutherland {
  .mu_0 = Mu {1.716e-5},
  .t_0 = Temp {273.15},
  .s = Temp {110.4},
};
constexpr gaslib::PowerLaw<double, si> power_law {
  .mu_0 = Mu {1.716e-5},
  .t_0 = Temp {273.15},
  .exponent = 0.76,
};
constexpr gaslib::ConstantViscosity<double, si> constant_mu {.mu = Mu {1.8e-5}};
constexpr gaslib::ConstantPrandtl<double> prandtl {.prandtl = gaslib::PrandtlNumber {0.71}};

using SutherlandAir =
  gaslib::Transport<gaslib::Sutherland<double, si>, gaslib::ConstantPrandtl<double>>;
using Air = gaslib::CaloricallyPerfect<double, si>;
using AirMixture = gaslib::CaloricallyPerfectMixture<double, si, 2>;

}  // namespace

static_assert(gaslib::ViscosityModel<gaslib::Sutherland<double, si>, double, si>);
static_assert(gaslib::ViscosityModel<gaslib::Sutherland<double, si>, D, si>);
static_assert(gaslib::ViscosityModel<gaslib::Sutherland<double, si>, float, si>);
static_assert(gaslib::ViscosityModel<gaslib::PowerLaw<float, scaled>, double, scaled>);
static_assert(gaslib::ViscosityModel<gaslib::ConstantViscosity<double, si>, D, si>);
static_assert(!gaslib::ViscosityModel<gaslib::Sutherland<double, si>, double, scaled>);
static_assert(!gaslib::ViscosityModel<gaslib::ConstantPrandtl<double>, double, si>);
static_assert(gaslib::ConductivityModel<gaslib::ConstantPrandtl<double>, double, si>);
static_assert(gaslib::ConductivityModel<gaslib::ConstantPrandtl<double>, D, scaled>);
static_assert(!gaslib::ConductivityModel<gaslib::Sutherland<double, si>, double, si>);
static_assert(gaslib::TransportModel<SutherlandAir, Air, double>);
static_assert(gaslib::TransportModel<SutherlandAir, Air, D>);
static_assert(gaslib::TransportModel<SutherlandAir, AirMixture, double>);
static_assert(
  !gaslib::TransportModel<SutherlandAir, gaslib::CaloricallyPerfect<double, scaled>, double>);
static_assert(std::is_trivially_copyable_v<SutherlandAir>);

void test_transport()
{
    using namespace boost::ut;

    "Sutherland's law at T_0 gives mu_0"_test = [] {
        expect(within_ulp(gaslib::viscosity(sutherland, Temp {273.15}).value(), 1.716e-5, 2.0));
    };

    "the power-law log-derivative d ln mu / d ln T equals the exponent"_test = [] {
        for (const double t : {100.0, 273.15, 1500.0}) {
            const auto mu =
              gaslib::viscosity(power_law, gaslib::Temperature<D, si> {D(t, 1.0)}).value();
            expect(within_ulp(mu.dual() * t / mu.real(), 0.76, 4.0)) << t;
        }
    };

    "the Sutherland derivative from dual numbers equals the analytic derivative"_test = [] {
        for (const double t : {100.0, 273.15, 1500.0}) {
            const auto mu =
              gaslib::viscosity(sutherland, gaslib::Temperature<D, si> {D(t, 1.0)}).value();
            const double analytic = mu.real() * ((1.5 / t) - (1.0 / (t + 110.4)));
            expect(within_ulp(mu.dual(), analytic, 8.0)) << t << mu.dual() << analytic;
        }
    };

    "constant viscosity ignores the temperature and has zero derivative"_test = [] {
        const auto mu =
          gaslib::viscosity(constant_mu, gaslib::Temperature<D, si> {D(300.0, 1.0)}).value();
        expect(mu.real() == 1.8e-5 && mu.dual() == 0.0);
    };

    "constant Prandtl gives k = mu cp / Pr"_test = [] {
        const gaslib::Properties<double, si> props {
          .cv = gaslib::SpecificHeat<double, si> {717.5},
          .cp = gaslib::SpecificHeat<double, si> {1004.5},
          .gamma = gaslib::HeatCapacityRatio {1.4},
          .sound_speed = gaslib::SoundSpeed<double, si> {347.0},
          .gas_constant = gaslib::SpecificGasConstant<double, si> {287.0},
        };
        const auto k = gaslib::conductivity(prandtl, Mu {1.8e-5}, props);
        expect(within_ulp(k.value(), 1.8e-5 * 1004.5 / 0.71, 2.0));
    };

    // The viscous models in cerberus Source/states/Eulerian/hydro/viscous/
    // MFP_hydro_viscous_impl.H (lines 15-83) at commit 9d45fc5.
    "values equal cerberus's viscous model formulas"_test = [] {
        const double mu_0 = 1.716e-5;
        const double t0 = 273.15;
        const double s = 110.4;
        const double n = 0.76;
        const double pr = 0.71;
        const double cp = 1004.5;
        for (const double t : {50.0, 273.15, 300.0, 2500.0}) {
            const double t_ = t / t0;
            const double mu_sutherland = mu_0 * t_ * std::sqrt(t_) * (t0 + s) / (t + s);
            const double mu_power = mu_0 * std::pow(t / t0, n);
            expect(within_ulp(gaslib::viscosity(sutherland, Temp {t}).value(), mu_sutherland, 2.0))
              << t;
            expect(within_ulp(gaslib::viscosity(power_law, Temp {t}).value(), mu_power, 2.0)) << t;
            const gaslib::Properties<double, si> props {
              .cv = gaslib::SpecificHeat<double, si> {cp / 1.4},
              .cp = gaslib::SpecificHeat<double, si> {cp},
              .gamma = gaslib::HeatCapacityRatio {1.4},
              .sound_speed = gaslib::SoundSpeed<double, si> {1.0},
              .gas_constant = gaslib::SpecificGasConstant<double, si> {cp - (cp / 1.4)},
            };
            expect(within_ulp(gaslib::conductivity(prandtl, Mu {mu_sutherland}, props).value(),
                              mu_sutherland * cp / pr,
                              2.0))
              << t;
        }
    };

    "transport at a state uses the state temperature and the gas cp"_test = [] {
        const auto air =
          gaslib::make_calorically_perfect(gaslib::SpecificGasConstant<double, si> {287.0},
                                           gaslib::HeatCapacityRatio {1.4});
        const SutherlandAir model {.viscosity = sutherland, .conductivity = prandtl};
        const auto c = gaslib::thermo_state(air, gaslib::Density<double, si> {1.2}, Temp {300.0});
        const auto coefficients = gaslib::transport(model, air, c.state);
        const double mu = gaslib::viscosity(sutherland, Temp {300.0}).value();
        expect(coefficients.viscosity.value() == mu);
        expect(within_ulp(coefficients.conductivity.value(), mu * air.cp.value() / 0.71, 2.0));
    };

    "transport for a mixture uses the mixture cp"_test = [] {
        using M = gaslib::ParticleMass<double, si>;
        using Q = gaslib::Charge<double, si>;
        using G = gaslib::HeatCapacityRatio<double>;
        const auto gas = gaslib::make_calorically_perfect_mixture<2>(
          std::array {M {4.65e-26}, M {6.63e-26}},
          std::array {Q {0.0}, Q {0.0}},
          gaslib::BoltzmannConstant<double, si> {1.380649e-23},
          std::array {G {1.4}, G {5.0 / 3.0}});
        const std::array y {0.7, 0.3};
        const auto mix = gaslib::mixture(
          gas.model,
          gaslib::as_kind<gaslib::kind::mass_fraction>(
            Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data(), 2)));
        const auto c =
          gaslib::thermo_state(gas.model, mix, gaslib::Density<double, si> {1.2}, Temp {300.0});
        const SutherlandAir model {.viscosity = sutherland, .conductivity = prandtl};
        const auto coefficients = gaslib::transport(model, gas.model, mix, c.state);
        expect(within_ulp(coefficients.conductivity.value(),
                          coefficients.viscosity.value() * mix.cp.value() / 0.71,
                          2.0));
    };

    "float states with double laws agree within float precision"_test = [] {
        const double reference = gaslib::viscosity(sutherland, Temp {300.0}).value();
        const auto mu = gaslib::viscosity(sutherland, gaslib::Temperature<float, si> {300.0F});
        static_assert(std::same_as<decltype(mu), const gaslib::DynamicViscosity<float, si>>);
        expect(std::abs(mu.value() - reference) <=
               4.0 * std::numeric_limits<float>::epsilon() * reference);
    };
}
