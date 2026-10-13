#include "dual.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <canary/models/calorically_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/thermo.hpp>
#include <canary/transport/constant.hpp>
#include <canary/transport/power_law.hpp>
#include <canary/transport/prandtl.hpp>
#include <canary/transport/sutherland.hpp>
#include <canary/transport/transport.hpp>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <type_traits>

namespace
{

using canary::scaled;
using canary::si;
using D = canary::test::Dual<double>;
using Mu = canary::DynamicViscosity<double, si>;
using Temp = canary::Temperature<double, si>;

constexpr double eps = std::numeric_limits<double>::epsilon();

bool within_ulp(double a, double b, double ulps)
{ return std::abs(a - b) <= ulps * eps * std::abs(b); }

// Air: mu_0 = 1.716e-5 Pa s at T_0 = 273.15 K, S = 110.4 K.
constexpr canary::Sutherland<double, si> sutherland {
  .mu_0 = Mu {1.716e-5},
  .t_0 = Temp {273.15},
  .s = Temp {110.4},
};
constexpr canary::PowerLaw<double, si> power_law {
  .mu_0 = Mu {1.716e-5},
  .t_0 = Temp {273.15},
  .exponent = 0.76,
};
constexpr canary::ConstantViscosity<double, si> constant_mu {.mu = Mu {1.8e-5}};
constexpr canary::ConstantPrandtl<double> prandtl {.prandtl = canary::PrandtlNumber {0.71}};

using SutherlandAir =
  canary::Transport<canary::Sutherland<double, si>, canary::ConstantPrandtl<double>>;
using Air = canary::CaloricallyPerfect<double, si>;
using AirMixture = canary::CaloricallyPerfectMixture<double, si, 2>;

}  // namespace

static_assert(canary::ViscosityModel<canary::Sutherland<double, si>, double, si>);
static_assert(canary::ViscosityModel<canary::Sutherland<double, si>, D, si>);
static_assert(canary::ViscosityModel<canary::Sutherland<double, si>, float, si>);
static_assert(canary::ViscosityModel<canary::PowerLaw<float, scaled>, double, scaled>);
static_assert(canary::ViscosityModel<canary::ConstantViscosity<double, si>, D, si>);
static_assert(!canary::ViscosityModel<canary::Sutherland<double, si>, double, scaled>);
static_assert(!canary::ViscosityModel<canary::ConstantPrandtl<double>, double, si>);
static_assert(canary::ConductivityModel<canary::ConstantPrandtl<double>, double, si>);
static_assert(canary::ConductivityModel<canary::ConstantPrandtl<double>, D, scaled>);
static_assert(!canary::ConductivityModel<canary::Sutherland<double, si>, double, si>);
static_assert(canary::TransportModel<SutherlandAir, Air, double>);
static_assert(canary::TransportModel<SutherlandAir, Air, D>);
static_assert(canary::TransportModel<SutherlandAir, AirMixture, double>);
static_assert(
  !canary::TransportModel<SutherlandAir, canary::CaloricallyPerfect<double, scaled>, double>);
static_assert(std::is_trivially_copyable_v<SutherlandAir>);

void test_transport()
{
    using namespace boost::ut;

    "Sutherland's law at T_0 gives mu_0"_test = [] {
        expect(within_ulp(canary::viscosity(sutherland, Temp {273.15}).value(), 1.716e-5, 2.0));
    };

    "the power-law log-derivative d ln mu / d ln T equals the exponent"_test = [] {
        for (const double t : {100.0, 273.15, 1500.0}) {
            const auto mu =
              canary::viscosity(power_law, canary::Temperature<D, si> {D(t, 1.0)}).value();
            expect(within_ulp(mu.dual() * t / mu.real(), 0.76, 4.0)) << t;
        }
    };

    "the Sutherland derivative from dual numbers equals the analytic derivative"_test = [] {
        for (const double t : {100.0, 273.15, 1500.0}) {
            const auto mu =
              canary::viscosity(sutherland, canary::Temperature<D, si> {D(t, 1.0)}).value();
            const double analytic = mu.real() * ((1.5 / t) - (1.0 / (t + 110.4)));
            expect(within_ulp(mu.dual(), analytic, 8.0)) << t << mu.dual() << analytic;
        }
    };

    "constant viscosity ignores the temperature and has zero derivative"_test = [] {
        const auto mu =
          canary::viscosity(constant_mu, canary::Temperature<D, si> {D(300.0, 1.0)}).value();
        expect(mu.real() == 1.8e-5 && mu.dual() == 0.0);
    };

    "constant Prandtl gives k = mu cp / Pr"_test = [] {
        const canary::Properties<double, si> props {
          .cv = canary::SpecificHeat<double, si> {717.5},
          .cp = canary::SpecificHeat<double, si> {1004.5},
          .gamma = canary::HeatCapacityRatio {1.4},
          .sound_speed = canary::SoundSpeed<double, si> {347.0},
          .gas_constant = canary::SpecificGasConstant<double, si> {287.0},
        };
        const auto k = canary::conductivity(prandtl, Mu {1.8e-5}, props);
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
            expect(within_ulp(canary::viscosity(sutherland, Temp {t}).value(), mu_sutherland, 2.0))
              << t;
            expect(within_ulp(canary::viscosity(power_law, Temp {t}).value(), mu_power, 2.0)) << t;
            const canary::Properties<double, si> props {
              .cv = canary::SpecificHeat<double, si> {cp / 1.4},
              .cp = canary::SpecificHeat<double, si> {cp},
              .gamma = canary::HeatCapacityRatio {1.4},
              .sound_speed = canary::SoundSpeed<double, si> {1.0},
              .gas_constant = canary::SpecificGasConstant<double, si> {cp - (cp / 1.4)},
            };
            expect(within_ulp(canary::conductivity(prandtl, Mu {mu_sutherland}, props).value(),
                              mu_sutherland * cp / pr,
                              2.0))
              << t;
        }
    };

    "transport at a state uses the state temperature and the gas cp"_test = [] {
        const auto air =
          canary::make_calorically_perfect(canary::SpecificGasConstant<double, si> {287.0},
                                           canary::HeatCapacityRatio {1.4});
        const SutherlandAir model {.viscosity = sutherland, .conductivity = prandtl};
        const auto c = canary::thermo_state(air, canary::Density<double, si> {1.2}, Temp {300.0});
        const auto coefficients = canary::transport(model, air, c.state);
        const double mu = canary::viscosity(sutherland, Temp {300.0}).value();
        expect(coefficients.viscosity.value() == mu);
        expect(within_ulp(coefficients.conductivity.value(), mu * air.cp.value() / 0.71, 2.0));
    };

    "transport for a mixture uses the mixture cp"_test = [] {
        using M = canary::ParticleMass<double, si>;
        using Q = canary::Charge<double, si>;
        using G = canary::HeatCapacityRatio<double>;
        const auto gas = canary::make_calorically_perfect_mixture<2>(
          std::array {M {4.65e-26}, M {6.63e-26}},
          std::array {Q {0.0}, Q {0.0}},
          canary::BoltzmannConstant<double, si> {1.380649e-23},
          std::array {G {1.4}, G {5.0 / 3.0}});
        const std::array y {0.7, 0.3};
        const auto mix = canary::mixture(
          gas.model,
          canary::as_kind<canary::kind::mass_fraction>(
            Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>(y.data(), 2)));
        const auto c =
          canary::thermo_state(gas.model, mix, canary::Density<double, si> {1.2}, Temp {300.0});
        const SutherlandAir model {.viscosity = sutherland, .conductivity = prandtl};
        const auto coefficients = canary::transport(model, gas.model, mix, c.state);
        expect(within_ulp(coefficients.conductivity.value(),
                          coefficients.viscosity.value() * mix.cp.value() / 0.71,
                          2.0));
    };

    "float states with double laws agree within float precision"_test = [] {
        const double reference = canary::viscosity(sutherland, Temp {300.0}).value();
        const auto mu = canary::viscosity(sutherland, canary::Temperature<float, si> {300.0F});
        static_assert(std::same_as<decltype(mu), const canary::DynamicViscosity<float, si>>);
        expect(std::abs(mu.value() - reference) <=
               4.0 * std::numeric_limits<float>::epsilon() * reference);
    };
}
