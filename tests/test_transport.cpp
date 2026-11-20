// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#include "dual.hpp"
#include "test_models.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
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
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
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
static_assert(canary::TransportModel<SutherlandAir, canary::test::IdealGas, D>);

template <typename Props>
concept PrandtlProperties =
  requires(const Props& props) { canary::conductivity(prandtl, Mu {1}, props); };
static_assert(PrandtlProperties<canary::test::IdealProperties<double>>);
static_assert(!PrandtlProperties<canary::MultiTemperatureProperties<double, si, 2>>);
static_assert(!PrandtlProperties<canary::SingleTemperatureProperties<double, scaled>>);

void test_transport()
{
    using namespace boost::ut;

    "checked transport factories validate float and double parameter boundaries"_test = [] {
        const auto check = []<typename P> {
            using Viscosity = canary::DynamicViscosity<P, si>;
            using Temperature = canary::Temperature<P, si>;
            using Prandtl = canary::PrandtlNumber<P>;
            const P nan = std::numeric_limits<P>::quiet_NaN();
            const P inf = std::numeric_limits<P>::infinity();
            for (const P value : {P(-1), -inf, inf, nan}) {
                expect(throws<canary::ParameterError>(
                  [&] { (void)canary::make_constant_viscosity(Viscosity {value}); }));
                expect(throws<canary::ParameterError>([&] {
                    (void)canary::make_power_law(Viscosity {value}, Temperature {P(1)}, P(0));
                }));
                expect(throws<canary::ParameterError>([&] {
                    (void)canary::make_sutherland(Viscosity {value},
                                                  Temperature {P(1)},
                                                  Temperature {P(0)});
                }));
            }
            for (const P value : {P(0), P(-0.0), P(-1), -inf, inf, nan}) {
                expect(throws<canary::ParameterError>([&] {
                    (void)canary::make_power_law(Viscosity {P(1)}, Temperature {value}, P(0));
                }));
                expect(throws<canary::ParameterError>([&] {
                    (void)canary::make_sutherland(Viscosity {P(1)},
                                                  Temperature {value},
                                                  Temperature {P(0)});
                }));
                expect(throws<canary::ParameterError>(
                  [&] { (void)canary::make_constant_prandtl(Prandtl {value}); }));
            }
            for (const P value : {-inf, inf, nan}) {
                expect(throws<canary::ParameterError>([&] {
                    (void)canary::make_power_law(Viscosity {P(1)}, Temperature {P(1)}, value);
                }));
                expect(throws<canary::ParameterError>([&] {
                    (void)canary::make_sutherland(Viscosity {P(1)},
                                                  Temperature {P(1)},
                                                  Temperature {value});
                }));
            }
            for (const P value : {P(-1), std::nextafter(P(-1), -inf)}) {
                expect(throws<canary::ParameterError>([&] {
                    (void)canary::make_sutherland(Viscosity {P(1)},
                                                  Temperature {P(1)},
                                                  Temperature {value});
                }));
            }
            const P above = std::nextafter(P(-1), P(0));
            expect(
              canary::make_sutherland(Viscosity {P(1)}, Temperature {P(1)}, Temperature {above})
                .s.value() == above);
            for (const P value : {P(0), P(-0.0)}) {
                expect(canary::make_constant_viscosity(Viscosity {value}).mu.value() == P(0));
                expect(canary::make_power_law(Viscosity {value}, Temperature {P(1)}, P(-1))
                         .mu_0.value() == P(0));
                expect(
                  canary::make_sutherland(Viscosity {value}, Temperature {P(1)}, Temperature {P(0)})
                    .mu_0.value() == P(0));
            }
            expect(canary::make_constant_prandtl(Prandtl {std::numeric_limits<P>::min()})
                     .prandtl.value() == std::numeric_limits<P>::min());
        };
        check.operator()<float>();
        check.operator()<double>();
    };

    "checked transport factories preserve float parameter storage and scaled results"_test = [] {
        const canary::UnitScales<double, scaled> scales {
          {.mass = 4.65e-26, .density = 1.2, .temperature = 300.0, .length = 0.1}};
        const auto check = [&]<typename P> {
            using Viscosity = canary::DynamicViscosity<P, si>;
            using Temperature = canary::Temperature<P, si>;
            const Viscosity mu {P(1.716e-5)};
            const Temperature t0 {P(273.15)};
            const Temperature offset {P(110.4)};
            const canary::DynamicViscosity<P, scaled> scaled_mu {
              canary::convert<scaled>(Mu {double(mu.value())}, scales)};
            const canary::Temperature<P, scaled> scaled_t0 {
              canary::convert<scaled>(Temp {double(t0.value())}, scales)};
            const canary::Temperature<P, scaled> scaled_offset {
              canary::convert<scaled>(Temp {double(offset.value())}, scales)};
            const auto constant = canary::make_constant_viscosity(mu);
            const auto power = canary::make_power_law(mu, t0, P(0.76));
            const auto sutherland_law = canary::make_sutherland(mu, t0, offset);
            const auto prandtl_law =
              canary::make_constant_prandtl(canary::PrandtlNumber<P> {P(0.71)});
            static_assert(std::same_as<decltype(constant), const canary::ConstantViscosity<P, si>>);
            static_assert(std::same_as<decltype(power), const canary::PowerLaw<P, si>>);
            static_assert(std::same_as<decltype(sutherland_law), const canary::Sutherland<P, si>>);
            static_assert(std::same_as<decltype(prandtl_law), const canary::ConstantPrandtl<P>>);
            static_assert(std::is_aggregate_v<canary::ConstantViscosity<P, si>>);
            static_assert(std::is_aggregate_v<canary::PowerLaw<P, si>>);
            static_assert(std::is_aggregate_v<canary::Sutherland<P, si>>);
            static_assert(std::is_aggregate_v<canary::ConstantPrandtl<P>>);
            const auto scaled_constant = canary::make_constant_viscosity(scaled_mu);
            const auto scaled_power = canary::make_power_law(scaled_mu, scaled_t0, P(0.76));
            const auto scaled_sutherland =
              canary::make_sutherland(scaled_mu, scaled_t0, scaled_offset);
            const auto close = [](double actual, double expected) {
                return std::abs(actual - expected) <=
                       16.0 * std::numeric_limits<P>::epsilon() * std::abs(expected);
            };
            for (const P t : {P(100), P(300), P(1500)}) {
                const Temperature temp {t};
                const canary::Temperature<P, scaled> scaled_temp {
                  canary::convert<scaled>(Temp {double(t)}, scales)};
                const auto compare = [&](const auto& law, const auto& scaled_law) {
                    const auto scaled_value = canary::viscosity(scaled_law, scaled_temp);
                    const double physical =
                      canary::convert<si>(
                        canary::DynamicViscosity<double, scaled> {double(scaled_value.value())},
                        scales)
                        .value();
                    expect(close(physical, canary::viscosity(law, temp).value()));
                };
                compare(constant, scaled_constant);
                compare(power, scaled_power);
                compare(sutherland_law, scaled_sutherland);
            }
            canary::SingleTemperatureProperties<P, si> properties {};
            properties.cp = canary::SpecificHeat<P, si> {P(1004.5)};
            canary::SingleTemperatureProperties<P, scaled> scaled_properties {};
            scaled_properties.cp = canary::SpecificHeat<P, scaled> {canary::convert<scaled>(
              canary::SpecificHeat<double, si> {double(properties.cp.value())},
              scales)};
            const auto scaled_k = canary::conductivity(prandtl_law, scaled_mu, scaled_properties);
            const auto physical_k = canary::convert<si>(
              canary::ThermalConductivity<double, scaled> {double(scaled_k.value())},
              scales);
            expect(
              close(physical_k.value(), canary::conductivity(prandtl_law, mu, properties).value()));
        };
        check.operator()<float>();
        check.operator()<double>();
    };

    "checked transport factories retain analytic temperature derivatives"_test = [] {
        const auto constant = canary::make_constant_viscosity(Mu {1.716e-5});
        const auto power = canary::make_power_law(Mu {1.716e-5}, Temp {273.15}, 0.76);
        const auto sutherland_law =
          canary::make_sutherland(Mu {1.716e-5}, Temp {273.15}, Temp {110.4});
        const canary::Temperature<D, si> temp {D(300.0, 1.0)};
        expect(canary::viscosity(constant, temp).value().dual() == 0.0);
        const auto power_mu = canary::viscosity(power, temp).value();
        expect(within_ulp(power_mu.dual(), power_mu.real() * 0.76 / 300.0, 4.0));
        const auto sutherland_mu = canary::viscosity(sutherland_law, temp).value();
        expect(within_ulp(sutherland_mu.dual(),
                          sutherland_mu.real() * (1.5 / 300.0 - 1.0 / 410.4),
                          8.0));
    };

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

    "transport accepts caller-owned single-temperature properties"_test = [] {
        const canary::test::IdealGas gas {{canary::SpecificGasConstant<double, si> {287.0},
                                           canary::SpecificHeat<double, si> {717.5}}};
        const auto state = canary::thermo_state(gas, canary::Density<double, si> {1.2}, Temp {300});
        const auto coefficients =
          canary::transport(SutherlandAir {sutherland, prandtl}, gas, state.state);
        expect(within_ulp(coefficients.conductivity.value(),
                          canary::viscosity(sutherland, Temp {300}).value() * 1004.5 / 0.71,
                          2.0));
    };

    "constant Prandtl gives k = mu cp / Pr"_test = [] {
        const canary::SingleTemperatureProperties<double, si> props {
          .cv = canary::SpecificHeat<double, si> {717.5},
          .cp = canary::SpecificHeat<double, si> {1004.5},
          .capacity_ratio = canary::HeatCapacityRatio {1.4},
          .gas_constant = canary::SpecificGasConstant<double, si> {287.0},
          .acoustics = {.sound_speed = canary::SoundSpeed<double, si> {347.0},
                        .acoustic_exponent = canary::AcousticExponent {1.4}},
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
            const canary::SingleTemperatureProperties<double, si> props {
              .cv = canary::SpecificHeat<double, si> {cp / 1.4},
              .cp = canary::SpecificHeat<double, si> {cp},
              .capacity_ratio = canary::HeatCapacityRatio {1.4},
              .gas_constant = canary::SpecificGasConstant<double, si> {cp - (cp / 1.4)},
              .acoustics = {.sound_speed = canary::SoundSpeed<double, si> {1.0},
                            .acoustic_exponent = canary::AcousticExponent {1.4}},
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
