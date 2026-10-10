#include "dual.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <concepts>
#include <gaslib/build.hpp>
#include <gaslib/models/calorically_perfect.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/state.hpp>
#include <gaslib/status.hpp>
#include <gaslib/thermo.hpp>
#include <gaslib/units.hpp>
#include <limits>
#include <type_traits>

namespace
{

using gaslib::scaled;
using gaslib::si;
using gaslib::Status;
using D = gaslib::test::Dual<double>;
using Rho = gaslib::Density<double, si>;
using E = gaslib::SpecificEnergy<double, si>;
using RhoE = gaslib::EnergyDensity<double, si>;
using P = gaslib::Pressure<double, si>;
using Temp = gaslib::Temperature<double, si>;

constexpr double eps = std::numeric_limits<double>::epsilon();

bool within_ulp(double a, double b, double ulps)
{ return std::abs(a - b) <= ulps * eps * std::abs(b); }
bool within_relative(double a, double b, double tol)
{ return std::abs(a - b) <= tol * std::abs(b); }

constexpr std::array densities {1.0e-6, 1.0, 1.0e3};
constexpr std::array energies {1.0e2, 2.5e5, 1.0e8};

}  // namespace

static_assert(gaslib::GasModel<gaslib::CaloricallyPerfect<double, si>, double>);
static_assert(gaslib::GasModel<gaslib::CaloricallyPerfect<double, scaled>, double>);
static_assert(gaslib::GasModel<gaslib::CaloricallyPerfect<double, si>, D>);
static_assert(gaslib::GasModel<gaslib::CaloricallyPerfect<double, si>, float>);
static_assert(gaslib::GasModel<gaslib::CaloricallyPerfect<float, si>, double>);
static_assert(gaslib::GasModel<gaslib::CaloricallyPerfect<float, si>, float>);
static_assert(!gaslib::MixtureGasModel<gaslib::CaloricallyPerfect<double, si>, double>);
static_assert(std::is_trivially_copyable_v<gaslib::CaloricallyPerfect<double, si>>);
static_assert(sizeof(gaslib::CaloricallyPerfect<double, si>) == 4 * sizeof(double));

void test_calorically_perfect()
{
    using namespace boost::ut;
    using gaslib::thermo_state;

    // SI air: R = 287, gamma = 1.4.
    const gaslib::CaloricallyPerfect<double, si> air =
      gaslib::make_calorically_perfect(gaslib::SpecificGasConstant<double, si> {287.0},
                                       gaslib::HeatCapacityRatio {1.4});

    "(rho, e) -> (p, T) -> (rho, e) round trips through every closure"_test = [&] {
        for (const double rho_v : densities) {
            for (const double e_v : energies) {
                const auto c = thermo_state(air, Rho {rho_v}, E {e_v});
                expect(c.status == Status::Valid);
                const auto from_rho_p = thermo_state(air, Rho {rho_v}, c.state.p);
                const auto from_rho_t = thermo_state(air, Rho {rho_v}, c.state.temp);
                const auto from_p_t = thermo_state(air, c.state.p, c.state.temp);
                expect(within_ulp(from_rho_p.state.e.value(), e_v, 4.0)) << rho_v << e_v;
                expect(within_ulp(from_rho_t.state.e.value(), e_v, 4.0)) << rho_v << e_v;
                expect(within_ulp(from_p_t.state.e.value(), e_v, 4.0)) << rho_v << e_v;
                expect(within_ulp(from_p_t.state.rho.value(), rho_v, 4.0)) << rho_v << e_v;
                expect(within_ulp(from_rho_t.state.p.value(), c.state.p.value(), 4.0))
                  << rho_v << e_v;
            }
        }
    };

    "the energy-density closure agrees with the specific-energy closure"_test = [&] {
        for (const double rho_v : densities) {
            for (const double e_v : energies) {
                const auto a = thermo_state(air, Rho {rho_v}, E {e_v});
                const auto b = thermo_state(air, Rho {rho_v}, RhoE {rho_v * e_v});
                expect(a.state.p == b.state.p && a.state.temp == b.state.temp);
                expect(within_ulp(b.state.e.value(), e_v, 2.0));
            }
        }
    };

    "cp - cv = R and cp / cv = gamma"_test = [&] {
        expect(within_ulp((air.cp - air.cv).value(), air.gas_constant.value(), 2.0))
          << (air.cp - air.cv).value() - air.gas_constant.value();
        expect(within_ulp(air.cp.value() / air.cv.value(), air.gamma.value(), 2.0));
    };

    "a^2 = gamma p / rho and equals dp/drho along an isentrope"_test = [&] {
        const auto c = thermo_state(air, Rho {1.2}, E {2.15e5});
        const auto props = gaslib::properties(air, c.state);
        const double a2 = props.sound_speed.value() * props.sound_speed.value();
        expect(within_ulp(a2, 1.4 * c.state.p.value() / 1.2, 4.0));
        // p = K rho^gamma along the isentrope through the state.
        const double k = c.state.p.value() / std::pow(1.2, 1.4);
        const double h = 1.0e-4 * 1.2;
        const double dp_drho =
          ((k * std::pow(1.2 + h, 1.4)) - (k * std::pow(1.2 - h, 1.4))) / (2.0 * h);
        expect(within_relative(a2, dp_drho, 1.0e-8)) << a2 << dp_drho;
    };

    "dual numbers give dp/drho at constant e = (gamma - 1) e and dT/de = 1/cv"_test = [&] {
        const double e_v = 2.15e5;
        const auto by_rho = thermo_state(air,
                                         gaslib::Density<D, si> {D(1.2, 1.0)},
                                         gaslib::SpecificEnergy<D, si> {D(e_v)});
        expect(within_ulp(by_rho.state.p.value().dual(), 0.4 * e_v, 2.0))
          << by_rho.state.p.value().dual();
        const auto by_e = thermo_state(air,
                                       gaslib::Density<D, si> {D(1.2)},
                                       gaslib::SpecificEnergy<D, si> {D(e_v, 1.0)});
        expect(within_ulp(by_e.state.temp.value().dual(), 1.0 / air.cv.value(), 2.0))
          << by_e.state.temp.value().dual();
        // Central differences.
        const double h_rho = 1.0e-5 * 1.2;
        const double dp_fd = (thermo_state(air, Rho {1.2 + h_rho}, E {e_v}).state.p.value() -
                              thermo_state(air, Rho {1.2 - h_rho}, E {e_v}).state.p.value()) /
                             (2.0 * h_rho);
        const double h_e = 1.0e-5 * e_v;
        const double dt_fd = (thermo_state(air, Rho {1.2}, E {e_v + h_e}).state.temp.value() -
                              thermo_state(air, Rho {1.2}, E {e_v - h_e}).state.temp.value()) /
                             (2.0 * h_e);
        expect(within_relative(by_rho.state.p.value().dual(), dp_fd, 1.0e-7));
        expect(within_relative(by_e.state.temp.value().dual(), dt_fd, 1.0e-7));
    };

    "non-positive or non-finite inputs give Invalid; valid inputs give positive results"_test =
      [&] {
          const double nan = std::numeric_limits<double>::quiet_NaN();
          expect(thermo_state(air, Rho {1.0}, E {0.0}).status == Status::Invalid);
          expect(thermo_state(air, Rho {1.0}, E {-1.0}).status == Status::Invalid);
          expect(thermo_state(air, Rho {0.0}, E {1.0}).status == Status::Invalid);
          expect(thermo_state(air, Rho {-1.0}, E {1.0}).status == Status::Invalid);
          expect(thermo_state(air, Rho {nan}, E {1.0}).status == Status::Invalid);
          expect(thermo_state(air, Rho {1.0}, P {nan}).status == Status::Invalid);
          expect(thermo_state(air, Rho {1.0}, Temp {-1.0}).status == Status::Invalid);
          expect(thermo_state(air, P {0.0}, Temp {300.0}).status == Status::Invalid);
          for (const double rho_v : densities) {
              for (const double e_v : energies) {
                  const auto c = thermo_state(air, Rho {rho_v}, E {e_v});
                  const auto props = gaslib::properties(air, c.state);
                  expect(c.status == Status::Valid);
                  expect(c.state.p.value() > 0.0 && c.state.temp.value() > 0.0);
                  expect(props.sound_speed.value() > 0.0 && props.cv.value() > 0.0 &&
                         props.cp.value() > 0.0);
              }
          }
      };

    "float and double parameter and state types agree within float precision"_test = [&] {
        const auto air_f =
          gaslib::make_calorically_perfect(gaslib::SpecificGasConstant<float, si> {287.0F},
                                           gaslib::HeatCapacityRatio {1.4F});
        const auto reference = thermo_state(air, Rho {1.2}, E {2.15e5});
        const auto f_params = thermo_state(air_f, Rho {1.2}, E {2.15e5});
        const auto f_state = thermo_state(air,
                                          gaslib::Density<float, si> {1.2F},
                                          gaslib::SpecificEnergy<float, si> {2.15e5F});
        static_assert(std::same_as<decltype(f_params), const gaslib::Closure<double, si>>);
        static_assert(std::same_as<decltype(f_state), const gaslib::Closure<float, si>>);
        const double tol = 4.0 * std::numeric_limits<float>::epsilon();
        expect(within_relative(f_params.state.p.value(), reference.state.p.value(), tol));
        expect(within_relative(f_params.state.temp.value(), reference.state.temp.value(), tol));
        expect(within_relative(f_state.state.p.value(), reference.state.p.value(), tol));
        expect(within_relative(f_state.state.temp.value(), reference.state.temp.value(), tol));
    };

    "a scaled model from cerberus-like inputs has R = k_B / m with k_B = 1"_test = [&] {
        const gaslib::UnitScales<double, scaled> scales {
          .mass = 1.6726e-27,
          .density = 1.0e-6,
          .temperature = 1.0e4,
          .length = 1.0e-3,
        };
        const auto built =
          gaslib::make_calorically_perfect(gaslib::ParticleMass<double, scaled> {2.0},
                                           gaslib::Charge<double, scaled> {1.0},
                                           gaslib::boltzmann_constant(scales),
                                           gaslib::HeatCapacityRatio {5.0 / 3.0});
        expect(built.model.gas_constant.value() == 0.5_d);
        expect(built.species.n == 1_i);
        expect(built.species.mass[0].value() == 2.0_d && built.species.charge[0].value() == 1.0_d);
        // p = rho T / m in scaled units.
        const auto c = thermo_state(built.model,
                                    gaslib::Density<double, scaled> {3.0},
                                    gaslib::Temperature<double, scaled> {0.25});
        expect(c.state.p.value() == 3.0 * 0.5 * 0.25);
    };

    "an SI model from particle mass and Boltzmann constant matches one from R"_test = [&] {
        const double m = 4.81e-26;
        const auto built = gaslib::make_calorically_perfect(
          gaslib::ParticleMass<double, si> {m},
          gaslib::Charge<double, si> {0.0},
          gaslib::boltzmann_constant(gaslib::UnitScales<double, si> {}),
          gaslib::HeatCapacityRatio {1.4});
        expect(built.model.gas_constant.value() == 1.380649e-23 / m);
    };

    "make_calorically_perfect rejects parameters out of range"_test = [&] {
        const auto from_r = [](double r, double gamma) {
            return gaslib::make_calorically_perfect(gaslib::SpecificGasConstant<double, si> {r},
                                                    gaslib::HeatCapacityRatio {gamma});
        };
        const auto from_m = [](double m, double q, double k) {
            return gaslib::make_calorically_perfect(gaslib::ParticleMass<double, si> {m},
                                                    gaslib::Charge<double, si> {q},
                                                    gaslib::BoltzmannConstant<double, si> {k},
                                                    gaslib::HeatCapacityRatio {1.4});
        };
        const double nan = std::numeric_limits<double>::quiet_NaN();
        expect(throws<gaslib::ParameterError>([&] { from_r(0.0, 1.4); }));
        expect(throws<gaslib::ParameterError>([&] { from_r(-287.0, 1.4); }));
        expect(throws<gaslib::ParameterError>([&] { from_r(287.0, 1.0); }));
        expect(throws<gaslib::ParameterError>([&] { from_r(287.0, nan); }));
        expect(throws<gaslib::ParameterError>([&] { from_m(0.0, 0.0, 1.0); }));
        expect(throws<gaslib::ParameterError>([&] { from_m(1.0, nan, 1.0); }));
        expect(throws<gaslib::ParameterError>([&] { from_m(1.0, 0.0, -1.0); }));
        expect(nothrow([&] { from_m(1.0, -1.0, 1.0); }));
    };
}
