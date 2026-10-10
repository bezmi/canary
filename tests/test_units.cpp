#include "dual.hpp"
#include "test_suites.hpp"

#include <boost/ut.hpp>
#include <cmath>
#include <concepts>
#include <gaslib/build.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/units.hpp>
#include <limits>

namespace
{

using gaslib::scaled;
using gaslib::si;
namespace kind = gaslib::kind;

constexpr double eps = std::numeric_limits<double>::epsilon();

// Cerberus-like references: a proton mass, a low density, 1e4 K and 1 mm.
constexpr gaslib::UnitScales<double, scaled> scales {
  .mass = 1.6726e-27,
  .density = 1.0e-6,
  .temperature = 1.0e4,
  .length = 1.0e-3,
};

bool within_ulp(double a, double b, double ulps)
{ return std::abs(a - b) <= ulps * eps * std::abs(b); }

// SI -> scaled -> SI for one dimensional kind, at several magnitudes.
template <gaslib::KindTag K> bool round_trips()
{
    bool ok = true;
    for (const double v : {1.0e-30, 3.7e-7, 1.0, 2.5e5, 6.1e23}) {
        const gaslib::Quantity<K, double, si> q {v};
        const auto back = gaslib::convert<si>(gaslib::convert<scaled>(q, scales), scales);
        static_assert(std::same_as<decltype(back), const gaslib::Quantity<K, double, si>>);
        ok = ok && within_ulp(back.value(), v, 2.0);
    }
    return ok;
}

template <typename Q, typename To>
concept ConvertsTo = requires(Q q) { gaslib::convert<To>(q, scales); };

}  // namespace

static_assert(ConvertsTo<gaslib::Density<double, si>, scaled>);
static_assert(ConvertsTo<gaslib::Density<double, scaled>, si>);
static_assert(ConvertsTo<gaslib::MolarMass<double>, si>);
static_assert(!ConvertsTo<gaslib::MolarMass<double>, scaled>);
static_assert(!ConvertsTo<gaslib::UniversalGasConstant<double>, scaled>);
static_assert(ConvertsTo<gaslib::MassFraction<double>, scaled>);

void test_units()
{
    using namespace boost::ut;
    using gaslib::convert;
    using gaslib::scale_factor;

    const double u_ref = std::sqrt(1.380649e-23 * scales.temperature / scales.mass);
    const double n0 =
      (scales.density / scales.mass) * scales.length * scales.length * scales.length;

    "the Boltzmann constant is exactly 1 in the scaled system"_test = [] {
        expect(gaslib::boltzmann_constant(scales).value() == 1.0_d);
        expect(gaslib::boltzmann_constant(gaslib::UnitScales<double, si> {}).value() ==
               1.380649e-23);
        constexpr gaslib::UnitScales<float, scaled> unit_float {
          .mass = 1.0F,
          .density = 1.0F,
          .temperature = 1.0F,
          .length = 1.0F,
        };
        expect(gaslib::boltzmann_constant(unit_float).value() == 1.0_f);
    };

    "converting the SI Boltzmann constant gives 1 within rounding"_test = [] {
        const auto k =
          convert<scaled>(gaslib::boltzmann_constant(gaslib::UnitScales<double, si> {}), scales);
        expect(within_ulp(k.value(), 1.0, 8.0)) << k.value() - 1.0;
    };

    "derived scales follow the cerberus reference definitions"_test = [u_ref, n0] {
        expect(gaslib::velocity_scale(scales) == u_ref);
        expect(gaslib::time_scale(scales) == scales.length / u_ref);
        expect(gaslib::count_scale(scales) == n0);
        expect(within_ulp(scale_factor<kind::density>(scales), scales.density, 4.0));
        expect(within_ulp(scale_factor<kind::particle_mass>(scales), scales.mass, 4.0));
        expect(within_ulp(scale_factor<kind::temperature>(scales), scales.temperature, 0.0));
        expect(within_ulp(scale_factor<kind::specific_energy>(scales), u_ref * u_ref, 4.0));
        expect(
          within_ulp(scale_factor<kind::pressure>(scales), scales.density * u_ref * u_ref, 4.0));
        expect(within_ulp(scale_factor<kind::energy_density>(scales),
                          scales.density * u_ref * u_ref,
                          4.0));
        expect(within_ulp(scale_factor<kind::dynamic_viscosity>(scales),
                          scales.density * u_ref * scales.length,
                          4.0));
        expect(within_ulp(scale_factor<kind::sound_speed>(scales), u_ref, 4.0));
        expect(within_ulp(scale_factor<kind::charge>(scales), 1.602176634e-19, 4.0));
        expect(scale_factor<kind::mass_fraction>(scales) == 1.0_d);
    };

    "SI to scaled to SI round trips within 2 ulp for every dimensional kind"_test = [] {
        expect(round_trips<kind::density>());
        expect(round_trips<kind::pressure>());
        expect(round_trips<kind::temperature>());
        expect(round_trips<kind::specific_energy>());
        expect(round_trips<kind::energy_density>());
        expect(round_trips<kind::specific_heat>());
        expect(round_trips<kind::specific_gas_constant>());
        expect(round_trips<kind::sound_speed>());
        expect(round_trips<kind::particle_mass>());
        expect(round_trips<kind::charge>());
        expect(round_trips<kind::boltzmann_constant>());
        expect(round_trips<kind::dynamic_viscosity>());
        expect(round_trips<kind::thermal_conductivity>());
        expect(round_trips<kind::mode_energy_density>());
        expect(round_trips<kind::mode_temperature>());
        expect(round_trips<kind::number_rate_density>());
        expect(round_trips<kind::mass_production_rate>());
    };

    "the scale factor of a dimension equals that of a kind of the same dimension"_test = [u_ref,
                                                                                          n0] {
        expect(scale_factor(kind::pressure::dimension, scales) ==
               scale_factor<kind::pressure>(scales));
        expect(scale_factor(kind::number_rate_density::dimension, scales) ==
               scale_factor<kind::number_rate_density>(scales));
        // Rate constant of a third-order reaction, (length^3 / count)^2 / time.
        const gaslib::Dimension third_order {.length = 6, .time = -1, .count = -2};
        const double x = scales.length;
        const double t_ref = x / u_ref;
        expect(within_ulp(scale_factor(third_order, scales),
                          (x * x * x * x * x * x) / (n0 * n0 * t_ref),
                          8.0));
        expect(scale_factor(third_order, gaslib::UnitScales<double, si> {}) == 1.0_d);
        expect(scale_factor(kind::molar_mass::dimension, gaslib::UnitScales<double, si> {}) ==
               1.0_d);
        expect(throws<gaslib::ParameterError>(
          [] { (void)scale_factor(kind::molar_mass::dimension, scales); }));
    };

    "SI and unitless quantities convert to themselves"_test = [] {
        const gaslib::MolarMass<double> molar {0.0289647};
        expect(convert<si>(molar, gaslib::UnitScales<double, si> {}).value() == 0.0289647);
        expect(convert<si>(molar, scales).value() == 0.0289647);
        const gaslib::MassFraction<double> y {0.3};
        expect(convert<scaled>(y, scales).value() == 0.3);
    };

    "the ideal-gas law holds with k_B = 1 after conversion"_test = [] {
        // p = rho k_B T / m in SI; in scaled units p = rho T / m.
        const double rho = 2.0e-6;
        const double temp = 2.3e4;
        const double m = 2.0 * 1.6726e-27;
        const double p = rho * 1.380649e-23 * temp / m;
        const double rho_s = convert<scaled>(gaslib::Density<double, si> {rho}, scales).value();
        const double temp_s =
          convert<scaled>(gaslib::Temperature<double, si> {temp}, scales).value();
        const double m_s = convert<scaled>(gaslib::ParticleMass<double, si> {m}, scales).value();
        const double p_s = convert<scaled>(gaslib::Pressure<double, si> {p}, scales).value();
        expect(within_ulp(p_s, rho_s * temp_s / m_s, 16.0)) << p_s << rho_s * temp_s / m_s;
    };

    "a dual-number quantity converts with its derivative"_test = [] {
        using D = gaslib::test::Dual<double>;
        const gaslib::Density<D, si> rho {D(2.0e-6, 1.0)};
        const auto rho_s = convert<scaled>(rho, scales);
        expect(within_ulp(rho_s.value().real(), 2.0, 4.0));
        expect(within_ulp(rho_s.value().dual(), 1.0 / scales.density, 4.0));
    };
}
