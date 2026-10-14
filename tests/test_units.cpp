#include "dual.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/quantity.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <limits>
#include <string>
#include <type_traits>

namespace
{

using canary::scaled;
using canary::si;
namespace kind = canary::kind;

constexpr double eps = std::numeric_limits<double>::epsilon();

// Cerberus-like references: a proton mass, a low density, 1e4 K and 1 mm.
constexpr canary::UnitScales<double, scaled> scales {{
  .mass = 1.6726e-27,
  .density = 1.0e-6,
  .temperature = 1.0e4,
  .length = 1.0e-3,
}};

bool within_ulp(double a, double b, double ulps)
{ return std::abs(a - b) <= ulps * eps * std::abs(b); }

// SI -> scaled -> SI for one dimensional kind, at several magnitudes.
template <canary::KindTag K> bool round_trips()
{
    bool ok = true;
    for (const double v : {1.0e-30, 3.7e-7, 1.0, 2.5e5, 6.1e23}) {
        const canary::Quantity<K, double, si> q {v};
        const auto back = canary::convert<si>(canary::convert<scaled>(q, scales), scales);
        static_assert(std::same_as<decltype(back), const canary::Quantity<K, double, si>>);
        ok = ok && within_ulp(back.value(), v, 2.0);
    }
    return ok;
}

template <typename Q, typename To>
concept ConvertsTo = requires(Q q) { canary::convert<To>(q, scales); };

}  // namespace

static_assert(ConvertsTo<canary::Density<double, si>, scaled>);
static_assert(ConvertsTo<canary::Density<double, scaled>, si>);
static_assert(ConvertsTo<canary::MolarMass<double>, si>);
static_assert(!ConvertsTo<canary::MolarMass<double>, scaled>);
static_assert(!ConvertsTo<canary::UniversalGasConstant<double>, scaled>);
static_assert(ConvertsTo<canary::MassFraction<double>, scaled>);

void test_units()
{
    using namespace boost::ut;
    using canary::convert;
    using canary::scale_factor;

    "unit scale construction rejects invalid or omitted base values"_test = [] {
        const auto check = []<typename P>() {
            using Scales = canary::UnitScales<P, scaled>;
            using Input = typename Scales::Input;
            const std::array members {&Input::mass,
                                      &Input::density,
                                      &Input::temperature,
                                      &Input::length};
            const std::array names {"mass", "density", "temperature", "length"};
            for (std::size_t i = 0; i < members.size(); ++i) {
                for (P value : {P(0),
                                -P(0),
                                P(-1),
                                std::numeric_limits<P>::infinity(),
                                -std::numeric_limits<P>::infinity(),
                                std::numeric_limits<P>::quiet_NaN()}) {
                    Input input {.mass = P(1),
                                 .density = P(1),
                                 .temperature = P(1),
                                 .length = P(1)};
                    input.*members[i] = value;
                    std::string message;
                    try {
                        (void)Scales {input};
                    } catch (const canary::ParameterError& error) {
                        message = error.what();
                    }
                    expect(message.starts_with(std::string(names[i]) +
                                               " scale must be finite and positive"))
                      << message;
                }
            }
            expect(throws<canary::ParameterError>(
              [] { (void)Scales {{.mass = P(1), .density = P(1), .temperature = P(1)}}; }));
        };
        check.template operator()<float>();
        check.template operator()<double>();
        check.template operator()<long double>();
    };

    "unit scales own their validated input and retain copy and assignment"_test = [] {
        using Scales = canary::UnitScales<double, scaled>;
        static_assert(!std::is_aggregate_v<Scales> && !std::is_default_constructible_v<Scales>);
        static_assert(std::is_trivially_copyable_v<Scales> && std::copyable<Scales>);
        static_assert(sizeof(Scales) == 4 * sizeof(double));
        static_assert(std::same_as<decltype(scales.mass()), double>);
        static_assert(std::same_as<decltype(scales.density()), double>);
        static_assert(std::same_as<decltype(scales.temperature()), double>);
        static_assert(std::same_as<decltype(scales.length()), double>);
        static_assert(scales.mass() == 1.6726e-27 && scales.density() == 1.0e-6);
        Scales::Input input {.mass = 2.0, .density = 3.0, .temperature = 4.0, .length = 5.0};
        const Scales original {input};
        input.mass = 9.0;
        const Scales copied = original;
        Scales assigned {input};
        assigned = copied;
        expect(original.mass() == 2.0_d && copied.mass() == 2.0_d && assigned.mass() == 2.0_d);
        expect(assigned.density() == 3.0_d && assigned.temperature() == 4.0_d &&
               assigned.length() == 5.0_d);
    };

    const double u_ref = std::sqrt(1.380649e-23 * scales.temperature() / scales.mass());
    const double n0 =
      (scales.density() / scales.mass()) * scales.length() * scales.length() * scales.length();

    "the Boltzmann constant is exactly 1 in the scaled system"_test = [] {
        expect(canary::boltzmann_constant(scales).value() == 1.0_d);
        expect(canary::boltzmann_constant(canary::UnitScales<double, si> {}).value() ==
               1.380649e-23);
        constexpr canary::UnitScales<float, scaled> unit_float {{
          .mass = 1.0F,
          .density = 1.0F,
          .temperature = 1.0F,
          .length = 1.0F,
        }};
        expect(canary::boltzmann_constant(unit_float).value() == 1.0_f);
    };

    "converting the SI Boltzmann constant gives 1 within rounding"_test = [] {
        const auto k =
          convert<scaled>(canary::boltzmann_constant(canary::UnitScales<double, si> {}), scales);
        expect(within_ulp(k.value(), 1.0, 8.0)) << k.value() - 1.0;
    };

    "derived scales follow the cerberus reference definitions"_test = [u_ref, n0] {
        expect(canary::velocity_scale(scales) == u_ref);
        expect(canary::time_scale(scales) == scales.length() / u_ref);
        expect(canary::count_scale(scales) == n0);
        expect(within_ulp(scale_factor<kind::density>(scales), scales.density(), 4.0));
        expect(within_ulp(scale_factor<kind::particle_mass>(scales), scales.mass(), 4.0));
        expect(within_ulp(scale_factor<kind::temperature>(scales), scales.temperature(), 0.0));
        expect(within_ulp(scale_factor<kind::specific_energy>(scales), u_ref * u_ref, 4.0));
        expect(
          within_ulp(scale_factor<kind::pressure>(scales), scales.density() * u_ref * u_ref, 4.0));
        expect(within_ulp(scale_factor<kind::energy_density>(scales),
                          scales.density() * u_ref * u_ref,
                          4.0));
        expect(within_ulp(scale_factor<kind::power_density>(scales),
                          scales.density() * u_ref * u_ref * u_ref / scales.length(),
                          8.0));
        expect(within_ulp(scale_factor<kind::dynamic_viscosity>(scales),
                          scales.density() * u_ref * scales.length(),
                          4.0));
        expect(within_ulp(scale_factor<kind::sound_speed>(scales), u_ref, 4.0));
        expect(within_ulp(scale_factor<kind::diffusion_coefficient>(scales),
                          scales.length() * u_ref,
                          4.0));
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
        expect(round_trips<kind::power_density>());
        expect(round_trips<kind::diffusion_coefficient>());
    };

    "the scale factor of a dimension equals that of a kind of the same dimension"_test = [u_ref,
                                                                                          n0] {
        expect(scale_factor(kind::pressure::dimension, scales) ==
               scale_factor<kind::pressure>(scales));
        expect(scale_factor(kind::number_rate_density::dimension, scales) ==
               scale_factor<kind::number_rate_density>(scales));
        // Rate constant of a third-order reaction, (length^3 / count)^2 / time.
        const canary::Dimension third_order {.length = 6, .time = -1, .count = -2};
        const double x = scales.length();
        const double t_ref = x / u_ref;
        expect(within_ulp(scale_factor(third_order, scales),
                          (x * x * x * x * x * x) / (n0 * n0 * t_ref),
                          8.0));
        expect(scale_factor(third_order, canary::UnitScales<double, si> {}) == 1.0_d);
        expect(scale_factor(kind::molar_mass::dimension, canary::UnitScales<double, si> {}) ==
               1.0_d);
        expect(throws<canary::ParameterError>(
          [] { (void)scale_factor(kind::molar_mass::dimension, scales); }));
    };

    "SI and unitless quantities convert to themselves"_test = [] {
        const canary::MolarMass<double> molar {0.0289647};
        expect(convert<si>(molar, canary::UnitScales<double, si> {}).value() == 0.0289647);
        expect(convert<si>(molar, scales).value() == 0.0289647);
        const canary::MassFraction<double> y {0.3};
        expect(convert<scaled>(y, scales).value() == 0.3);
    };

    "the ideal-gas law holds with k_B = 1 after conversion"_test = [] {
        // p = rho k_B T / m in SI; in scaled units p = rho T / m.
        const double rho = 2.0e-6;
        const double temp = 2.3e4;
        const double m = 2.0 * 1.6726e-27;
        const double p = rho * 1.380649e-23 * temp / m;
        const double rho_s = convert<scaled>(canary::Density<double, si> {rho}, scales).value();
        const double temp_s =
          convert<scaled>(canary::Temperature<double, si> {temp}, scales).value();
        const double m_s = convert<scaled>(canary::ParticleMass<double, si> {m}, scales).value();
        const double p_s = convert<scaled>(canary::Pressure<double, si> {p}, scales).value();
        expect(within_ulp(p_s, rho_s * temp_s / m_s, 16.0)) << p_s << rho_s * temp_s / m_s;
    };

    "unit conversion retains each input scalar type across scale precisions"_test = [] {
        const auto check = []<typename V, typename P>() {
            const canary::UnitScales<P, scaled> reference {
              {.mass = P(1), .density = P(8), .temperature = P(4), .length = P(1)}};
            const auto converted = convert<scaled>(canary::Density<V, si> {V(16)}, reference);
            static_assert(std::same_as<decltype(converted), const canary::Density<V, scaled>>);
            expect(converted.value() == V(2));
            const auto back = convert<si>(converted, reference);
            static_assert(std::same_as<decltype(back), const canary::Density<V, si>>);
            expect(back.value() == V(16));

            using D = canary::test::Dual<V>;
            const auto dual = convert<scaled>(canary::Density<D, si> {D(V(16), V(-8))}, reference);
            static_assert(std::same_as<decltype(dual), const canary::Density<D, scaled>>);
            expect(dual.value().real() == V(2) && dual.value().dual() == V(-1));
            const auto dual_back = convert<si>(dual, reference);
            static_assert(std::same_as<decltype(dual_back), const canary::Density<D, si>>);
            expect(dual_back.value().real() == V(16) && dual_back.value().dual() == V(-8));

            const auto identity = convert<si>(canary::Density<D, si> {D(V(16), V(-8))}, reference);
            expect(identity.value().real() == V(16) && identity.value().dual() == V(-8));
            const auto unitless =
              convert<scaled>(canary::MassFraction<D> {D(V(0.5), V(-1))}, reference);
            static_assert(std::same_as<decltype(unitless), const canary::MassFraction<D>>);
            expect(unitless.value().real() == V(0.5) && unitless.value().dual() == V(-1));
        };
        const auto check_scales = [&]<typename V>() {
            check.template operator()<V, float>();
            check.template operator()<V, double>();
            check.template operator()<V, long double>();
        };
        check_scales.template operator()<float>();
        check_scales.template operator()<double>();
        check_scales.template operator()<long double>();
    };

    "float conversion widens intermediate particle counts"_test = [] {
        const canary::UnitScales<float, scaled> reference {
          {.mass = 1.0F, .density = 1.0e20F, .temperature = 1.0F, .length = 1.0e10F}};
        const auto rho =
          convert<scaled>(canary::Density<float, si> {reference.density()}, reference);
        expect(rho.value() == 1.0_f);
        using D = canary::test::Dual<float>;
        const auto dual =
          convert<scaled>(canary::Density<D, si> {D(reference.density(), -reference.density())},
                          reference);
        expect(dual.value().real() == 1.0_f && dual.value().dual() == -1.0_f);
    };

    "dual float conversion preserves components with factors outside float range"_test = [] {
        using D = canary::test::Dual<float>;
        const canary::UnitScales<double, scaled> reference {
          {.mass = 1.0, .density = 1.0e40, .temperature = 1.0, .length = 1.0}};
        const auto small = convert<scaled>(canary::Density<D, si> {D(1.0F, -2.0F)}, reference);
        expect(small.value().real() == static_cast<float>(1.0 / 1.0e40));
        expect(small.value().dual() == static_cast<float>(-2.0 / 1.0e40));
        const auto large =
          convert<si>(canary::Density<D, scaled> {D(1.0e-30F, -2.0e-30F)}, reference);
        expect(large.value().real() == static_cast<float>(double(1.0e-30F) * 1.0e40));
        expect(large.value().dual() == static_cast<float>(double(-2.0e-30F) * 1.0e40));
    };

    "unit conversion retains long double precision and range"_test = [] {
        const long double adjacent = std::nextafter(1.0L, 2.0L);
        const canary::UnitScales<long double, scaled> reference {
          {.mass = 1.0L, .density = 1.0L, .temperature = adjacent, .length = 1.0L}};
        const auto scaled_temperature =
          convert<scaled>(canary::Temperature<long double, si> {adjacent}, reference);
        expect(scaled_temperature.value() == 1.0L);
        const auto temperature =
          convert<si>(canary::Temperature<long double, scaled> {1.0L}, reference);
        expect(temperature.value() == adjacent);
        const long double largest = std::numeric_limits<long double>::max();
        const canary::UnitScales<double, scaled> unit {
          {.mass = 1.0, .density = 1.0, .temperature = 1.0, .length = 1.0}};
        using D = canary::test::Dual<long double>;
        const auto dual = convert<scaled>(canary::Temperature<D, si> {D(largest, -largest)}, unit);
        expect(dual.value().real() == largest && dual.value().dual() == -largest);
    };

    "a dual-number quantity converts with its derivative"_test = [] {
        using D = canary::test::Dual<double>;
        const canary::Density<D, si> rho {D(2.0e-6, 1.0)};
        const auto rho_s = convert<scaled>(rho, scales);
        expect(within_ulp(rho_s.value().real(), 2.0, 4.0));
        expect(within_ulp(rho_s.value().dual(), 1.0 / scales.density(), 4.0));
    };
}
