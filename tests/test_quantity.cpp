#include "dual.hpp"
#include "test_suites.hpp"

#include <boost/ut.hpp>
#include <cmath>
#include <concepts>
#include <gaslib/quantity.hpp>
#include <gaslib/scalar.hpp>
#include <type_traits>

namespace
{

using gaslib::scaled;
using gaslib::si;
using D = gaslib::test::Dual<double>;

template <typename Q>
concept Addable = requires(Q a, Q b) { a + b; };

template <typename A, typename B>
concept AddableWith = requires(A a, B b) { a + b; };

template <typename A, typename B>
concept MultipliableWith = requires(A a, B b) { a * b; };

}  // namespace

static_assert(sizeof(gaslib::Density<double, si>) == sizeof(double));
static_assert(sizeof(gaslib::Density<float, scaled>) == sizeof(float));
static_assert(sizeof(gaslib::Density<D, si>) == sizeof(D));
static_assert(std::is_trivially_copyable_v<gaslib::Density<double, si>>);
static_assert(std::is_trivially_copyable_v<gaslib::Density<D, si>>);
static_assert(std::is_standard_layout_v<gaslib::Pressure<double, scaled>>);

// Construction is explicit only.
static_assert(!std::convertible_to<double, gaslib::Density<double, si>>);
static_assert(std::constructible_from<gaslib::Density<double, si>, double>);

static_assert(gaslib::KindOf<gaslib::Temperature<double, si>, gaslib::kind::temperature>);
static_assert(!gaslib::KindOf<gaslib::ModeTemperature<double, si>, gaslib::kind::temperature>);
static_assert(gaslib::QuantityKind<const gaslib::Pressure<float, scaled>&>);
static_assert(!gaslib::QuantityKind<double>);
static_assert(!gaslib::Scalar<gaslib::Density<double, si>>);

// Dimensionless kinds have no unit system; molar kinds exist only in SI.
static_assert(gaslib::KindUnits<gaslib::kind::mass_fraction, gaslib::unitless>);
static_assert(!gaslib::KindUnits<gaslib::kind::mass_fraction, si>);
static_assert(!gaslib::KindUnits<gaslib::kind::density, gaslib::unitless>);
static_assert(gaslib::KindUnits<gaslib::kind::molar_mass, si>);
static_assert(!gaslib::KindUnits<gaslib::kind::molar_mass, scaled>);
static_assert(!gaslib::KindUnits<gaslib::kind::universal_gas_constant, scaled>);

// Same-kind arithmetic only, in one unit system.
static_assert(Addable<gaslib::Pressure<double, si>>);
static_assert(!AddableWith<gaslib::Pressure<double, si>, gaslib::Temperature<double, si>>);
static_assert(!AddableWith<gaslib::Pressure<double, si>, gaslib::EnergyDensity<double, si>>);
static_assert(!AddableWith<gaslib::Density<double, si>, gaslib::Density<double, scaled>>);
static_assert(!AddableWith<gaslib::Density<double, si>, double>);
static_assert(AddableWith<gaslib::Density<D, si>, gaslib::Density<double, si>>);
static_assert(!MultipliableWith<gaslib::Density<double, si>, gaslib::Temperature<double, si>>);
static_assert(MultipliableWith<gaslib::Density<double, si>, double>);

void test_quantity()
{
    using namespace boost::ut;
    using gaslib::Density;
    using gaslib::EnergyDensity;
    using gaslib::HeatCapacityRatio;
    using gaslib::KindOf;
    using gaslib::Pressure;
    using gaslib::SpecificEnergy;
    using gaslib::SpecificGasConstant;
    using gaslib::SpecificHeat;
    using gaslib::Temperature;
    namespace kind = gaslib::kind;

    "unitless kinds deduce the value type from the constructor argument"_test = [] {
        const gaslib::MassFraction y {0.25};
        const gaslib::HeatCapacityRatio gamma {1.4F};
        static_assert(std::same_as<decltype(y), const gaslib::MassFraction<double>>);
        static_assert(std::same_as<decltype(gamma), const gaslib::HeatCapacityRatio<float>>);
        expect(y.value() == 0.25_d);
    };

    "same-kind addition, subtraction and comparison"_test = [] {
        const Pressure<double, si> a {3.0};
        const Pressure<double, si> b {1.0};
        expect((a + b).value() == 4.0_d);
        expect((a - b).value() == 2.0_d);
        expect(b < a);
        expect(!(a < b));
        expect(a == Pressure<double, si> {3.0});
        expect(a != b);
    };

    "scaling by a plain value keeps the kind"_test = [] {
        const Density<double, si> rho {2.0};
        const auto twice = rho * 2.0;
        const auto also_twice = 2.0F * rho;
        const auto half = rho / 4.0;
        static_assert(std::same_as<decltype(twice), const Density<double, si>>);
        static_assert(std::same_as<decltype(also_twice), const Density<double, si>>);
        expect(twice.value() == 4.0_d && also_twice.value() == 4.0_d && half.value() == 0.5_d);
    };

    "mixed value types give the type of the raw operation"_test = [] {
        const Density<float, si> rho {1.5F};
        const Density<double, si> rho_ref {0.5};
        static_assert(std::same_as<decltype(rho + rho_ref), Density<double, si>>);
        static_assert(std::same_as<decltype(Density<D, si> {D(1.0)} + rho_ref), Density<D, si>>);
        expect((rho + rho_ref).value() == 2.0_d);
    };

    "pressure from density, gas constant and temperature"_test = [] {
        const Density<double, si> rho {1.2};
        const SpecificGasConstant<double, si> gas_constant {287.0};
        const Temperature<double, si> temp {300.0};
        const KindOf<kind::pressure> auto p = gaslib::pressure(rho, gas_constant, temp);
        static_assert(std::same_as<decltype(p), const Pressure<double, si>>);
        expect(p.value() == 1.2 * 287.0 * 300.0);
    };

    "temperature from pressure, density and gas constant inverts the ideal-gas relation"_test = [] {
        const Pressure<double, si> p {101325.0};
        const Density<double, si> rho {1.2};
        const SpecificGasConstant<double, si> gas_constant {287.0};
        const KindOf<kind::temperature> auto temp = gaslib::temperature(p, rho, gas_constant);
        expect(temp.value() == 101325.0 / (1.2 * 287.0));
        const auto p_again = gaslib::pressure(rho, gas_constant, temp);
        expect(std::abs(p_again.value() - p.value()) <= 2.0 * 0x1p-52 * p.value());
    };

    "pressure from heat capacity ratio and energy density"_test = [] {
        const KindOf<kind::pressure> auto p =
          gaslib::pressure(HeatCapacityRatio {1.4}, EnergyDensity<double, si> {2.5e5});
        expect(p.value() == (1.4 - 1.0) * 2.5e5);
    };

    "temperature from specific energy and heat capacity"_test = [] {
        const KindOf<kind::temperature> auto temp =
          gaslib::temperature(SpecificEnergy<double, si> {2.15e5},
                              SpecificHeat<double, si> {717.5});
        expect(temp.value() == 2.15e5 / 717.5);
    };

    "sound speed from heat capacity ratio, pressure and density"_test = [] {
        const KindOf<kind::sound_speed> auto a = gaslib::sound_speed(HeatCapacityRatio {1.4},
                                                                     Pressure<double, si> {1.0e5},
                                                                     Density<double, si> {1.2});
        expect(a.value() == std::sqrt(1.4 * 1.0e5 / 1.2));
    };

    "specific energy from energy density and density"_test = [] {
        const KindOf<kind::specific_energy> auto e =
          gaslib::specific_energy(EnergyDensity<double, scaled> {3.0},
                                  Density<double, scaled> {2.0});
        static_assert(std::same_as<decltype(e), const SpecificEnergy<double, scaled>>);
        expect(e.value() == 1.5_d);
    };

    "density from pressure, gas constant and temperature"_test = [] {
        const KindOf<kind::density> auto rho =
          gaslib::density(Pressure<double, si> {101325.0},
                          SpecificGasConstant<double, si> {287.0},
                          Temperature<double, si> {300.0});
        expect(rho.value() == 101325.0 / (287.0 * 300.0));
    };

    "specific energy from heat capacity and temperature inverts temperature(e, cv)"_test = [] {
        const SpecificHeat<double, si> cv {717.5};
        const KindOf<kind::specific_energy> auto e =
          gaslib::specific_energy(cv, Temperature<double, si> {300.0});
        expect(e.value() == 717.5 * 300.0);
        expect(gaslib::temperature(e, cv).value() == 300.0_d);
    };

    "conductivity from viscosity, heat capacity and Prandtl number"_test = [] {
        const KindOf<kind::thermal_conductivity> auto k =
          gaslib::conductivity(gaslib::DynamicViscosity<double, si> {1.8e-5},
                               SpecificHeat<double, si> {1004.5},
                               gaslib::PrandtlNumber {0.71});
        expect(k.value() == 1.8e-5 * 1004.5 / 0.71);
    };

    "relations accept a parameter value type that differs from the state value type"_test = [] {
        const Density<float, si> rho {1.2F};
        const SpecificGasConstant<double, si> gas_constant {287.0};
        const Temperature<float, si> temp {300.0F};
        const auto p = gaslib::pressure(rho, gas_constant, temp);
        static_assert(std::same_as<decltype(p), const Pressure<double, si>>);
        expect(p.value() == static_cast<double>(1.2F) * 287.0 * static_cast<double>(300.0F));
    };

    "dual-number quantities carry the derivative through arithmetic and relations"_test = [] {
        // Derivatives with respect to rho at fixed rho_e = 3, gamma = 1.4.
        const Density<D, si> rho {D(2.0, 1.0)};
        const EnergyDensity<D, si> rho_e {D(3.0)};
        const auto e = gaslib::specific_energy(rho_e, rho);  // de/drho = -rho_e/rho^2
        expect(e.value().real() == 1.5_d);
        expect(e.value().dual() == -0.75_d);
        const auto p = gaslib::pressure(HeatCapacityRatio {1.4}, rho_e);
        const auto a =
          gaslib::sound_speed(HeatCapacityRatio {1.4}, p, rho);  // da/drho = -a/(2 rho)
        expect(std::abs(a.value().dual() + (a.value().real() / 4.0)) <=
               4.0 * 0x1p-52 * a.value().real());
        const auto rho_sum = rho + Density<double, si> {1.0};
        expect(rho_sum.value().real() == 3.0_d && rho_sum.value().dual() == 1.0_d);
        const auto rho_scaled = rho * 3.0;
        expect(rho_scaled.value().dual() == 3.0_d);
    };
}
