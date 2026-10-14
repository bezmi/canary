#include "dual.hpp"
#include "test_suites.hpp"

#include <boost/ut.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <cmath>
#include <concepts>
#include <type_traits>

namespace
{

using canary::scaled;
using canary::si;
using D = canary::test::Dual<double>;

template <typename Q>
concept Addable = requires(Q a, Q b) { a + b; };

template <typename A, typename B>
concept AddableWith = requires(A a, B b) { a + b; };

template <typename A, typename B>
concept MultipliableWith = requires(A a, B b) { a * b; };

template <typename Ratio>
concept HasCaloricallyPerfectRelation =
  requires(Ratio ratio) { canary::pressure(ratio, canary::EnergyDensity<double, si> {1}); } ||
  requires(Ratio ratio) { canary::energy_density(ratio, canary::Pressure<double, si> {1}); } ||
  requires(Ratio ratio) {
      canary::isochoric_specific_heat(canary::SpecificGasConstant<double, si> {1}, ratio);
  } || requires(Ratio ratio) {
      canary::isobaric_specific_heat(canary::SpecificHeat<double, si> {1}, ratio);
  };

}  // namespace

static_assert(sizeof(canary::Density<double, si>) == sizeof(double));
static_assert(sizeof(canary::Density<float, scaled>) == sizeof(float));
static_assert(sizeof(canary::Density<D, si>) == sizeof(D));
static_assert(std::is_trivially_copyable_v<canary::Density<double, si>>);
static_assert(std::is_trivially_copyable_v<canary::Density<D, si>>);
static_assert(std::is_standard_layout_v<canary::Pressure<double, scaled>>);

static_assert(sizeof(canary::AcousticExponent<double>) == sizeof(double));
static_assert(sizeof(canary::AcousticExponent<D>) == sizeof(D));
static_assert(std::is_trivially_copyable_v<canary::AcousticExponent<float>>);
static_assert(canary::KindOf<canary::AcousticExponent<double>, canary::kind::acoustic_exponent>);
static_assert(
  !std::constructible_from<canary::AcousticExponent<double>, canary::HeatCapacityRatio<double>>);
static_assert(
  !std::constructible_from<canary::HeatCapacityRatio<double>, canary::AcousticExponent<double>>);
static_assert(!AddableWith<canary::AcousticExponent<double>, canary::HeatCapacityRatio<double>>);
static_assert(HasCaloricallyPerfectRelation<canary::HeatCapacityRatio<double>>);
static_assert(!HasCaloricallyPerfectRelation<canary::AcousticExponent<double>>);

// Construction is explicit only.
static_assert(!std::convertible_to<double, canary::Density<double, si>>);
static_assert(std::constructible_from<canary::Density<double, si>, double>);

// Value-type conversion within one kind and unit system is explicit only, and
// does not drop the derivative of a dual number.
static_assert(std::constructible_from<canary::Pressure<float, si>, canary::Pressure<double, si>>);
static_assert(!std::convertible_to<canary::Pressure<double, si>, canary::Pressure<float, si>>);
static_assert(
  !std::constructible_from<canary::Pressure<float, si>, canary::Temperature<double, si>>);
static_assert(
  !std::constructible_from<canary::Pressure<double, si>, canary::Pressure<double, scaled>>);
static_assert(std::constructible_from<canary::Density<D, si>, canary::Density<double, si>>);
static_assert(!std::constructible_from<canary::Density<double, si>, canary::Density<D, si>>);

static_assert(canary::KindOf<canary::Temperature<double, si>, canary::kind::temperature>);
static_assert(!canary::KindOf<canary::ModeTemperature<double, si>, canary::kind::temperature>);
static_assert(canary::QuantityKind<const canary::Pressure<float, scaled>&>);
static_assert(!canary::QuantityKind<double>);
static_assert(!canary::Scalar<canary::Density<double, si>>);

// Dimensionless kinds have no unit system; molar kinds exist only in SI.
static_assert(canary::KindUnits<canary::kind::mass_fraction, canary::unitless>);
static_assert(!canary::KindUnits<canary::kind::mass_fraction, si>);
static_assert(!canary::KindUnits<canary::kind::density, canary::unitless>);
static_assert(canary::KindUnits<canary::kind::molar_mass, si>);
static_assert(!canary::KindUnits<canary::kind::molar_mass, scaled>);
static_assert(!canary::KindUnits<canary::kind::universal_gas_constant, scaled>);

// Same-kind arithmetic only, in one unit system.
static_assert(Addable<canary::Pressure<double, si>>);
static_assert(!AddableWith<canary::Pressure<double, si>, canary::Temperature<double, si>>);
static_assert(!AddableWith<canary::Pressure<double, si>, canary::EnergyDensity<double, si>>);
static_assert(!AddableWith<canary::Density<double, si>, canary::Density<double, scaled>>);
static_assert(!AddableWith<canary::Density<double, si>, double>);
static_assert(AddableWith<canary::Density<D, si>, canary::Density<double, si>>);
static_assert(!MultipliableWith<canary::Density<double, si>, canary::Temperature<double, si>>);
static_assert(MultipliableWith<canary::Density<double, si>, double>);

void test_quantity()
{
    using namespace boost::ut;
    using canary::Density;
    using canary::EnergyDensity;
    using canary::HeatCapacityRatio;
    using canary::KindOf;
    using canary::Pressure;
    using canary::SpecificEnergy;
    using canary::SpecificGasConstant;
    using canary::SpecificHeat;
    using canary::Temperature;
    namespace kind = canary::kind;

    "unitless kinds deduce the value type from the constructor argument"_test = [] {
        const canary::MassFraction y {0.25};
        const canary::HeatCapacityRatio gamma {1.4F};
        static_assert(std::same_as<decltype(y), const canary::MassFraction<double>>);
        static_assert(std::same_as<decltype(gamma), const canary::HeatCapacityRatio<float>>);
        expect(y.value() == 0.25_d);
    };

    "explicit conversion between value types rounds to the target type"_test = [] {
        const canary::Pressure<float, si> p {canary::Pressure<double, si> {0.1}};
        expect(p.value() == 0.1_f);
        const canary::Density<D, si> rho {canary::Density<double, si> {2.0}};
        expect(rho.value().real() == 2.0_d && rho.value().dual() == 0.0_d);
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
        const KindOf<kind::pressure> auto p = canary::pressure(rho, gas_constant, temp);
        static_assert(std::same_as<decltype(p), const Pressure<double, si>>);
        expect(p.value() == 1.2 * 287.0 * 300.0);
    };

    "temperature from pressure, density and gas constant inverts the ideal-gas relation"_test = [] {
        const Pressure<double, si> p {101325.0};
        const Density<double, si> rho {1.2};
        const SpecificGasConstant<double, si> gas_constant {287.0};
        const KindOf<kind::temperature> auto temp = canary::temperature(p, rho, gas_constant);
        expect(temp.value() == 101325.0 / (1.2 * 287.0));
        const auto p_again = canary::pressure(rho, gas_constant, temp);
        expect(std::abs(p_again.value() - p.value()) <= 2.0 * 0x1p-52 * p.value());
    };

    "pressure from heat capacity ratio and energy density"_test = [] {
        const KindOf<kind::pressure> auto p =
          canary::pressure(HeatCapacityRatio {1.4}, EnergyDensity<double, si> {2.5e5});
        expect(p.value() == (1.4 - 1.0) * 2.5e5);
    };

    "temperature from specific energy and heat capacity"_test = [] {
        const KindOf<kind::temperature> auto temp =
          canary::temperature(SpecificEnergy<double, si> {2.15e5},
                              SpecificHeat<double, si> {717.5});
        expect(temp.value() == 2.15e5 / 717.5);
    };

    "sound speed from heat capacity ratio, pressure and density"_test = [] {
        const KindOf<kind::sound_speed> auto a = canary::sound_speed(HeatCapacityRatio {1.4},
                                                                     Pressure<double, si> {1.0e5},
                                                                     Density<double, si> {1.2});
        expect(a.value() == std::sqrt(1.4 * 1.0e5 / 1.2));
    };

    "specific energy from energy density and density"_test = [] {
        const KindOf<kind::specific_energy> auto e =
          canary::specific_energy(EnergyDensity<double, scaled> {3.0},
                                  Density<double, scaled> {2.0});
        static_assert(std::same_as<decltype(e), const SpecificEnergy<double, scaled>>);
        expect(e.value() == 1.5_d);
    };

    "density from pressure, gas constant and temperature"_test = [] {
        const KindOf<kind::density> auto rho =
          canary::density(Pressure<double, si> {101325.0},
                          SpecificGasConstant<double, si> {287.0},
                          Temperature<double, si> {300.0});
        expect(rho.value() == 101325.0 / (287.0 * 300.0));
    };

    "specific energy from heat capacity and temperature inverts temperature(e, cv)"_test = [] {
        const SpecificHeat<double, si> cv {717.5};
        const KindOf<kind::specific_energy> auto e =
          canary::specific_energy(cv, Temperature<double, si> {300.0});
        expect(e.value() == 717.5 * 300.0);
        expect(canary::temperature(e, cv).value() == 300.0_d);
    };

    "energy density from density and specific energy, and from pressure and gamma"_test = [] {
        const KindOf<kind::energy_density> auto rho_e =
          canary::energy_density(Density<double, si> {1.2}, SpecificEnergy<double, si> {2.15e5});
        expect(rho_e.value() == 1.2 * 2.15e5);
        const KindOf<kind::energy_density> auto from_p =
          canary::energy_density(HeatCapacityRatio {1.4}, Pressure<double, si> {1.0e5});
        expect(from_p.value() == 1.0e5 / (1.4 - 1.0));
        expect(canary::pressure(
                 HeatCapacityRatio {1.5},
                 canary::energy_density(HeatCapacityRatio {1.5}, Pressure<double, si> {3.0}))
                 .value() == 3.0_d);
    };

    "gas constant and heat capacities from particle mass, Boltzmann constant and gamma"_test = [] {
        const KindOf<kind::specific_gas_constant> auto r =
          canary::specific_gas_constant(canary::BoltzmannConstant<double, si> {1.380649e-23},
                                        canary::ParticleMass<double, si> {4.65e-26});
        expect(r.value() == 1.380649e-23 / 4.65e-26);
        const KindOf<kind::specific_heat> auto cv =
          canary::isochoric_specific_heat(r, HeatCapacityRatio {1.4});
        const KindOf<kind::specific_heat> auto cp =
          canary::isobaric_specific_heat(cv, HeatCapacityRatio {1.4});
        expect(cv.value() == r.value() / (1.4 - 1.0));
        expect(cp.value() == 1.4 * cv.value());
        // cp - cv = R and cp / cv = gamma within rounding.
        expect(std::abs((cp - cv).value() - r.value()) <= 4.0 * 0x1p-52 * r.value());
        expect(std::abs((cp.value() / cv.value()) - 1.4) <= 2.0 * 0x1p-52 * 1.4);
    };

    "heat capacity ratio from cp and cv"_test = [] {
        const KindOf<kind::heat_capacity_ratio> auto gamma =
          canary::heat_capacity_ratio(SpecificHeat<double, si> {1004.5},
                                      SpecificHeat<double, si> {717.5});
        expect(gamma.value() == 1004.5 / 717.5);
    };

    "conductivity from viscosity, heat capacity and Prandtl number"_test = [] {
        const KindOf<kind::thermal_conductivity> auto k =
          canary::conductivity(canary::DynamicViscosity<double, si> {1.8e-5},
                               SpecificHeat<double, si> {1004.5},
                               canary::PrandtlNumber {0.71});
        expect(k.value() == 1.8e-5 * 1004.5 / 0.71);
    };

    "relations accept a parameter value type that differs from the state value type"_test = [] {
        const Density<float, si> rho {1.2F};
        const SpecificGasConstant<double, si> gas_constant {287.0};
        const Temperature<float, si> temp {300.0F};
        const auto p = canary::pressure(rho, gas_constant, temp);
        static_assert(std::same_as<decltype(p), const Pressure<double, si>>);
        expect(p.value() == static_cast<double>(1.2F) * 287.0 * static_cast<double>(300.0F));
    };

    "dual-number quantities carry the derivative through arithmetic and relations"_test = [] {
        // Derivatives with respect to rho at fixed rho_e = 3, gamma = 1.4.
        const Density<D, si> rho {D(2.0, 1.0)};
        const EnergyDensity<D, si> rho_e {D(3.0)};
        const auto e = canary::specific_energy(rho_e, rho);  // de/drho = -rho_e/rho^2
        expect(e.value().real() == 1.5_d);
        expect(e.value().dual() == -0.75_d);
        const auto p = canary::pressure(HeatCapacityRatio {1.4}, rho_e);
        const auto a =
          canary::sound_speed(HeatCapacityRatio {1.4}, p, rho);  // da/drho = -a/(2 rho)
        expect(std::abs(a.value().dual() + (a.value().real() / 4.0)) <=
               4.0 * 0x1p-52 * a.value().real());
        const auto rho_sum = rho + Density<double, si> {1.0};
        expect(rho_sum.value().real() == 3.0_d && rho_sum.value().dual() == 1.0_d);
        const auto rho_scaled = rho * 3.0;
        expect(rho_scaled.value().dual() == 3.0_d);
    };
}
