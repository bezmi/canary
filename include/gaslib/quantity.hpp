#ifndef GASLIB_QUANTITY_HPP
#define GASLIB_QUANTITY_HPP

#include <concepts>
#include <gaslib/config.hpp>
#include <gaslib/scalar.hpp>
#include <type_traits>

namespace gaslib
{

// Unit systems. Values of a dimensional quantity are in SI units (si) or
// in cerberus's reference-scaled units (scaled). Dimensionless kinds carry
// the tag unitless, since their values are the same in every system.
struct si {};
struct scaled {};
struct unitless {};

template <typename U>
concept UnitSystem = std::same_as<U, si> || std::same_as<U, scaled>;

template <typename U>
concept UnitTag = UnitSystem<U> || std::same_as<U, unitless>;

// Exponents of the base dimensions. count is the number of particles, so
// a per-particle quantity has count -1; it lets the scaled system use a
// density scale independent of the particle mass and length scales.
// amount is the amount of substance (mole) and is used only in SI.
struct Dimension {
    int mass = 0;
    int length = 0;
    int time = 0;
    int temperature = 0;
    int charge = 0;
    int count = 0;
    int amount = 0;

    friend constexpr bool operator==(const Dimension&, const Dimension&) = default;
};

inline constexpr Dimension dimensionless {};

// A kind fixes the physical meaning of a value and its dimension, not its
// unit. Each kind is an empty tag with a static member dimension.
template <typename K>
concept KindTag = std::is_empty_v<K> && requires {
    { K::dimension } -> std::convertible_to<Dimension>;
};

namespace kind
{

// clang-format off
struct density                { static constexpr Dimension dimension {.mass = 1, .length = -3}; };
struct pressure               { static constexpr Dimension dimension {.mass = 1, .length = -1, .time = -2}; };
struct temperature            { static constexpr Dimension dimension {.temperature = 1}; };
// Internal energy per unit mass.
struct specific_energy        { static constexpr Dimension dimension {.length = 2, .time = -2}; };
// Internal energy per unit volume.
struct energy_density         { static constexpr Dimension dimension {.mass = 1, .length = -1, .time = -2}; };
// cv or cp, per unit mass.
struct specific_heat          { static constexpr Dimension dimension {.length = 2, .time = -2, .temperature = -1}; };
struct specific_gas_constant  { static constexpr Dimension dimension {.length = 2, .time = -2, .temperature = -1}; };
struct sound_speed            { static constexpr Dimension dimension {.length = 1, .time = -1}; };
struct mass_fraction          { static constexpr Dimension dimension {}; };
struct particle_mass          { static constexpr Dimension dimension {.mass = 1, .count = -1}; };
struct molar_mass             { static constexpr Dimension dimension {.mass = 1, .amount = -1}; };
// Charge of one particle.
struct charge                 { static constexpr Dimension dimension {.charge = 1, .count = -1}; };
struct heat_capacity_ratio    { static constexpr Dimension dimension {}; };
struct prandtl_number         { static constexpr Dimension dimension {}; };
struct boltzmann_constant     { static constexpr Dimension dimension {.mass = 1, .length = 2, .time = -2, .temperature = -1, .count = -1}; };
struct universal_gas_constant { static constexpr Dimension dimension {.mass = 1, .length = 2, .time = -2, .temperature = -1, .amount = -1}; };
struct dynamic_viscosity      { static constexpr Dimension dimension {.mass = 1, .length = -1, .time = -1}; };
struct thermal_conductivity   { static constexpr Dimension dimension {.mass = 1, .length = 1, .time = -3, .temperature = -1}; };
// Energy per unit volume and temperature of one internal energy mode.
struct mode_energy_density    { static constexpr Dimension dimension {.mass = 1, .length = -1, .time = -2}; };
struct mode_temperature       { static constexpr Dimension dimension {.temperature = 1}; };
// clang-format on

}  // namespace kind

// A dimensionless kind has the tag unitless and a dimensional kind has a unit
// system. Kinds with an amount of substance exist only in SI.
template <typename K, typename U>
concept KindUnits =
  KindTag<K> && UnitTag<U> && ((K::dimension == dimensionless) == std::same_as<U, unitless>) &&
  (K::dimension.amount == 0 || std::same_as<U, si>);

// A value of kind K, value type T and unit system U. It has the size of T,
// is trivially copyable when T is, and is constructed only explicitly.
template <KindTag K, Scalar T, UnitTag U>
    requires KindUnits<K, U>
class Quantity
{
  public:
    using kind = K;
    using value_type = T;
    using units = U;

    constexpr Quantity() = default;
    GASLIB_HOST_DEVICE constexpr explicit Quantity(T value) : value_(value) {}

    // Conversion from another value type of the same kind and unit system,
    // for example a double result stored in a float state, or a constant
    // stored as a dual number. Explicit, so every precision change is
    // visible.
    template <Scalar V>
        requires(!std::same_as<V, T>) && std::constructible_from<T, V>
    GASLIB_HOST_DEVICE constexpr explicit Quantity(const Quantity<K, V, U>& other) :
        value_(T(other.value()))
    {
    }

    [[nodiscard]] GASLIB_HOST_DEVICE constexpr T value() const { return value_; }

  private:
    T value_ {};
};

namespace detail
{
template <typename Q> struct is_quantity : std::false_type {};
template <KindTag K, Scalar T, UnitTag U> struct is_quantity<Quantity<K, T, U>> : std::true_type {};
}  // namespace detail

namespace detail
{
// Wraps the result of a relation formula; the value type is that of the
// formula.
template <KindTag K, UnitTag U, Scalar T>
GASLIB_HOST_DEVICE constexpr Quantity<K, T, U> make_quantity(T value)
{ return Quantity<K, T, U> {value}; }
}  // namespace detail

template <typename Q>
concept QuantityKind = detail::is_quantity<std::remove_cvref_t<Q>>::value;

template <typename Q, typename K>
concept KindOf = QuantityKind<Q> && std::same_as<typename std::remove_cvref_t<Q>::kind, K>;

// clang-format off
template <Scalar T, UnitSystem U> using Density              = Quantity<kind::density, T, U>;
template <Scalar T, UnitSystem U> using Pressure             = Quantity<kind::pressure, T, U>;
template <Scalar T, UnitSystem U> using Temperature          = Quantity<kind::temperature, T, U>;
template <Scalar T, UnitSystem U> using SpecificEnergy       = Quantity<kind::specific_energy, T, U>;
template <Scalar T, UnitSystem U> using EnergyDensity        = Quantity<kind::energy_density, T, U>;
template <Scalar T, UnitSystem U> using SpecificHeat         = Quantity<kind::specific_heat, T, U>;
template <Scalar T, UnitSystem U> using SpecificGasConstant  = Quantity<kind::specific_gas_constant, T, U>;
template <Scalar T, UnitSystem U> using SoundSpeed           = Quantity<kind::sound_speed, T, U>;
template <Scalar T>               using MassFraction         = Quantity<kind::mass_fraction, T, unitless>;
template <Scalar T, UnitSystem U> using ParticleMass         = Quantity<kind::particle_mass, T, U>;
template <Scalar T>               using MolarMass            = Quantity<kind::molar_mass, T, si>;
template <Scalar T, UnitSystem U> using Charge               = Quantity<kind::charge, T, U>;
template <Scalar T>               using HeatCapacityRatio    = Quantity<kind::heat_capacity_ratio, T, unitless>;
template <Scalar T>               using PrandtlNumber        = Quantity<kind::prandtl_number, T, unitless>;
template <Scalar T, UnitSystem U> using BoltzmannConstant    = Quantity<kind::boltzmann_constant, T, U>;
template <Scalar T>               using UniversalGasConstant = Quantity<kind::universal_gas_constant, T, si>;
template <Scalar T, UnitSystem U> using DynamicViscosity     = Quantity<kind::dynamic_viscosity, T, U>;
template <Scalar T, UnitSystem U> using ThermalConductivity  = Quantity<kind::thermal_conductivity, T, U>;
template <Scalar T, UnitSystem U> using ModeEnergyDensity    = Quantity<kind::mode_energy_density, T, U>;
template <Scalar T, UnitSystem U> using ModeTemperature      = Quantity<kind::mode_temperature, T, U>;
// clang-format on

// Same-kind arithmetic and comparison. The operands may have different value
// types (a state value and a model parameter); the result has the type of
// the raw operation.
template <KindTag K, Scalar A, Scalar B, UnitTag U>
    requires ScalarPair<A, B>
GASLIB_HOST_DEVICE constexpr auto operator+(Quantity<K, A, U> a, Quantity<K, B, U> b)
{ return Quantity<K, decltype(a.value() + b.value()), U> {a.value() + b.value()}; }

template <KindTag K, Scalar A, Scalar B, UnitTag U>
    requires ScalarPair<A, B>
GASLIB_HOST_DEVICE constexpr auto operator-(Quantity<K, A, U> a, Quantity<K, B, U> b)
{ return Quantity<K, decltype(a.value() - b.value()), U> {a.value() - b.value()}; }

template <KindTag K, Scalar T, UnitTag U>
GASLIB_HOST_DEVICE constexpr bool operator<(Quantity<K, T, U> a, Quantity<K, T, U> b)
{ return a.value() < b.value(); }

template <KindTag K, Scalar T, UnitTag U>
GASLIB_HOST_DEVICE constexpr bool operator==(Quantity<K, T, U> a, Quantity<K, T, U> b)
{ return a.value() == b.value(); }

// Scaling by a dimensionless plain value keeps the kind.
template <KindTag K, Scalar A, UnitTag U, Scalar S>
    requires ScalarPair<A, S>
GASLIB_HOST_DEVICE constexpr auto operator*(Quantity<K, A, U> q, S s)
{ return Quantity<K, decltype(q.value() * s), U> {q.value() * s}; }

template <KindTag K, Scalar A, UnitTag U, Scalar S>
    requires ScalarPair<A, S>
GASLIB_HOST_DEVICE constexpr auto operator*(S s, Quantity<K, A, U> q)
{ return Quantity<K, decltype(s * q.value()), U> {s * q.value()}; }

template <KindTag K, Scalar A, UnitTag U, Scalar S>
    requires ScalarPair<A, S>
GASLIB_HOST_DEVICE constexpr auto operator/(Quantity<K, A, U> q, S s)
{ return Quantity<K, decltype(q.value() / s), U> {q.value() / s}; }

// Relations between different kinds. Each is one formula with a fixed order
// of operations; there is no other product or quotient of different kinds.
// Argument value types may differ (state values and model parameters).

namespace detail
{
template <Scalar A, Scalar B>
    requires ScalarPair<A, B>
using product_t = decltype(std::declval<A>() * std::declval<B>());

template <Scalar A, Scalar B>
    requires ScalarPair<A, B>
using quotient_t = decltype(std::declval<A>() / std::declval<B>());

template <Scalar A, Scalar B>
    requires ScalarPair<A, B>
using difference_t = decltype(std::declval<A>() - std::declval<B>());
}  // namespace detail

// p = (rho * R) * T
template <Scalar A, Scalar B, Scalar C, UnitSystem U>
    requires ScalarPair<A, B> && ScalarPair<detail::product_t<A, B>, C>
GASLIB_HOST_DEVICE constexpr auto
  pressure(Density<A, U> rho, SpecificGasConstant<B, U> gas_constant, Temperature<C, U> temp)
{
    return detail::make_quantity<kind::pressure, U>(rho.value() * gas_constant.value() *
                                                    temp.value());
}

// T = p / (rho * R)
template <Scalar A, Scalar B, Scalar C, UnitSystem U>
    requires ScalarPair<B, C> && ScalarPair<A, detail::product_t<B, C>>
GASLIB_HOST_DEVICE constexpr auto
  temperature(Pressure<A, U> p, Density<B, U> rho, SpecificGasConstant<C, U> gas_constant)
{
    return detail::make_quantity<kind::temperature, U>(p.value() /
                                                       (rho.value() * gas_constant.value()));
}

// p = (gamma - 1) * rho_e
template <Scalar G, Scalar A, UnitSystem U>
    requires ScalarPair<G, typename ScalarTraits<G>::value_type> &&
             ScalarPair<detail::difference_t<G, typename ScalarTraits<G>::value_type>, A>
GASLIB_HOST_DEVICE constexpr auto pressure(HeatCapacityRatio<G> gamma, EnergyDensity<A, U> rho_e)
{
    using V = ScalarTraits<G>::value_type;
    return detail::make_quantity<kind::pressure, U>((gamma.value() - V(1)) * rho_e.value());
}

// T = e / cv
template <Scalar A, Scalar B, UnitSystem U>
    requires ScalarPair<A, B>
GASLIB_HOST_DEVICE constexpr auto temperature(SpecificEnergy<A, U> e, SpecificHeat<B, U> cv)
{ return detail::make_quantity<kind::temperature, U>(e.value() / cv.value()); }

// a = sqrt((gamma * p) / rho)
template <Scalar G, Scalar A, Scalar B, UnitSystem U>
    requires ScalarPair<G, A> && ScalarPair<detail::product_t<G, A>, B> &&
             HasSqrt<detail::quotient_t<detail::product_t<G, A>, B>>
GASLIB_HOST_DEVICE auto sound_speed(HeatCapacityRatio<G> gamma, Pressure<A, U> p, Density<B, U> rho)
{
    return detail::make_quantity<kind::sound_speed, U>(
      math::sqrt(gamma.value() * p.value() / rho.value()));
}

// e = rho_e / rho
template <Scalar A, Scalar B, UnitSystem U>
    requires ScalarPair<A, B>
GASLIB_HOST_DEVICE constexpr auto specific_energy(EnergyDensity<A, U> rho_e, Density<B, U> rho)
{ return detail::make_quantity<kind::specific_energy, U>(rho_e.value() / rho.value()); }

// rho = p / (R * T)
template <Scalar A, Scalar B, Scalar C, UnitSystem U>
    requires ScalarPair<B, C> && ScalarPair<A, detail::product_t<B, C>>
GASLIB_HOST_DEVICE constexpr auto
  density(Pressure<A, U> p, SpecificGasConstant<B, U> gas_constant, Temperature<C, U> temp)
{
    return detail::make_quantity<kind::density, U>(p.value() /
                                                   (gas_constant.value() * temp.value()));
}

// e = cv * T
template <Scalar A, Scalar B, UnitSystem U>
    requires ScalarPair<A, B>
GASLIB_HOST_DEVICE constexpr auto specific_energy(SpecificHeat<A, U> cv, Temperature<B, U> temp)
{ return detail::make_quantity<kind::specific_energy, U>(cv.value() * temp.value()); }

// rho_e = rho * e
template <Scalar A, Scalar B, UnitSystem U>
    requires ScalarPair<A, B>
GASLIB_HOST_DEVICE constexpr auto energy_density(Density<A, U> rho, SpecificEnergy<B, U> e)
{ return detail::make_quantity<kind::energy_density, U>(rho.value() * e.value()); }

// rho_e = p / (gamma - 1)
template <Scalar G, Scalar A, UnitSystem U>
    requires ScalarPair<G, typename ScalarTraits<G>::value_type> &&
             ScalarPair<A, detail::difference_t<G, typename ScalarTraits<G>::value_type>>
GASLIB_HOST_DEVICE constexpr auto energy_density(HeatCapacityRatio<G> gamma, Pressure<A, U> p)
{
    using V = ScalarTraits<G>::value_type;
    return detail::make_quantity<kind::energy_density, U>(p.value() / (gamma.value() - V(1)));
}

// R = k_B / m
template <Scalar A, Scalar B, UnitSystem U>
    requires ScalarPair<A, B>
GASLIB_HOST_DEVICE constexpr auto specific_gas_constant(BoltzmannConstant<A, U> boltzmann,
                                                        ParticleMass<B, U> mass)
{ return detail::make_quantity<kind::specific_gas_constant, U>(boltzmann.value() / mass.value()); }

// cv = R / (gamma - 1)
template <Scalar A, Scalar G, UnitSystem U>
    requires ScalarPair<G, typename ScalarTraits<G>::value_type> &&
             ScalarPair<A, detail::difference_t<G, typename ScalarTraits<G>::value_type>>
GASLIB_HOST_DEVICE constexpr auto isochoric_specific_heat(SpecificGasConstant<A, U> gas_constant,
                                                          HeatCapacityRatio<G> gamma)
{
    using V = ScalarTraits<G>::value_type;
    return detail::make_quantity<kind::specific_heat, U>(gas_constant.value() /
                                                         (gamma.value() - V(1)));
}

// cp = gamma * cv
template <Scalar A, Scalar G, UnitSystem U>
    requires ScalarPair<G, A>
GASLIB_HOST_DEVICE constexpr auto isobaric_specific_heat(SpecificHeat<A, U> cv,
                                                         HeatCapacityRatio<G> gamma)
{ return detail::make_quantity<kind::specific_heat, U>(gamma.value() * cv.value()); }

// gamma = cp / cv
template <Scalar A, Scalar B, UnitSystem U>
    requires ScalarPair<A, B>
GASLIB_HOST_DEVICE constexpr auto heat_capacity_ratio(SpecificHeat<A, U> cp, SpecificHeat<B, U> cv)
{ return detail::make_quantity<kind::heat_capacity_ratio, unitless>(cp.value() / cv.value()); }

// k = (mu * cp) / Pr
template <Scalar A, Scalar B, Scalar C, UnitSystem U>
    requires ScalarPair<A, B> && ScalarPair<detail::product_t<A, B>, C>
GASLIB_HOST_DEVICE constexpr auto
  conductivity(DynamicViscosity<A, U> mu, SpecificHeat<B, U> cp, PrandtlNumber<C> prandtl)
{
    return detail::make_quantity<kind::thermal_conductivity, U>(mu.value() * cp.value() /
                                                                prandtl.value());
}

}  // namespace gaslib

#endif  // GASLIB_QUANTITY_HPP
