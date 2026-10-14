#ifndef CANARY_UNITS_HPP
#define CANARY_UNITS_HPP

#include <canary/build.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <format>
#include <limits>
#include <type_traits>
#include <utility>

namespace canary
{

namespace constants
{
// Exact SI values (2019 definition of the SI base units).
template <std::floating_point P> inline constexpr P boltzmann = P(1.380649e-23);             // J/K
template <std::floating_point P> inline constexpr P elementary_charge = P(1.602176634e-19);  // C
template <std::floating_point P> inline constexpr P avogadro = P(6.02214076e23);  // 1/mol
// CODATA 2018 recommended value of the electric constant [F/m].
template <std::floating_point P> inline constexpr P vacuum_permittivity = P(8.8541878128e-12);
// Universal gas constant N_A k_B [J/(mol K)], the product of two exact
// values rounded to P.
template <std::floating_point P>
inline constexpr P universal_gas_constant = avogadro<P> * boltzmann<P>;
}  // namespace constants

namespace detail
{
// Positive host-side scale arithmetic with a separate binary exponent. Products
// and quotients do not overflow before the final result is converted to W.
template <std::floating_point W> struct BinaryScale {
    W fraction;
    std::int64_t exponent;

    explicit BinaryScale(W value)
    {
        int e = 0;
        fraction = std::frexp(value, &e);
        exponent = e;
    }

    BinaryScale(W value, std::int64_t e) : BinaryScale(value) { exponent += e; }

    BinaryScale operator*(const BinaryScale& other) const
    { return {fraction * other.fraction, exponent + other.exponent}; }

    BinaryScale operator/(const BinaryScale& other) const
    { return {fraction / other.fraction, exponent - other.exponent}; }

    BinaryScale square_root() const
    {
        const auto even_exponent = exponent - (exponent % 2 != 0 ? 1 : 0);
        const W argument = exponent == even_exponent ? fraction : W(2) * fraction;
        return {std::sqrt(argument), even_exponent / 2};
    }

    W value() const
    {
        if (exponent > std::numeric_limits<int>::max()) {
            return std::numeric_limits<W>::infinity();
        }
        if (exponent < std::numeric_limits<int>::min()) { return W(0); }
        return std::ldexp(fraction, static_cast<int>(exponent));
    }
};

template <std::floating_point W> void check_scale(const char* name, W value)
{
    if (!(value > W(0)) || !std::isfinite(value)) {
        throw ParameterError(
          std::format("{} scale must be finite and positive, got {}", name, value));
    }
}

template <std::floating_point W> W reference_velocity(W mass, W temperature)
{
    const W energy = constants::boltzmann<W> * temperature;
    const W squared = energy / mass;
    if (std::isnormal(energy) && std::isnormal(squared)) { return std::sqrt(squared); }
    return (BinaryScale<W>(constants::boltzmann<W>) * BinaryScale<W>(temperature) /
            BinaryScale<W>(mass))
      .square_root()
      .value();
}

template <std::floating_point W> W reference_count(W mass, W density, W length)
{
    const W number_density = density / mass;
    const W line_count = number_density * length;
    const W area_count = line_count * length;
    const W count = area_count * length;
    if (std::isnormal(number_density) && std::isnormal(line_count) && std::isnormal(area_count) &&
        std::isnormal(count)) {
        return count;
    }
    const BinaryScale<W> x(length);
    return (BinaryScale<W>(density) / BinaryScale<W>(mass) * x * x * x).value();
}
}  // namespace detail

// Reference scales of a unit system. A value in the system is the SI value
// divided by the scale of its kind.
template <std::floating_point P, UnitSystem U> struct UnitScales;

// SI is the identity.
template <std::floating_point P> struct UnitScales<P, si> {};

// Scales of cerberus's nondimensional system, in SI units: the reference
// particle mass m_ref [kg], mass density rho_ref [kg/m^3], temperature
// T_ref [K] and length x_ref [m]. The velocity scale is sqrt(k_B T_ref /
// m_ref), so the Boltzmann constant is 1 in this system. Charges are in
// units of the elementary charge.
// Host construction requires finite, positive base scales and derived
// velocity, time and particle-count scales in common_type_t<double, P>. Input
// is copied; subsequent changes do not affect the validated value.
template <std::floating_point P> struct UnitScales<P, scaled>
{
  public:
    struct Input {
        P mass {};
        P density {};
        P temperature {};
        P length {};
    };

    // Throws ParameterError naming an invalid base or derived scale. An
    // omitted input field is zero and is rejected.
    explicit UnitScales(Input input) : values_(input)
    {
        detail::check_scale("mass", input.mass);
        detail::check_scale("density", input.density);
        detail::check_scale("temperature", input.temperature);
        detail::check_scale("length", input.length);
        using W = std::common_type_t<double, P>;
        const W velocity = detail::reference_velocity(W(input.mass), W(input.temperature));
        detail::check_scale("velocity", velocity);
        detail::check_scale("time", W(input.length) / velocity);
        detail::check_scale(
          "particle count",
          detail::reference_count(W(input.mass), W(input.density), W(input.length)));
    }

    [[nodiscard]] constexpr P mass() const { return values_.mass; }
    [[nodiscard]] constexpr P density() const { return values_.density; }
    [[nodiscard]] constexpr P temperature() const { return values_.temperature; }
    [[nodiscard]] constexpr P length() const { return values_.length; }

  private:
    Input values_;
};

// Host-side derived scales use at least double precision and retain long
// double precision. The base scales remain stored in P.

// u_ref = sqrt(k_B T_ref / m_ref) [m/s]
template <std::floating_point P>
std::common_type_t<double, P> velocity_scale(const UnitScales<P, scaled>& s)
{
    using W = std::common_type_t<double, P>;
    return detail::reference_velocity(W(s.mass()), W(s.temperature()));
}

// t_ref = x_ref / u_ref [s]
template <std::floating_point P>
std::common_type_t<double, P> time_scale(const UnitScales<P, scaled>& s)
{ return s.length() / velocity_scale(s); }

// Number of particles in a reference volume, n0 = (rho_ref / m_ref) x_ref^3.
// It is the scale of the particle-count dimension, so that the density
// scale is rho_ref and the particle-mass scale is m_ref.
template <std::floating_point P>
std::common_type_t<double, P> count_scale(const UnitScales<P, scaled>& s)
{
    using W = std::common_type_t<double, P>;
    return detail::reference_count(W(s.mass()), W(s.density()), W(s.length()));
}

namespace detail
{
// Multiplies numerator by base^exponent when exponent > 0, or denominator
// by base^-exponent when exponent < 0. The factor is numerator/denominator,
// so it has a single division.
template <std::floating_point P>
constexpr void apply_power(P& numerator, P& denominator, P base, int exponent)
{
    for (int i = 0; i < exponent; ++i) { numerator *= base; }
    for (int i = 0; i > exponent; --i) { denominator *= base; }
}
}  // namespace detail

// SI value of one unit of dimension d in the system of the scales: the
// product of the base scales raised to the exponents of d. Throws
// ParameterError when d has an amount of substance, which has no scaled
// unit.
template <std::floating_point P>
std::common_type_t<double, P> scale_factor(const Dimension& d, const UnitScales<P, scaled>& s)
{
    using W = std::common_type_t<double, P>;
    if (d.amount != 0) {
        throw ParameterError("a dimension with an amount of substance has no scaled unit");
    }
    const W count = count_scale(s);
    W numerator = W(1);
    W denominator = W(1);
    detail::apply_power(numerator, denominator, s.mass() * count, d.mass);
    detail::apply_power(numerator, denominator, W(s.length()), d.length);
    detail::apply_power(numerator, denominator, time_scale(s), d.time);
    detail::apply_power(numerator, denominator, W(s.temperature()), d.temperature);
    detail::apply_power(numerator, denominator, constants::elementary_charge<W> * count, d.charge);
    detail::apply_power(numerator, denominator, count, d.count);
    return numerator / denominator;
}

template <std::floating_point P>
constexpr std::common_type_t<double, P> scale_factor(const Dimension& /*d*/,
                                                     const UnitScales<P, si>& /*scales*/)
{ return std::common_type_t<double, P>(1); }

// SI value of one unit of kind K in the system of the scales.
template <KindTag K, std::floating_point P>
    requires(K::dimension.amount == 0)
std::common_type_t<double, P> scale_factor(const UnitScales<P, scaled>& s)
{ return scale_factor(K::dimension, s); }

template <KindTag K, std::floating_point P>
constexpr std::common_type_t<double, P> scale_factor(const UnitScales<P, si>& /*scales*/)
{ return std::common_type_t<double, P>(1); }

// A conversion of kind K from unit system From to To with scales of system
// S: between SI and S, within one system, or of a dimensionless kind.
template <typename K, typename From, typename To, typename S>
concept ConvertibleUnits =
  KindTag<K> && UnitTag<From> && UnitSystem<To> && UnitSystem<S> &&
  (std::same_as<From, unitless> ||
   (KindUnits<K, To> &&
    (std::same_as<From, To> || (std::same_as<From, si> && std::same_as<To, S>) ||
     (std::same_as<From, S> && std::same_as<To, si>))));

namespace detail
{
template <std::floating_point W, Scalar T>
    requires(
      std::same_as<W, typename ScalarTraits<T>::value_type> ||
      requires(const T& x) {
          { ScalarTraits<T>::template cast<W>(x) } -> Scalar;
          requires std::same_as<
            typename ScalarTraits<decltype(ScalarTraits<T>::template cast<W>(x))>::value_type,
            W>;
      })
constexpr auto scalar_cast(const T& x)
{
    if constexpr (std::same_as<W, typename ScalarTraits<T>::value_type>) {
        return x;
    } else {
        return ScalarTraits<T>::template cast<W>(x);
    }
}

template <Scalar T, std::floating_point P>
using conversion_precision_t = std::common_type_t<double, P, typename ScalarTraits<T>::value_type>;

template <Scalar T, std::floating_point P>
using conversion_scalar_t = decltype(scalar_cast<conversion_precision_t<T, P>>(std::declval<T>()));

template <typename T, typename P>
concept ConvertibleScalar = Scalar<T> && std::floating_point<P> && requires(const T& x) {
    { scalar_cast<conversion_precision_t<T, P>>(x) } -> Scalar;
    requires ScalarWith<conversion_scalar_t<T, P>, conversion_precision_t<T, P>>;
    {
        scalar_cast<typename ScalarTraits<T>::value_type>(std::declval<conversion_scalar_t<T, P>>())
    } -> std::same_as<T>;
};

template <std::floating_point W, std::floating_point P>
UnitScales<W, scaled> scales_cast(const UnitScales<P, scaled>& s)
{
    return UnitScales<W, scaled> {{.mass = static_cast<W>(s.mass()),
                                   .density = static_cast<W>(s.density()),
                                   .temperature = static_cast<W>(s.temperature()),
                                   .length = static_cast<W>(s.length())}};
}
}  // namespace detail

// The quantity q in unit system To. The kind does not change; a
// dimensionless quantity is returned unchanged. Host-side: called when
// models are built, not in per-cell code. The result retains scalar type T.
// Arithmetic uses the common precision of double, P and T's value_type;
// ScalarTraits::cast converts all scalar components when precision changes.
template <UnitSystem To, KindTag K, Scalar T, UnitTag From, std::floating_point P, UnitSystem S>
    requires ConvertibleUnits<K, From, To, S> && ScalarPair<T, P> &&
             (std::same_as<From, unitless> || std::same_as<From, To> ||
              detail::ConvertibleScalar<T, P>)
auto convert(Quantity<K, T, From> q, const UnitScales<P, S>& scales)
{
    if constexpr (std::same_as<From, unitless> || std::same_as<From, To>) {
        return q;
    } else {
        using W = detail::conversion_precision_t<T, P>;
        using WideScalar = detail::conversion_scalar_t<T, P>;
        const W factor = scale_factor<K>(detail::scales_cast<W>(scales));
        const WideScalar value = detail::scalar_cast<W>(q.value());
        const WideScalar converted = [&] {
            if constexpr (std::same_as<From, si>) {
                return WideScalar(value / factor);
            } else {
                return WideScalar(value * factor);
            }
        }();
        return Quantity<K, T, To> {
          detail::scalar_cast<typename ScalarTraits<T>::value_type>(converted)};
    }
}

namespace detail
{
// The scales with double values, for converting constants whose scale
// factors may lie outside the range of P.
template <std::floating_point P> UnitScales<double, si> double_scales(const UnitScales<P, si>&)
{ return {}; }

template <std::floating_point P>
UnitScales<double, scaled> double_scales(const UnitScales<P, scaled>& s)
{
    return UnitScales<double, scaled> {{
      .mass = static_cast<double>(s.mass()),
      .density = static_cast<double>(s.density()),
      .temperature = static_cast<double>(s.temperature()),
      .length = static_cast<double>(s.length()),
    }};
}
}  // namespace detail

// The Boltzmann constant in the system of the scales: the exact SI value, or
// exactly 1 in the scaled system, where the velocity scale is defined by it.
template <std::floating_point P>
constexpr BoltzmannConstant<P, si> boltzmann_constant(const UnitScales<P, si>& /*scales*/)
{ return BoltzmannConstant<P, si> {constants::boltzmann<P>}; }

template <std::floating_point P>
constexpr BoltzmannConstant<P, scaled> boltzmann_constant(const UnitScales<P, scaled>& /*scales*/)
{ return BoltzmannConstant<P, scaled> {P(1)}; }

}  // namespace canary

#endif  // CANARY_UNITS_HPP
