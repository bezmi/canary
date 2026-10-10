#ifndef GASLIB_UNITS_HPP
#define GASLIB_UNITS_HPP

#include <cmath>
#include <concepts>
#include <gaslib/build.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/scalar.hpp>
#include <type_traits>

namespace gaslib
{

namespace constants
{
// Exact SI values (2019 definition of the SI base units).
template <std::floating_point P> inline constexpr P boltzmann = P(1.380649e-23);             // J/K
template <std::floating_point P> inline constexpr P elementary_charge = P(1.602176634e-19);  // C
template <std::floating_point P> inline constexpr P avogadro = P(6.02214076e23);  // 1/mol
// Universal gas constant N_A k_B [J/(mol K)], the product of two exact
// values rounded to P.
template <std::floating_point P>
inline constexpr P universal_gas_constant = avogadro<P> * boltzmann<P>;
}  // namespace constants

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
template <std::floating_point P> struct UnitScales<P, scaled> {
    P mass;
    P density;
    P temperature;
    P length;
};

// Derived scales, computed in the order of cerberus's MFP::update_ref.

// u_ref = sqrt(k_B T_ref / m_ref) [m/s]
template <std::floating_point P> P velocity_scale(const UnitScales<P, scaled>& s)
{ return std::sqrt(constants::boltzmann<P> * s.temperature / s.mass); }

// t_ref = x_ref / u_ref [s]
template <std::floating_point P> P time_scale(const UnitScales<P, scaled>& s)
{ return s.length / velocity_scale(s); }

// Number of particles in a reference volume, n0 = (rho_ref / m_ref) x_ref^3.
// It is the scale of the particle-count dimension, so that the density
// scale is rho_ref and the particle-mass scale is m_ref.
template <std::floating_point P> P count_scale(const UnitScales<P, scaled>& s)
{
    const P number_density = s.density / s.mass;
    return number_density * s.length * s.length * s.length;
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
template <std::floating_point P> P scale_factor(const Dimension& d, const UnitScales<P, scaled>& s)
{
    if (d.amount != 0) {
        throw ParameterError("a dimension with an amount of substance has no scaled unit");
    }
    const P count = count_scale(s);
    P numerator = P(1);
    P denominator = P(1);
    detail::apply_power(numerator, denominator, s.mass * count, d.mass);
    detail::apply_power(numerator, denominator, s.length, d.length);
    detail::apply_power(numerator, denominator, time_scale(s), d.time);
    detail::apply_power(numerator, denominator, s.temperature, d.temperature);
    detail::apply_power(numerator, denominator, constants::elementary_charge<P> * count, d.charge);
    detail::apply_power(numerator, denominator, count, d.count);
    return numerator / denominator;
}

template <std::floating_point P>
constexpr P scale_factor(const Dimension& /*d*/, const UnitScales<P, si>& /*scales*/)
{ return P(1); }

// SI value of one unit of kind K in the system of the scales.
template <KindTag K, std::floating_point P>
    requires(K::dimension.amount == 0)
P scale_factor(const UnitScales<P, scaled>& s)
{ return scale_factor(K::dimension, s); }

template <KindTag K, std::floating_point P>
constexpr P scale_factor(const UnitScales<P, si>& /*scales*/)
{ return P(1); }

// A conversion of kind K from unit system From to To with scales of system
// S: between SI and S, within one system, or of a dimensionless kind.
template <typename K, typename From, typename To, typename S>
concept ConvertibleUnits =
  KindTag<K> && UnitTag<From> && UnitSystem<To> && UnitSystem<S> &&
  (std::same_as<From, unitless> ||
   (KindUnits<K, To> &&
    (std::same_as<From, To> || (std::same_as<From, si> && std::same_as<To, S>) ||
     (std::same_as<From, S> && std::same_as<To, si>))));

// The quantity q in unit system To. The kind does not change; a
// dimensionless quantity is returned unchanged. Host-side: called when
// models are built, not in per-cell code.
template <UnitSystem To, KindTag K, Scalar T, UnitTag From, std::floating_point P, UnitSystem S>
    requires ConvertibleUnits<K, From, To, S> && ScalarPair<T, P>
auto convert(Quantity<K, T, From> q, const UnitScales<P, S>& scales)
{
    if constexpr (std::same_as<From, unitless> || std::same_as<From, To>) {
        return q;
    } else if constexpr (std::same_as<From, si>) {
        const auto value = q.value() / scale_factor<K>(scales);
        return Quantity<K, std::remove_const_t<decltype(value)>, To> {value};
    } else {
        const auto value = q.value() * scale_factor<K>(scales);
        return Quantity<K, std::remove_const_t<decltype(value)>, To> {value};
    }
}

// The Boltzmann constant in the system of the scales: the exact SI value, or
// exactly 1 in the scaled system, where the velocity scale is defined by it.
template <std::floating_point P>
constexpr BoltzmannConstant<P, si> boltzmann_constant(const UnitScales<P, si>& /*scales*/)
{ return BoltzmannConstant<P, si> {constants::boltzmann<P>}; }

template <std::floating_point P>
constexpr BoltzmannConstant<P, scaled> boltzmann_constant(const UnitScales<P, scaled>& /*scales*/)
{ return BoltzmannConstant<P, scaled> {P(1)}; }

}  // namespace gaslib

#endif  // GASLIB_UNITS_HPP
