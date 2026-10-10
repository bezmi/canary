#ifndef GASLIB_TRANSPORT_POWER_LAW_HPP
#define GASLIB_TRANSPORT_POWER_LAW_HPP

#include <concepts>
#include <gaslib/config.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/scalar.hpp>

namespace gaslib
{

// Power law: mu = mu_0 (T / T_0)^n.
template <std::floating_point P, UnitSystem U> struct PowerLaw {
    DynamicViscosity<P, U> mu_0;
    Temperature<P, U> t_0;
    P exponent;
};

// The law is an empirical fit, evaluated on value parts.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P> && HasPow<T>
GASLIB_HOST_DEVICE DynamicViscosity<T, U> viscosity(const PowerLaw<P, U>& law,
                                                    Temperature<T, U> temp)
{
    const T ratio = temp.value() / law.t_0.value();
    return DynamicViscosity<T, U> {T(law.mu_0.value() * math::pow(ratio, T(law.exponent)))};
}

}  // namespace gaslib

#endif  // GASLIB_TRANSPORT_POWER_LAW_HPP
