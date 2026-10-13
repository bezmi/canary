#ifndef CANARY_TRANSPORT_POWER_LAW_HPP
#define CANARY_TRANSPORT_POWER_LAW_HPP

#include <concepts>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>

namespace canary
{

// Power law: mu = mu_0 (T / T_0)^n. mu_0 must be finite and non-negative,
// t_0 finite and positive, and exponent finite. Evaluation requires finite
// positive temperature within the fit range and does not validate inputs
// or detect arithmetic overflow.
template <std::floating_point P, UnitSystem U> struct PowerLaw {
    DynamicViscosity<P, U> mu_0;
    Temperature<P, U> t_0;
    P exponent;
};

// The empirical fit uses the full scalar and propagates temperature derivatives.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P> && HasPow<T>
CANARY_HOST_DEVICE DynamicViscosity<T, U> viscosity(const PowerLaw<P, U>& law,
                                                    Temperature<T, U> temp)
{
    const T ratio = temp.value() / law.t_0.value();
    return DynamicViscosity<T, U> {T(law.mu_0.value() * math::pow(ratio, T(law.exponent)))};
}

}  // namespace canary

#endif  // CANARY_TRANSPORT_POWER_LAW_HPP
