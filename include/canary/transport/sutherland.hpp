#ifndef CANARY_TRANSPORT_SUTHERLAND_HPP
#define CANARY_TRANSPORT_SUTHERLAND_HPP

#include <concepts>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>

namespace canary
{

// Sutherland's law: mu = mu_0 (T / T_0)^(3/2) (T_0 + S) / (T + S).
template <std::floating_point P, UnitSystem U> struct Sutherland {
    DynamicViscosity<P, U> mu_0;
    Temperature<P, U> t_0;
    Temperature<P, U> s;
};

// The law is an empirical fit, evaluated on value parts in the order
// mu_0 * r * sqrt(r) * (T_0 + S) / (T + S) with r = T / T_0.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P> && HasSqrt<T>
CANARY_HOST_DEVICE DynamicViscosity<T, U> viscosity(const Sutherland<P, U>& law,
                                                    Temperature<T, U> temp)
{
    const T ratio = temp.value() / law.t_0.value();
    return DynamicViscosity<T, U> {
      T(law.mu_0.value() * ratio * math::sqrt(ratio) * (law.t_0.value() + law.s.value()) /
        (temp.value() + law.s.value())),
    };
}

}  // namespace canary

#endif  // CANARY_TRANSPORT_SUTHERLAND_HPP
