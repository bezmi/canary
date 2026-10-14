#ifndef CANARY_TRANSPORT_SUTHERLAND_HPP
#define CANARY_TRANSPORT_SUTHERLAND_HPP

#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <cmath>
#include <concepts>
#include <format>

namespace canary
{

// Sutherland's law: mu = mu_0 (T / T_0)^(3/2) (T_0 + S) / (T + S).
// mu_0 must be finite and non-negative, t_0 finite and positive, and s finite
// with t_0 + s > 0. Evaluation requires finite positive temperature within
// the fit range and temp + s > 0. Inputs and arithmetic overflow are unchecked.
template <std::floating_point P, UnitSystem U> struct Sutherland {
    DynamicViscosity<P, U> mu_0;
    Temperature<P, U> t_0;
    Temperature<P, U> s;
};

// Host-side construction in unit system U. Throws ParameterError when mu_0
// is nonfinite or negative, t_0 is nonfinite or nonpositive, s is nonfinite,
// or t_0 + s is nonpositive. Negative s is permitted under this constraint;
// callers must still satisfy temp + s > 0 during unchecked evaluation.
template <std::floating_point P, UnitSystem U>
Sutherland<P, U>
  make_sutherland(DynamicViscosity<P, U> mu_0, Temperature<P, U> t_0, Temperature<P, U> s)
{
    if (!std::isfinite(mu_0.value()) || mu_0.value() < P(0)) {
        throw ParameterError(
          std::format("mu_0 must be finite and non-negative, got {}", mu_0.value()));
    }
    if (!std::isfinite(t_0.value()) || !(t_0.value() > P(0))) {
        throw ParameterError(std::format("t_0 must be finite and positive, got {}", t_0.value()));
    }
    if (!std::isfinite(s.value())) {
        throw ParameterError(std::format("s must be finite, got {}", s.value()));
    }
    // For finite t_0 > 0, this tests the sum's sign without forming an overflowing sum.
    if (!(s.value() > -t_0.value())) {
        throw ParameterError(
          std::format("t_0 + s must be positive, got t_0={} and s={}", t_0.value(), s.value()));
    }
    return {.mu_0 = mu_0, .t_0 = t_0, .s = s};
}

// The empirical fit uses the full scalar and propagates temperature derivatives.
// Evaluation order: mu_0 * r * sqrt(r) * (T_0 + S) / (T + S), with r = T / T_0.
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
