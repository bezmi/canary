// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_TRANSPORT_POWER_LAW_HPP
#define CANARY_TRANSPORT_POWER_LAW_HPP

#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <cmath>
#include <concepts>
#include <format>

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

// Host-side construction in unit system U. Throws ParameterError when mu_0
// is nonfinite or negative, t_0 is nonfinite or nonpositive, or exponent is
// nonfinite. The returned law does not check evaluation temperatures.
template <std::floating_point P, UnitSystem U>
PowerLaw<P, U> make_power_law(DynamicViscosity<P, U> mu_0, Temperature<P, U> t_0, P exponent)
{
    if (!std::isfinite(mu_0.value()) || mu_0.value() < P(0)) {
        throw ParameterError(
          std::format("mu_0 must be finite and non-negative, got {}", mu_0.value()));
    }
    if (!std::isfinite(t_0.value()) || !(t_0.value() > P(0))) {
        throw ParameterError(std::format("t_0 must be finite and positive, got {}", t_0.value()));
    }
    if (!std::isfinite(exponent)) {
        throw ParameterError(std::format("exponent must be finite, got {}", exponent));
    }
    return {.mu_0 = mu_0, .t_0 = t_0, .exponent = exponent};
}

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
