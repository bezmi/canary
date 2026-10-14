#ifndef CANARY_TRANSPORT_CONSTANT_HPP
#define CANARY_TRANSPORT_CONSTANT_HPP

#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <cmath>
#include <concepts>
#include <format>

namespace canary
{

// Viscosity independent of temperature. mu must be finite and non-negative;
// evaluation does not validate it. The temperature argument is unused.
template <std::floating_point P, UnitSystem U> struct ConstantViscosity {
    DynamicViscosity<P, U> mu;
};

// Host-side construction in unit system U. Throws ParameterError when mu
// is nonfinite or negative. Zero viscosity is permitted.
template <std::floating_point P, UnitSystem U>
ConstantViscosity<P, U> make_constant_viscosity(DynamicViscosity<P, U> mu)
{
    if (!std::isfinite(mu.value()) || mu.value() < P(0)) {
        throw ParameterError(std::format("mu must be finite and non-negative, got {}", mu.value()));
    }
    return {.mu = mu};
}

template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE DynamicViscosity<T, U> viscosity(const ConstantViscosity<P, U>& law,
                                                    Temperature<T, U> /*temp*/)
{ return DynamicViscosity<T, U> {law.mu}; }

}  // namespace canary

#endif  // CANARY_TRANSPORT_CONSTANT_HPP
