#ifndef CANARY_TRANSPORT_CONSTANT_HPP
#define CANARY_TRANSPORT_CONSTANT_HPP

#include <concepts>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>

namespace canary
{

// Viscosity independent of temperature. mu must be finite and non-negative;
// evaluation does not validate it. The temperature argument is unused.
template <std::floating_point P, UnitSystem U> struct ConstantViscosity {
    DynamicViscosity<P, U> mu;
};

template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE DynamicViscosity<T, U> viscosity(const ConstantViscosity<P, U>& law,
                                                    Temperature<T, U> /*temp*/)
{ return DynamicViscosity<T, U> {law.mu}; }

}  // namespace canary

#endif  // CANARY_TRANSPORT_CONSTANT_HPP
