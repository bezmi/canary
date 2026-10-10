#ifndef GASLIB_TRANSPORT_CONSTANT_HPP
#define GASLIB_TRANSPORT_CONSTANT_HPP

#include <concepts>
#include <gaslib/config.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/scalar.hpp>

namespace gaslib
{

// Viscosity independent of temperature.
template <std::floating_point P, UnitSystem U> struct ConstantViscosity {
    DynamicViscosity<P, U> mu;
};

template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE DynamicViscosity<T, U> viscosity(const ConstantViscosity<P, U>& law,
                                                    Temperature<T, U> /*temp*/)
{ return DynamicViscosity<T, U> {law.mu}; }

}  // namespace gaslib

#endif  // GASLIB_TRANSPORT_CONSTANT_HPP
