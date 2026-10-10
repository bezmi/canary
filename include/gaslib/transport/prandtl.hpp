#ifndef GASLIB_TRANSPORT_PRANDTL_HPP
#define GASLIB_TRANSPORT_PRANDTL_HPP

#include <concepts>
#include <gaslib/config.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/scalar.hpp>
#include <gaslib/state.hpp>

namespace gaslib
{

// Thermal conductivity at a constant Prandtl number: k = mu cp / Pr.
template <std::floating_point P> struct ConstantPrandtl {
    PrandtlNumber<P> prandtl;
};

template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE ThermalConductivity<T, U> conductivity(const ConstantPrandtl<P>& law,
                                                          DynamicViscosity<T, U> mu,
                                                          const Properties<T, U>& properties)
{ return ThermalConductivity<T, U> {conductivity(mu, properties.cp, law.prandtl)}; }

}  // namespace gaslib

#endif  // GASLIB_TRANSPORT_PRANDTL_HPP
