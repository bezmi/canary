#ifndef CANARY_TRANSPORT_PRANDTL_HPP
#define CANARY_TRANSPORT_PRANDTL_HPP

#include <concepts>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/state.hpp>

namespace canary
{

// Thermal conductivity at a constant Prandtl number: k = mu cp / Pr.
template <std::floating_point P> struct ConstantPrandtl {
    PrandtlNumber<P> prandtl;
};

template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE ThermalConductivity<T, U> conductivity(const ConstantPrandtl<P>& law,
                                                          DynamicViscosity<T, U> mu,
                                                          const Properties<T, U>& properties)
{ return ThermalConductivity<T, U> {conductivity(mu, properties.cp, law.prandtl)}; }

}  // namespace canary

#endif  // CANARY_TRANSPORT_PRANDTL_HPP
