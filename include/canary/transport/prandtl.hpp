#ifndef CANARY_TRANSPORT_PRANDTL_HPP
#define CANARY_TRANSPORT_PRANDTL_HPP

#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/state.hpp>
#include <cmath>
#include <concepts>
#include <format>

namespace canary
{

// Thermal conductivity at a constant Prandtl number: k = mu cp / Pr.
// prandtl must be finite and positive. Evaluation requires finite
// non-negative viscosity and finite positive properties.cp; inputs and
// arithmetic overflow are unchecked.
template <std::floating_point P> struct ConstantPrandtl {
    PrandtlNumber<P> prandtl;
};

// Host-side construction. Throws ParameterError when prandtl is nonfinite
// or nonpositive. The returned law uses the stored scalar type P unchanged.
template <std::floating_point P> ConstantPrandtl<P> make_constant_prandtl(PrandtlNumber<P> prandtl)
{
    if (!std::isfinite(prandtl.value()) || !(prandtl.value() > P(0))) {
        throw ParameterError(
          std::format("prandtl must be finite and positive, got {}", prandtl.value()));
    }
    return {.prandtl = prandtl};
}

template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE ThermalConductivity<T, U> conductivity(const ConstantPrandtl<P>& law,
                                                          DynamicViscosity<T, U> mu,
                                                          const Properties<T, U>& properties)
{ return ThermalConductivity<T, U> {conductivity(mu, properties.cp, law.prandtl)}; }

}  // namespace canary

#endif  // CANARY_TRANSPORT_PRANDTL_HPP
