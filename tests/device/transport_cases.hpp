#ifndef CANARY_TESTS_DEVICE_TRANSPORT_CASES_HPP
#define CANARY_TESTS_DEVICE_TRANSPORT_CASES_HPP

// Transport checks for the host and native compilations (see cell_inputs.hpp).
// Transport operations take owning states, which the caller supplies.

#include "cell_inputs.hpp"

#include <canary/models/calorically_perfect.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/state.hpp>
#include <canary/transport/constant.hpp>
#include <canary/transport/mixture.hpp>
#include <canary/transport/power_law.hpp>
#include <canary/transport/prandtl.hpp>
#include <canary/transport/sutherland.hpp>
#include <canary/transport/transport.hpp>
#include <concepts>
#include <cstddef>
#include <utility>

namespace canary::test
{

// The including file defines CANARY_CHECK_ENTRY (see thermodynamics_cases.hpp).

template <typename L, typename T, UnitSystem U>
CANARY_CHECK_ENTRY void
  viscosity_check(Output<DynamicViscosity<T, U>> out, const L* law, Temperature<T, U> temp)
{ out[0] = viscosity(*law, temp); }

// Transport of a model without composition input at an owning state.
template <typename X, typename G, typename State>
CANARY_CHECK_ENTRY void transport_check(
  Output<decltype(transport(
    std::declval<const X&>(), std::declval<const G&>(), std::declval<const State&>()))> out,
  const X* model,
  const G* gas,
  State state)
{ out[0] = transport(*model, *gas, state); }

// Transport and diffusion of a mixture model at an owning state, with the
// mixture prepared from the mass fractions of the selected cell.
template <typename G, typename State> using CellMixture = mixture_t<G, typename State::value_type>;

template <Storage S, std::size_t N, typename X, typename G, typename State>
CANARY_CHECK_ENTRY void
  mixture_transport_check(Output<decltype(transport(std::declval<const X&>(),
                                                    std::declval<const G&>(),
                                                    std::declval<CellMixture<G, State>>(),
                                                    std::declval<const State&>()))> out,
                          const X* model,
                          const G* gas,
                          const typename State::value_type* composition,
                          State state)
{
    out[0] = transport(*model,
                       *gas,
                       cell_composition<Composition::Mixture, S, N>(*gas, composition),
                       state);
}

template <std::size_t N, typename X, typename G, typename State>
CANARY_CHECK_ENTRY void binary_diffusion_check(
  Output<decltype(binary_diffusion_coefficients(std::declval<const X&>(),
                                                std::declval<const G&>(),
                                                std::declval<CellMixture<G, State>>(),
                                                std::declval<const State&>()))> out,
  const X* model,
  const G* gas,
  const typename State::value_type* composition,
  State state)
{
    out[0] = binary_diffusion_coefficients(
      *model,
      *gas,
      cell_composition<Composition::Mixture, Storage::Contiguous, N>(*gas, composition),
      state);
}

template <std::size_t N, typename X, typename G, typename State>
CANARY_CHECK_ENTRY void mixture_diffusion_check(
  Output<decltype(mixture_diffusion_coefficients(std::declval<const X&>(),
                                                 std::declval<const G&>(),
                                                 std::declval<CellMixture<G, State>>(),
                                                 std::declval<const State&>()))> out,
  const X* model,
  const G* gas,
  const typename State::value_type* composition,
  State state)
{
    out[0] = mixture_diffusion_coefficients(
      *model,
      *gas,
      cell_composition<Composition::Mixture, Storage::Contiguous, N>(*gas, composition),
      state);
}

// Required instances; see thermodynamics_cases.hpp for the convention.

// Transport, full binary and mixture-averaged diffusion of one mixture model on
// one gas model.
template <std::size_t N, typename X, typename G, typename T> void require_mixture_transport()
{
    using State = state_t<G, T>;
    (void)&mixture_transport_check<Storage::Contiguous, N, X, G, State>;
    (void)&binary_diffusion_check<N, X, G, State>;
    (void)&mixture_diffusion_check<N, X, G, State>;
}

// Wilke and Gupta-Yos on a thermally perfect gas of N species (structural
// species counts repeat air5 records).
template <std::floating_point P, UnitSystem U, std::size_t N, typename T>
void require_structural_transport()
{
    using Gas = ThermallyPerfectMixture<P, U, N, 4>;
    require_mixture_transport<N, WilkeTransport<P, U, N, 3>, Gas, T>();
    require_mixture_transport<N, GuptaYosTransport<P, U, N>, Gas, T>();
}

// Structural species counts 1, 2 and 16, for every P/T pair.
template <std::floating_point P, typename T, UnitSystem U> void transport_structural_checks()
{
    require_structural_transport<P, U, 1, T>();
    require_structural_transport<P, U, 2, T>();
    require_structural_transport<P, U, 16, T>();
}

template <std::floating_point P, typename T, UnitSystem U> void transport_checks()
{
    // Simple laws: viscosity, and transport with a constant Prandtl number on a
    // perfect gas and a five-species perfect mixture.
    const auto simple = [](auto* law) {
        using L = std::remove_pointer_t<decltype(law)>;
        using X = Transport<L, ConstantPrandtl<P>>;
        (void)&viscosity_check<L, T, U>;
        (void)&transport_check<X, CaloricallyPerfect<P, U>, ThermoState<T, U>>;
        (void)&mixture_transport_check<Storage::Contiguous,
                                       5,
                                       X,
                                       CaloricallyPerfectMixture<P, U, 5>,
                                       ThermoState<T, U>>;
    };
    simple(static_cast<ConstantViscosity<P, U>*>(nullptr));
    simple(static_cast<Sutherland<P, U>*>(nullptr));
    simple(static_cast<PowerLaw<P, U>*>(nullptr));
    // Named air5 Wilke and Gupta-Yos on the one- and two-temperature models.
    using Wilke = WilkeTransport<P, U, 5, 3>;
    using GuptaYos = GuptaYosTransport<P, U, 5>;
    using OneTemperature = ThermallyPerfectMixture<P, U, 5, 4>;
    using TwoTemperature = TwoTemperatureMixture<P, U, 5, 4>;
    require_mixture_transport<5, Wilke, OneTemperature, T>();
    require_mixture_transport<5, GuptaYos, OneTemperature, T>();
    require_mixture_transport<5, Wilke, TwoTemperature, T>();
    require_mixture_transport<5, GuptaYos, TwoTemperature, T>();
    // Each grid storage of the composition once.
    (void)&mixture_transport_check<Storage::LayoutRight,
                                   5,
                                   Wilke,
                                   OneTemperature,
                                   ThermoState<T, U>>;
    (void)&mixture_transport_check<Storage::LayoutLeft,
                                   5,
                                   GuptaYos,
                                   OneTemperature,
                                   ThermoState<T, U>>;
    transport_structural_checks<P, T, U>();
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_TRANSPORT_CASES_HPP
