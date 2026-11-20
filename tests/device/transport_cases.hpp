// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_TESTS_DEVICE_TRANSPORT_CASES_HPP
#define CANARY_TESTS_DEVICE_TRANSPORT_CASES_HPP

// Transport checks for the host and native compilations (see cell_inputs.hpp).
// Transport operations take owning states, which the caller supplies.

#include "cell_inputs.hpp"
#include "check_groups.hpp"
#include "derived_inputs.hpp"
#include "host_inputs.hpp"
#include "native_runtime.hpp"
#include "transport_inputs.hpp"

#include <algorithm>
#include <array>
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
#include <string>
#include <utility>
#include <vector>

namespace canary::test
{

// Entries: see thermodynamics_cases.hpp.

template <typename L, typename T, UnitSystem U>
CANARY_HOST_DEVICE void
  viscosity_check(Output<DynamicViscosity<T, U>> out, const L* law, Temperature<T, U> temp)
{ out[0] = viscosity(*law, temp); }

// Transport of a model without composition input at an owning state.
template <typename X, typename G, typename State>
CANARY_HOST_DEVICE void transport_check(
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
CANARY_HOST_DEVICE void
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
CANARY_HOST_DEVICE void binary_diffusion_check(
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
CANARY_HOST_DEVICE void mixture_diffusion_check(
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

// Runners. Each evaluates the entries of one transport model at every state
// of a fixture and every seed (see host_inputs.hpp).

// Transport, full binary and mixture-averaged diffusion of model on gas at
// every state, with the composition in storage S; diffusion only with
// contiguous storage.
template <Storage S,
          bool Diffusion,
          std::size_t N,
          typename X,
          typename G,
          typename T,
          typename Declared,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void mixture_transport(Report& report,
                       const std::string& name,
                       const X& model,
                       const G& gas,
                       const std::array<Declared, C>& states,
                       const std::array<DerivedState<M>, C>& derived,
                       const UnitScales<P, U>& scales)
{
    using Owning = state_t<G, T>;
    for (std::size_t s = 0; s < C; ++s) {
        for_each_seed<T>([&](Seeds<T>& seeds) {
            const auto composition = cell_storage<S>(seeds.fractions("y", states[s].y));
            const Owning state = owning_state<Owning>(seeds, states[s], derived[s], scales);
            const std::string label = state_label(name, s, seeds.label());
            evaluate<&mixture_transport_check<S, N, X, G, Owning>>(report,
                                                                   label,
                                                                   model,
                                                                   gas,
                                                                   composition,
                                                                   state);
            if constexpr (Diffusion) {
                evaluate<&binary_diffusion_check<N, X, G, Owning>>(report,
                                                                   label,
                                                                   model,
                                                                   gas,
                                                                   composition,
                                                                   state);
                evaluate<&mixture_diffusion_check<N, X, G, Owning>>(report,
                                                                    label,
                                                                    model,
                                                                    gas,
                                                                    composition,
                                                                    state);
            }
        });
    }
}

// Wilke and Gupta-Yos transport and diffusion on N structural species, for
// every order of the species.
template <std::size_t N, typename T, std::size_t D, std::floating_point P, UnitSystem U>
void structural_transport(Report& report,
                          const std::array<DerivedState<0>, D>& derived,
                          const UnitScales<P, U>& scales)
{
    using Gas = ThermallyPerfectMixture<P, U, N, 4>;
    const auto orders = structural_orders<N>();
    static_assert(D % 6 == 0);
    if (orders.size() * 6 != D) {
        report.error("structural fixture of " + std::to_string(N) +
                     " species: derived inputs do not match the orders");
        return;
    }
    for (std::size_t o = 0; o < orders.size(); ++o) {
        const auto models = make_structural_models<N>(orders[o], scales);
        const auto states = structural_states(orders[o]);
        std::array<DerivedState<0>, 6> order_derived {};
        std::copy_n(derived.begin() + static_cast<std::ptrdiff_t>(o * 6), 6, order_derived.begin());
        const std::string name =
          std::to_string(N) + " structural species, order " + std::to_string(o);
        mixture_transport<Storage::Contiguous, true, N, WilkeTransport<P, U, N, 3>, Gas, T>(
          report,
          name + ", Wilke",
          models.wilke,
          models.gas,
          states,
          order_derived,
          scales);
        mixture_transport<Storage::Contiguous, true, N, GuptaYosTransport<P, U, N>, Gas, T>(
          report,
          name + ", Gupta-Yos",
          models.gupta_yos,
          models.gas,
          states,
          order_derived,
          scales);
    }
}

// Structural species counts 1, 2 and 16, for every P/T pair.
template <std::floating_point P, typename T, UnitSystem U> Report transport_structural_checks()
{
    Report report;
    const auto scales = fixture_scales<P, U>();
    structural_transport<1, T>(report, structural1_derived, scales);
    structural_transport<2, T>(report, structural2_derived, scales);
    structural_transport<16, T>(report, structural16_derived, scales);
    return report;
}

// A simple viscosity law alone and with a constant Prandtl number, on the
// perfect gas and the five-species perfect mixture.
template <typename T, typename L, std::floating_point P, UnitSystem U>
void simple_law(Report& report,
                const std::string& name,
                const L& law,
                const UnitScales<P, U>& scales)
{
    using X = Transport<L, ConstantPrandtl<P>>;
    using Perfect = CaloricallyPerfect<P, U>;
    using Mixture5 = CaloricallyPerfectMixture<P, U, 5>;
    const X model {law, make_prandtl<P>()};
    const Perfect perfect = make_perfect_gas(scales);
    const Mixture5 mixture5 = make_perfect_mixture<5>(5, scales);
    const auto mixture_states = simple_law_mixture_states();
    for (std::size_t s = 0; s < simple_law_states.size(); ++s) {
        for_each_seed<T>([&](Seeds<T>& seeds) {
            const auto temp =
              seeds.template quantity<kind::temperature>("temp", simple_law_states[s].temp, scales);
            evaluate<&viscosity_check<L, T, U>>(report,
                                                state_label(name, s, seeds.label()),
                                                law,
                                                temp);
        });
        for_each_seed<T>([&](Seeds<T>& seeds) {
            const auto state = owning_state<ThermoState<T, U>>(seeds,
                                                               simple_law_states[s],
                                                               simple_law_perfect_derived[s],
                                                               scales);
            evaluate<&transport_check<X, Perfect, ThermoState<T, U>>>(
              report,
              state_label(name + " on the perfect gas", s, seeds.label()),
              model,
              perfect,
              state);
        });
    }
    mixture_transport<Storage::Contiguous, false, 5, X, Mixture5, T>(report,
                                                                     name +
                                                                       " on the perfect mixture",
                                                                     model,
                                                                     mixture5,
                                                                     mixture_states,
                                                                     simple_law_mixture_derived,
                                                                     scales);
}

template <std::floating_point P, typename T, UnitSystem U> Report transport_checks()
{
    Report report;
    const auto scales = fixture_scales<P, U>();
    // Simple laws: viscosity, and transport with a constant Prandtl number on a
    // perfect gas and a five-species perfect mixture.
    simple_law<T>(report, "constant viscosity", make_constant_law(scales), scales);
    simple_law<T>(report, "Sutherland", make_sutherland_law(scales), scales);
    simple_law<T>(report, "power law", make_power_law_law(scales), scales);
    // Named air5 Wilke and Gupta-Yos on the one- and two-temperature models.
    using Wilke = WilkeTransport<P, U, 5, 3>;
    using GuptaYos = GuptaYosTransport<P, U, 5>;
    using OneTemperature = ThermallyPerfectMixture<P, U, 5, 4>;
    using TwoTemperature = TwoTemperatureMixture<P, U, 5, 4>;
    {
        const auto built = make_air5(scales);
        const Wilke wilke = make_air5_wilke_transport(built.species, scales);
        const GuptaYos gupta_yos = make_air5_gupta_yos_transport(built.species, scales);
        const auto& states = named_transport_states;
        const auto& derived = named_transport_derived;
        mixture_transport<Storage::Contiguous, true, 5, Wilke, OneTemperature, T>(report,
                                                                                  "air5 Wilke",
                                                                                  wilke,
                                                                                  built.model,
                                                                                  states,
                                                                                  derived,
                                                                                  scales);
        mixture_transport<Storage::Contiguous, true, 5, GuptaYos, OneTemperature, T>(
          report,
          "air5 Gupta-Yos",
          gupta_yos,
          built.model,
          states,
          derived,
          scales);
        // Each grid storage of the composition once.
        mixture_transport<Storage::LayoutRight, false, 5, Wilke, OneTemperature, T>(report,
                                                                                    "air5 Wilke",
                                                                                    wilke,
                                                                                    built.model,
                                                                                    states,
                                                                                    derived,
                                                                                    scales);
        mixture_transport<Storage::LayoutLeft, false, 5, GuptaYos, OneTemperature, T>(
          report,
          "air5 Gupta-Yos",
          gupta_yos,
          built.model,
          states,
          derived,
          scales);
    }
    {
        const auto built = make_air5_two_temperature(scales);
        const auto& states = named_transport_states;
        const auto& derived = named_transport_two_temperature_derived;
        mixture_transport<Storage::Contiguous, true, 5, Wilke, TwoTemperature, T>(
          report,
          "two-temperature air5 Wilke",
          make_air5_wilke_transport(built.species, scales),
          built.model,
          states,
          derived,
          scales);
        mixture_transport<Storage::Contiguous, true, 5, GuptaYos, TwoTemperature, T>(
          report,
          "two-temperature air5 Gupta-Yos",
          make_air5_gupta_yos_transport(built.species, scales),
          built.model,
          states,
          derived,
          scales);
    }
    report.merge(transport_structural_checks<P, T, U>());
    return report;
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_TRANSPORT_CASES_HPP
