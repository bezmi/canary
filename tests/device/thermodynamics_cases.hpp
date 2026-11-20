// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_TESTS_DEVICE_THERMODYNAMICS_CASES_HPP
#define CANARY_TESTS_DEVICE_THERMODYNAMICS_CASES_HPP

// Thermodynamic checks for the host and native compilations (see cell_inputs.hpp).

#include "cell_inputs.hpp"
#include "check_groups.hpp"
#include "derived_inputs.hpp"
#include "host_inputs.hpp"
#include "native_runtime.hpp"
#include "thermodynamics_inputs.hpp"

#include <array>
#include <canary/models/calorically_perfect.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/nasa9.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/state.hpp>
#include <canary/thermo.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace canary::test
{

// Checks. Each entry is a host-and-device function that writes the complete
// public result to out[0]. A backend calls it on the host or from a kernel on
// the device (see native_runtime.hpp).

// Closure of a model without composition input from two state variables.
template <typename G, typename A, typename B>
CANARY_HOST_DEVICE void
  closure_check(Output<ClosureOf<state_t<G, typename A::value_type>>> out, const G* gas, A a, B b)
{ out[0] = thermo_state(*gas, a, b); }

// Closure of a mixture model from the composition of the selected cell, two
// state variables and optional further inputs (a temperature guess).
template <Composition F,
          Storage S,
          std::size_t N,
          typename G,
          typename A,
          typename B,
          typename... Extra>
CANARY_HOST_DEVICE void
  mixture_closure_check(Output<ClosureOf<state_t<G, typename A::value_type>>> out,
                        const G* gas,
                        const typename A::value_type* composition,
                        A a,
                        B b,
                        Extra... extra)
{ out[0] = thermo_state(*gas, cell_composition<F, S, N>(*gas, composition), a, b, extra...); }

// Closure of a modal model with the mode inputs of kind K of the selected cell.
template <Composition F,
          Storage S,
          ModeForm MF,
          KindTag K,
          std::size_t N,
          typename G,
          typename A,
          typename B>
CANARY_HOST_DEVICE void
  modal_closure_check(Output<ClosureOf<state_t<G, typename A::value_type>>> out,
                      const G* gas,
                      const typename A::value_type* composition,
                      A a,
                      B b,
                      const typename A::value_type* modes)
{
    out[0] = thermo_state(*gas,
                          cell_composition<F, S, N>(*gas, composition),
                          a,
                          b,
                          cell_modes<MF, S, K, typename G::units, G::mode_count>(modes));
}

// Mixture of the selected cell from mass fractions.
template <Storage S, std::size_t N, typename G, typename T>
CANARY_HOST_DEVICE void
  fractions_mixture_check(Output<mixture_t<G, T>> out, const G* gas, const T* composition)
{ out[0] = cell_composition<Composition::Mixture, S, N>(*gas, composition); }

// Mixture of the selected cell from partial densities and the density.
template <Storage S, std::size_t N, typename G, typename T>
CANARY_HOST_DEVICE void densities_mixture_check(Output<mixture_t<G, T>> out,
                                                const G* gas,
                                                const T* composition,
                                                Density<T, typename G::units> rho)
{
    out[0] =
      mixture(*gas, cell_composition<Composition::PartialDensities, S, N>(*gas, composition), rho);
}

// Properties of an owning state of a model without composition input.
template <typename G, typename State>
CANARY_HOST_DEVICE void properties_check(
  Output<decltype(properties(std::declval<const G&>(), std::declval<State>()))> out,
  const G* gas,
  State state)
{ out[0] = properties(*gas, state); }

// Properties of an owning state of a mixture model, with the mixture prepared
// from the mass fractions of the selected cell.
template <std::size_t N, typename G, typename State>
CANARY_HOST_DEVICE void mixture_properties_check(
  Output<decltype(properties(std::declval<const G&>(),
                             std::declval<mixture_t<G, typename State::value_type>>(),
                             std::declval<State>()))> out,
  const G* gas,
  const typename State::value_type* composition,
  State state)
{
    out[0] =
      properties(*gas,
                 cell_composition<Composition::Mixture, Storage::Contiguous, N>(*gas, composition),
                 state);
}

// Electron pressure of an owning state of a multi-temperature model.
template <std::size_t N, typename G, typename State>
CANARY_HOST_DEVICE void electron_pressure_check(
  Output<decltype(electron_pressure(std::declval<const G&>(),
                                    std::declval<mixture_t<G, typename State::value_type>>(),
                                    std::declval<State>()))> out,
  const G* gas,
  const typename State::value_type* composition,
  State state)
{
    out[0] = electron_pressure(
      *gas,
      cell_composition<Composition::Mixture, Storage::Contiguous, N>(*gas, composition),
      state);
}

// Per-species properties at one temperature.
template <typename G, typename T>
CANARY_HOST_DEVICE void species_properties_check(
  Output<decltype(species_properties(std::declval<const G&>(),
                                     std::declval<Temperature<T, typename G::units>>()))> out,
  const G* gas,
  Temperature<T, typename G::units> temp)
{ out[0] = species_properties(*gas, temp); }

template <typename G, typename T>
CANARY_HOST_DEVICE void species_isobaric_specific_heat_check(
  Output<decltype(species_isobaric_specific_heat(
    std::declval<const G&>(), std::declval<Temperature<T, typename G::units>>()))> out,
  const G* gas,
  Temperature<T, typename G::units> temp)
{ out[0] = species_isobaric_specific_heat(*gas, temp); }

template <typename G, typename T>
CANARY_HOST_DEVICE void species_vibrational_electronic_energy_check(
  Output<decltype(species_vibrational_electronic_energy(
    std::declval<const G&>(), std::declval<Temperature<T, typename G::units>>()))> out,
  const G* gas,
  Temperature<T, typename G::units> temp)
{ out[0] = species_vibrational_electronic_energy(*gas, temp); }

template <typename G, typename T>
CANARY_HOST_DEVICE void species_mode_energy_check(
  Output<decltype(species_mode_energy(std::declval<const G&>(),
                                      ModeIndex {},
                                      std::declval<Temperature<T, typename G::units>>()))> out,
  const G* gas,
  ModeIndex m,
  Temperature<T, typename G::units> temp)
{ out[0] = species_mode_energy(*gas, m, temp); }

// NASA9 cp/R and h/R, and s/R, of a curve at temperature t.
template <typename Curve, typename T>
CANARY_HOST_DEVICE void nasa9_cp_h_check(Output<Nasa9CpH<T>> out, const Curve* curve, T t)
{ out[0] = nasa9_cp_h(*curve, t, math::log(t)); }

template <typename Curve, typename T>
CANARY_HOST_DEVICE void nasa9_s_check(Output<T> out, const Curve* curve, T t)
{ out[0] = nasa9_s(*curve, t, math::log(t)); }

// Runners. Each evaluates the entries of one model at every state of a
// fixture and every seed (see host_inputs.hpp): every closure pair with
// contiguous mass fractions, and each storage, composition form, mode form
// and mode kind at least once on (rho, T).

// A thermodynamic fixture: declared states and the inputs derived from them,
// for closures without mode arguments (equilibrium) and for closures with mode
// arguments and owning states (modal; the same as equilibrium for models
// without modes).
template <std::size_t N, std::size_t K, std::size_t M, std::size_t C> struct ThermoFixture {
    std::string name;
    std::array<FixtureState<N, K>, C> states;
    std::array<DerivedState<M>, C> equilibrium;
    std::array<DerivedState<M>, C> modal;
};

template <std::size_t N, std::size_t M, std::size_t C>
ThermoFixture<N, 0, M, C> single_fixture(std::string name,
                                         const std::array<FixtureState<N, 0>, C>& states,
                                         const std::array<DerivedState<M>, C>& derived)
{ return {std::move(name), states, derived, derived}; }

template <typename T, UnitSystem U>
using ClosurePairs = std::tuple<std::pair<Density<T, U>, SpecificEnergy<T, U>>,
                                std::pair<Density<T, U>, EnergyDensity<T, U>>,
                                std::pair<Density<T, U>, Pressure<T, U>>,
                                std::pair<Density<T, U>, Temperature<T, U>>,
                                std::pair<Pressure<T, U>, Temperature<T, U>>>;
template <typename T, UnitSystem U>
using DensityTemperature = std::tuple<std::pair<Density<T, U>, Temperature<T, U>>>;
template <typename T, UnitSystem U>
using EnergyPairs = std::tuple<std::pair<Density<T, U>, SpecificEnergy<T, U>>,
                               std::pair<Density<T, U>, EnergyDensity<T, U>>>;

// The closure input of quantity type Q at a state of temperature temp [K]
// with derived inputs d.
template <typename Q, typename T, std::size_t M, typename P, typename U>
Q pair_input(Seeds<T>& seeds, double temp, const DerivedState<M>& d, const UnitScales<P, U>& scales)
{
    using K = typename Q::kind;
    if constexpr (std::same_as<K, kind::density>) {
        return seeds.template quantity<K>("rho", d.rho, scales);
    } else if constexpr (std::same_as<K, kind::specific_energy>) {
        return seeds.template quantity<K>("e", d.e, scales);
    } else if constexpr (std::same_as<K, kind::energy_density>) {
        return seeds.template quantity<K>("rho_e", d.rho * d.e, scales);
    } else if constexpr (std::same_as<K, kind::pressure>) {
        return seeds.template quantity<K>("p", d.p, scales);
    } else {
        static_assert(std::same_as<K, kind::temperature>);
        return seeds.template quantity<K>("temp", temp, scales);
    }
}

// The composition input of form F in storage S: raw mass fractions, or the
// partial densities rho y_i.
template <Composition F,
          Storage S,
          typename T,
          std::size_t N,
          std::size_t M,
          typename P,
          typename U>
std::vector<T> composition_input(Seeds<T>& seeds,
                                 const std::array<double, N>& y,
                                 const DerivedState<M>& d,
                                 const UnitScales<P, U>& scales)
{
    if constexpr (F == Composition::PartialDensities) {
        std::array<double, N> rho_y {};
        for (std::size_t i = 0; i < N; ++i) { rho_y[i] = d.rho * y[i]; }
        return cell_storage<S>(seeds.template values<kind::density>("rho_y", rho_y, scales));
    } else {
        return cell_storage<S>(seeds.fractions("y", y));
    }
}

// The mode inputs of kind MK in storage S: the mode temperatures of the state
// or the mode energy densities rho e_m.
template <Storage S,
          KindTag MK,
          typename T,
          std::size_t N,
          std::size_t K,
          std::size_t M,
          typename P,
          typename U>
std::vector<T> mode_input(Seeds<T>& seeds,
                          const FixtureState<N, K>& state,
                          const DerivedState<M>& d,
                          const UnitScales<P, U>& scales)
{
    if constexpr (std::same_as<MK, kind::mode_temperature>) {
        return cell_storage<S>(seeds.template values<MK>("mode_temp", state.modes, scales));
    } else {
        static_assert(std::same_as<MK, kind::mode_energy_density>);
        std::array<double, M> rho_em {};
        for (std::size_t m = 0; m < M; ++m) { rho_em[m] = d.rho * d.mode_e[m]; }
        return cell_storage<S>(seeds.template values<MK>("rho_em", rho_em, scales));
    }
}

// Closures of a model without composition input for each pair.
template <typename G,
          typename T,
          typename... Pair,
          std::size_t N,
          std::size_t K,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void closures(Report& report,
              const G& gas,
              const ThermoFixture<N, K, M, C>& fixture,
              const UnitScales<P, U>& scales,
              std::tuple<Pair...>* /*pairs*/)
{
    const auto pair = [&]<typename A, typename B>(std::pair<A, B>*) {
        for (std::size_t s = 0; s < C; ++s) {
            for_each_seed<T>([&](Seeds<T>& seeds) {
                const auto& d = fixture.equilibrium[s];
                const A a = pair_input<A>(seeds, fixture.states[s].temp, d, scales);
                const B b = pair_input<B>(seeds, fixture.states[s].temp, d, scales);
                evaluate_expecting<&closure_check<G, A, B>>(
                  report,
                  state_label(fixture.name, s, seeds.label()),
                  ClosureExpectation {},
                  gas,
                  a,
                  b);
            });
        }
    };
    (pair(static_cast<Pair*>(nullptr)), ...);
}

// Closures of a mixture model from composition form F in storage S for each
// pair; with Guess, once for each guess factor of the generating temperature.
template <Composition F,
          Storage S,
          bool Guess,
          std::size_t N,
          typename G,
          typename T,
          typename... Pair,
          std::size_t K,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void mixture_closures(Report& report,
                      const G& gas,
                      const ThermoFixture<N, K, M, C>& fixture,
                      const UnitScales<P, U>& scales,
                      std::tuple<Pair...>* /*pairs*/)
{
    using Guesses = std::conditional_t<Guess, decltype(guess_factors), std::array<double, 1>>;
    const Guesses factors = [] {
        if constexpr (Guess) {
            return guess_factors;
        } else {
            return std::array<double, 1> {0.0};
        }
    }();
    const auto pair = [&]<typename A, typename B>(std::pair<A, B>*) {
        for (std::size_t s = 0; s < C; ++s) {
            for (const double factor : factors) {
                for_each_seed<T>([&](Seeds<T>& seeds) {
                    const auto& d = fixture.equilibrium[s];
                    const double temp = fixture.states[s].temp;
                    const auto composition =
                      composition_input<F, S>(seeds, fixture.states[s].y, d, scales);
                    const A a = pair_input<A>(seeds, temp, d, scales);
                    const B b = pair_input<B>(seeds, temp, d, scales);
                    std::string label = state_label(fixture.name, s, seeds.label());
                    if constexpr (Guess) {
                        const TemperatureGuess<T, U> guess {
                          fixed_quantity<kind::temperature, T>(factor * temp, scales)};
                        label += ", guess " + std::to_string(factor) + " T";
                        evaluate_expecting<
                          &mixture_closure_check<F, S, N, G, A, B, TemperatureGuess<T, U>>>(
                          report,
                          label,
                          ClosureExpectation {},
                          gas,
                          composition,
                          a,
                          b,
                          guess);
                    } else {
                        evaluate_expecting<&mixture_closure_check<F, S, N, G, A, B>>(
                          report,
                          label,
                          ClosureExpectation {},
                          gas,
                          composition,
                          a,
                          b);
                    }
                });
            }
        }
    };
    (pair(static_cast<Pair*>(nullptr)), ...);
}

// Closures of a modal model with mode inputs of kind MK in form MF, from
// composition form F in storage S, for each pair.
template <Composition F,
          Storage S,
          ModeForm MF,
          KindTag MK,
          std::size_t N,
          typename G,
          typename T,
          typename... Pair,
          std::size_t K,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void modal_closures(Report& report,
                    const G& gas,
                    const ThermoFixture<N, K, M, C>& fixture,
                    const UnitScales<P, U>& scales,
                    std::tuple<Pair...>* /*pairs*/)
{
    const auto pair = [&]<typename A, typename B>(std::pair<A, B>*) {
        for (std::size_t s = 0; s < C; ++s) {
            for_each_seed<T>([&](Seeds<T>& seeds) {
                const auto& state = fixture.states[s];
                const auto& d = fixture.modal[s];
                const auto composition = composition_input<F, S>(seeds, state.y, d, scales);
                const A a = pair_input<A>(seeds, state.temp, d, scales);
                const B b = pair_input<B>(seeds, state.temp, d, scales);
                const auto modes = mode_input<S, MK>(seeds, state, d, scales);
                evaluate_expecting<&modal_closure_check<F, S, MF, MK, N, G, A, B>>(
                  report,
                  state_label(fixture.name, s, seeds.label()),
                  ClosureExpectation {},
                  gas,
                  composition,
                  a,
                  b,
                  modes);
            });
        }
    };
    (pair(static_cast<Pair*>(nullptr)), ...);
}

// Mixtures of the selected cell from mass fractions and, with Densities, from
// partial densities and the density.
template <std::size_t N,
          typename G,
          typename T,
          bool Densities,
          std::size_t K,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void mixtures(Report& report,
              const G& gas,
              const ThermoFixture<N, K, M, C>& fixture,
              const UnitScales<P, U>& scales)
{
    for (std::size_t s = 0; s < C; ++s) {
        for_each_seed<T>([&](Seeds<T>& seeds) {
            const auto composition =
              composition_input<Composition::Fractions, Storage::Contiguous>(seeds,
                                                                             fixture.states[s].y,
                                                                             fixture.equilibrium[s],
                                                                             scales);
            evaluate_expecting<&fractions_mixture_check<Storage::Contiguous, N, G, T>>(
              report,
              state_label(fixture.name, s, seeds.label()),
              CompositionExpectation {},
              gas,
              composition);
        });
        if constexpr (Densities) {
            for_each_seed<T>([&](Seeds<T>& seeds) {
                const auto& d = fixture.equilibrium[s];
                const auto composition =
                  composition_input<Composition::PartialDensities, Storage::Contiguous>(
                    seeds,
                    fixture.states[s].y,
                    d,
                    scales);
                const auto rho = seeds.template quantity<kind::density>("rho", d.rho, scales);
                evaluate_expecting<&densities_mixture_check<Storage::Contiguous, N, G, T>>(
                  report,
                  state_label(fixture.name, s, seeds.label()),
                  CompositionExpectation {},
                  gas,
                  composition,
                  rho);
            });
        }
    }
}

// Closures and mixtures of a mixture model of N species, with partial-density
// input when Densities.
template <std::size_t N,
          typename G,
          typename T,
          bool Densities,
          std::size_t K,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void mixture_model(Report& report,
                   const G& gas,
                   const ThermoFixture<N, K, M, C>& fixture,
                   const UnitScales<P, U>& scales)
{
    DensityTemperature<T, U>* const rho_t = nullptr;
    mixture_closures<Composition::Fractions, Storage::Contiguous, false, N, G, T>(
      report,
      gas,
      fixture,
      scales,
      static_cast<ClosurePairs<T, U>*>(nullptr));
    mixture_closures<Composition::Fractions, Storage::LayoutRight, false, N, G, T>(report,
                                                                                   gas,
                                                                                   fixture,
                                                                                   scales,
                                                                                   rho_t);
    mixture_closures<Composition::Fractions, Storage::LayoutLeft, false, N, G, T>(report,
                                                                                  gas,
                                                                                  fixture,
                                                                                  scales,
                                                                                  rho_t);
    mixture_closures<Composition::Mixture, Storage::Contiguous, false, N, G, T>(report,
                                                                                gas,
                                                                                fixture,
                                                                                scales,
                                                                                rho_t);
    if constexpr (Densities) {
        mixture_closures<Composition::PartialDensities, Storage::Contiguous, false, N, G, T>(
          report,
          gas,
          fixture,
          scales,
          rho_t);
    }
    mixtures<N, G, T, Densities>(report, gas, fixture, scales);
}

// Modal closures with mode inputs of kind MK.
template <KindTag MK,
          std::size_t N,
          typename G,
          typename T,
          bool Densities,
          std::size_t K,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void modal_kind(Report& report,
                const G& gas,
                const ThermoFixture<N, K, M, C>& fixture,
                const UnitScales<P, U>& scales)
{
    DensityTemperature<T, U>* const rho_t = nullptr;
    modal_closures<Composition::Fractions, Storage::Contiguous, ModeForm::Quantities, MK, N, G, T>(
      report,
      gas,
      fixture,
      scales,
      static_cast<ClosurePairs<T, U>*>(nullptr));
    modal_closures<Composition::Fractions, Storage::LayoutRight, ModeForm::Indexed, MK, N, G, T>(
      report,
      gas,
      fixture,
      scales,
      rho_t);
    modal_closures<Composition::Fractions, Storage::LayoutLeft, ModeForm::Quantities, MK, N, G, T>(
      report,
      gas,
      fixture,
      scales,
      rho_t);
    modal_closures<Composition::Mixture, Storage::Contiguous, ModeForm::Indexed, MK, N, G, T>(
      report,
      gas,
      fixture,
      scales,
      rho_t);
    if constexpr (Densities) {
        modal_closures<Composition::PartialDensities,
                       Storage::Contiguous,
                       ModeForm::Quantities,
                       MK,
                       N,
                       G,
                       T>(report, gas, fixture, scales, rho_t);
    }
}

// Properties of the owning state of every fixture state, with the mixture of
// its mass fractions.
template <std::size_t N,
          typename G,
          typename T,
          std::size_t K,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void mixture_properties(Report& report,
                        const G& gas,
                        const ThermoFixture<N, K, M, C>& fixture,
                        const UnitScales<P, U>& scales)
{
    using State = state_t<G, T>;
    for (std::size_t s = 0; s < C; ++s) {
        for_each_seed<T>([&](Seeds<T>& seeds) {
            const auto composition =
              cell_storage<Storage::Contiguous>(seeds.fractions("y", fixture.states[s].y));
            const State state =
              owning_state<State>(seeds, fixture.states[s], fixture.modal[s], scales);
            evaluate<&mixture_properties_check<N, G, State>>(
              report,
              state_label(fixture.name, s, seeds.label()),
              gas,
              composition,
              state);
        });
    }
}

// Species operations at the temperature of every fixture state.
template <typename G,
          typename T,
          std::size_t N,
          std::size_t K,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void species_operations(Report& report,
                        const G& gas,
                        const ThermoFixture<N, K, M, C>& fixture,
                        const UnitScales<P, U>& scales)
{
    for (std::size_t s = 0; s < C; ++s) {
        for_each_seed<T>([&](Seeds<T>& seeds) {
            const auto temp =
              seeds.template quantity<kind::temperature>("temp", fixture.states[s].temp, scales);
            const std::string label = state_label(fixture.name, s, seeds.label());
            evaluate<&species_properties_check<G, T>>(report, label, gas, temp);
            evaluate<&species_isobaric_specific_heat_check<G, T>>(report, label, gas, temp);
        });
    }
}

// Species energies of every mode at the mode temperatures of every fixture
// state.
template <typename G,
          typename T,
          std::size_t N,
          std::size_t K,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void species_mode_energies(Report& report,
                           const G& gas,
                           const ThermoFixture<N, K, M, C>& fixture,
                           const UnitScales<P, U>& scales)
{
    for (std::size_t s = 0; s < C; ++s) {
        for (const ModeIndex m : mode_indices(0, K)) {
            for_each_seed<T>([&](Seeds<T>& seeds) {
                const auto temp =
                  seeds.template quantity<kind::temperature>("mode_temp",
                                                             fixture.states[s].modes[m.value()],
                                                             scales);
                evaluate<&species_mode_energy_check<G, T>>(
                  report,
                  state_label(fixture.name, s, seeds.label()) + ", mode " +
                    std::to_string(m.value()),
                  gas,
                  m,
                  temp);
            });
        }
    }
}

template <std::size_t N,
          typename G,
          typename T,
          bool Densities,
          std::size_t K,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void modal_model(Report& report,
                 const G& gas,
                 const ThermoFixture<N, K, M, C>& fixture,
                 const UnitScales<P, U>& scales)
{
    mixture_model<N, G, T, Densities>(report, gas, fixture, scales);
    modal_kind<kind::mode_temperature, N, G, T, Densities>(report, gas, fixture, scales);
    modal_kind<kind::mode_energy_density, N, G, T, Densities>(report, gas, fixture, scales);
    mixture_properties<N, G, T>(report, gas, fixture, scales);
    species_operations<G, T>(report, gas, fixture, scales);
    species_mode_energies<G, T>(report, gas, fixture, scales);
}

template <std::size_t N,
          typename G,
          typename T,
          bool Densities,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void thermally_perfect(Report& report,
                       const G& gas,
                       const ThermoFixture<N, 0, M, C>& fixture,
                       const UnitScales<P, U>& scales)
{
    mixture_model<N, G, T, Densities>(report, gas, fixture, scales);
    mixture_closures<Composition::Fractions, Storage::Contiguous, true, N, G, T>(
      report,
      gas,
      fixture,
      scales,
      static_cast<EnergyPairs<T, U>*>(nullptr));
    mixture_properties<N, G, T>(report, gas, fixture, scales);
    species_operations<G, T>(report, gas, fixture, scales);
}

// Energy closures of a thermally perfect model, with and without a guess, for
// a parameter type P that differs in precision from the state scalar T.
template <std::size_t N,
          typename G,
          typename T,
          std::size_t M,
          std::size_t C,
          typename P,
          typename U>
void energy_closures(Report& report,
                     const G& gas,
                     const ThermoFixture<N, 0, M, C>& fixture,
                     const UnitScales<P, U>& scales)
{
    EnergyPairs<T, U>* const pairs = nullptr;
    mixture_closures<Composition::Fractions, Storage::Contiguous, false, N, G, T>(report,
                                                                                  gas,
                                                                                  fixture,
                                                                                  scales,
                                                                                  pairs);
    mixture_closures<Composition::Fractions, Storage::Contiguous, true, N, G, T>(report,
                                                                                 gas,
                                                                                 fixture,
                                                                                 scales,
                                                                                 pairs);
}

// The thermally perfect fixtures: controlled one-segment and named nitrogen,
// and five-species air at the air and mixed compositions.
inline ThermoFixture<1, 0, 0, 3> controlled_nitrogen_fixture()
{
    return single_fixture("controlled nitrogen",
                          controlled_nitrogen_states,
                          controlled_nitrogen_derived);
}
inline ThermoFixture<1, 0, 0, 10> named_nitrogen_fixture()
{ return single_fixture("named nitrogen", named_states<1>({1.0}), named_nitrogen_derived); }
inline ThermoFixture<5, 0, 0, 10> air5_air_fixture()
{ return single_fixture("air5 air", named_states(air5_air), air5_air_derived); }
inline ThermoFixture<5, 0, 0, 10> air5_mixed_fixture()
{ return single_fixture("air5 mixed", named_states(air5_mixed), air5_mixed_derived); }

template <std::floating_point P, typename T, UnitSystem U> Report thermodynamic_mixed_checks()
{
    Report report;
    const auto scales = fixture_scales<P, U>();
    energy_closures<1, ThermallyPerfectMixture<P, U, 1, 1>, T>(
      report,
      make_thermally_perfect_mixture(controlled_nitrogen, scales).model,
      controlled_nitrogen_fixture(),
      scales);
    energy_closures<1, ThermallyPerfectMixture<P, U, 1, 4>, T>(
      report,
      make_thermally_perfect_mixture(named_nitrogen, scales).model,
      named_nitrogen_fixture(),
      scales);
    const auto air5 = make_air5(scales).model;
    energy_closures<5, ThermallyPerfectMixture<P, U, 5, 4>, T>(report,
                                                               air5,
                                                               air5_air_fixture(),
                                                               scales);
    energy_closures<5, ThermallyPerfectMixture<P, U, 5, 4>, T>(report,
                                                               air5,
                                                               air5_mixed_fixture(),
                                                               scales);
    return report;
}

// NASA9 cp/R, h/R and s/R of curve at each temperature value t of the
// evaluation units, seeded by t. Within 1e-6 t of 298.15 K the
// limits of agreement of cp/R and h/R gain the absolute part 1e-9 t: there
// h/R of a species with zero enthalpy of formation is a rounding residue of
// terms of magnitude t.
template <typename T, typename Curve, typename P, typename U>
void nasa9_curve(Report& report,
                 const std::string& name,
                 const Curve& curve,
                 const std::vector<typename ScalarTraits<T>::value_type>& temperatures,
                 const UnitScales<P, U>& scales)
{
    using V = ScalarTraits<T>::value_type;
    const V reference = fixed_quantity<kind::temperature, V>(298.15, scales).value();
    for (const V value : temperatures) {
        Agreement limits {};
        if (std::abs(value - reference) <= V(1e-6) * value) {
            limits.absolute = 1e-9 * static_cast<double>(value);
        }
        for_each_seed<T>([&](Seeds<T>& seeds) {
            const T t = seeds.converted("t", value);
            const std::string label =
              name + " at t = " + std::to_string(value) + ", " + seeds.label();
            evaluate_within<&nasa9_cp_h_check<Curve, T>>(report, label, limits, curve, t);
            evaluate<&nasa9_s_check<Curve, T>>(report, label, curve, t);
        });
    }
}

// The temperatures temps_k [K] in the evaluation units.
template <typename V, typename P, typename U, std::size_t C>
std::vector<V> converted_temperatures(const std::array<double, C>& temps_k,
                                      const UnitScales<P, U>& scales)
{
    std::vector<V> out;
    for (const double t : temps_k) {
        out.push_back(fixed_quantity<kind::temperature, V>(t, scales).value());
    }
    return out;
}

// The named-curve temperatures of fit: half the first break, the first and
// last breaks, 1.1 times the last break, each segment midpoint, and the
// lower, middle and upper temperature of each transition with the adjacent
// representable values of each, in the evaluation units without duplicates.
template <typename V, std::size_t S, typename P, typename U>
std::vector<V> named_curve_temperatures(const Nasa9Data<S>& fit, const UnitScales<P, U>& scales)
{
    const auto n = static_cast<std::size_t>(fit.segments);
    std::vector<double> kelvin {0.5 * fit.breaks[0],
                                fit.breaks[0],
                                fit.breaks[n],
                                1.1 * fit.breaks[n]};
    for (std::size_t j = 0; j < n; ++j) {
        kelvin.push_back(0.5 * (fit.breaks[j] + fit.breaks[j + 1]));
    }
    std::vector<V> out;
    for (const double t : kelvin) {
        out.push_back(fixed_quantity<kind::temperature, V>(t, scales).value());
    }
    for (std::size_t i = 1; i < n; ++i) {
        const double width = fit.blends[i - 1];
        if (width <= 0.0) { continue; }
        for (const double t :
             {fit.breaks[i] - (0.5 * width), fit.breaks[i], fit.breaks[i] + (0.5 * width)}) {
            const V value = fixed_quantity<kind::temperature, V>(t, scales).value();
            out.push_back(std::nextafter(value, -std::numeric_limits<V>::infinity()));
            out.push_back(value);
            out.push_back(std::nextafter(value, std::numeric_limits<V>::infinity()));
        }
    }
    std::ranges::sort(out);
    const auto [first, last] = std::ranges::unique(out);
    out.erase(first, last);
    return out;
}

template <std::size_t S> bool same_fit(const Nasa9Data<S>& a, const Nasa9Data<S>& b)
{ return a.segments == b.segments && a.breaks == b.breaks && a.blends == b.blends && a.a == b.a; }

// The temperature unit of U in kelvin, the temperature scale of a direct curve.
template <std::floating_point P, UnitSystem U> P curve_temperature_scale()
{
    if constexpr (std::same_as<U, si>) {
        return P(1);
    } else {
        return P(1000);
    }
}

// The thermodynamic checks for parameter type P, state scalar T and units U, in
// two groups of similar size so that each compiles as its own unit: models with
// one temperature and NASA9 curves, then two- and multi-temperature models.
template <std::floating_point P, typename T, UnitSystem U> Report thermodynamic_single_checks()
{
    using V = ScalarTraits<T>::value_type;
    Report report;
    const auto scales = fixture_scales<P, U>();
    // Perfect gas.
    {
        using Perfect = CaloricallyPerfect<P, U>;
        const Perfect gas = make_perfect_gas(scales);
        const auto fixture = single_fixture("perfect gas", perfect_states, perfect_derived);
        closures<Perfect, T>(report,
                             gas,
                             fixture,
                             scales,
                             static_cast<ClosurePairs<T, U>*>(nullptr));
        for (std::size_t s = 0; s < fixture.states.size(); ++s) {
            for_each_seed<T>([&](Seeds<T>& seeds) {
                const auto state = owning_state<ThermoState<T, U>>(seeds,
                                                                   fixture.states[s],
                                                                   fixture.modal[s],
                                                                   scales);
                evaluate<&properties_check<Perfect, ThermoState<T, U>>>(
                  report,
                  state_label(fixture.name, s, seeds.label()),
                  gas,
                  state);
            });
        }
    }
    // Perfect mixtures: 5 and 11 species, at and below capacity.
    {
        const auto five = single_fixture("perfect mixture 5",
                                         perfect_mixture_states<5>(),
                                         perfect_mixture5_derived);
        const auto eleven = single_fixture("perfect mixture 11",
                                           perfect_mixture_states<11>(),
                                           perfect_mixture11_derived);
        const auto at5 = make_perfect_mixture<5>(5, scales);
        const auto below5 = make_perfect_mixture<16>(5, scales);
        const auto at11 = make_perfect_mixture<11>(11, scales);
        const auto below11 = make_perfect_mixture<16>(11, scales);
        mixture_model<5, CaloricallyPerfectMixture<P, U, 5>, T, true>(report, at5, five, scales);
        mixture_model<5, CaloricallyPerfectMixture<P, U, 16>, T, true>(report,
                                                                       below5,
                                                                       five,
                                                                       scales);
        mixture_model<11, CaloricallyPerfectMixture<P, U, 11>, T, true>(report,
                                                                        at11,
                                                                        eleven,
                                                                        scales);
        mixture_model<11, CaloricallyPerfectMixture<P, U, 16>, T, true>(report,
                                                                        below11,
                                                                        eleven,
                                                                        scales);
        mixture_properties<5, CaloricallyPerfectMixture<P, U, 5>, T>(report, at5, five, scales);
        mixture_properties<11, CaloricallyPerfectMixture<P, U, 16>, T>(report,
                                                                       below11,
                                                                       eleven,
                                                                       scales);
    }
    // Thermally perfect: controlled one-segment, named N2 and five-species air.
    thermally_perfect<1, ThermallyPerfectMixture<P, U, 1, 1>, T, true>(
      report,
      make_thermally_perfect_mixture(controlled_nitrogen, scales).model,
      controlled_nitrogen_fixture(),
      scales);
    thermally_perfect<1, ThermallyPerfectMixture<P, U, 1, 4>, T, false>(
      report,
      make_thermally_perfect_mixture(named_nitrogen, scales).model,
      named_nitrogen_fixture(),
      scales);
    {
        const auto air5 = make_air5(scales).model;
        thermally_perfect<5, ThermallyPerfectMixture<P, U, 5, 4>, T, true>(report,
                                                                           air5,
                                                                           air5_air_fixture(),
                                                                           scales);
        thermally_perfect<5, ThermallyPerfectMixture<P, U, 5, 4>, T, true>(report,
                                                                           air5,
                                                                           air5_mixed_fixture(),
                                                                           scales);
    }
    // NASA9 curves of one, two and four segments: the constant curve, the
    // seven term-isolated curves and every air11 curve (and each air5 curve
    // that differs from the air11 curve of its name).
    const P scale = curve_temperature_scale<P, U>();
    nasa9_curve<T>(report,
                   "constant curve",
                   make_nasa9_curve(constant_fit(4.0), scale),
                   converted_temperatures<V>(constant_curve_temperatures, scales),
                   scales);
    for (std::size_t k = 0; k < 7; ++k) {
        nasa9_curve<T>(report,
                       "term-isolated curve " + std::to_string(k),
                       make_nasa9_curve(term_isolated_fit(k), scale),
                       converted_temperatures<V>(term_isolated_temperatures, scales),
                       scales);
    }
    for (const auto& species : air11_species) {
        nasa9_curve<T>(report,
                       std::string("air11 ") + species.name,
                       make_nasa9_curve(species.fit, scale),
                       named_curve_temperatures<V>(species.fit, scales),
                       scales);
    }
    for (const auto& species : air5_species) {
        for (const auto& named : air11_species) {
            if (std::string_view(named.name) == species.name && !same_fit(named.fit, species.fit)) {
                nasa9_curve<T>(report,
                               std::string("air5 ") + species.name,
                               make_nasa9_curve(species.fit, scale),
                               named_curve_temperatures<V>(species.fit, scales),
                               scales);
            }
        }
    }
    return report;
}

template <std::floating_point P, typename T, UnitSystem U> Report thermodynamic_modal_checks()
{
    Report report;
    const auto scales = fixture_scales<P, U>();
    // Two-temperature air.
    {
        using TwoTemperature = TwoTemperatureMixture<P, U, 5, 4>;
        const auto gas = make_air5_two_temperature(scales).model;
        const ThermoFixture<5, 1, 1, 6> fixture {"two-temperature air5",
                                                 two_temperature_states,
                                                 two_temperature_equilibrium_derived,
                                                 two_temperature_modal_derived};
        modal_model<5, TwoTemperature, T, true>(report, gas, fixture, scales);
        for (std::size_t s = 0; s < fixture.states.size(); ++s) {
            for_each_seed<T>([&](Seeds<T>& seeds) {
                const auto temp_ve =
                  seeds.template quantity<kind::temperature>("temp_ve",
                                                             fixture.states[s].modes[0],
                                                             scales);
                evaluate<&species_vibrational_electronic_energy_check<TwoTemperature, T>>(
                  report,
                  state_label(fixture.name, s, seeds.label()),
                  gas,
                  temp_ve);
            });
        }
    }
    // Multi-temperature: four-temperature air5, three-temperature air11 and the
    // controlled two-species model.
    modal_model<5, MultiTemperatureMixture<P, U, 5, 3, 4, 8>, T, false>(
      report,
      make_air5_four_temperature(scales).model,
      ThermoFixture<5, 3, 3, 4> {"four-temperature air5",
                                 four_temperature_states,
                                 four_temperature_equilibrium_derived,
                                 four_temperature_modal_derived},
      scales);
    {
        using ThreeTemperature = MultiTemperatureMixture<P, U, 11, 2, 4, 8>;
        const auto gas = make_air11_three_temperature(scales).model;
        const ThermoFixture<11, 2, 2, 6> fixture {"three-temperature air11",
                                                  three_temperature_states,
                                                  three_temperature_equilibrium_derived,
                                                  three_temperature_modal_derived};
        modal_model<11, ThreeTemperature, T, true>(report, gas, fixture, scales);
        using State = state_t<ThreeTemperature, T>;
        for (std::size_t s = 0; s < fixture.states.size(); ++s) {
            for_each_seed<T>([&](Seeds<T>& seeds) {
                const auto composition =
                  cell_storage<Storage::Contiguous>(seeds.fractions("y", fixture.states[s].y));
                const State state =
                  owning_state<State>(seeds, fixture.states[s], fixture.modal[s], scales);
                evaluate<&electron_pressure_check<11, ThreeTemperature, State>>(
                  report,
                  state_label(fixture.name, s, seeds.label()),
                  gas,
                  composition,
                  state);
            });
        }
    }
    {
        using Controlled = MultiTemperatureMixture<P, U, 2, 2, 1, 2>;
        modal_model<2, Controlled, T, false>(
          report,
          make_controlled_modal(controlled_two_levels, scales),
          ThermoFixture<2, 2, 2, 4> {"controlled modal",
                                     controlled_modal_states,
                                     controlled_modal_equilibrium_derived,
                                     controlled_modal_derived},
          scales);
        // The N ground level alone: forward closure, mode energies and properties.
        const auto ground = make_controlled_modal(controlled_ground_level, scales);
        const ThermoFixture<2, 2, 2, 2> fixture {"controlled ground level",
                                                 controlled_ground_level_states,
                                                 controlled_ground_level_derived,
                                                 controlled_ground_level_derived};
        modal_closures<Composition::Fractions,
                       Storage::Contiguous,
                       ModeForm::Quantities,
                       kind::mode_temperature,
                       2,
                       Controlled,
                       T>(report,
                          ground,
                          fixture,
                          scales,
                          static_cast<DensityTemperature<T, U>*>(nullptr));
        species_mode_energies<Controlled, T>(report, ground, fixture, scales);
        mixture_properties<2, Controlled, T>(report, ground, fixture, scales);
    }
    return report;
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_THERMODYNAMICS_CASES_HPP
