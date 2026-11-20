// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_TESTS_DEVICE_EXCHANGE_CASES_HPP
#define CANARY_TESTS_DEVICE_EXCHANGE_CASES_HPP

// Energy-exchange checks for the host and native compilations (see
// cell_inputs.hpp). Standalone chemistry-coupled terms receive their rates of
// progress as input; the combined sources evaluate their own.

#include "cell_inputs.hpp"
#include "check_groups.hpp"
#include "derived_inputs.hpp"
#include "exchange_inputs.hpp"
#include "host_inputs.hpp"
#include "native_runtime.hpp"

#include <array>
#include <canary/energy_exchange.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/state.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace canary::test
{

// Entries: see thermodynamics_cases.hpp.

// Two-temperature terms at (rho, T, T_ve) of the selected cell.
template <std::size_t N, typename X, typename G, typename T>
CANARY_HOST_DEVICE void
  two_temperature_vibration_translation_check(Output<PowerDensity<T, typename G::units>> out,
                                              const X* exchange,
                                              const G* gas,
                                              const T* composition,
                                              Density<T, typename G::units> rho,
                                              Temperature<T, typename G::units> temp,
                                              Temperature<T, typename G::units> temp_ve)
{
    out[0] = vibration_translation_source(
      *exchange,
      *gas,
      cell_composition<Composition::Mixture, Storage::Contiguous, N>(*gas, composition),
      rho,
      temp,
      temp_ve);
}

template <typename X, typename M, typename T>
CANARY_HOST_DEVICE void two_temperature_chemistry_vibration_check(
  Output<PowerDensity<T, typename X::units>> out,
  const X* exchange,
  const M* mechanism,
  const RatesOfProgress<T, typename X::units, M::reaction_count>* rates,
  Temperature<T, typename X::units> temp,
  Temperature<T, typename X::units> temp_ve)
{ out[0] = chemistry_vibration_source(*exchange, *mechanism, *rates, temp, temp_ve); }

template <std::size_t N, typename M, typename X, typename G, typename T>
CANARY_HOST_DEVICE void
  two_temperature_sources_check(Output<TwoTemperatureSources<T, typename G::units, N>> out,
                                const M* mechanism,
                                const X* exchange,
                                const G* gas,
                                const T* composition,
                                Density<T, typename G::units> rho,
                                Temperature<T, typename G::units> temp,
                                Temperature<T, typename G::units> temp_ve)
{
    out[0] = two_temperature_sources(
      *mechanism,
      *exchange,
      *gas,
      cell_composition<Composition::Mixture, Storage::Contiguous, N>(*gas, composition),
      rho,
      temp,
      temp_ve);
}

// Multi-temperature terms of the selected cell, with the state in form F (see
// cell_state). The mode-energy and mode-temperature grids are read only by the
// borrowed form.
template <typename G, typename T>
using ModeSources = std::array<PowerDensity<T, typename G::units>, G::mode_count>;
template <typename G, typename T>
using CellState = MultiTemperatureState<T, typename G::units, G::mode_count>;

template <StateForm F, std::size_t N, typename X, typename G, typename T>
CANARY_HOST_DEVICE void vibration_translation_check(Output<ModeSources<G, T>> out,
                                                    const X* exchange,
                                                    const G* gas,
                                                    const T* composition,
                                                    CellState<G, T> state,
                                                    const T* mode_energies,
                                                    const T* mode_temperatures)
{
    out[0] = vibration_translation_source(
      *exchange,
      *gas,
      cell_composition<Composition::Mixture, Storage::Contiguous, N>(*gas, composition),
      cell_state<F>(state, mode_energies, mode_temperatures));
}

template <StateForm F, typename X, typename M, typename G, typename T>
CANARY_HOST_DEVICE void
  chemistry_vibration_check(Output<ModeSources<G, T>> out,
                            const X* exchange,
                            const M* mechanism,
                            const RatesOfProgress<T, typename G::units, M::reaction_count>* rates,
                            CellState<G, T> state,
                            const T* mode_energies,
                            const T* mode_temperatures)
{
    out[0] = chemistry_vibration_source(*exchange,
                                        *mechanism,
                                        *rates,
                                        cell_state<F>(state, mode_energies, mode_temperatures));
}

template <StateForm F, std::size_t N, typename X, typename G, typename T>
CANARY_HOST_DEVICE void electron_translation_check(Output<ModeSources<G, T>> out,
                                                   const X* exchange,
                                                   const G* gas,
                                                   const T* composition,
                                                   CellState<G, T> state,
                                                   const T* mode_energies,
                                                   const T* mode_temperatures)
{
    out[0] = electron_translation_source(
      *exchange,
      cell_composition<Composition::Mixture, Storage::Contiguous, N>(*gas, composition),
      cell_state<F>(state, mode_energies, mode_temperatures));
}

template <StateForm F, std::size_t N, typename X, typename G, typename T>
CANARY_HOST_DEVICE void vibration_electron_check(Output<ModeSources<G, T>> out,
                                                 const X* exchange,
                                                 const G* gas,
                                                 const T* composition,
                                                 CellState<G, T> state,
                                                 const T* mode_energies,
                                                 const T* mode_temperatures)
{
    out[0] = vibration_electron_source(
      *exchange,
      *gas,
      cell_composition<Composition::Mixture, Storage::Contiguous, N>(*gas, composition),
      cell_state<F>(state, mode_energies, mode_temperatures));
}

template <typename X, typename M, typename G, typename T>
CANARY_HOST_DEVICE void electron_impact_ionisation_check(
  Output<ModeSources<G, T>> out,
  const X* exchange,
  const RatesOfProgress<T, typename G::units, M::reaction_count>* rates)
{ out[0] = electron_impact_ionisation_source(*exchange, *rates); }

template <StateForm F, std::size_t N, typename M, typename X, typename G, typename T>
CANARY_HOST_DEVICE void multi_temperature_sources_check(
  Output<MultiTemperatureSources<T, typename G::units, N, G::mode_count>> out,
  const M* mechanism,
  const X* exchange,
  const G* gas,
  const T* composition,
  CellState<G, T> state,
  const T* mode_energies,
  const T* mode_temperatures)
{
    out[0] = multi_temperature_sources(
      *mechanism,
      *exchange,
      *gas,
      cell_composition<Composition::Mixture, Storage::Contiguous, N>(*gas, composition),
      cell_state<F>(state, mode_energies, mode_temperatures));
}

// Runners. Each evaluates the entries of one exchange at every state of a
// fixture and every seed (see host_inputs.hpp). A standalone term of a
// process the exchange does not configure is exactly zero, value and
// derivative; every result is compared with the host result.

// The process counts of an exchange type.
template <typename X> struct ExchangeCounts;
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t V, std::size_t C>
struct ExchangeCounts<EnergyExchange<P, U, N, V, C>> {
    static constexpr std::size_t vibration_translation = V;
    static constexpr std::size_t chemistry_vibration = C;
};
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t M,
          std::size_t V,
          std::size_t C,
          std::size_t E,
          std::size_t VE,
          std::size_t I>
struct ExchangeCounts<MultiTemperatureExchange<P, U, N, M, V, C, E, VE, I>> {
    static constexpr std::size_t vibration_translation = V;
    static constexpr std::size_t chemistry_vibration = C;
    static constexpr std::size_t electron_translation = E;
    static constexpr std::size_t vibration_electron = VE;
    static constexpr std::size_t electron_impact_ionisation = I;
};

// Evaluates Entry and checks its result for exact zeros when Zero.
template <auto Entry, bool Zero, typename... Args>
void evaluate_term(Report& report, const std::string& fixture, const Args&... args)
{
    if constexpr (Zero) {
        evaluate_expecting<Entry>(report, fixture, ZeroExpectation {}, args...);
    } else {
        evaluate<Entry>(report, fixture, args...);
    }
}

// The rates of progress rates as an input, every rate one input.
template <typename T, std::size_t R, typename P, typename U>
RatesOfProgress<T, U, R>
  rates_input(Seeds<T>& seeds, const DerivedRates<R>& rates, const UnitScales<P, U>& scales)
{
    const auto forward =
      seeds.template values<kind::number_rate_density>("q_f", rates.forward, scales);
    const auto backward =
      seeds.template values<kind::number_rate_density>("q_b", rates.backward, scales);
    RatesOfProgress<T, U, R> out {};
    for (const ReactionIndex k : reaction_indices(0, R)) {
        out.forward()[k] = NumberRateDensity<T, U> {forward[k.value()]};
        out.backward()[k] = NumberRateDensity<T, U> {backward[k.value()]};
    }
    return out;
}

// An exchange fixture: declared states, the inputs derived from them, and the
// rates of progress supplied to the standalone chemistry-coupled terms at each
// state.
template <std::size_t N, std::size_t K, std::size_t R, std::size_t C> struct ExchangeFixture {
    std::string name;
    std::array<FixtureState<N, K>, C> states;
    std::array<DerivedState<K>, C> derived;
    std::vector<std::vector<DerivedRates<R>>> rates;
};

// The rates of every state: the derived rates of each state, or the same
// independent rate cases at every state.
template <std::size_t R, std::size_t C>
std::vector<std::vector<DerivedRates<R>>> state_rates(const std::array<DerivedRates<R>, C>& derived)
{
    std::vector<std::vector<DerivedRates<R>>> out;
    for (const DerivedRates<R>& rates : derived) { out.push_back({rates}); }
    return out;
}

template <std::size_t R, std::size_t K>
std::vector<std::vector<DerivedRates<R>>> independent_rates(
  std::size_t states, std::size_t reaction, const std::array<std::array<double, 2>, K>& cases)
{
    std::vector<DerivedRates<R>> rates;
    for (const auto& c : cases) { rates.push_back(isolated_rates<R>(reaction, c)); }
    return std::vector<std::vector<DerivedRates<R>>>(states, rates);
}

// The two-temperature terms of exchange at every state of fixture.
template <std::size_t N,
          typename M,
          typename X,
          typename G,
          typename T,
          std::size_t R,
          std::size_t C,
          typename P,
          typename U>
void two_temperature_terms(Report& report,
                           const M& mechanism,
                           const X& exchange,
                           const G& gas,
                           const ExchangeFixture<N, 1, R, C>& fixture,
                           const UnitScales<P, U>& scales)
{
    using Counts = ExchangeCounts<X>;
    for (std::size_t s = 0; s < C; ++s) {
        const auto& state = fixture.states[s];
        // The vibration-translation source is proportional to the difference of
        // the vibrational energies at T and T_ve. When the temperatures differ
        // by at most 1e-6 T, rounding differences are amplified by at least
        // 1e6, and the relative limit of agreement is 1e-5.
        Agreement vibration_translation_limits {};
        if (std::abs(state.temp - state.modes[0]) <= 1e-6 * state.temp) {
            vibration_translation_limits.relative = 1e-5;
        }
        for_each_seed<T>([&](Seeds<T>& seeds) {
            const auto composition =
              cell_storage<Storage::Contiguous>(seeds.fractions("y", state.y));
            const auto rho = seeds.template quantity<kind::density>("rho", state.rho, scales);
            const auto temp =
              seeds.template quantity<kind::temperature>("temp", state.temp, scales);
            const auto temp_ve =
              seeds.template quantity<kind::temperature>("temp_ve", state.modes[0], scales);
            const std::string label = state_label(fixture.name, s, seeds.label());
            if constexpr (Counts::vibration_translation == 0) {
                evaluate_expecting<&two_temperature_vibration_translation_check<N, X, G, T>>(
                  report,
                  label,
                  ZeroExpectation {},
                  exchange,
                  gas,
                  composition,
                  rho,
                  temp,
                  temp_ve);
            } else {
                evaluate_within<&two_temperature_vibration_translation_check<N, X, G, T>>(
                  report,
                  label,
                  vibration_translation_limits,
                  exchange,
                  gas,
                  composition,
                  rho,
                  temp,
                  temp_ve);
            }
            evaluate<&two_temperature_sources_check<N, M, X, G, T>>(report,
                                                                    label,
                                                                    mechanism,
                                                                    exchange,
                                                                    gas,
                                                                    composition,
                                                                    rho,
                                                                    temp,
                                                                    temp_ve);
        });
        for (std::size_t r = 0; r < fixture.rates[s].size(); ++r) {
            for_each_seed<T>([&](Seeds<T>& seeds) {
                const auto rates = rates_input<T>(seeds, fixture.rates[s][r], scales);
                const auto temp =
                  seeds.template quantity<kind::temperature>("temp", state.temp, scales);
                const auto temp_ve =
                  seeds.template quantity<kind::temperature>("temp_ve", state.modes[0], scales);
                evaluate_term<&two_temperature_chemistry_vibration_check<X, M, T>,
                              Counts::chemistry_vibration == 0>(
                  report,
                  state_label(fixture.name, s, seeds.label()) + ", rates " + std::to_string(r),
                  exchange,
                  mechanism,
                  rates,
                  temp,
                  temp_ve);
            });
        }
    }
}

// The multi-temperature terms of exchange at every state of fixture, with the
// state in form F; electron-impact ionisation with the owning form only.
template <StateForm F,
          std::size_t N,
          typename M,
          typename X,
          typename G,
          typename T,
          std::size_t K,
          std::size_t R,
          std::size_t C,
          typename P,
          typename U>
void multi_temperature_terms(Report& report,
                             const M& mechanism,
                             const X& exchange,
                             const G& gas,
                             const ExchangeFixture<N, K, R, C>& fixture,
                             const UnitScales<P, U>& scales)
{
    using Counts = ExchangeCounts<X>;
    using State = CellState<G, T>;
    // The mode energies and temperatures of state in the grids the borrowed
    // form reads.
    const auto grids = [](const State& state) {
        std::array<T, K> mode_e {};
        std::array<T, K> mode_temp {};
        for (std::size_t m = 0; m < K; ++m) {
            mode_e[m] = state.mode_e[m].value();
            mode_temp[m] = state.mode_temp[m].value();
        }
        return std::pair {cell_storage<Storage::LayoutRight>(mode_e),
                          cell_storage<Storage::LayoutLeft>(mode_temp)};
    };
    for (std::size_t s = 0; s < C; ++s) {
        const auto& declared = fixture.states[s];
        for_each_seed<T>([&](Seeds<T>& seeds) {
            const auto composition =
              cell_storage<Storage::Contiguous>(seeds.fractions("y", declared.y));
            const State state = owning_state<State>(seeds, declared, fixture.derived[s], scales);
            const auto [mode_e, mode_temp] = grids(state);
            const std::string label = state_label(fixture.name, s, seeds.label());
            evaluate_term<&vibration_translation_check<F, N, X, G, T>,
                          Counts::vibration_translation ==
                            0>(report, label, exchange, gas, composition, state, mode_e, mode_temp);
            evaluate_term<&electron_translation_check<F, N, X, G, T>,
                          Counts::electron_translation ==
                            0>(report, label, exchange, gas, composition, state, mode_e, mode_temp);
            evaluate_term<&vibration_electron_check<F, N, X, G, T>,
                          Counts::vibration_electron ==
                            0>(report, label, exchange, gas, composition, state, mode_e, mode_temp);
            evaluate<&multi_temperature_sources_check<F, N, M, X, G, T>>(report,
                                                                         label,
                                                                         mechanism,
                                                                         exchange,
                                                                         gas,
                                                                         composition,
                                                                         state,
                                                                         mode_e,
                                                                         mode_temp);
        });
        for (std::size_t r = 0; r < fixture.rates[s].size(); ++r) {
            const auto& rate_case = fixture.rates[s][r];
            for_each_seed<T>([&](Seeds<T>& seeds) {
                const auto rates = rates_input<T>(seeds, rate_case, scales);
                const State state =
                  owning_state<State>(seeds, declared, fixture.derived[s], scales);
                const auto [mode_e, mode_temp] = grids(state);
                evaluate_term<&chemistry_vibration_check<F, X, M, G, T>,
                              Counts::chemistry_vibration == 0>(
                  report,
                  state_label(fixture.name, s, seeds.label()) + ", rates " + std::to_string(r),
                  exchange,
                  mechanism,
                  rates,
                  state,
                  mode_e,
                  mode_temp);
            });
            if constexpr (F == StateForm::Owning) {
                for_each_seed<T>([&](Seeds<T>& seeds) {
                    const auto rates = rates_input<T>(seeds, rate_case, scales);
                    evaluate_term<&electron_impact_ionisation_check<X, M, G, T>,
                                  Counts::electron_impact_ionisation == 0>(
                      report,
                      state_label(fixture.name, s, seeds.label()) + ", rates " + std::to_string(r),
                      exchange,
                      rates);
                });
            }
        }
    }
}

// Energy exchange of two-temperature air5 (full, empty, each correction of
// the first vibration-translation pair alone, the first chemistry-vibration
// coupling alone), three-temperature air11 (full with owning and borrowed
// states, empty, and each process alone) and four-temperature air5 without
// free electrons.
template <std::floating_point P, typename T, UnitSystem U> Report exchange_checks()
{
    Report report;
    const auto scales = fixture_scales<P, U>();
    const auto p0 = standard_pressure(scales);
    {
        const auto built = make_air5_two_temperature(scales);
        const auto& gas = built.model;
        const auto mechanism = make_air5_park(scales, p0, gas);
        using M = std::remove_const_t<decltype(mechanism)>;
        using G = std::remove_cvref_t<decltype(gas)>;
        constexpr std::size_t R = M::reaction_count;
        const auto names = species_names(air5_species);
        const std::size_t coupled =
          static_cast<std::size_t>(air5_park_exchange.chemistry_vibration[0].reaction);
        const ExchangeFixture<5, 1, R, 8> named {"two-temperature air5 Park",
                                                 air5_park_states,
                                                 {},
                                                 state_rates(air5_park_rates_derived)};
        two_temperature_terms<5, M, EnergyExchange<P, U, 5, 15, 3>, G, T>(
          report,
          mechanism,
          make_air5_park_exchange(built.species, mechanism, scales),
          gas,
          named,
          scales);
        const auto process = [&](std::string name) {
            return ExchangeFixture<5, 1, R, 2> {std::move(name),
                                                two_temperature_process_states,
                                                {},
                                                independent_rates<R>(2, coupled, coupling_rates)};
        };
        two_temperature_terms<5, M, EnergyExchange<P, U, 5, 0, 0>, G, T>(
          report,
          mechanism,
          make_energy_exchange(no_exchange, names, built.species, mechanism, scales),
          gas,
          process("two-temperature air5, no process"),
          scales);
        for (const RelaxationCorrection correction : corrections) {
            two_temperature_terms<5, M, EnergyExchange<P, U, 5, 1, 0>, G, T>(
              report,
              mechanism,
              make_energy_exchange(first_vibration_translation(air5_park_exchange, correction),
                                   names,
                                   built.species,
                                   mechanism,
                                   scales),
              gas,
              process("two-temperature air5, vibration-translation only, correction " +
                      std::to_string(static_cast<int>(correction))),
              scales);
        }
        two_temperature_terms<5, M, EnergyExchange<P, U, 5, 0, 1>, G, T>(
          report,
          mechanism,
          make_energy_exchange(first_chemistry_vibration(air5_park_exchange),
                               names,
                               built.species,
                               mechanism,
                               scales),
          gas,
          process("two-temperature air5, chemistry-vibration only"),
          scales);
    }
    {
        const auto built = make_air11_three_temperature(scales);
        const auto& gas = built.model;
        const auto mechanism = make_air11_kim_jo(scales, p0, gas);
        using M = std::remove_const_t<decltype(mechanism)>;
        using G = std::remove_cvref_t<decltype(gas)>;
        constexpr std::size_t R = M::reaction_count;
        const auto names = species_names(air11_species);
        const auto exchange = [&](const auto& input, const auto& electrons) {
            return make_multi_temperature_exchange(input,
                                                   electrons,
                                                   names,
                                                   built.species,
                                                   mechanism,
                                                   gas,
                                                   scales);
        };
        const ExchangeFixture<11, 2, R, 4> named {"three-temperature air11 Kim-Jo",
                                                  air11_kim_jo_states,
                                                  air11_kim_jo_derived,
                                                  state_rates(air11_kim_jo_rates_derived)};
        using Full = MultiTemperatureExchange<P, U, 11, 2, 30, 21, 10, 1, 2>;
        const Full full = make_air11_kim_jo_exchange(built.species, mechanism, gas, scales);
        multi_temperature_terms<StateForm::Owning, 11, M, Full, G, T>(report,
                                                                      mechanism,
                                                                      full,
                                                                      gas,
                                                                      named,
                                                                      scales);
        multi_temperature_terms<StateForm::Borrowed, 11, M, Full, G, T>(report,
                                                                        mechanism,
                                                                        full,
                                                                        gas,
                                                                        named,
                                                                        scales);
        const std::size_t coupled =
          static_cast<std::size_t>(air11_kim_jo_exchange.chemistry_vibration[0].reaction);
        const std::size_t ionised = static_cast<std::size_t>(
          air11_kim_jo_exchange_electrons.electron_impact_ionisation[0].reaction);
        const auto fixture = [&]<std::size_t C>(std::string name,
                                                const std::array<FixtureState<11, 2>, C>& states,
                                                const std::array<DerivedState<2>, C>& derived,
                                                std::size_t reaction,
                                                const auto& cases) {
            return ExchangeFixture<11, 2, R, C> {std::move(name),
                                                 states,
                                                 derived,
                                                 independent_rates<R>(C, reaction, cases)};
        };
        const auto process = [&](std::string name) {
            return fixture(std::move(name),
                           multi_temperature_process_states,
                           multi_temperature_process_derived,
                           coupled,
                           coupling_rates);
        };
        multi_temperature_terms<StateForm::Owning,
                                11,
                                M,
                                MultiTemperatureExchange<P, U, 11, 2, 0, 0, 0, 0, 0>,
                                G,
                                T>(report,
                                   mechanism,
                                   exchange(no_exchange, no_electron_exchange),
                                   gas,
                                   process("three-temperature air11, no process"),
                                   scales);
        for (const RelaxationCorrection correction : corrections) {
            multi_temperature_terms<StateForm::Owning,
                                    11,
                                    M,
                                    MultiTemperatureExchange<P, U, 11, 2, 1, 0, 0, 0, 0>,
                                    G,
                                    T>(
              report,
              mechanism,
              exchange(first_vibration_translation(air11_kim_jo_exchange, correction),
                       no_electron_exchange),
              gas,
              process("three-temperature air11, vibration-translation only, correction " +
                      std::to_string(static_cast<int>(correction))),
              scales);
        }
        multi_temperature_terms<StateForm::Owning,
                                11,
                                M,
                                MultiTemperatureExchange<P, U, 11, 2, 0, 1, 0, 0, 0>,
                                G,
                                T>(
          report,
          mechanism,
          exchange(first_chemistry_vibration(air11_kim_jo_exchange), no_electron_exchange),
          gas,
          process("three-temperature air11, chemistry-vibration only"),
          scales);
        using ElectronTranslationOnly = MultiTemperatureExchange<P, U, 11, 2, 0, 0, 1, 0, 0>;
        multi_temperature_terms<StateForm::Owning, 11, M, ElectronTranslationOnly, G, T>(
          report,
          mechanism,
          exchange(no_exchange, air11_electron_translation(neutral_pair)),
          gas,
          fixture("three-temperature air11, neutral electron-translation only",
                  neutral_electron_states,
                  neutral_electron_derived,
                  coupled,
                  coupling_rates),
          scales);
        multi_temperature_terms<StateForm::Owning, 11, M, ElectronTranslationOnly, G, T>(
          report,
          mechanism,
          exchange(no_exchange, air11_electron_translation(charged_pair)),
          gas,
          fixture("three-temperature air11, charged electron-translation only",
                  charged_electron_states,
                  charged_electron_derived,
                  coupled,
                  coupling_rates),
          scales);
        multi_temperature_terms<StateForm::Owning,
                                11,
                                M,
                                MultiTemperatureExchange<P, U, 11, 2, 0, 0, 0, 1, 0>,
                                G,
                                T>(report,
                                   mechanism,
                                   exchange(no_exchange, air11_vibration_electron),
                                   gas,
                                   fixture("three-temperature air11, vibration-electron only",
                                           vibration_electron_states,
                                           vibration_electron_derived,
                                           coupled,
                                           coupling_rates),
                                   scales);
        multi_temperature_terms<StateForm::Owning,
                                11,
                                M,
                                MultiTemperatureExchange<P, U, 11, 2, 0, 0, 0, 0, 1>,
                                G,
                                T>(
          report,
          mechanism,
          exchange(no_exchange, air11_electron_impact_ionisation),
          gas,
          fixture("three-temperature air11, electron-impact ionisation only",
                  multi_temperature_process_states,
                  multi_temperature_process_derived,
                  ionised,
                  ionisation_rates),
          scales);
    }
    // Four-temperature air5 without a free-electron species.
    {
        const auto built = make_air5_four_temperature(scales);
        const auto& gas = built.model;
        const auto mechanism = make_air5_park(scales, p0, gas);
        using M = std::remove_const_t<decltype(mechanism)>;
        using G = std::remove_cvref_t<decltype(gas)>;
        constexpr std::size_t R = M::reaction_count;
        const std::size_t coupled =
          static_cast<std::size_t>(air5_park_exchange.chemistry_vibration[0].reaction);
        multi_temperature_terms<StateForm::Owning,
                                5,
                                M,
                                MultiTemperatureExchange<P, U, 5, 3, 15, 3, 0, 0, 0>,
                                G,
                                T>(
          report,
          mechanism,
          make_multi_temperature_exchange(air5_park_exchange,
                                          no_electron_exchange,
                                          species_names(air5_species),
                                          built.species,
                                          mechanism,
                                          gas,
                                          scales),
          gas,
          ExchangeFixture<5, 3, R, 4> {"four-temperature air5 Park",
                                       four_temperature_states,
                                       four_temperature_modal_derived,
                                       independent_rates<R>(4, coupled, coupling_rates)},
          scales);
    }
    return report;
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_EXCHANGE_CASES_HPP
