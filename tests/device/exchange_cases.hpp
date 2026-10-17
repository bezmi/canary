#ifndef CANARY_TESTS_DEVICE_EXCHANGE_CASES_HPP
#define CANARY_TESTS_DEVICE_EXCHANGE_CASES_HPP

// Energy-exchange checks for the host and native compilations (see
// cell_inputs.hpp). Standalone chemistry-coupled terms receive their rates of
// progress as input; the combined sources evaluate their own.

#include "cell_inputs.hpp"

#include <array>
#include <canary/energy_exchange.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/state.hpp>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace canary::test
{

// The including file defines CANARY_CHECK_ENTRY (see thermodynamics_cases.hpp).

// Two-temperature terms at (rho, T, T_ve) of the selected cell.
template <std::size_t N, typename X, typename G, typename T>
CANARY_CHECK_ENTRY void
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
CANARY_CHECK_ENTRY void two_temperature_chemistry_vibration_check(
  Output<PowerDensity<T, typename X::units>> out,
  const X* exchange,
  const M* mechanism,
  const RatesOfProgress<T, typename X::units, M::reaction_count>* rates,
  Temperature<T, typename X::units> temp,
  Temperature<T, typename X::units> temp_ve)
{ out[0] = chemistry_vibration_source(*exchange, *mechanism, *rates, temp, temp_ve); }

template <std::size_t N, typename M, typename X, typename G, typename T>
CANARY_CHECK_ENTRY void
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
CANARY_CHECK_ENTRY void vibration_translation_check(Output<ModeSources<G, T>> out,
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
CANARY_CHECK_ENTRY void
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
CANARY_CHECK_ENTRY void electron_translation_check(Output<ModeSources<G, T>> out,
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
CANARY_CHECK_ENTRY void vibration_electron_check(Output<ModeSources<G, T>> out,
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
CANARY_CHECK_ENTRY void electron_impact_ionisation_check(
  Output<ModeSources<G, T>> out,
  const X* exchange,
  const RatesOfProgress<T, typename G::units, M::reaction_count>* rates)
{ out[0] = electron_impact_ionisation_source(*exchange, *rates); }

template <StateForm F, std::size_t N, typename M, typename X, typename G, typename T>
CANARY_CHECK_ENTRY void multi_temperature_sources_check(
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

// Required instances; see thermodynamics_cases.hpp for the convention.

template <std::size_t N, typename M, typename X, typename G, typename T>
void require_two_temperature_exchange()
{
    (void)&two_temperature_vibration_translation_check<N, X, G, T>;
    (void)&two_temperature_chemistry_vibration_check<X, M, T>;
    (void)&two_temperature_sources_check<N, M, X, G, T>;
}

template <StateForm F, std::size_t N, typename M, typename X, typename G, typename T>
void require_state_terms()
{
    (void)&vibration_translation_check<F, N, X, G, T>;
    (void)&chemistry_vibration_check<F, X, M, G, T>;
    (void)&electron_translation_check<F, N, X, G, T>;
    (void)&vibration_electron_check<F, N, X, G, T>;
    (void)&multi_temperature_sources_check<F, N, M, X, G, T>;
}

template <std::size_t N, typename M, typename X, typename G, typename T>
void require_multi_temperature_exchange()
{
    require_state_terms<StateForm::Owning, N, M, X, G, T>();
    (void)&electron_impact_ionisation_check<X, M, G, T>;
}

// Exchange configurations: the full named exchange, every process configured
// empty, and each process alone, so that zero and nonzero counts of every
// process are compiled.
template <std::floating_point P, typename T, UnitSystem U> void exchange_checks()
{
    using Park = ReactionMechanism<P, U, 5, 5>;
    using TwoTemperature = TwoTemperatureMixture<P, U, 5, 4>;
    require_two_temperature_exchange<5, Park, EnergyExchange<P, U, 5, 15, 3>, TwoTemperature, T>();
    require_two_temperature_exchange<5, Park, EnergyExchange<P, U, 5, 0, 0>, TwoTemperature, T>();
    require_two_temperature_exchange<5, Park, EnergyExchange<P, U, 5, 1, 0>, TwoTemperature, T>();
    require_two_temperature_exchange<5, Park, EnergyExchange<P, U, 5, 0, 1>, TwoTemperature, T>();

    using KimJo = ReactionMechanism<P, U, 11, 40>;
    using ThreeTemperature = MultiTemperatureMixture<P, U, 11, 2, 4, 8>;
    using Full = MultiTemperatureExchange<P, U, 11, 2, 30, 21, 10, 1, 2>;
    require_multi_temperature_exchange<11, KimJo, Full, ThreeTemperature, T>();
    require_state_terms<StateForm::Borrowed, 11, KimJo, Full, ThreeTemperature, T>();
    const auto only = [](auto v, auto c, auto e, auto ve, auto i) {
        using X = MultiTemperatureExchange<P,
                                           U,
                                           11,
                                           2,
                                           decltype(v)::value,
                                           decltype(c)::value,
                                           decltype(e)::value,
                                           decltype(ve)::value,
                                           decltype(i)::value>;
        require_multi_temperature_exchange<11, KimJo, X, ThreeTemperature, T>();
    };
    using Zero = std::integral_constant<std::size_t, 0>;
    using One = std::integral_constant<std::size_t, 1>;
    only(Zero {}, Zero {}, Zero {}, Zero {}, Zero {});
    only(One {}, Zero {}, Zero {}, Zero {}, Zero {});
    only(Zero {}, One {}, Zero {}, Zero {}, Zero {});
    only(Zero {}, Zero {}, One {}, Zero {}, Zero {});
    only(Zero {}, Zero {}, Zero {}, One {}, Zero {});
    only(Zero {}, Zero {}, Zero {}, Zero {}, One {});

    // Four-temperature air5 without a free-electron species.
    require_multi_temperature_exchange<5,
                                       Park,
                                       MultiTemperatureExchange<P, U, 5, 3, 15, 3, 0, 0, 0>,
                                       MultiTemperatureMixture<P, U, 5, 3, 4, 8>,
                                       T>();
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_EXCHANGE_CASES_HPP
