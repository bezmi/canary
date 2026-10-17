#ifndef CANARY_TESTS_DEVICE_REACTIONS_CASES_HPP
#define CANARY_TESTS_DEVICE_REACTIONS_CASES_HPP

// Reaction checks for the host and native compilations (see cell_inputs.hpp).
// rates_of_progress and production_rates are separate checks: production uses
// its own net-rate evaluation.

#include "cell_inputs.hpp"

#include <array>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/state.hpp>
#include <concepts>
#include <cstddef>

namespace canary::test
{

// The including file defines CANARY_CHECK_ENTRY (see thermodynamics_cases.hpp).

// Rates of progress of the selected cell at one temperature, or with the scalar
// T_ve of a one-mode model as the extra argument.
template <Storage S, std::size_t N, typename M, typename G, typename T, typename... Extra>
CANARY_CHECK_ENTRY void
  rates_check(Output<RatesOfProgress<T, typename G::units, M::reaction_count>> out,
              const M* mechanism,
              const G* gas,
              const T* composition,
              Density<T, typename G::units> rho,
              Temperature<T, typename G::units> temp,
              Extra... extra)
{
    out[0] = rates_of_progress(*mechanism,
                               *gas,
                               cell_composition<Composition::Mixture, S, N>(*gas, composition),
                               rho,
                               temp,
                               extra...);
}

template <Storage S, std::size_t N, typename M, typename G, typename T, typename... Extra>
CANARY_CHECK_ENTRY void
  production_check(Output<std::array<MassProductionRate<T, typename G::units>, N>> out,
                   const M* mechanism,
                   const G* gas,
                   const T* composition,
                   Density<T, typename G::units> rho,
                   Temperature<T, typename G::units> temp,
                   Extra... extra)
{
    out[0] = production_rates(*mechanism,
                              *gas,
                              cell_composition<Composition::Mixture, S, N>(*gas, composition),
                              rho,
                              temp,
                              extra...);
}

// As above with the mode temperatures of the selected cell in form F.
template <ModeForm F, Storage S, std::size_t N, typename M, typename G, typename T>
CANARY_CHECK_ENTRY void
  modal_rates_check(Output<RatesOfProgress<T, typename G::units, M::reaction_count>> out,
                    const M* mechanism,
                    const G* gas,
                    const T* composition,
                    Density<T, typename G::units> rho,
                    Temperature<T, typename G::units> temp,
                    const T* modes)
{
    out[0] = rates_of_progress(
      *mechanism,
      *gas,
      cell_composition<Composition::Mixture, S, N>(*gas, composition),
      rho,
      temp,
      cell_modes<F, S, kind::mode_temperature, typename G::units, G::mode_count>(modes));
}

template <ModeForm F, Storage S, std::size_t N, typename M, typename G, typename T>
CANARY_CHECK_ENTRY void
  modal_production_check(Output<std::array<MassProductionRate<T, typename G::units>, N>> out,
                         const M* mechanism,
                         const G* gas,
                         const T* composition,
                         Density<T, typename G::units> rho,
                         Temperature<T, typename G::units> temp,
                         const T* modes)
{
    out[0] = production_rates(
      *mechanism,
      *gas,
      cell_composition<Composition::Mixture, S, N>(*gas, composition),
      rho,
      temp,
      cell_modes<F, S, kind::mode_temperature, typename G::units, G::mode_count>(modes));
}

// Required instances; see thermodynamics_cases.hpp for the convention.

// Rates and production at one temperature (no mode argument).
template <std::size_t N, typename M, typename G, typename T> void require_single_rates()
{
    (void)&rates_check<Storage::Contiguous, N, M, G, T>;
    (void)&production_check<Storage::Contiguous, N, M, G, T>;
}

// Rates and production with owning mode temperatures.
template <std::size_t N, typename M, typename G, typename T> void require_owning_rates()
{
    (void)&modal_rates_check<ModeForm::Owning, Storage::Contiguous, N, M, G, T>;
    (void)&modal_production_check<ModeForm::Owning, Storage::Contiguous, N, M, G, T>;
}

// Rates and production with each borrowed mode form and grid storage once.
template <std::size_t N, typename M, typename G, typename T> void require_borrowed_rates()
{
    (void)&modal_rates_check<ModeForm::Quantities, Storage::LayoutRight, N, M, G, T>;
    (void)&modal_production_check<ModeForm::Quantities, Storage::LayoutRight, N, M, G, T>;
    (void)&modal_rates_check<ModeForm::Indexed, Storage::LayoutLeft, N, M, G, T>;
    (void)&modal_production_check<ModeForm::Indexed, Storage::LayoutLeft, N, M, G, T>;
}

template <std::floating_point P, typename T, UnitSystem U> void reaction_checks()
{
    // Single temperature: air5 with the 13-reaction mechanism.
    require_single_rates<5,
                         ReactionMechanism<P, U, 5, 13>,
                         ThermallyPerfectMixture<P, U, 5, 4>,
                         T>();
    // Two temperature: Park's 5 reactions with scalar T_ve, owning modes and
    // explicit single temperature.
    using Park = ReactionMechanism<P, U, 5, 5>;
    using TwoTemperature = TwoTemperatureMixture<P, U, 5, 4>;
    (void)&rates_check<Storage::Contiguous, 5, Park, TwoTemperature, T, Temperature<T, U>>;
    (void)&production_check<Storage::Contiguous, 5, Park, TwoTemperature, T, Temperature<T, U>>;
    require_owning_rates<5, Park, TwoTemperature, T>();
    require_single_rates<5, Park, TwoTemperature, T>();
    // Three temperature: air11 with 40 reactions, every mode form.
    using KimJo = ReactionMechanism<P, U, 11, 40>;
    using ThreeTemperature = MultiTemperatureMixture<P, U, 11, 2, 4, 8>;
    require_owning_rates<11, KimJo, ThreeTemperature, T>();
    require_borrowed_rates<11, KimJo, ThreeTemperature, T>();
    require_single_rates<11, KimJo, ThreeTemperature, T>();
    // Four internal modes: the 12-reaction interleaved mechanism.
    using Interleaved = ReactionMechanism<P, U, 11, 12>;
    using FourModes = MultiTemperatureMixture<P, U, 11, 4, 4, 8>;
    require_owning_rates<11, Interleaved, FourModes, T>();
    require_single_rates<11, Interleaved, FourModes, T>();
    // Park's reactions against four-temperature air5 (used without electrons).
    require_owning_rates<5, Park, MultiTemperatureMixture<P, U, 5, 3, 4, 8>, T>();
}

// Multi-temperature rates and production with P and T of different precision.
template <std::floating_point P, typename T, UnitSystem U> void reaction_mixed_checks()
{
    require_owning_rates<11,
                         ReactionMechanism<P, U, 11, 40>,
                         MultiTemperatureMixture<P, U, 11, 2, 4, 8>,
                         T>();
    require_owning_rates<11,
                         ReactionMechanism<P, U, 11, 12>,
                         MultiTemperatureMixture<P, U, 11, 4, 4, 8>,
                         T>();
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_REACTIONS_CASES_HPP
