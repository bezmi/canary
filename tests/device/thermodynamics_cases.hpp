#ifndef CANARY_TESTS_DEVICE_THERMODYNAMICS_CASES_HPP
#define CANARY_TESTS_DEVICE_THERMODYNAMICS_CASES_HPP

// Thermodynamic checks for the host and native compilations (see cell_inputs.hpp).

#include "cell_inputs.hpp"
#include "check_groups.hpp"

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
#include <concepts>
#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>

namespace canary::test
{

// Checks. The including file defines CANARY_CHECK_ENTRY as the function
// specifier of an entry: __global__ for a native kernel, empty for the host
// companion. Each entry writes the complete public result to out[0].

// Closure of a model without composition input from two state variables.
template <typename G, typename A, typename B>
CANARY_CHECK_ENTRY void
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
CANARY_CHECK_ENTRY void
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
CANARY_CHECK_ENTRY void
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
CANARY_CHECK_ENTRY void
  fractions_mixture_check(Output<mixture_t<G, T>> out, const G* gas, const T* composition)
{ out[0] = cell_composition<Composition::Mixture, S, N>(*gas, composition); }

// Mixture of the selected cell from partial densities and the density.
template <Storage S, std::size_t N, typename G, typename T>
CANARY_CHECK_ENTRY void densities_mixture_check(Output<mixture_t<G, T>> out,
                                                const G* gas,
                                                const T* composition,
                                                Density<T, typename G::units> rho)
{
    out[0] =
      mixture(*gas, cell_composition<Composition::PartialDensities, S, N>(*gas, composition), rho);
}

// Properties of an owning state of a model without composition input.
template <typename G, typename State>
CANARY_CHECK_ENTRY void properties_check(
  Output<decltype(properties(std::declval<const G&>(), std::declval<State>()))> out,
  const G* gas,
  State state)
{ out[0] = properties(*gas, state); }

// Properties of an owning state of a mixture model, with the mixture prepared
// from the mass fractions of the selected cell.
template <std::size_t N, typename G, typename State>
CANARY_CHECK_ENTRY void mixture_properties_check(
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
CANARY_CHECK_ENTRY void electron_pressure_check(
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
CANARY_CHECK_ENTRY void species_properties_check(
  Output<decltype(species_properties(std::declval<const G&>(),
                                     std::declval<Temperature<T, typename G::units>>()))> out,
  const G* gas,
  Temperature<T, typename G::units> temp)
{ out[0] = species_properties(*gas, temp); }

template <typename G, typename T>
CANARY_CHECK_ENTRY void species_isobaric_specific_heat_check(
  Output<decltype(species_isobaric_specific_heat(
    std::declval<const G&>(), std::declval<Temperature<T, typename G::units>>()))> out,
  const G* gas,
  Temperature<T, typename G::units> temp)
{ out[0] = species_isobaric_specific_heat(*gas, temp); }

template <typename G, typename T>
CANARY_CHECK_ENTRY void species_vibrational_electronic_energy_check(
  Output<decltype(species_vibrational_electronic_energy(
    std::declval<const G&>(), std::declval<Temperature<T, typename G::units>>()))> out,
  const G* gas,
  Temperature<T, typename G::units> temp)
{ out[0] = species_vibrational_electronic_energy(*gas, temp); }

template <typename G, typename T>
CANARY_CHECK_ENTRY void species_mode_energy_check(
  Output<decltype(species_mode_energy(std::declval<const G&>(),
                                      ModeIndex {},
                                      std::declval<Temperature<T, typename G::units>>()))> out,
  const G* gas,
  ModeIndex m,
  Temperature<T, typename G::units> temp)
{ out[0] = species_mode_energy(*gas, m, temp); }

// NASA9 cp/R and h/R, and s/R, of a curve at temperature t.
template <typename Curve, typename T>
CANARY_CHECK_ENTRY void nasa9_cp_h_check(Output<Nasa9CpH<T>> out, const Curve* curve, T t)
{ out[0] = nasa9_cp_h(*curve, t, math::log(t)); }

template <typename Curve, typename T>
CANARY_CHECK_ENTRY void nasa9_s_check(Output<T> out, const Curve* curve, T t)
{ out[0] = nasa9_s(*curve, t, math::log(t)); }

// Required instances. Each list names the entries compiled for one model, one
// parameter type P, one state scalar T and one unit system U: every closure pair
// with contiguous mass fractions, and each storage form, composition form, mode
// form and mode kind at least once on (rho, T).

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

template <typename G, typename... Pair>
void require_closures(Report& report, std::tuple<Pair...>* /*pairs*/)
{
    (unevaluated<&closure_check<G, typename Pair::first_type, typename Pair::second_type>>(report),
     ...);
}

template <Composition F, Storage S, std::size_t N, typename G, typename... Extra, typename... Pair>
void require_mixture_closures(Report& report, std::tuple<Pair...>* /*pairs*/)
{
    (unevaluated<&mixture_closure_check<F,
                                        S,
                                        N,
                                        G,
                                        typename Pair::first_type,
                                        typename Pair::second_type,
                                        Extra...>>(report),
     ...);
}

template <Composition F,
          Storage S,
          ModeForm MF,
          KindTag K,
          std::size_t N,
          typename G,
          typename... Pair>
void require_modal_closures(Report& report, std::tuple<Pair...>* /*pairs*/)
{
    (unevaluated<&modal_closure_check<F,
                                      S,
                                      MF,
                                      K,
                                      N,
                                      G,
                                      typename Pair::first_type,
                                      typename Pair::second_type>>(report),
     ...);
}

// Closures and mixtures of a mixture model of N species, with partial-density
// input when Densities.
template <std::size_t N, typename G, typename T, bool Densities>
void require_mixture_model(Report& report)
{
    using U = typename G::units;
    require_mixture_closures<Composition::Fractions, Storage::Contiguous, N, G>(
      report,
      static_cast<ClosurePairs<T, U>*>(nullptr));
    DensityTemperature<T, U>* const rho_t = nullptr;
    require_mixture_closures<Composition::Fractions, Storage::LayoutRight, N, G>(report, rho_t);
    require_mixture_closures<Composition::Fractions, Storage::LayoutLeft, N, G>(report, rho_t);
    require_mixture_closures<Composition::Mixture, Storage::Contiguous, N, G>(report, rho_t);
    unevaluated<&fractions_mixture_check<Storage::Contiguous, N, G, T>>(report);
    if constexpr (Densities) {
        require_mixture_closures<Composition::PartialDensities, Storage::Contiguous, N, G>(report,
                                                                                           rho_t);
        unevaluated<&densities_mixture_check<Storage::Contiguous, N, G, T>>(report);
    }
}

// Modal closures of a two- or multi-temperature model of N species, with mode
// inputs of kind K.
template <KindTag K, std::size_t N, typename G, typename T, bool Densities>
void require_modal_kind(Report& report)
{
    using U = typename G::units;
    require_modal_closures<Composition::Fractions,
                           Storage::Contiguous,
                           ModeForm::Quantities,
                           K,
                           N,
                           G>(report, static_cast<ClosurePairs<T, U>*>(nullptr));
    DensityTemperature<T, U>* const rho_t = nullptr;
    require_modal_closures<Composition::Fractions,
                           Storage::LayoutRight,
                           ModeForm::Indexed,
                           K,
                           N,
                           G>(report, rho_t);
    require_modal_closures<Composition::Fractions,
                           Storage::LayoutLeft,
                           ModeForm::Quantities,
                           K,
                           N,
                           G>(report, rho_t);
    require_modal_closures<Composition::Mixture, Storage::Contiguous, ModeForm::Indexed, K, N, G>(
      report,
      rho_t);
    if constexpr (Densities) {
        require_modal_closures<Composition::PartialDensities,
                               Storage::Contiguous,
                               ModeForm::Quantities,
                               K,
                               N,
                               G>(report, rho_t);
    }
}

template <std::size_t N, typename G, typename T, bool Densities>
void require_modal_model(Report& report)
{
    require_mixture_model<N, G, T, Densities>(report);
    require_modal_kind<kind::mode_temperature, N, G, T, Densities>(report);
    require_modal_kind<kind::mode_energy_density, N, G, T, Densities>(report);
    unevaluated<&mixture_properties_check<N, G, state_t<G, T>>>(report);
    unevaluated<&species_properties_check<G, T>>(report);
    unevaluated<&species_isobaric_specific_heat_check<G, T>>(report);
    unevaluated<&species_mode_energy_check<G, T>>(report);
}

template <std::size_t N, typename G, typename T, bool Densities>
void require_thermally_perfect(Report& report)
{
    require_mixture_model<N, G, T, Densities>(report);
    require_mixture_closures<Composition::Fractions,
                             Storage::Contiguous,
                             N,
                             G,
                             TemperatureGuess<T, typename G::units>>(
      report,
      static_cast<EnergyPairs<T, typename G::units>*>(nullptr));
    unevaluated<&mixture_properties_check<N, G, state_t<G, T>>>(report);
    unevaluated<&species_properties_check<G, T>>(report);
    unevaluated<&species_isobaric_specific_heat_check<G, T>>(report);
}

// Energy closures of the thermally perfect models, with and without a guess, for
// a parameter type P that differs in precision from the state scalar T.
template <std::size_t N, typename G, typename T> void require_energy_closures(Report& report)
{
    EnergyPairs<T, typename G::units>* const pairs = nullptr;
    require_mixture_closures<Composition::Fractions, Storage::Contiguous, N, G>(report, pairs);
    require_mixture_closures<Composition::Fractions,
                             Storage::Contiguous,
                             N,
                             G,
                             TemperatureGuess<T, typename G::units>>(report, pairs);
}

template <std::floating_point P, typename T, UnitSystem U> Report thermodynamic_mixed_checks()
{
    Report report;
    require_energy_closures<1, ThermallyPerfectMixture<P, U, 1, 1>, T>(report);
    require_energy_closures<1, ThermallyPerfectMixture<P, U, 1, 4>, T>(report);
    require_energy_closures<5, ThermallyPerfectMixture<P, U, 5, 4>, T>(report);
    return report;
}

// The thermodynamic checks for parameter type P, state scalar T and units U, in
// two groups of similar size so that each compiles as its own unit: models with
// one temperature and NASA9 curves, then two- and multi-temperature models.
template <std::floating_point P, typename T, UnitSystem U> Report thermodynamic_single_checks()
{
    Report report;
    // Perfect gas.
    using Perfect = CaloricallyPerfect<P, U>;
    require_closures<Perfect>(report, static_cast<ClosurePairs<T, U>*>(nullptr));
    unevaluated<&properties_check<Perfect, ThermoState<T, U>>>(report);
    // Perfect mixtures: 5 and 11 species, at and below capacity.
    require_mixture_model<5, CaloricallyPerfectMixture<P, U, 5>, T, true>(report);
    require_mixture_model<5, CaloricallyPerfectMixture<P, U, 16>, T, true>(report);
    require_mixture_model<11, CaloricallyPerfectMixture<P, U, 11>, T, true>(report);
    require_mixture_model<11, CaloricallyPerfectMixture<P, U, 16>, T, true>(report);
    unevaluated<
      &mixture_properties_check<5, CaloricallyPerfectMixture<P, U, 5>, ThermoState<T, U>>>(report);
    unevaluated<
      &mixture_properties_check<11, CaloricallyPerfectMixture<P, U, 16>, ThermoState<T, U>>>(
      report);
    // Thermally perfect: controlled one-segment, named N2 and five-species air.
    require_thermally_perfect<1, ThermallyPerfectMixture<P, U, 1, 1>, T, true>(report);
    require_thermally_perfect<1, ThermallyPerfectMixture<P, U, 1, 4>, T, false>(report);
    require_thermally_perfect<5, ThermallyPerfectMixture<P, U, 5, 4>, T, true>(report);
    // NASA9 curves of one, two and four segments.
    unevaluated<&nasa9_cp_h_check<Nasa9Curve<P, 1>, T>>(report);
    unevaluated<&nasa9_s_check<Nasa9Curve<P, 1>, T>>(report);
    unevaluated<&nasa9_cp_h_check<Nasa9Curve<P, 2>, T>>(report);
    unevaluated<&nasa9_s_check<Nasa9Curve<P, 2>, T>>(report);
    unevaluated<&nasa9_cp_h_check<Nasa9Curve<P, 4>, T>>(report);
    unevaluated<&nasa9_s_check<Nasa9Curve<P, 4>, T>>(report);
    return report;
}

template <std::floating_point P, typename T, UnitSystem U> Report thermodynamic_modal_checks()
{
    Report report;
    // Two-temperature air.
    using TwoTemperature = TwoTemperatureMixture<P, U, 5, 4>;
    require_modal_model<5, TwoTemperature, T, true>(report);
    unevaluated<&species_vibrational_electronic_energy_check<TwoTemperature, T>>(report);
    // Multi-temperature: four-temperature air5, three-temperature air11 and the
    // controlled two-species model.
    require_modal_model<5, MultiTemperatureMixture<P, U, 5, 3, 4, 8>, T, false>(report);
    using ThreeTemperature = MultiTemperatureMixture<P, U, 11, 2, 4, 8>;
    require_modal_model<11, ThreeTemperature, T, true>(report);
    unevaluated<&electron_pressure_check<11, ThreeTemperature, state_t<ThreeTemperature, T>>>(
      report);
    require_modal_model<2, MultiTemperatureMixture<P, U, 2, 2, 1, 2>, T, false>(report);
    return report;
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_THERMODYNAMICS_CASES_HPP
