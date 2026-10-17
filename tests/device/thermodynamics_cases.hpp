#ifndef CANARY_TESTS_DEVICE_THERMODYNAMICS_CASES_HPP
#define CANARY_TESTS_DEVICE_THERMODYNAMICS_CASES_HPP

// Input preparation shared by the host and native thermodynamic checks. Each
// check reads one cell from caller storage, builds the borrowed views the public
// operation accepts, and calls that operation in the entry file. Storage stays
// with the caller; nothing here owns or copies a result. Define
// CANARY_HOST_DEVICE before including this header.

#include "../dual.hpp"

#include <array>
#include <canary/indexed_view.hpp>
#include <canary/models/calorically_perfect.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/nasa9.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/thermo.hpp>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <mdspan/mdspan.hpp>
#include <tuple>
#include <type_traits>
#include <utility>

namespace canary::test
{

// The single complete result of one check, written by the native kernel.
template <typename Result> using Output = Kokkos::mdspan<Result, Kokkos::extents<std::size_t, 1>>;

// Storage of the C components (species or modes) of the selected cell:
// - Contiguous: C consecutive values at the address;
// - LayoutRight / LayoutLeft: the middle row of a [3, C] grid in that layout,
//   whose first element follows one padding element at the address.
enum class Storage : std::uint8_t { Contiguous, LayoutRight, LayoutLeft };

// The components of the selected cell: a contiguous rank-one mdspan, or a
// strided slice whose offset and stride come from the grid mapping.
template <Storage S, std::size_t C, typename V>
CANARY_HOST_DEVICE auto cell_components(const V* data)
{
    using Extents = Kokkos::extents<std::size_t, C>;
    if constexpr (S == Storage::Contiguous) {
        return Kokkos::mdspan<const V, Extents>(data);
    } else {
        using Layout =
          std::conditional_t<S == Storage::LayoutRight, Kokkos::layout_right, Kokkos::layout_left>;
        const Kokkos::mdspan<const V, Kokkos::extents<std::size_t, 3, C>, Layout> grid(data + 1);
        using Slice = Kokkos::mdspan<const V, Extents, Kokkos::layout_stride>;
        return Slice(
          data + 1 + grid.mapping()(1, 0),
          typename Slice::mapping_type(Extents {}, std::array<std::size_t, 1> {grid.stride(1)}));
    }
}

// The composition of the selected cell as a public closure accepts it: raw mass
// fractions, a mixture prepared from them, or partial densities (with the
// density supplied separately to the closure).
enum class Composition : std::uint8_t { Fractions, Mixture, PartialDensities };

template <Composition F, Storage S, std::size_t N, typename G, typename T>
CANARY_HOST_DEVICE auto cell_composition(const G& gas, const T* data)
{
    const auto values = cell_components<S, N>(data);
    if constexpr (F == Composition::Fractions) {
        return as_kind<kind::mass_fraction>(values);
    } else if constexpr (F == Composition::Mixture) {
        return mixture(gas, as_kind<kind::mass_fraction>(values));
    } else {
        return as_kind<kind::density, typename G::units>(values);
    }
}

// The M mode inputs (mode temperatures or mode energy densities, kind K) of the
// selected cell: an owning array, a quantity mdspan, or a ModeIndex view.
enum class ModeForm : std::uint8_t { Owning, Quantities, Indexed };

template <ModeForm F, Storage S, KindTag K, UnitSystem U, std::size_t M, typename T>
    requires(F != ModeForm::Owning || S == Storage::Contiguous)
CANARY_HOST_DEVICE auto cell_modes(const T* data)
{
    if constexpr (F == ModeForm::Owning) {
        std::array<Quantity<K, T, U>, M> modes {};
        for (const ModeIndex m : mode_indices(0, M)) {
            modes[m.value()] = Quantity<K, T, U> {data[m.value()]};
        }
        return modes;
    } else if constexpr (F == ModeForm::Quantities) {
        return as_kind<K, U>(cell_components<S, M>(data));
    } else {
        return as_indexed<ModeIndex>(as_kind<K, U>(cell_components<S, M>(data)));
    }
}

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

template <typename G, typename... Pair> void require_closures(std::tuple<Pair...>* /*pairs*/)
{ ((void)&closure_check<G, typename Pair::first_type, typename Pair::second_type>, ...); }

template <Composition F, Storage S, std::size_t N, typename G, typename... Extra, typename... Pair>
void require_mixture_closures(std::tuple<Pair...>* /*pairs*/)
{
    ((void)&mixture_closure_check<F,
                                  S,
                                  N,
                                  G,
                                  typename Pair::first_type,
                                  typename Pair::second_type,
                                  Extra...>,
     ...);
}

template <Composition F,
          Storage S,
          ModeForm MF,
          KindTag K,
          std::size_t N,
          typename G,
          typename... Pair>
void require_modal_closures(std::tuple<Pair...>* /*pairs*/)
{
    ((void)&modal_closure_check<F,
                                S,
                                MF,
                                K,
                                N,
                                G,
                                typename Pair::first_type,
                                typename Pair::second_type>,
     ...);
}

// Closures and mixtures of a mixture model of N species, with partial-density
// input when Densities.
template <std::size_t N, typename G, typename T, bool Densities> void require_mixture_model()
{
    using U = typename G::units;
    require_mixture_closures<Composition::Fractions, Storage::Contiguous, N, G>(
      static_cast<ClosurePairs<T, U>*>(nullptr));
    DensityTemperature<T, U>* const rho_t = nullptr;
    require_mixture_closures<Composition::Fractions, Storage::LayoutRight, N, G>(rho_t);
    require_mixture_closures<Composition::Fractions, Storage::LayoutLeft, N, G>(rho_t);
    require_mixture_closures<Composition::Mixture, Storage::Contiguous, N, G>(rho_t);
    (void)&fractions_mixture_check<Storage::Contiguous, N, G, T>;
    if constexpr (Densities) {
        require_mixture_closures<Composition::PartialDensities, Storage::Contiguous, N, G>(rho_t);
        (void)&densities_mixture_check<Storage::Contiguous, N, G, T>;
    }
}

// Modal closures of a two- or multi-temperature model of N species, with mode
// inputs of kind K.
template <KindTag K, std::size_t N, typename G, typename T, bool Densities>
void require_modal_kind()
{
    using U = typename G::units;
    require_modal_closures<Composition::Fractions,
                           Storage::Contiguous,
                           ModeForm::Quantities,
                           K,
                           N,
                           G>(static_cast<ClosurePairs<T, U>*>(nullptr));
    DensityTemperature<T, U>* const rho_t = nullptr;
    require_modal_closures<Composition::Fractions,
                           Storage::LayoutRight,
                           ModeForm::Indexed,
                           K,
                           N,
                           G>(rho_t);
    require_modal_closures<Composition::Fractions,
                           Storage::LayoutLeft,
                           ModeForm::Quantities,
                           K,
                           N,
                           G>(rho_t);
    require_modal_closures<Composition::Mixture, Storage::Contiguous, ModeForm::Indexed, K, N, G>(
      rho_t);
    if constexpr (Densities) {
        require_modal_closures<Composition::PartialDensities,
                               Storage::Contiguous,
                               ModeForm::Quantities,
                               K,
                               N,
                               G>(rho_t);
    }
}

template <std::size_t N, typename G, typename T, bool Densities> void require_modal_model()
{
    require_mixture_model<N, G, T, Densities>();
    require_modal_kind<kind::mode_temperature, N, G, T, Densities>();
    require_modal_kind<kind::mode_energy_density, N, G, T, Densities>();
    (void)&mixture_properties_check<N, G, state_t<G, T>>;
    (void)&species_properties_check<G, T>;
    (void)&species_isobaric_specific_heat_check<G, T>;
    (void)&species_mode_energy_check<G, T>;
}

template <std::size_t N, typename G, typename T, bool Densities> void require_thermally_perfect()
{
    require_mixture_model<N, G, T, Densities>();
    require_mixture_closures<Composition::Fractions,
                             Storage::Contiguous,
                             N,
                             G,
                             TemperatureGuess<T, typename G::units>>(
      static_cast<EnergyPairs<T, typename G::units>*>(nullptr));
    (void)&mixture_properties_check<N, G, state_t<G, T>>;
    (void)&species_properties_check<G, T>;
    (void)&species_isobaric_specific_heat_check<G, T>;
}

// Energy closures of the thermally perfect models, with and without a guess, for
// a parameter type P that differs in precision from the state scalar T.
template <std::size_t N, typename G, typename T> void require_energy_closures()
{
    EnergyPairs<T, typename G::units>* const pairs = nullptr;
    require_mixture_closures<Composition::Fractions, Storage::Contiguous, N, G>(pairs);
    require_mixture_closures<Composition::Fractions,
                             Storage::Contiguous,
                             N,
                             G,
                             TemperatureGuess<T, typename G::units>>(pairs);
}

template <std::floating_point P, typename T, UnitSystem U> void thermodynamic_mixed_checks()
{
    require_energy_closures<1, ThermallyPerfectMixture<P, U, 1, 1>, T>();
    require_energy_closures<1, ThermallyPerfectMixture<P, U, 1, 4>, T>();
    require_energy_closures<5, ThermallyPerfectMixture<P, U, 5, 4>, T>();
}

// The thermodynamic checks for parameter type P, state scalar T and units U, in
// two groups of similar size so that each compiles as its own unit: models with
// one temperature and NASA9 curves, then two- and multi-temperature models.
template <std::floating_point P, typename T, UnitSystem U> void thermodynamic_single_checks()
{
    // Perfect gas.
    using Perfect = CaloricallyPerfect<P, U>;
    require_closures<Perfect>(static_cast<ClosurePairs<T, U>*>(nullptr));
    (void)&properties_check<Perfect, ThermoState<T, U>>;
    // Perfect mixtures: 5 and 11 species, at and below capacity.
    require_mixture_model<5, CaloricallyPerfectMixture<P, U, 5>, T, true>();
    require_mixture_model<5, CaloricallyPerfectMixture<P, U, 16>, T, true>();
    require_mixture_model<11, CaloricallyPerfectMixture<P, U, 11>, T, true>();
    require_mixture_model<11, CaloricallyPerfectMixture<P, U, 16>, T, true>();
    (void)&mixture_properties_check<5, CaloricallyPerfectMixture<P, U, 5>, ThermoState<T, U>>;
    (void)&mixture_properties_check<11, CaloricallyPerfectMixture<P, U, 16>, ThermoState<T, U>>;
    // Thermally perfect: controlled one-segment, named N2 and five-species air.
    require_thermally_perfect<1, ThermallyPerfectMixture<P, U, 1, 1>, T, true>();
    require_thermally_perfect<1, ThermallyPerfectMixture<P, U, 1, 4>, T, false>();
    require_thermally_perfect<5, ThermallyPerfectMixture<P, U, 5, 4>, T, true>();
    // NASA9 curves of one, two and four segments.
    (void)&nasa9_cp_h_check<Nasa9Curve<P, 1>, T>;
    (void)&nasa9_s_check<Nasa9Curve<P, 1>, T>;
    (void)&nasa9_cp_h_check<Nasa9Curve<P, 2>, T>;
    (void)&nasa9_s_check<Nasa9Curve<P, 2>, T>;
    (void)&nasa9_cp_h_check<Nasa9Curve<P, 4>, T>;
    (void)&nasa9_s_check<Nasa9Curve<P, 4>, T>;
}

template <std::floating_point P, typename T, UnitSystem U> void thermodynamic_modal_checks()
{
    // Two-temperature air.
    using TwoTemperature = TwoTemperatureMixture<P, U, 5, 4>;
    require_modal_model<5, TwoTemperature, T, true>();
    (void)&species_vibrational_electronic_energy_check<TwoTemperature, T>;
    // Multi-temperature: four-temperature air5, three-temperature air11 and the
    // controlled two-species model.
    require_modal_model<5, MultiTemperatureMixture<P, U, 5, 3, 4, 8>, T, false>();
    using ThreeTemperature = MultiTemperatureMixture<P, U, 11, 2, 4, 8>;
    require_modal_model<11, ThreeTemperature, T, true>();
    (void)&electron_pressure_check<11, ThreeTemperature, state_t<ThreeTemperature, T>>;
    require_modal_model<2, MultiTemperatureMixture<P, U, 2, 2, 1, 2>, T, false>();
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_THERMODYNAMICS_CASES_HPP
