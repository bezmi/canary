// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_TESTS_DEVICE_CELL_INPUTS_HPP
#define CANARY_TESTS_DEVICE_CELL_INPUTS_HPP

// Input preparation shared by the host and native numerical checks. A check reads
// one cell from caller storage, builds the borrowed views the public operation
// accepts, and calls that operation. Storage stays with the caller; nothing here
// owns or copies a result. Define CANARY_HOST_DEVICE before including this header.

#include "../dual.hpp"

#include <array>
#include <canary/indexed_view.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/thermo.hpp>
#include <cstddef>
#include <cstdint>
#include <mdspan/mdspan.hpp>
#include <type_traits>

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

// The composition of the selected cell as a public operation accepts it: raw
// mass fractions, a mixture prepared from them, or partial densities (with the
// density supplied separately).
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

// The M mode inputs (of kind K) of the selected cell: an owning array, a
// quantity mdspan, or a ModeIndex view.
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

// A multi-temperature state as a public source operation accepts it: the owning
// state itself, or a borrowed view with the scalars of the owning state and its
// mode energies and temperatures read from grids of opposite layouts.
enum class StateForm : std::uint8_t { Owning, Borrowed };

template <StateForm F, Scalar T, UnitSystem U, std::size_t M>
CANARY_HOST_DEVICE auto cell_state(const MultiTemperatureState<T, U, M>& state,
                                   const T* mode_energies,
                                   const T* mode_temperatures)
{
    if constexpr (F == StateForm::Owning) {
        return state;
    } else {
        const auto mode_e =
          cell_modes<ModeForm::Indexed, Storage::LayoutRight, kind::specific_energy, U, M>(
            mode_energies);
        const auto mode_temp =
          cell_modes<ModeForm::Indexed, Storage::LayoutLeft, kind::mode_temperature, U, M>(
            mode_temperatures);
        return MultiTemperatureStateView<T, U, decltype(mode_e), decltype(mode_temp)> {state.rho,
                                                                                       state.e,
                                                                                       state.p,
                                                                                       state.temp,
                                                                                       mode_e,
                                                                                       mode_temp};
    }
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_CELL_INPUTS_HPP
