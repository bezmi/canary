// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#include "dual.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <canary/indexed_view.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <concepts>
#include <cstddef>
#include <mdspan/mdspan.hpp>
#include <type_traits>
#include <utility>

namespace
{
using Extents = Kokkos::dextents<std::size_t, 1>;
template <typename V, typename I>
concept ReadableAt = requires(V v, I i) { v[i]; };
template <typename V, typename I>
concept WritableAt = requires(V v, I i) { v[i] = 1.0; };
template <typename I, typename V>
concept CanIndex = requires(V v) { canary::as_indexed<I>(v); };

// Accessor state must survive wrapping and explicit mdspan interoperability.
struct OffsetAccessor {
    using offset_policy = OffsetAccessor;
    using element_type = double;
    using reference = double&;
    using data_handle_type = double*;
    std::size_t offset_value;
    constexpr reference access(data_handle_type p, std::size_t i) const noexcept
    { return p[offset_value + i]; }
    constexpr data_handle_type offset(data_handle_type p, std::size_t i) const noexcept
    { return p + i; }
};

using Raw = Kokkos::mdspan<double, Extents>;
using View = canary::IndexedView<canary::ModeIndex, Raw>;
using ConstView = canary::IndexedView<canary::ModeIndex, Kokkos::mdspan<const double, Extents>>;
static_assert(ReadableAt<View, canary::ModeIndex>);
static_assert(!ReadableAt<View, canary::SpeciesIndex>);
static_assert(!ReadableAt<View, std::size_t>);
static_assert(!ReadableAt<View, int>);
static_assert(WritableAt<const View, canary::ModeIndex>);
static_assert(!WritableAt<ConstView, canary::ModeIndex>);
static_assert(!std::is_convertible_v<View, Raw>);
static_assert(!std::is_convertible_v<View, canary::IndexedView<canary::SpeciesIndex, Raw>>);
static_assert(
  !CanIndex<canary::ModeIndex, Kokkos::mdspan<double, Kokkos::dextents<std::size_t, 2>>>);
static_assert(!CanIndex<std::size_t, Raw>);
static_assert(!CanIndex<canary::ModeIndex, std::array<double, 3>>);
static_assert(sizeof(View) == sizeof(Raw));
static_assert(alignof(View) == alignof(Raw));
static_assert(std::is_trivially_copyable_v<View> == std::is_trivially_copyable_v<Raw>);
static_assert(noexcept(canary::as_indexed<canary::ModeIndex>(std::declval<Raw>())) ==
              std::is_nothrow_copy_constructible_v<Raw>);
static_assert(noexcept(std::declval<View>()[canary::ModeIndex {0}]) ==
              noexcept(std::declval<Raw>()[std::size_t {0}]));
}  // namespace

void test_indexed_view()
{
    using namespace boost::ut;
    "indexed views preserve layouts offsets and mutation"_test = [] {
        const auto check = []<typename Layout>() {
            using Mapping = typename Layout::template mapping<Extents>;
            std::array<double, 12> data {};
            const Mapping mapping = [] {
                if constexpr (std::same_as<Layout, Kokkos::layout_stride>) {
                    return Mapping(Extents(3), std::array<std::size_t, 1> {3});
                } else {
                    return Mapping(Extents(3));
                }
            }();
            const Kokkos::mdspan<double, Extents, Layout> raw(data.data() + 1, mapping);
            const auto typed = canary::as_indexed<canary::ModeIndex>(raw);
            static_assert(sizeof(typed) == sizeof(raw));
            expect(typed.rank() == 1_u && typed.extent(0) == 3_u);
            for (const auto m : canary::mode_indices(0, 3)) {
                typed[m] = double(m.value() + 10);
                expect(&typed[m] == &raw[m.value()]);
                expect(data[1 + mapping(m.value())] == double(m.value() + 10));
            }
            expect(data[0] == 0.0);
            expect(typed.base().data_handle() == raw.data_handle());
            expect(typed.base().mapping() == mapping);
            const Kokkos::mdspan<const double, Extents, Layout> read_only(raw);
            const auto borrowed = canary::as_indexed<canary::SpeciesIndex>(read_only);
            expect(borrowed[canary::SpeciesIndex {2}] == 12.0);
        };
        check.template operator()<Kokkos::layout_left>();
        check.template operator()<Kokkos::layout_right>();
        check.template operator()<Kokkos::layout_stride>();
    };

    "indexed views preserve accessor state and empty extents"_test = [] {
        std::array<double, 5> data {};
        using RawOffset = Kokkos::mdspan<double, Extents, Kokkos::layout_right, OffsetAccessor>;
        const RawOffset raw(data.data(), Extents(2), OffsetAccessor {2});
        const auto view = canary::as_indexed<canary::SpeciesIndex>(raw);
        static_assert(sizeof(view) == sizeof(raw));
        view[canary::SpeciesIndex {1}] = 7.0;
        expect(data[3] == 7.0);
        expect(view.base().accessor().offset_value == 2_u);
        expect(&view[canary::SpeciesIndex {0}] == &data[2]);
        const auto empty = canary::as_indexed<canary::ModeIndex>(Raw(nullptr, 0));
        expect(empty.extent(0) == 0_u);
        expect(empty.base().data_handle() == nullptr);
    };

    "indexed views compose quantity access with full scalar derivatives"_test = [] {
        const auto check = []<typename T>() {
            const std::array<T, 3> data {T(2, 3), T(99), T(5, 7)};
            using Mapping = Kokkos::layout_stride::mapping<Extents>;
            const Kokkos::mdspan<const T, Extents, Kokkos::layout_stride> raw(
              data.data(),
              Mapping(Extents(2), std::array<std::size_t, 1> {2}));
            const auto quantities =
              canary::as_kind<canary::kind::mode_temperature, canary::si>(raw);
            const auto view = canary::as_indexed<canary::ModeIndex>(quantities);
            static_assert(sizeof(view) == sizeof(quantities));
            static_assert(std::same_as<decltype(view[canary::ModeIndex {0}]),
                                       canary::ModeTemperature<T, canary::si>>);
            expect(view[canary::ModeIndex {1}].value().real() == 5);
            expect(view[canary::ModeIndex {1}].value().dual() == 7);
            expect(view.base()[0].value().dual() == 3);
        };
        check.template operator()<canary::test::Dual<float>>();
        check.template operator()<canary::test::Dual<double>>();
    };
}
