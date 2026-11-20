// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_INDEXED_VIEW_HPP
#define CANARY_INDEXED_VIEW_HPP

#include <canary/config.hpp>
#include <concepts>
#include <cstddef>
#include <mdspan/mdspan.hpp>
#include <type_traits>
#include <utility>

namespace canary
{
namespace detail
{
template <typename M> struct is_mdspan : std::false_type {};
template <typename E, typename X, typename L, typename A>
struct is_mdspan<Kokkos::mdspan<E, X, L, A>> : std::true_type {};

template <typename I>
concept ViewIndex = std::copyable<I> && requires(const I index) {
    { index.value() } -> std::unsigned_integral;
};

template <typename M>
concept RankOneMdspan = is_mdspan<M>::value && (M::rank() == 1);
}  // namespace detail

// Non-owning, rank-one mdspan indexed only by Index. Layout, accessor and
// element constness are retained; copying the view copies only its descriptor.
// Referenced storage and any accessor resources must outlive all uses of the
// view. Access requires index.value() < extent(0), representable by the mdspan's
// index_type. The consuming model checks its required extent at the operation
// boundary. A const descriptor can still write mutable elements, as with mdspan.
// Exception specifications follow the wrapped descriptor and access operations.
template <detail::ViewIndex Index, detail::RankOneMdspan Mdspan> class IndexedView
{
  public:
    using index_type = Index;
    using element_type = typename Mdspan::element_type;
    using value_type = typename Mdspan::value_type;
    using reference = typename Mdspan::reference;
    using extents_type = typename Mdspan::extents_type;
    using layout_type = typename Mdspan::layout_type;
    using accessor_type = typename Mdspan::accessor_type;

    CANARY_HOST_DEVICE constexpr explicit IndexedView(const Mdspan& view) noexcept(
      std::is_nothrow_copy_constructible_v<Mdspan>) : view_(view)
    {
    }

    template <std::same_as<Index> I>
    [[nodiscard]] CANARY_HOST_DEVICE constexpr reference operator[](I index) const
      noexcept(noexcept(std::declval<const Mdspan&>()[index.value()]))
    { return view_[index.value()]; }

    [[nodiscard]] CANARY_HOST_DEVICE static constexpr std::size_t rank() noexcept { return 1; }

    // r must be zero, as for the underlying rank-one mdspan.
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto extent(std::size_t r) const
      noexcept(noexcept(std::declval<const Mdspan&>().extent(r)))
    { return view_.extent(r); }

    // Explicit interoperability with interfaces that use native mdspan indices.
    // The returned descriptor borrows the same storage and accessor resources.
    [[nodiscard]] CANARY_HOST_DEVICE constexpr Mdspan base() const
      noexcept(std::is_nothrow_copy_constructible_v<Mdspan>)
    { return view_; }

  private:
    Mdspan view_;
};

// Apply an index domain to an existing rank-one mdspan without copying elements.
// Construct an mdspan over an array to borrow its storage. No implicit conversion
// to the underlying mdspan or to another index domain is provided.
template <detail::ViewIndex Index, detail::RankOneMdspan Mdspan>
[[nodiscard]] CANARY_HOST_DEVICE constexpr auto
  as_indexed(const Mdspan& view) noexcept(std::is_nothrow_copy_constructible_v<Mdspan>)
{ return IndexedView<Index, Mdspan> {view}; }

}  // namespace canary

#endif  // CANARY_INDEXED_VIEW_HPP
