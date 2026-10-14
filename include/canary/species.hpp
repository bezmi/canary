#ifndef CANARY_SPECIES_HPP
#define CANARY_SPECIES_HPP

#include <algorithm>
#include <array>
#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/status.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <optional>
#include <span>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace canary
{

// A species identity, separate from array sizes, reaction indices and mode
// indices. Construction requires value < UINT32_MAX; the containing model
// determines the species-count bound. Default construction denotes species 0.
class SpeciesIndex
{
  public:
    CANARY_HOST_DEVICE constexpr SpeciesIndex() noexcept {}
    CANARY_HOST_DEVICE constexpr explicit SpeciesIndex(std::size_t value) noexcept
    {
        CANARY_ASSERT(value < std::numeric_limits<std::uint32_t>::max());
        value_ = static_cast<std::uint32_t>(value);
    }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr std::uint32_t value() const noexcept
    { return value_; }
    friend CANARY_HOST_DEVICE constexpr bool operator==(SpeciesIndex,
                                                        SpeciesIndex) noexcept = default;

  private:
    std::uint32_t value_ {0};
};

// An optional species identity in one word. Default construction denotes
// absence. Dereferencing requires a present index; it does not check the
// containing model's species-count bound. No allocation or throwing access.
class OptionalSpeciesIndex
{
  public:
    CANARY_HOST_DEVICE constexpr OptionalSpeciesIndex() noexcept {}
    CANARY_HOST_DEVICE constexpr explicit OptionalSpeciesIndex(SpeciesIndex index) noexcept :
        value_(index.value())
    {
    }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr bool has_value() const noexcept
    { return value_ != std::numeric_limits<std::uint32_t>::max(); }
    CANARY_HOST_DEVICE constexpr explicit operator bool() const noexcept { return has_value(); }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr SpeciesIndex operator*() const noexcept
    {
        CANARY_ASSERT(has_value());
        return SpeciesIndex {value_};
    }
    friend CANARY_HOST_DEVICE constexpr bool operator==(OptionalSpeciesIndex,
                                                        OptionalSpeciesIndex) noexcept = default;

  private:
    std::uint32_t value_ {std::numeric_limits<std::uint32_t>::max()};
};

// Consecutive species identities in [first, last), in storage order. The
// iterator holds a native-size position and yields SpeciesIndex values.
// Bounds and the end sentinel hold uint32_t values. The endpoint may be
// UINT32_MAX; it is never yielded as a species identity.
class SpeciesIndexRange
{
  public:
    struct sentinel {
        std::uint32_t last {0};
    };

    class iterator
    {
      public:
        using value_type = SpeciesIndex;
        using difference_type = std::ptrdiff_t;
        using iterator_concept = std::forward_iterator_tag;

        CANARY_HOST_DEVICE constexpr iterator() noexcept {}
        [[nodiscard]] CANARY_HOST_DEVICE constexpr SpeciesIndex operator*() const noexcept
        { return SpeciesIndex {position_}; }
        CANARY_HOST_DEVICE constexpr iterator& operator++() noexcept
        {
            CANARY_ASSERT(position_ < std::numeric_limits<std::uint32_t>::max());
            ++position_;
            return *this;
        }
        CANARY_HOST_DEVICE constexpr iterator operator++(int) noexcept
        {
            const iterator previous = *this;
            ++*this;
            return previous;
        }
        friend CANARY_HOST_DEVICE constexpr bool operator==(iterator a, iterator b) noexcept
        { return a.position_ == b.position_; }
        friend CANARY_HOST_DEVICE constexpr bool operator==(iterator a, sentinel b) noexcept
        { return a.position_ >= b.last; }

      private:
        CANARY_HOST_DEVICE constexpr explicit iterator(std::uint32_t position) noexcept :
            position_(position)
        {
        }
        friend class SpeciesIndexRange;
        std::size_t position_ {0};
    };

    [[nodiscard]] CANARY_HOST_DEVICE constexpr iterator begin() const noexcept
    { return iterator {first_}; }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr sentinel end() const noexcept
    { return sentinel {last_}; }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr std::uint32_t size() const noexcept
    { return last_ - first_; }

  private:
    CANARY_HOST_DEVICE constexpr SpeciesIndexRange(std::uint32_t first, std::uint32_t last) noexcept
        : first_(first), last_(last)
    {
    }
    friend CANARY_HOST_DEVICE constexpr SpeciesIndexRange
      species_indices(std::size_t first, std::size_t last) noexcept;
    std::uint32_t first_;
    std::uint32_t last_;
};

// Species identities from first up to, excluding, last. Requires
// first <= last <= UINT32_MAX. Equal bounds produce an empty range.
// Bounds are converted once; iteration does not allocate or throw.
[[nodiscard]] CANARY_HOST_DEVICE constexpr SpeciesIndexRange
  species_indices(std::size_t first, std::size_t last) noexcept
{
    CANARY_ASSERT(first <= last && last <= std::numeric_limits<std::uint32_t>::max());
    return {static_cast<std::uint32_t>(first), static_cast<std::uint32_t>(last)};
}

// Host construction requires an explicit role, independent of species names
// and ordering. The numerical mixture supports at most one free-electron species.
enum class SpeciesRole { unspecified, heavy, electron };

// Nuclear composition of one species. Host construction requires nonempty,
// distinct element identifiers and positive counts. Free electrons have no
// nuclear elements; their charge is specified separately.
struct ElementCount {
    std::string_view element;
    unsigned count;
};

// Borrowed host metadata. Element entries and their identifiers must remain
// alive and unmodified until construction finishes. Missing elements indicate
// unknown composition; a present empty span specifies no nuclear elements.
// Reaction construction requires known composition for every participant.
struct SpeciesMetadataView {
    SpeciesRole role {SpeciesRole::unspecified};
    std::optional<std::span<const ElementCount>> elements {};
};

// Particle mass and charge of each species, in the model's unit system. Only
// the first n entries are used. Arithmetic requires 1 <= n <= Capacity
// and finite positive particle masses in those entries.
template <std::floating_point P, UnitSystem U, std::size_t Capacity>
    requires(Capacity > 0)
struct SpeciesData {
    static_assert(Capacity <= std::numeric_limits<std::uint32_t>::max(),
                  "species capacity exceeds the index range");
    int n;
    std::array<ParticleMass<P, U>, Capacity> mass;
    std::array<Charge<P, U>, Capacity> charge;
};

namespace detail
{
template <typename M> struct is_mdspan : std::false_type {};
template <typename E, typename X, typename L, typename A>
struct is_mdspan<Kokkos::mdspan<E, X, L, A>> : std::true_type {};
}  // namespace detail

// An mdspan accessor policy: the member types and functions that
// Kokkos::mdspan requires of one.
template <typename A>
concept MdspanAccessor =
  std::copyable<A> && requires(const A a, A::data_handle_type p, std::size_t i) {
      typename A::element_type;
      typename A::offset_policy;
      { a.access(p, i) } -> std::convertible_to<typename A::reference>;
      { a.offset(p, i) } -> std::convertible_to<typename A::offset_policy::data_handle_type>;
  };

// A Kokkos::mdspan of plain scalar values.
template <typename M>
concept ScalarMdspan = detail::is_mdspan<std::remove_cvref_t<M>>::value &&
                       Scalar<typename std::remove_cvref_t<M>::value_type>;

// Accessor that reads each element of an mdspan of plain values as a
// Quantity of kind K in unit system U. Elements are read-only; storage is
// not reinterpreted.
template <KindTag K, UnitTag U, MdspanAccessor Base>
    requires KindUnits<K, U>
struct KindAccessor {
    using value_type = std::remove_const_t<typename Base::element_type>;
    using element_type = const Quantity<K, value_type, U>;
    using reference = Quantity<K, value_type, U>;
    using data_handle_type = Base::data_handle_type;
    using offset_policy = KindAccessor<K, U, typename Base::offset_policy>;

    constexpr KindAccessor() = default;
    CANARY_HOST_DEVICE constexpr explicit KindAccessor(const Base& base) : base(base) {}

    // Conversion from the accessor of another base, as mdspan's offset_policy
    // requires.
    template <MdspanAccessor Other>
        requires std::constructible_from<Base, const Other&>
    CANARY_HOST_DEVICE constexpr explicit KindAccessor(const KindAccessor<K, U, Other>& other) :
        base(other.base)
    {
    }

    [[nodiscard]] CANARY_HOST_DEVICE constexpr reference access(data_handle_type p,
                                                                std::size_t i) const
    { return reference {base.access(p, i)}; }

    [[nodiscard]] CANARY_HOST_DEVICE constexpr offset_policy::data_handle_type
      offset(data_handle_type p, std::size_t i) const
    { return base.offset(p, i); }

    Base base;
};

// The same elements as v, read as quantities of kind K in unit system U
// (unitless for dimensionless kinds such as mass fractions). The extents
// and the layout of v are kept.
template <KindTag K, UnitTag U = unitless, ScalarMdspan M>
    requires KindUnits<K, U>
CANARY_HOST_DEVICE constexpr auto as_kind(const M& v)
{
    using accessor = KindAccessor<K, U, typename M::accessor_type>;
    return Kokkos::mdspan<typename accessor::element_type,
                          typename M::extents_type,
                          typename M::layout_type,
                          accessor>(v.data_handle(), v.mapping(), accessor {v.accessor()});
}

// A rank-1 range whose elements read as quantities.
template <typename V>
concept SpeciesRange = requires(const V v, std::size_t i) {
    requires std::remove_cvref_t<V>::rank() == 1;
    { v.extent(0) } -> std::convertible_to<std::size_t>;
    { v[i] } -> QuantityKind;
};

// A rank-1 range whose elements read as quantities of kind K.
template <typename V, typename K>
concept SpeciesRangeOf = SpeciesRange<V> && KindTag<K> && requires(const V v, std::size_t i) {
    { v[i] } -> KindOf<K>;
};

template <typename V>
using range_quantity_t = std::remove_cvref_t<decltype(std::declval<const V&>()[std::size_t {}])>;

// A species range with unit system U and a value type that combines with
// the parameter type P and can be clamped.
template <typename V, typename K, typename U, typename P>
concept SpeciesWeights =
  SpeciesRangeOf<V, K> && std::same_as<typename range_quantity_t<V>::units, U> &&
  ScalarPair<typename range_quantity_t<V>::value_type, P> &&
  HasMin<typename range_quantity_t<V>::value_type> &&
  HasMax<typename range_quantity_t<V>::value_type>;

namespace detail
{

// Value of weight i clamped to [0, total]. The weights are mass fractions
// or partial densities; only their ratios to the sum of the clamped weights
// are used, so their value parts enter the mixing rules directly.
template <SpeciesRange V, Scalar T>
    requires std::same_as<typename range_quantity_t<V>::value_type, T> && HasMin<T> && HasMax<T>
CANARY_HOST_DEVICE constexpr T clamped_weight(const V& weights, SpeciesIndex index, T total)
{
    using Value = ScalarTraits<T>::value_type;
    return math::min(math::max(weights[index.value()].value(), T(Value(0))), total);
}

// Sums over the species of the clamped weights w, of w / m and of
// (w / m) q, where w is each weight clamped to [0, total].
template <Scalar W, Scalar N, Scalar C> struct WeightSums {
    W weight;
    N number;
    C charge;
};

template <std::floating_point P, UnitSystem U, std::size_t Capacity, SpeciesRange V, Scalar T>
CANARY_HOST_DEVICE constexpr auto
  weight_sums(const SpeciesData<P, U, Capacity>& species, const V& weights, T total)
{
    // NOLINTNEXTLINE(misc-static-assert): run-time species count and extent
    CANARY_ASSERT(species.n >= 1 && static_cast<std::size_t>(species.n) <= Capacity &&
                  static_cast<std::size_t>(species.n) ==
                    static_cast<std::size_t>(weights.extent(0)));
    // Number weight w / m and charge weight (w / m) q of one species.
    using Number = decltype(std::declval<T>() / std::declval<P>());
    using ChargeWeight = decltype(std::declval<Number>() * std::declval<P>());
    const T w0 = clamped_weight(weights, SpeciesIndex {0}, total);
    const Number n0 = w0 / species.mass[0].value();
    WeightSums<T, Number, ChargeWeight> sums {
      .weight = w0,
      .number = n0,
      .charge = n0 * species.charge[0].value(),
    };
    for (const SpeciesIndex i : species_indices(1, static_cast<std::size_t>(species.n))) {
        const T w = clamped_weight(weights, i, total);
        const Number n = w / species.mass[i.value()].value();
        sums.weight = sums.weight + w;
        sums.number = sums.number + n;
        sums.charge = sums.charge + n * species.charge[i.value()].value();
    }
    return sums;
}

// Index of name in names; absent when it is not there. Host-side.
template <std::size_t N>
OptionalSpeciesIndex species_index(const std::array<std::string_view, N>& names,
                                   std::string_view name)
{
    static_assert(N <= std::numeric_limits<std::uint32_t>::max(),
                  "species capacity exceeds the index range");
    for (const SpeciesIndex i : species_indices(0, N)) {
        if (names[i.value()] == name) { return OptionalSpeciesIndex {i}; }
    }
    return {};
}

template <SpeciesRange V, Scalar T>
CANARY_HOST_DEVICE constexpr Status composition_status(const V& weights, T total)
{
    using Value = ScalarTraits<T>::value_type;
    const auto n = static_cast<std::size_t>(weights.extent(0));
    if (n > std::numeric_limits<std::uint32_t>::max()) { return Status::Invalid; }
    const Value total_value = ScalarTraits<T>::value(total);
    if (!std::isfinite(total_value) || !(total_value > Value(0))) { return Status::Invalid; }
    const Value bound =
      Value(128) * std::numeric_limits<Value>::epsilon() * static_cast<Value>(n) * total_value;
    Value sum = 0;
    Value clamped_sum = 0;
    bool finite = true;
    bool bounded = true;
    for (const SpeciesIndex i : species_indices(0, n)) {
        const Value w = ScalarTraits<T>::value(weights[i.value()].value());
        finite = finite && std::isfinite(w);
        bounded = bounded && w >= -bound && w <= total_value + bound;
        sum += w;
        clamped_sum += std::clamp(w, Value(0), total_value);
    }
    if (!finite || !(clamped_sum > Value(0))) { return Status::Invalid; }
    if (!bounded || std::abs(sum - total_value) > bound) { return Status::CompositionViolation; }
    return Status::Valid;
}

}  // namespace detail

// Mean particle mass of a mixture, sum(w) / sum(w / m_i), with each weight
// clamped to [0, total] (total 1 for mass fractions, rho for partial
// densities). Normalising by the sum of the clamped weights keeps the
// result a weighted mean when the weights do not sum to the total. Check
// composition_status(fractions) before calling and reject Status::Invalid.
// For example, all-zero fractions or a nonfinite fraction are Invalid and
// must not reach this division. The arithmetic does not report a status.
template <std::floating_point P, UnitSystem U, std::size_t Capacity, typename V>
    requires SpeciesWeights<V, kind::mass_fraction, unitless, P>
CANARY_HOST_DEVICE constexpr auto mean_particle_mass(const SpeciesData<P, U, Capacity>& species,
                                                     const V& fractions)
{
    using T = range_quantity_t<V>::value_type;
    const auto sums =
      detail::weight_sums(species, fractions, T(typename ScalarTraits<T>::value_type(1)));
    const auto m = sums.weight / sums.number;
    return ParticleMass<std::remove_const_t<decltype(m)>, U> {m};
}

// Check composition_status(partial_densities, rho) before calling and reject
// Status::Invalid, including zero total density and all-zero partial densities.
template <std::floating_point P, UnitSystem U, std::size_t Capacity, typename V, Scalar T>
    requires SpeciesWeights<V, kind::density, U, P> &&
             std::same_as<typename range_quantity_t<V>::value_type, T>
CANARY_HOST_DEVICE constexpr auto mean_particle_mass(const SpeciesData<P, U, Capacity>& species,
                                                     const V& partial_densities,
                                                     Density<T, U> rho)
{
    const auto sums = detail::weight_sums(species, partial_densities, rho.value());
    const auto m = sums.weight / sums.number;
    return ParticleMass<std::remove_const_t<decltype(m)>, U> {m};
}

// Mean particle charge, sum((w / m_i) q_i) / sum(w / m_i), with the same
// clamped weights. Check composition_status(fractions) before calling and
// reject Status::Invalid. All-zero or nonfinite fractions must not reach
// this division. The arithmetic does not report a status.
template <std::floating_point P, UnitSystem U, std::size_t Capacity, typename V>
    requires SpeciesWeights<V, kind::mass_fraction, unitless, P>
CANARY_HOST_DEVICE constexpr auto mean_charge(const SpeciesData<P, U, Capacity>& species,
                                              const V& fractions)
{
    using T = range_quantity_t<V>::value_type;
    const auto sums =
      detail::weight_sums(species, fractions, T(typename ScalarTraits<T>::value_type(1)));
    const auto q = sums.charge / sums.number;
    return Charge<std::remove_const_t<decltype(q)>, U> {q};
}

// Check composition_status(partial_densities, rho) before calling and reject
// Status::Invalid, including zero total density and all-zero partial densities.
template <std::floating_point P, UnitSystem U, std::size_t Capacity, typename V, Scalar T>
    requires SpeciesWeights<V, kind::density, U, P> &&
             std::same_as<typename range_quantity_t<V>::value_type, T>
CANARY_HOST_DEVICE constexpr auto mean_charge(const SpeciesData<P, U, Capacity>& species,
                                              const V& partial_densities,
                                              Density<T, U> rho)
{
    const auto sums = detail::weight_sums(species, partial_densities, rho.value());
    const auto q = sums.charge / sums.number;
    return Charge<std::remove_const_t<decltype(q)>, U> {q};
}

// Validity of a composition, from the value parts:
// - Invalid: the species count exceeds UINT32_MAX, a weight or the total
//   is not finite, the total is not positive, or every clamped weight is zero.
// - CompositionViolation: a weight lies outside [-b, total + b], or the
//   weights sum to the total only within more than b, where
//   b = 128 * epsilon * n * total.
// - Valid otherwise.
template <typename V>
    requires SpeciesRangeOf<V, kind::mass_fraction>
CANARY_HOST_DEVICE constexpr Status composition_status(const V& fractions)
{
    using T = range_quantity_t<V>::value_type;
    return detail::composition_status(fractions, T(typename ScalarTraits<T>::value_type(1)));
}

template <typename V, Scalar T, UnitSystem U>
    requires SpeciesRangeOf<V, kind::density> &&
             std::same_as<typename range_quantity_t<V>::value_type, T> &&
             std::same_as<typename range_quantity_t<V>::units, U>
CANARY_HOST_DEVICE constexpr Status composition_status(const V& partial_densities,
                                                       Density<T, U> rho)
{ return detail::composition_status(partial_densities, rho.value()); }

namespace detail
{

struct OwnedElementCount {
    OwnedElementName element;
    unsigned count;
};

template <> struct InputStorage<SpeciesMetadataView> {
    SpeciesRole role {};
    std::optional<std::vector<OwnedElementCount>> elements;
    using Scratch = std::vector<ElementCount>;
    InputStorage() = default;
    explicit InputStorage(const SpeciesMetadataView& input) : role(input.role)
    {
        if (input.elements) {
            elements.emplace();
            elements->reserve(input.elements->size());
            for (const ElementCount& item : *input.elements) {
                elements->push_back({OwnedElementName {item.element}, item.count});
            }
        }
    }
    SpeciesMetadataView borrow(Scratch& scratch) const
    {
        SpeciesMetadataView input {.role = role};
        if (elements) {
            scratch.reserve(elements->size());
            for (const OwnedElementCount& item : *elements) {
                scratch.push_back({*item.element.value, item.count});
            }
            input.elements = std::span<const ElementCount> {scratch};
        }
        return input;
    }
};

}  // namespace detail

}  // namespace canary

#endif  // CANARY_SPECIES_HPP
