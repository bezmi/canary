#ifndef CANARY_THERMO_HPP
#define CANARY_THERMO_HPP

#include <array>
#include <concepts>
#include <cstddef>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <mdspan/mdspan.hpp>
#include <tuple>
#include <type_traits>
#include <utility>

namespace canary
{

// A gas model: a trivially copyable struct that names its parameter type,
// its unit system and its state type, and (for mixtures) its per-cell
// mixture type:
//   using parameter_type = P;                        floating point
//   using units = U;                                 si or scaled
//   template <Scalar T> using State = ...;           ThermoStateOf<State<T>, T, U>
//   template <Scalar T> using Mixture = ...;         mixtures only
// Closures are free functions close(model, [mixture,] a, b, ...) found by
// argument-dependent lookup.
template <typename G>
concept GasModelType = std::is_trivially_copyable_v<G> && requires {
    typename G::parameter_type;
    typename G::units;
    requires std::floating_point<typename G::parameter_type>;
    requires UnitSystem<typename G::units>;
};

template <GasModelType G, Scalar T> using state_t = G::template State<T>;

// A model that reduces a composition to a per-cell mixture value before its
// closures.
template <typename G, typename T>
concept HasMixture = GasModelType<G> && Scalar<T> && requires { typename G::template Mixture<T>; };

template <GasModelType G, Scalar T>
    requires HasMixture<G, T>
using mixture_t = G::template Mixture<T>;

// An argument passed through thermo_state to a closure: a state variable,
// a mixture value, a view or an extra input, all passed by value.
template <typename A>
concept ClosureArgument = std::copyable<A>;

namespace detail
{

// Rank of a state variable, or -1 for any other argument. thermo_state passes
// state variables to close in increasing rank, so callers may give them in
// any order.
template <typename A> constexpr int state_rank()
{
    if constexpr (QuantityKind<A>) {
        using K = A::kind;
        if constexpr (std::same_as<K, kind::density>) {
            return 0;
        } else if constexpr (std::same_as<K, kind::pressure>) {
            return 1;
        } else if constexpr (std::same_as<K, kind::temperature>) {
            return 2;
        } else if constexpr (std::same_as<K, kind::specific_energy> ||
                             std::same_as<K, kind::energy_density>) {
            return 3;
        } else {
            return -1;
        }
    } else if constexpr (SpeciesRangeOf<A, kind::mode_energy_density> ||
                         SpeciesRangeOf<A, kind::mode_temperature>) {
        return 4;
    } else {
        return -1;
    }
}

}  // namespace detail

// A thermo_state argument that thermo_state reorders by rank.
template <typename A>
concept StateVariable = detail::state_rank<std::remove_cvref_t<A>>() >= 0;

// close(model, args...) can be called.
template <typename G, typename... Args>
concept CloseCallable =
  GasModelType<G> && requires(const G& gas, const Args&... args) { close(gas, args...); };

// Closure from density and energy per unit volume through the model's
// closure from density and energy per unit mass. A model may provide its
// own overload, which overload resolution prefers.
template <GasModelType G, Scalar T, UnitSystem U, ClosureArgument... Rest>
    requires CloseCallable<G, Density<T, U>, SpecificEnergy<T, U>, Rest...>
CANARY_HOST_DEVICE auto
  close(const G& gas, Density<T, U> rho, EnergyDensity<T, U> rho_e, const Rest&... rest)
{ return close(gas, rho, specific_energy(rho_e, rho), rest...); }

template <GasModelType G, Scalar T, UnitSystem U, ClosureArgument... Rest>
    requires HasMixture<G, T> &&
             CloseCallable<G, mixture_t<G, T>, Density<T, U>, SpecificEnergy<T, U>, Rest...>
CANARY_HOST_DEVICE auto close(const G& gas,
                              const mixture_t<G, T>& mix,
                              Density<T, U> rho,
                              EnergyDensity<T, U> rho_e,
                              const Rest&... rest)
{ return close(gas, mix, rho, specific_energy(rho_e, rho), rest...); }

namespace detail
{

// Order in which thermo_state passes its arguments to close: arguments before
// the first state variable keep their place, then the state variables in
// increasing rank (stable), then the other arguments in their given order.
template <ClosureArgument... Args>
constexpr std::array<std::size_t, sizeof...(Args)> argument_order()
{
    constexpr std::size_t n = sizeof...(Args);
    constexpr std::array<int, n> ranks {state_rank<Args>()...};
    std::array<std::size_t, n> order {};
    std::size_t first = n;
    for (std::size_t i = 0; i < n; ++i) {
        if (ranks[i] >= 0) {
            first = i;
            break;
        }
    }
    std::size_t k = 0;
    for (std::size_t i = 0; i < first; ++i) { order[k++] = i; }
    for (int rank = 0; rank <= 4; ++rank) {
        for (std::size_t i = first; i < n; ++i) {
            if (ranks[i] == rank) { order[k++] = i; }
        }
    }
    for (std::size_t i = first; i < n; ++i) {
        if (ranks[i] < 0) { order[k++] = i; }
    }
    return order;
}

template <std::size_t I, ClosureArgument... Args>
using nth_t = std::tuple_element_t<I, std::tuple<Args...>>;

template <std::size_t I, ClosureArgument First, ClosureArgument... Rest>
CANARY_HOST_DEVICE constexpr const auto& nth(const First& first, const Rest&... rest)
{
    if constexpr (I == 0) {
        return first;
    } else {
        return nth<I - 1>(rest...);
    }
}

template <typename G, typename Index, typename... Args> struct closes_sorted;

template <typename G, std::size_t... I, typename... Args>
struct closes_sorted<G, std::index_sequence<I...>, Args...> {
    static constexpr std::array<std::size_t, sizeof...(Args)> order = argument_order<Args...>();
    static constexpr bool value = CloseCallable<G, nth_t<order[I], Args...>...>;
};

template <GasModelType G, std::size_t... I, ClosureArgument... Args>
CANARY_HOST_DEVICE auto
  close_sorted(const G& gas, std::index_sequence<I...> /*indices*/, const Args&... args)
{
    constexpr std::array<std::size_t, sizeof...(Args)> order = argument_order<Args...>();
    return close(gas, nth<order[I]>(args...)...);
}

}  // namespace detail

// close can be called with the arguments in the order thermo_state uses.
template <typename G, typename... Args>
concept ClosesInAnyOrder =
  GasModelType<G> && (ClosureArgument<Args> && ...) &&
  detail::closes_sorted<G, std::index_sequence_for<Args...>, Args...>::value;

// State of a model from two state variables, given in either order, with the
// model's other per-cell inputs: thermo_state(gas, [mixture,] a, b[, modes][,
// extra]). Returns what the model's close returns, a ClosureOf its state.
template <GasModelType G, ClosureArgument... Args>
    requires ClosesInAnyOrder<G, Args...>
CANARY_HOST_DEVICE auto thermo_state(const G& gas, const Args&... args)
{ return detail::close_sorted(gas, std::index_sequence_for<Args...> {}, args...); }

// mixture(gas, fractions) can be called with the mass-fraction view V.
template <typename G, typename V>
concept MixesFromFractions =
  GasModelType<G> && SpeciesRangeOf<V, kind::mass_fraction> &&
  requires(const G& gas, const V& fractions) { mixture(gas, fractions); };

// mixture(gas, partial_densities, rho) can be called with the
// partial-density view V and the density Rho.
template <typename G, typename V, typename Rho>
concept MixesFromDensities =
  GasModelType<G> && SpeciesRangeOf<V, kind::density> && KindOf<Rho, kind::density> &&
  requires(const G& gas, const V& partial_densities, const Rho& rho) {
      mixture(gas, partial_densities, rho);
  };

namespace detail
{

template <typename K, ClosureArgument... Args> constexpr std::size_t count_of_kind()
{ return ((KindOf<Args, K> ? std::size_t {1} : std::size_t {0}) + ... + std::size_t {0}); }

// Index of the first argument of kind K; sizeof...(Args) if there is none.
template <typename K, ClosureArgument... Args> constexpr std::size_t index_of_kind()
{
    constexpr std::array<bool, sizeof...(Args)> is_kind {KindOf<Args, K>...};
    for (std::size_t i = 0; i < sizeof...(Args); ++i) {
        if (is_kind[i]) { return i; }
    }
    return sizeof...(Args);
}

template <GasModelType G, SpeciesRange V>
    requires MixesFromFractions<G, V>
using mixture_from_fractions_t =
  decltype(mixture(std::declval<const G&>(), std::declval<const V&>()));

template <GasModelType G, SpeciesRange V, typename Rho>
    requires MixesFromDensities<G, V, Rho>
using mixture_from_densities_t =
  decltype(mixture(std::declval<const G&>(), std::declval<const V&>(), std::declval<const Rho&>()));

}  // namespace detail

// thermo_state of a mixture model from mass fractions: the mixture value is
// computed first, then thermo_state(gas, mixture, args...).
template <GasModelType G, SpeciesRange V, ClosureArgument... Args>
    requires MixesFromFractions<G, V> &&
             ClosesInAnyOrder<G, detail::mixture_from_fractions_t<G, V>, Args...>
CANARY_HOST_DEVICE auto thermo_state(const G& gas, const V& fractions, const Args&... args)
{ return thermo_state(gas, mixture(gas, fractions), args...); }

// thermo_state of a mixture model from partial densities. The density must
// be one of the state variables.
template <GasModelType G, SpeciesRange V, ClosureArgument... Args>
    requires(detail::count_of_kind<kind::density, Args...>() == 1) &&
            MixesFromDensities<
              G,
              V,
              detail::nth_t<detail::index_of_kind<kind::density, Args...>(), Args...>> &&
            ClosesInAnyOrder<
              G,
              detail::mixture_from_densities_t<
                G,
                V,
                detail::nth_t<detail::index_of_kind<kind::density, Args...>(), Args...>>,
              Args...>
CANARY_HOST_DEVICE auto thermo_state(const G& gas, const V& partial_densities, const Args&... args)
{
    const auto& rho = detail::nth<detail::index_of_kind<kind::density, Args...>()>(args...);
    return thermo_state(gas, mixture(gas, partial_densities, rho), args...);
}

namespace detail
{
template <typename G, typename T, typename A, typename B>
concept ClosesFromAlone = requires(const G& gas, A a, B b) {
    { thermo_state(gas, a, b) } -> std::same_as<ClosureOf<state_t<G, T>>>;
};

template <typename G, typename T, typename A, typename B>
concept ClosesFromWithMixture = requires(const G& gas, const mixture_t<G, T>& mix, A a, B b) {
    { thermo_state(gas, mix, a, b) } -> std::same_as<ClosureOf<state_t<G, T>>>;
};

template <typename G, typename T>
concept PropertiesAlone = requires(const G& gas, const state_t<G, T>& s) {
    { properties(gas, s) } -> std::same_as<Properties<T, typename G::units>>;
};

template <typename G, typename T>
concept PropertiesWithMixture =
  requires(const G& gas, const mixture_t<G, T>& mix, const state_t<G, T>& s) {
      { properties(gas, mix, s) } -> std::same_as<Properties<T, typename G::units>>;
  };
}  // namespace detail

// thermo_state of model G with value type T accepts the state variables A
// and B (in either order), after the mixture value for a mixture model.
template <typename G, typename T, typename A, typename B>
concept ClosesFrom = GasModelType<G> && Scalar<T> && StateVariable<A> && StateVariable<B> &&
                     ((HasMixture<G, T> && detail::ClosesFromWithMixture<G, T, A, B>) ||
                      (!HasMixture<G, T> && detail::ClosesFromAlone<G, T, A, B>));

// A gas model usable with state value type T: closures from (rho, e),
// (rho, p), (rho, T) and (p, T), and properties at a state.
template <typename G, typename T>
concept GasModel =
  GasModelType<G> && ScalarWith<T, typename G::parameter_type> &&
  requires { typename G::template State<T>; } &&
  ThermoStateOf<state_t<G, T>, T, typename G::units> &&
  ClosesFrom<G, T, Density<T, typename G::units>, SpecificEnergy<T, typename G::units>> &&
  ClosesFrom<G, T, Density<T, typename G::units>, Pressure<T, typename G::units>> &&
  ClosesFrom<G, T, Density<T, typename G::units>, Temperature<T, typename G::units>> &&
  ClosesFrom<G, T, Pressure<T, typename G::units>, Temperature<T, typename G::units>> &&
  ((HasMixture<G, T> && detail::PropertiesWithMixture<G, T>) ||
   (!HasMixture<G, T> && detail::PropertiesAlone<G, T>));

// View of mass fractions with value type T, as passed to mixture().
template <Scalar T>
using MassFractionView = decltype(as_kind<kind::mass_fraction>(
  std::declval<Kokkos::mdspan<const T, Kokkos::dextents<std::size_t, 1>>>()));

// A gas model whose closures take a per-cell mixture value, which
// mixture(gas, fractions) computes from the composition.
template <typename G, typename T>
concept MixtureGasModel =
  GasModel<G, T> && HasMixture<G, T> && requires(const G& gas, MassFractionView<T> y) {
      { mixture(gas, y) } -> std::same_as<mixture_t<G, T>>;
  };

// View of mode energy densities with value type T, as passed to
// thermo_state.
template <Scalar T, UnitSystem U>
using ModeEnergyDensityView = decltype(as_kind<kind::mode_energy_density, U>(
  std::declval<Kokkos::mdspan<const T, Kokkos::dextents<std::size_t, 1>>>()));

// View of mode temperatures with value type T, as passed to thermo_state.
template <Scalar T, UnitSystem U>
using ModeTemperatureView = decltype(as_kind<kind::mode_temperature, U>(
  std::declval<Kokkos::mdspan<const T, Kokkos::dextents<std::size_t, 1>>>()));

// A rank-1 range of quantities of kind K with value type T in unit system U,
// of any layout and accessor; for example a ModeEnergyDensityView or a
// ModeTemperatureView. Models with internal energy modes take their mode
// inputs as such ranges.
template <typename V, typename K, typename T, typename U>
concept ModeRangeOf = SpeciesRangeOf<V, K> && Scalar<T> && UnitSystem<U> &&
                      std::same_as<typename range_quantity_t<V>::value_type, T> &&
                      std::same_as<typename range_quantity_t<V>::units, U>;

namespace detail
{
template <typename G, typename T>
concept ClosesWithModes = requires(const G& gas,
                                   Density<T, typename G::units> rho,
                                   EnergyDensity<T, typename G::units> rho_e,
                                   ModeEnergyDensityView<T, typename G::units> modes) {
    { thermo_state(gas, rho, rho_e, modes) } -> std::same_as<ClosureOf<state_t<G, T>>>;
};

template <typename G, typename T>
concept ClosesWithModesAndMixture = requires(const G& gas,
                                             const mixture_t<G, T>& mix,
                                             Density<T, typename G::units> rho,
                                             EnergyDensity<T, typename G::units> rho_e,
                                             ModeEnergyDensityView<T, typename G::units> modes) {
    { thermo_state(gas, mix, rho, rho_e, modes) } -> std::same_as<ClosureOf<state_t<G, T>>>;
};
}  // namespace detail

// A gas model with internal energy modes: its state has mode energies and
// mode temperatures, and it closes from density, energy density and the
// mode energy densities.
template <typename G, typename T>
concept MultiTemperatureGasModel =
  GasModel<G, T> &&
  requires(const state_t<G, T>& s) {
      { s.mode_e[0] } -> std::convertible_to<SpecificEnergy<T, typename G::units>>;
      { s.mode_temp[0] } -> std::convertible_to<ModeTemperature<T, typename G::units>>;
  } &&
  ((HasMixture<G, T> && detail::ClosesWithModesAndMixture<G, T>) ||
   (!HasMixture<G, T> && detail::ClosesWithModes<G, T>));

// Model G has a closure from (rho, e) that takes the extra per-cell input X
// after the state variables. Checked with the parameter type as value type.
template <typename G, typename X>
concept AcceptsInput =
  GasModelType<G> && ClosureArgument<X> &&
  ((HasMixture<G, typename G::parameter_type> &&
    CloseCallable<G,
                  mixture_t<G, typename G::parameter_type>,
                  Density<typename G::parameter_type, typename G::units>,
                  SpecificEnergy<typename G::parameter_type, typename G::units>,
                  X>) ||
   (!HasMixture<G, typename G::parameter_type> &&
    CloseCallable<G,
                  Density<typename G::parameter_type, typename G::units>,
                  SpecificEnergy<typename G::parameter_type, typename G::units>,
                  X>));

}  // namespace canary

#endif  // CANARY_THERMO_HPP
