#ifndef CANARY_REACTIONS_HPP
#define CANARY_REACTIONS_HPP

// Finite-rate chemistry of a mixture: a mechanism of reversible reactions
// with Arrhenius forward rates, backward rates from the equilibrium constant
// or from their own Arrhenius fits, and optional third bodies; per cell, the
// rates of progress of every reaction and the mass production rate of every
// species, at one temperature, at two (the translational temperature T and
// the vibrational-electronic temperature T_ve of a two-temperature model) or
// at T and the M mode temperatures T_m of a multi-temperature model.
//
// Concentrations are number densities n_i = rho y_i / m_i (particles per
// unit volume), so that the same functions work in SI and in scaled units.
// T_l is T clamped to the mechanism's temperature limits. Each rate has a
// second temperature T_o: the primary T_l or the internal temperature selected
// by its physical requirement. Internal temperatures are not clamped. Explicit
// single-temperature evaluation collapses every requirement to T_l.
// An Arrhenius rate constant with
// translational exponent s is evaluated at the rate-controlling temperature
// T_c = T_l^s T_o^(1 - s):
//   k(T_c) = A T_c^n exp(-theta / T_c).
// For a reaction with reactant coefficients nu'_i and product coefficients
// nu''_i:
//   k_f = k(T_c) with the forward A, n, theta and s,
//   [M] = sum_i eps_i n_i with a third body, 1 without,
//   q_f = k_f [M] prod_i n_i^nu'_i,   q_b = k_b [M] prod_i n_i^nu''_i,
//   omega_i = m_i sum_r (nu''_ir - nu'_ir) (q_f,r - q_b,r).
// The backward rate constant is either an Arrhenius rate of its own, at its
// own T_c, or k_b = k_f(T_eq) / K_c(T_eq) at the equilibrium temperature
// T_eq selected by the independent equilibrium_temperature requirement. k_f(T_eq) is the forward
// Arrhenius expression at T_c = T_eq, and K_c is the equilibrium constant in concentration form
//   K_c = exp(-sum_i dnu_i g_i m_i / (k_B T_eq)) (p0 / (k_B T_eq))^(sum_i dnu_i),
// where dnu_i = nu''_i - nu'_i, g_i is the standard-state Gibbs energy per
// unit mass of species i at T_eq from the species_properties of the thermo
// model, and p0 the standard pressure. With T_eq = T_l, k_b does not depend
// on T_o.

#include <algorithm>
#include <array>
#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/indexed_view.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/reaction_species.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/thermo.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iterator>
#include <limits>
#include <map>
#include <mdspan/mdspan.hpp>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <variant>

namespace canary
{

// A reaction identity, separate from species and mode identities and
// array extents. The containing mechanism determines the reaction-count bound.
// Default construction denotes reaction 0; value() is the storage-access boundary.
class ReactionIndex
{
  public:
    CANARY_HOST_DEVICE constexpr ReactionIndex() noexcept {}
    CANARY_HOST_DEVICE constexpr explicit ReactionIndex(std::size_t value) noexcept : value_(value)
    {
    }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr std::size_t value() const noexcept { return value_; }
    friend CANARY_HOST_DEVICE constexpr bool operator==(ReactionIndex,
                                                        ReactionIndex) noexcept = default;

  private:
    std::size_t value_ {0};
};

// Consecutive reaction identities in [first, last). The endpoint may be SIZE_MAX;
// it is never yielded. Iteration retains a native-size position.
class ReactionIndexRange
{
  public:
    struct sentinel {
        std::size_t last {0};
    };
    class iterator
    {
      public:
        using value_type = ReactionIndex;
        using difference_type = std::ptrdiff_t;
        using iterator_concept = std::forward_iterator_tag;
        CANARY_HOST_DEVICE constexpr iterator() noexcept {}
        [[nodiscard]] CANARY_HOST_DEVICE constexpr ReactionIndex operator*() const noexcept
        { return ReactionIndex {position_}; }
        CANARY_HOST_DEVICE constexpr iterator& operator++() noexcept
        {
            CANARY_ASSERT(position_ < std::numeric_limits<std::size_t>::max());
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
        CANARY_HOST_DEVICE constexpr explicit iterator(std::size_t position) noexcept :
            position_(position)
        {
        }
        friend class ReactionIndexRange;
        std::size_t position_ {0};
    };
    [[nodiscard]] CANARY_HOST_DEVICE constexpr iterator begin() const noexcept
    { return iterator {first_}; }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr sentinel end() const noexcept { return {last_}; }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr std::size_t size() const noexcept
    { return last_ - first_; }

  private:
    CANARY_HOST_DEVICE constexpr ReactionIndexRange(std::size_t first, std::size_t last) noexcept :
        first_(first), last_(last)
    {
    }
    friend CANARY_HOST_DEVICE constexpr ReactionIndexRange
      reaction_indices(std::size_t first, std::size_t last) noexcept;
    std::size_t first_;
    std::size_t last_;
};

// Reaction identities from first up to, excluding, last. Requires first <= last.
// Equal bounds produce an empty range. No allocation or throwing access.
[[nodiscard]] CANARY_HOST_DEVICE constexpr ReactionIndexRange
  reaction_indices(std::size_t first, std::size_t last) noexcept
{
    CANARY_ASSERT(first <= last);
    return {first, last};
}

// Physical temperature requirements are host inputs. Construction resolves them
// against the model's component assignments; no variant enters numerical data.
struct TranslationTemperatureInput {};
struct RotationTemperatureInput {};
struct VibrationTemperatureInput {
    const char* species;
};
struct BoundElectronicTemperatureInput {
    const char* species;
};
struct FreeElectronTemperatureInput {};
using PhysicalTemperatureInput = std::variant<TranslationTemperatureInput,
                                              RotationTemperatureInput,
                                              VibrationTemperatureInput,
                                              BoundElectronicTemperatureInput,
                                              FreeElectronTemperatureInput>;

// SI per mole: k = A T_c^n exp(-theta/T_c), T_c = T^s T_o^(1-s).
// A has units (m^3/mol)^(order-1)/s, theta and T_c are in K, n and s
// are dimensionless. A model-aware builder requires other_temperature when
// s != 1, and validates every supplied requirement, including unused ones.
struct ArrheniusInput {
    double A;
    double n;
    double activation_temperature;
    double translational_exponent = 1.0;
    std::optional<PhysicalTemperatureInput> other_temperature {};
};

// Host input of a species name and a number: a stoichiometric coefficient
// or a third-body efficiency.
struct SpeciesCoefficient {
    const char* name;
    double value;
};

// Host input of up to K species coefficients; only the first count entries
// are used.
template <std::size_t K>
    requires(K > 0)
struct SpeciesCoefficients {
    int count;
    std::array<SpeciesCoefficient, K> entries;
};

// How the backward rate constant of a reaction is found.
enum class BackwardRateModel : std::uint8_t {
    // k_b = k_f / K_c.
    equilibrium,
    // An Arrhenius rate constant of its own.
    arrhenius,
};

// Host input of one reversible reaction: its equation (for error messages),
// the reactant and product coefficients (positive integers), the forward
// rate, the backward rate model and, for BackwardRateModel::arrhenius, the
// backward rate. Equilibrium rates use the independent equilibrium_temperature;
// explicit backward Arrhenius rates forbid it. Also specifies whether a third
// body takes part, with the efficiencies
// of the listed species (species not listed have efficiency 1). With no
// third body, efficiencies must be empty.
template <std::size_t K>
    requires(K > 0)
struct ReactionInput {
    const char* equation;
    SpeciesCoefficients<K> reactants;
    SpeciesCoefficients<K> products;
    ArrheniusInput forward;
    BackwardRateModel backward_model;
    ArrheniusInput backward;
    std::optional<PhysicalTemperatureInput> equilibrium_temperature {};
    bool third_body;
    SpeciesCoefficients<K> efficiencies;
};

// Host input of a mechanism of R reactions: the temperature limits of rate
// evaluation in K and the reactions, each with at most K species per list.
template <std::size_t R, std::size_t K>
    requires(R > 0 && K > 0)
struct MechanismInput {
    double temperature_lower;
    double temperature_upper;
    std::array<ReactionInput<K>, R> reactions;
};

// Arrhenius rate constant k = A T_c^n exp(-theta / T_c) in model units: T_c
// and theta in the temperature unit of U, and A per particle in
// (length^3 / count)^(order - 1) / time of U, times the temperature scale
// raised to n. A and n are plain values, since the dimension of A depends
// on the reaction order. translational_exponent is s of
// T_c = T^s T_o^(1 - s); other_mode is the index of the energy mode whose
// temperature is T_o. Absence denotes the primary temperature, including
// every requirement collapsed by explicit single-temperature construction.
template <std::floating_point P, UnitSystem U> struct ArrheniusRate {
    P A;
    P n;
    Temperature<P, U> activation_temperature;
    P translational_exponent;
    OptionalModeIndex other_mode;
};

// One reaction of a mechanism over N species, indexed as the species of the
// thermo model: reactant and product coefficients, the sum of the
// coefficient changes sum_i (products[i] - reactants[i]), the forward rate,
// the backward rate model and rate (used with BackwardRateModel::arrhenius),
// whether an equilibrium backward rate is evaluated at the T_o of
// backward.other_mode instead of T_l, whether a third body takes part and
// the efficiency of every species (used with a third body).
template <std::floating_point P, UnitSystem U, std::size_t N>
    requires(N > 0)
struct Reaction {
  private:
    std::array<int, N> reactants_;
    std::array<int, N> products_;

  public:
    int coefficient_change;
    ArrheniusRate<P, U> forward;
    BackwardRateModel backward_model;
    ArrheniusRate<P, U> backward;
    bool equilibrium_at_other_temperature;
    bool third_body;

  private:
    std::array<P, N> efficiencies_;

  public:
    // Views borrow this record. Reacquire them after moving or replacing the owner.
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto reactants() & noexcept
    {
        return as_indexed<SpeciesIndex>(
          Kokkos::mdspan<int, Kokkos::extents<std::size_t, N>>(reactants_.data()));
    }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto reactants() const& noexcept
    {
        return as_indexed<SpeciesIndex>(
          Kokkos::mdspan<const int, Kokkos::extents<std::size_t, N>>(reactants_.data()));
    }
    auto reactants() && = delete;
    auto reactants() const&& = delete;
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto products() & noexcept
    {
        return as_indexed<SpeciesIndex>(
          Kokkos::mdspan<int, Kokkos::extents<std::size_t, N>>(products_.data()));
    }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto products() const& noexcept
    {
        return as_indexed<SpeciesIndex>(
          Kokkos::mdspan<const int, Kokkos::extents<std::size_t, N>>(products_.data()));
    }
    auto products() && = delete;
    auto products() const&& = delete;
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto efficiencies() & noexcept
    {
        return as_indexed<SpeciesIndex>(
          Kokkos::mdspan<P, Kokkos::extents<std::size_t, N>>(efficiencies_.data()));
    }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto efficiencies() const& noexcept
    {
        return as_indexed<SpeciesIndex>(
          Kokkos::mdspan<const P, Kokkos::extents<std::size_t, N>>(efficiencies_.data()));
    }
    auto efficiencies() && = delete;
    auto efficiencies() const&& = delete;
};

template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
    requires(N > 0 && R > 0)
struct ReactionMechanism;

namespace detail
{
// Behavior is fixed after parameter conversion. Temperature identities are
// separate from the arithmetic selected for a homogeneous reaction block.
enum class RateBehavior : std::uint8_t {
    ArrheniusTranslationTranslation,
    ArrheniusTranslationTranslationThirdBody,
    ArrheniusTranslationOther,
    ArrheniusTranslationOtherThirdBody,
    ArrheniusOtherTranslation,
    ArrheniusOtherTranslationThirdBody,
    ArrheniusOtherOther,
    ArrheniusOtherOtherThirdBody,
    EquilibriumForwardTranslation,
    EquilibriumForwardTranslationThirdBody,
    EquilibriumForwardOther,
    EquilibriumForwardOtherThirdBody,
};
struct RateBlock {
    std::size_t end {};
    OptionalModeIndex equilibrium_mode {};
    RateBehavior behavior {};
};
template <std::size_t R> struct RateSchedule {
    std::array<ReactionIndex, R> order {};
    std::array<RateBlock, R> blocks {};
    std::size_t block_count {};
    std::size_t explicit_end {};
    std::size_t translation_end {};
};
template <std::size_t R> struct RateModeReferences {
    std::array<ModeIndex, R> forward {};
    std::array<ModeIndex, R> backward_only {};
    std::size_t forward_count {};
    std::size_t backward_only_count {};
};
template <std::size_t M, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
CANARY_HOST_DEVICE constexpr const RateSchedule<R>&
  reaction_schedule(const ReactionMechanism<P, U, N, R>& mechanism) noexcept;
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
CANARY_HOST_DEVICE constexpr const RateModeReferences<R>&
  reaction_mode_references(const ReactionMechanism<P, U, N, R>& mechanism) noexcept;

template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          std::size_t K,
          typename Source,
          typename OtherMode>
ReactionMechanism<P, U, N, R> build_mechanism(const MechanismInput<R, K>& input,
                                              const ReactionSpeciesView<Source, N>& species,
                                              const UnitScales<P, U>& scales,
                                              Pressure<P, U> standard_pressure,
                                              const OtherMode& other_mode);
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
CANARY_HOST_DEVICE constexpr std::size_t
  reaction_mode_count(const ReactionMechanism<P, U, N, R>& mechanism) noexcept;
}  // namespace detail

// R reactions over the N species of a thermo model, in unit system U, as
// derived by make_reaction_mechanism: the reactions, the particle mass of
// every species, the temperature limits of rate evaluation, the standard
// pressure p0 of the equilibrium constants and the Boltzmann constant.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
    requires(N > 0 && R > 0)
struct ReactionMechanism {
    static constexpr std::size_t reaction_count = R;
    using parameter_type = P;
    using units = U;

  private:
    std::array<Reaction<P, U, N>, R> reactions_;
    std::array<ParticleMass<P, U>, N> mass_;
    std::size_t configured_mode_count_ {};
    detail::RateSchedule<R> schedule_ {};
    detail::RateModeReferences<R> controlling_modes_ {};

  public:
    Temperature<P, U> temperature_lower;
    Temperature<P, U> temperature_upper;
    Pressure<P, U> standard_pressure;
    BoltzmannConstant<P, U> boltzmann;

    // Mechanism data are read-only through borrowed views. Reacquire views after
    // moving or replacing the owner; construction validates all stored reactions.
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto reactions() const& noexcept
    {
        return as_indexed<ReactionIndex>(
          Kokkos::mdspan<const Reaction<P, U, N>, Kokkos::extents<std::size_t, R>>(
            reactions_.data()));
    }
    auto reactions() && = delete;
    auto reactions() const&& = delete;
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto mass() const& noexcept
    {
        return as_indexed<SpeciesIndex>(
          Kokkos::mdspan<const ParticleMass<P, U>, Kokkos::extents<std::size_t, N>>(mass_.data()));
    }
    auto mass() && = delete;
    auto mass() const&& = delete;

  private:
    ReactionMechanism() = default;
    template <std::size_t M, std::floating_point Q, UnitSystem V, std::size_t Ns, std::size_t Rs>
    friend CANARY_HOST_DEVICE constexpr const detail::RateSchedule<Rs>&
      detail::reaction_schedule(const ReactionMechanism<Q, V, Ns, Rs>&) noexcept;
    template <std::floating_point Q, UnitSystem V, std::size_t Ns, std::size_t Rs>
    friend CANARY_HOST_DEVICE constexpr const detail::RateModeReferences<Rs>&
      detail::reaction_mode_references(const ReactionMechanism<Q, V, Ns, Rs>&) noexcept;
    template <std::floating_point Q, UnitSystem V, std::size_t Ns, std::size_t Rs>
    friend CANARY_HOST_DEVICE constexpr std::size_t
      detail::reaction_mode_count(const ReactionMechanism<Q, V, Ns, Rs>&) noexcept;
    template <std::floating_point Q,
              UnitSystem V,
              std::size_t Ns,
              std::size_t Rs,
              std::size_t K,
              typename Source,
              typename OtherMode>
    friend ReactionMechanism<Q, V, Ns, Rs>
      detail::build_mechanism(const MechanismInput<Rs, K>& input,
                              const ReactionSpeciesView<Source, Ns>& species,
                              const UnitScales<Q, V>& scales,
                              Pressure<Q, V> standard_pressure,
                              const OtherMode& other_mode);
};

namespace detail
{

template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
CANARY_HOST_DEVICE constexpr std::size_t
  reaction_mode_count(const ReactionMechanism<P, U, N, R>& mechanism) noexcept
{ return mechanism.configured_mode_count_; }

template <std::size_t M, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
CANARY_HOST_DEVICE constexpr const RateSchedule<R>&
  reaction_schedule(const ReactionMechanism<P, U, N, R>& mechanism) noexcept
{ return mechanism.schedule_; }
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
CANARY_HOST_DEVICE constexpr const RateModeReferences<R>&
  reaction_mode_references(const ReactionMechanism<P, U, N, R>& mechanism) noexcept
{ return mechanism.controlling_modes_; }

template <std::floating_point P, UnitSystem U>
bool internal_rate_control(const ArrheniusRate<P, U>& rate) noexcept
{ return rate.other_mode.has_value() && rate.translational_exponent != P(1); }

template <std::floating_point P, UnitSystem U, std::size_t N>
RateBehavior rate_behavior(const Reaction<P, U, N>& r) noexcept
{
    const bool forward = internal_rate_control(r.forward);
    if (r.backward_model == BackwardRateModel::equilibrium) {
        if (forward) {
            return r.third_body ? RateBehavior::EquilibriumForwardOtherThirdBody
                                : RateBehavior::EquilibriumForwardOther;
        }
        return r.third_body ? RateBehavior::EquilibriumForwardTranslationThirdBody
                            : RateBehavior::EquilibriumForwardTranslation;
    }
    const bool backward = internal_rate_control(r.backward);
    if (forward && backward) {
        return r.third_body ? RateBehavior::ArrheniusOtherOtherThirdBody
                            : RateBehavior::ArrheniusOtherOther;
    }
    if (forward) {
        return r.third_body ? RateBehavior::ArrheniusOtherTranslationThirdBody
                            : RateBehavior::ArrheniusOtherTranslation;
    }
    if (backward) {
        return r.third_body ? RateBehavior::ArrheniusTranslationOtherThirdBody
                            : RateBehavior::ArrheniusTranslationOther;
    }
    return r.third_body ? RateBehavior::ArrheniusTranslationTranslationThirdBody
                        : RateBehavior::ArrheniusTranslationTranslation;
}

// The total sort key preserves original order within equal source/behavior
// groups. Only indices are permuted; coefficients remain in public record order.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
RateSchedule<R> prepare_rate_schedule(const std::array<Reaction<P, U, N>, R>& records)
{
    RateSchedule<R> out {};
    const auto key = [&](ReactionIndex i) {
        const auto& r = records[i.value()];
        const bool equilibrium = r.backward_model == BackwardRateModel::equilibrium;
        const auto mode = equilibrium ? r.backward.other_mode : OptionalModeIndex {};
        return std::tuple {equilibrium ? (mode ? 2 : 1) : 0,
                           mode ? (*mode).value() : 0U,
                           rate_behavior(r),
                           i.value()};
    };
    for (const ReactionIndex i : reaction_indices(0, R)) { out.order[i.value()] = i; }
    std::sort(out.order.begin(), out.order.end(), [&](ReactionIndex a, ReactionIndex b) {
        return key(a) < key(b);
    });
    for (std::size_t position = 0; position < R; ++position) {
        const auto& r = records[out.order[position].value()];
        const bool equilibrium = r.backward_model == BackwardRateModel::equilibrium;
        const auto mode = equilibrium ? r.backward.other_mode : OptionalModeIndex {};
        const auto behavior = rate_behavior(r);
        if (out.block_count == 0 || out.blocks[out.block_count - 1].behavior != behavior ||
            out.blocks[out.block_count - 1].equilibrium_mode != mode) {
            out.blocks[out.block_count++] = {.end = position + 1,
                                             .equilibrium_mode = mode,
                                             .behavior = behavior};
        } else {
            out.blocks[out.block_count - 1].end = position + 1;
        }
        if (!equilibrium) { out.explicit_end = out.block_count; }
        if (!equilibrium || !mode) { out.translation_end = out.block_count; }
    }
    return out;
}

// Modes of the M configured modes, in ascending order, that control a forward
// Arrhenius rate, and those that control only an explicit backward rate. Each
// listed mode controls at least one of the R reactions, so neither count
// exceeds R.
template <std::size_t M, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
RateModeReferences<R> prepare_rate_mode_references(const std::array<Reaction<P, U, N>, R>& records)
{
    RateModeReferences<R> out {};
    for (const ModeIndex m : mode_indices(0, M)) {
        const auto controls = [&](const ArrheniusRate<P, U>& rate) {
            return internal_rate_control(rate) && *rate.other_mode == m;
        };
        if (std::ranges::any_of(records, [&](const auto& r) { return controls(r.forward); })) {
            out.forward[out.forward_count++] = m;
        } else if (std::ranges::any_of(records, [&](const auto& r) {
                       return r.backward_model == BackwardRateModel::arrhenius &&
                              controls(r.backward);
                   })) {
            out.backward_only[out.backward_only_count++] = m;
        }
    }
    return out;
}

// Checks the entry count of a list and the names of its first count
// entries: each known and listed once. Returns the species indices.
template <std::size_t N, std::size_t K>
std::array<SpeciesIndex, K> checked_indices(const SpeciesCoefficients<K>& list,
                                            const std::array<std::string_view, N>& names,
                                            std::string_view role,
                                            int min_count)
{
    if (list.count < min_count || list.count > static_cast<int>(K)) {
        throw ParameterError(std::format("number of {} must be between {} and {}, got {}",
                                         role,
                                         min_count,
                                         K,
                                         list.count));
    }
    std::array<SpeciesIndex, K> index {};
    std::array<bool, N> seen {};
    for (std::size_t e = 0; e < static_cast<std::size_t>(list.count); ++e) {
        const char* name = list.entries[e].name;
        if (name == nullptr) {
            throw ParameterError(std::format("{} entry {} has no name", role, e));
        }
        const auto found = species_index(names, name);
        if (!found) {
            throw ParameterError(std::format("unknown species \"{}\" in the {}", name, role));
        }
        index[e] = *found;
        if (seen[index[e].value()]) {
            throw ParameterError(
              std::format("species \"{}\" is listed twice in the {}", name, role));
        }
        seen[index[e].value()] = true;
    }
    return index;
}

// Stoichiometric coefficients per species from a list of at least one
// species. Each must be a positive integer.
template <std::size_t N, std::size_t K>
std::array<int, N> stoichiometric_coefficients(const SpeciesCoefficients<K>& list,
                                               const std::array<std::string_view, N>& names,
                                               std::string_view role)
{
    const std::array<SpeciesIndex, K> index = checked_indices(list, names, role, 1);
    std::array<int, N> out {};
    for (std::size_t e = 0; e < static_cast<std::size_t>(list.count); ++e) {
        const double nu = list.entries[e].value;
        if (!(nu >= 1.0) || !(nu < std::ldexp(1.0, std::numeric_limits<int>::digits)) ||
            nu != std::floor(nu)) {
            throw ParameterError(
              std::format("coefficient of \"{}\" in the {} must be a positive integer, got {}",
                          list.entries[e].name,
                          role,
                          nu));
        }
        out[index[e].value()] = static_cast<int>(nu);
    }
    return out;
}

// Check the finite storage range before narrowing dimensionless rate parameters.
template <std::floating_point P>
P stored_reaction_parameter(std::common_type_t<double, P> value, std::string_view name)
{
    if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<P>::max()) {
        throw ParameterError(std::format("{} is outside the finite parameter range", name));
    }
    return static_cast<P>(value);
}

// Third-body efficiencies per species: 1 unless listed. Each listed value
// must be finite and non-negative.
template <std::floating_point P, std::size_t N, std::size_t K>
std::array<P, N> third_body_efficiencies(const SpeciesCoefficients<K>& list,
                                         const std::array<std::string_view, N>& names)
{
    const std::array<SpeciesIndex, K> index =
      checked_indices(list, names, "third-body efficiencies", 0);
    std::array<P, N> out {};
    out.fill(P(1));
    for (std::size_t e = 0; e < static_cast<std::size_t>(list.count); ++e) {
        const double eff = list.entries[e].value;
        if (!std::isfinite(eff) || eff < 0.0) {
            throw ParameterError(
              std::format("third-body efficiency of \"{}\" must be finite and non-negative, got {}",
                          list.entries[e].name,
                          eff));
        }
        out[index[e].value()] = stored_reaction_parameter<P>(
          eff, std::format("third-body efficiency of \"{}\"", list.entries[e].name));
    }
    return out;
}

// Host-side reaction orders include an optional third body. Accumulation
// is checked before addition; stored integer coefficients remain unchanged.
template <std::size_t N>
std::intmax_t
  reaction_order(const std::array<int, N>& coefficients, bool third_body, std::string_view role)
{
    std::intmax_t order = third_body ? 1 : 0;
    for (const int coefficient : coefficients) {
        if (coefficient < 0 || order > std::numeric_limits<std::intmax_t>::max() - coefficient) {
            throw ParameterError(std::format("{} reaction order is not representable", role));
        }
        order += coefficient;
    }
    return order;
}

// The integer exponents must be representable before constructing Dimension.
inline Dimension rate_dimension(std::intmax_t order, std::string_view role)
{
    constexpr std::intmax_t maximum_order = std::numeric_limits<int>::max() / 3 + 1;
    if (order < 1 || order > maximum_order) {
        throw ParameterError(std::format(
          "{} rate: reaction order {} is outside the representable dimension range [1, {}]",
          role,
          order,
          maximum_order));
    }
    const auto power = static_cast<int>(order - 1);
    return {.length = 3 * power, .time = -1, .count = -power};
}

// Rate constant of the given order in model units. A must be finite and
// positive, n finite, theta finite and non-negative and s finite in [0, 1].
// A per mole in SI
// becomes A / N_A^(order - 1) per particle; with T_SI = tau T, the model A
// is that times tau^n divided by the scale of (length^3 / count)^(order - 1)
// / time. Conversion uses at least double precision and retains long double.
template <std::floating_point P, UnitSystem U>
ArrheniusRate<P, U> arrhenius_rate(const ArrheniusInput& input,
                                   std::intmax_t order,
                                   const UnitScales<P, U>& scales,
                                   std::string_view role)
{
    using W = std::common_type_t<double, P>;
    if (!std::isfinite(input.A) || !(input.A > 0.0)) {
        throw ParameterError(
          std::format("{} rate: A must be finite and positive, got {}", role, input.A));
    }
    if (!std::isfinite(input.n)) {
        throw ParameterError(std::format("{} rate: n must be finite, got {}", role, input.n));
    }
    if (!std::isfinite(input.activation_temperature) || input.activation_temperature < 0.0) {
        throw ParameterError(
          std::format("{} rate: activation temperature must be finite and non-negative, got {}",
                      role,
                      input.activation_temperature));
    }
    const double s_t = input.translational_exponent;
    if (!std::isfinite(s_t) || s_t < 0.0 || s_t > 1.0) {
        throw ParameterError(
          std::format("{} rate: translational exponent must be finite and between 0 and 1, got {}",
                      role,
                      s_t));
    }
    const P n = stored_reaction_parameter<P>(input.n, std::format("{} rate: n", role));
    const P exponent = stored_reaction_parameter<P>(
      s_t, std::format("{} rate: translational exponent", role));
    const auto theta_name = std::format("{} rate: activation temperature", role);
    const auto theta = thermodynamic_parameter<P>(
      convert<U>(Temperature<W, si> {static_cast<W>(input.activation_temperature)}, scales),
      theta_name.c_str());
    const auto& s = scales;
    const Dimension dimension = rate_dimension(order, role);
    const auto tau = scale_factor<kind::temperature>(s);
    const auto per_particle =
      input.A / std::pow(constants::avogadro<double>, static_cast<W>(order - 1));
    const W converted_a = per_particle * std::pow(tau, input.n) / scale_factor(dimension, s);
    if (!std::isfinite(converted_a) || !(converted_a > W(0)) ||
        converted_a > std::numeric_limits<P>::max()) {
        throw ParameterError(
          std::format("{} rate: A = {} is not representable in the model's units", role, input.A));
    }
    const P a = static_cast<P>(converted_a);
    if (!(a > P(0))) {
        throw ParameterError(
          std::format("{} rate: A = {} is not representable in the model's units", role, input.A));
    }
    return {
      .A = a,
      .n = n,
      .activation_temperature = theta,
      .translational_exponent = exponent,
      .other_mode = {},
    };
}

// Validate species-bearing requirements before resolving model components.
// The host definitions and optional model reference outlive construction only.
template <typename Source, std::size_t N> struct TemperatureSpecies {
    const ReactionSpeciesView<Source, N>& species;

    SpeciesMetadataView metadata(SpeciesIndex i) const
    {
        const auto value = species.metadata(i);
        try {
            check_species_metadata(value, species.charge_number(i).value());
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("species \"{}\": {}", species.name(i), error.what()));
        }
        return value;
    }
    SpeciesIndex named(const char* name) const
    {
        if (name == nullptr || std::string_view {name}.empty()) {
            throw ParameterError("temperature requirement needs a species");
        }
        for (const SpeciesIndex i : species_indices(0, N)) {
            if (species.name(i) == name) { return i; }
        }
        throw ParameterError(std::format("unknown temperature species \"{}\"", name));
    }
    SpeciesIndex vibration(const char* name) const
    {
        const auto i = named(name);
        const auto m = metadata(i);
        const bool molecular =
          m.elements &&
          (m.elements->size() > 1 || (m.elements->size() == 1 && (*m.elements)[0].count > 1));
        if (m.role != SpeciesRole::heavy || !molecular) {
            throw ParameterError(
              std::format("vibration temperature requires a molecular heavy species: \"{}\"",
                          name));
        }
        return i;
    }
    SpeciesIndex electronic(const char* name) const
    {
        const auto i = named(name);
        if (metadata(i).role != SpeciesRole::heavy) {
            throw ParameterError(
              std::format("bound electronic temperature requires a heavy species: \"{}\"", name));
        }
        return i;
    }
    SpeciesIndex electron() const
    {
        for (const SpeciesIndex i : species_indices(0, N)) {
            if (species.metadata(i).role == SpeciesRole::electron) {
                (void)metadata(i);
                return i;
            }
        }
        throw ParameterError("free-electron temperature requires a free-electron species");
    }
};

template <typename Source, std::size_t N> struct SingleOtherTemperature {
    static constexpr std::size_t mode_count = 0;
    TemperatureSpecies<Source, N> species;
    OptionalModeIndex operator()(TranslationTemperatureInput) const { return {}; }
    OptionalModeIndex operator()(RotationTemperatureInput) const { return {}; }
    OptionalModeIndex operator()(VibrationTemperatureInput input) const
    {
        (void)species.vibration(input.species);
        return {};
    }
    OptionalModeIndex operator()(BoundElectronicTemperatureInput input) const
    {
        (void)species.electronic(input.species);
        return {};
    }
    OptionalModeIndex operator()(FreeElectronTemperatureInput) const
    {
        (void)species.electron();
        return {};
    }
};

template <typename Resolver>
OptionalModeIndex resolve_temperature(const std::optional<PhysicalTemperatureInput>& input,
                                      const Resolver& resolver,
                                      bool required,
                                      std::string_view role,
                                      std::string_view missing)
{
    try {
        if (input) { return std::visit(resolver, *input); }
        if (required && Resolver::mode_count > 0) { throw ParameterError(std::string(missing)); }
        return {};
    } catch (const ParameterError& error) {
        throw ParameterError(std::format("{} rate: {}", role, error.what()));
    }
}

// Reaction of input with the other mode of each rate from
// the host resolver, which rejects unavailable physical components.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t K, typename OtherMode>
Reaction<P, U, N> make_reaction(const ReactionInput<K>& input,
                                const std::array<std::string_view, N>& names,
                                const UnitScales<P, U>& scales,
                                const OtherMode& other_mode)
{
    Reaction<P, U, N> r {};
    const auto reactants = stoichiometric_coefficients(input.reactants, names, "reactants");
    const auto products = stoichiometric_coefficients(input.products, names, "products");
    for (const SpeciesIndex i : species_indices(0, N)) {
        r.reactants()[i] = reactants[i.value()];
        r.products()[i] = products[i.value()];
    }
    r.third_body = input.third_body;
    if (input.third_body) {
        const auto efficiencies = third_body_efficiencies<P>(input.efficiencies, names);
        for (const SpeciesIndex i : species_indices(0, N)) {
            r.efficiencies()[i] = efficiencies[i.value()];
        }
    } else if (input.efficiencies.count != 0) {
        throw ParameterError("third-body efficiencies are given without a third body");
    }
    const std::intmax_t reactant_order = reaction_order(reactants, input.third_body, "reactant");
    const std::intmax_t product_order = reaction_order(products, input.third_body, "product");
    const std::intmax_t change = product_order - reactant_order;
    if (change < std::numeric_limits<int>::min() || change > std::numeric_limits<int>::max()) {
        throw ParameterError(
          std::format("net stoichiometric coefficient {} is not representable", change));
    }
    r.coefficient_change = static_cast<int>(change);
    r.forward = arrhenius_rate(input.forward, reactant_order, scales, "forward");
    r.forward.other_mode = resolve_temperature(input.forward.other_temperature,
                                               other_mode,
                                               input.forward.translational_exponent != 1.0,
                                               "forward",
                                               "other temperature requirement is missing");
    r.backward_model = input.backward_model;
    r.equilibrium_at_other_temperature = false;
    if (input.backward_model == BackwardRateModel::arrhenius) {
        if (input.equilibrium_temperature) {
            throw ParameterError("backward rate: equilibrium temperature is only valid for an "
                                 "equilibrium backward rate");
        }
        r.backward = arrhenius_rate(input.backward, product_order, scales, "backward");
        r.backward.other_mode = resolve_temperature(input.backward.other_temperature,
                                                    other_mode,
                                                    input.backward.translational_exponent != 1.0,
                                                    "backward",
                                                    "other temperature requirement is missing");
    } else if (input.backward_model == BackwardRateModel::equilibrium) {
        if (input.backward.other_temperature) {
            throw ParameterError(
              "backward rate: an equilibrium backward rate uses equilibrium_temperature");
        }
        r.backward.other_mode =
          resolve_temperature(input.equilibrium_temperature,
                              other_mode,
                              true,
                              "backward",
                              "equilibrium temperature requirement is missing");
        r.equilibrium_at_other_temperature = r.backward.other_mode.has_value();
    } else {
        throw ParameterError("unknown backward rate model");
    }
    return r;
}

// Products fit in int64_t because coefficients and charge numbers use signed
// 32-bit storage and element counts use unsigned 32-bit storage. Check sums
// before addition so large mechanisms cannot overflow conservation arithmetic.
inline void add_conservation_term(std::int64_t& sum, std::int64_t term)
{
    if ((term > 0 && sum > std::numeric_limits<std::int64_t>::max() - term) ||
        (term < 0 && sum < std::numeric_limits<std::int64_t>::min() - term)) {
        throw ParameterError("conservation sum exceeds the signed 64-bit range");
    }
    sum += term;
}

template <std::floating_point P, UnitSystem U, typename Source, std::size_t N>
void check_reaction_conservation(const Reaction<P, U, N>& reaction,
                                 const ReactionSpeciesView<Source, N>& species)
{
    std::int64_t charge_balance = 0;
    std::map<std::string_view, std::int64_t> element_balance;
    for (const SpeciesIndex i : species_indices(0, N)) {
        const int reactant = reaction.reactants()[i];
        const int product = reaction.products()[i];
        if (reactant == 0 && product == 0) { continue; }
        const ChargeNumber charge = species.charge_number(i);
        const auto metadata = species.metadata(i);
        try {
            check_species_metadata(metadata, charge.value());
            if (!metadata.elements) {
                throw ParameterError("elemental composition is required for reaction participants");
            }
        } catch (const ParameterError& error) {
            throw ParameterError(std::format(R"(species "{}": {})", species.name(i), error.what()));
        }
        const std::int64_t net = std::int64_t {product} - reactant;
        add_conservation_term(charge_balance, net * charge.value());
        for (const auto& element : *metadata.elements) {
            add_conservation_term(element_balance[element.element],
                                  net * std::int64_t {element.count});
        }
    }
    if (charge_balance != 0) {
        throw ParameterError(
          std::format("charge is not conserved (net elementary charges {})", charge_balance));
    }
    for (const auto& [element, balance] : element_balance) {
        if (balance != 0) {
            throw ParameterError(
              std::format("element {} is not conserved (net nuclei {})", element, balance));
        }
    }
}

// Builds the reactions of input into mech; see make_reaction_mechanism.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          std::size_t K,
          typename Source,
          typename OtherMode>
ReactionMechanism<P, U, N, R> build_mechanism(const MechanismInput<R, K>& input,
                                              const ReactionSpeciesView<Source, N>& species,
                                              const UnitScales<P, U>& scales,
                                              Pressure<P, U> standard_pressure,
                                              const OtherMode& other_mode)
{
    std::array<std::string_view, N> names {};
    std::array<ParticleMass<P, U>, N> mass {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        try {
            names[i.value()] = species.name(i);
            mass[i.value()] = particle_mass_parameter(species.molar_mass(i), scales);
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("species {}: {}", i.value(), error.what()));
        }
    }
    for (const SpeciesIndex i : species_indices(0, N)) {
        if (species_index(names, names[i.value()]) != OptionalSpeciesIndex {i}) {
            throw ParameterError(std::format(R"(species "{}" is named twice)", names[i.value()]));
        }
    }
    const double lower = input.temperature_lower;
    const double upper = input.temperature_upper;
    if (!std::isfinite(lower) || !std::isfinite(upper) || !(lower > 0.0) || !(lower < upper)) {
        throw ParameterError(
          std::format("temperature limits must satisfy 0 < lower < upper, got {} and {}",
                      lower,
                      upper));
    }
    const P p0 = standard_pressure.value();
    if (!std::isfinite(p0) || !(p0 > P(0))) {
        throw ParameterError(
          std::format("standard pressure must be finite and positive, got {}", p0));
    }
    ReactionMechanism<P, U, N, R> mech {};
    for (const ReactionIndex k : reaction_indices(0, R)) {
        const ReactionInput<K>& reaction = input.reactions[k.value()];
        try {
            mech.reactions_[k.value()] = make_reaction(reaction, names, scales, other_mode);
            check_reaction_conservation(mech.reactions_[k.value()], species);
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("reaction {} ({}): {}",
                                             k.value(),
                                             reaction.equation != nullptr ? reaction.equation : "",
                                             error.what()));
        }
    }
    mech.mass_ = mass;
    mech.configured_mode_count_ = OtherMode::mode_count;
    using W = std::common_type_t<double, P>;
    mech.temperature_lower = positive_thermodynamic_parameter<P>(
      convert<U>(Temperature<W, si> {static_cast<W>(lower)}, scales), "temperature lower limit");
    mech.temperature_upper = positive_thermodynamic_parameter<P>(
      convert<U>(Temperature<W, si> {static_cast<W>(upper)}, scales), "temperature upper limit");
    if (!(mech.temperature_lower < mech.temperature_upper)) {
        throw ParameterError("temperature limits must remain ordered after conversion");
    }
    mech.standard_pressure = standard_pressure;
    mech.boltzmann = boltzmann_constant(scales);
    mech.schedule_ = prepare_rate_schedule(mech.reactions_);
    mech.controlling_modes_ = prepare_rate_mode_references<OtherMode::mode_count>(mech.reactions_);
    return mech;
}

}  // namespace detail

// Mechanism of the reactions of input over the N species of a thermo model,
// from one host species view in model order. Names, particle masses, charges
// and elemental composition come from those definitions, in the unit system of the
// scales, with standard pressure p0. Rate constants are converted to
// per-particle model units, activation temperatures and the temperature
// limits to the temperature unit of U. Throws ParameterError when:
// - a name is null, empty or listed twice, or a molar/stored particle mass
//   is not finite and positive;
// - the temperature limits are not finite with 0 < lower < upper, or p0 is
//   not finite and positive; converted temperature limits must also remain
//   finite, positive and strictly ordered;
// and, with the reaction index and equation in the message, when:
// - the reactants or products are empty, name an unknown species or one
//   species twice, or a coefficient is not a positive integer;
// - a third-body efficiency names an unknown species or is negative or not
//   finite, or efficiencies are given without a third body;
// - an Arrhenius A is not finite and positive (also after conversion), n is
//   not finite, an activation temperature is negative or not finite, or a
//   translational exponent is not finite or outside [0, 1] (the backward
//   rate is checked only with BackwardRateModel::arrhenius).
// Stored efficiencies, n, activation temperatures and translational exponents
// must remain finite after conversion to the parameter type.
// Every reaction participant must have valid explicit metadata, known elemental
// composition and an integral signed 32-bit charge number. Elemental and charge
// balances are exact checked integer sums; nonzero or overflowing balances are
// rejected with reaction context. Unchanged explicit participants are validated.
// The host species and thermo model used for rate evaluation must share ordering.
// Without a gas model this constructs an explicitly single-temperature
// mechanism. Supplied physical requirements are validated and collapsed to T_l;
// evaluating this mechanism with internal temperatures violates its contract.
template <std::floating_point P,
          UnitSystem U,
          typename Source,
          std::size_t N,
          std::size_t R,
          std::size_t K>
    requires(N > 0 && R > 0 && K > 0)
ReactionMechanism<P, U, N, R> make_reaction_mechanism(const MechanismInput<R, K>& input,
                                                      const ReactionSpeciesView<Source, N>& species,
                                                      const UnitScales<P, U>& scales,
                                                      Pressure<P, U> standard_pressure)
{
    return detail::build_mechanism(input,
                                   species,
                                   scales,
                                   standard_pressure,
                                   detail::SingleOtherTemperature<Source, N> {{species}});
}

// As above, with the standard pressure 1 bar (1e5 Pa).
template <std::floating_point P,
          UnitSystem U,
          typename Source,
          std::size_t N,
          std::size_t R,
          std::size_t K>
    requires(N > 0 && R > 0 && K > 0)
ReactionMechanism<P, U, N, R> make_reaction_mechanism(const MechanismInput<R, K>& input,
                                                      const ReactionSpeciesView<Source, N>& species,
                                                      const UnitScales<P, U>& scales)
{
    return make_reaction_mechanism(input,
                                   species,
                                   scales,
                                   convert<U>(Pressure<P, si> {P(1.0e5)}, scales));
}

// A thermo model with energy modes: G::mode_count modes, mode m with its
// own temperature T_m. species_vibration_mode(gas, i) is the mode that holds
// the vibration of species i (in the model's species order) and
// electron_energy_mode(gas) the mode that holds the energy of the free
// electrons; each returns absence when there is none.
template <typename G>
concept HasEnergyModes = GasModelType<G> && requires(const G& gas, SpeciesIndex i) {
    { G::mode_count } -> std::convertible_to<std::size_t>;
    { species_vibration_mode(gas, i) } -> std::same_as<OptionalModeIndex>;
    { electron_energy_mode(gas) } -> std::same_as<OptionalModeIndex>;
};

// Reaction construction also resolves bound electronic components.
template <typename G>
concept HasReactionTemperatureModes =
  HasEnergyModes<G> && detail::HasModeCount<G> && requires(const G& gas, SpeciesIndex i) {
      { species_electronic_mode(gas, i) } -> std::same_as<OptionalModeIndex>;
  };

namespace detail
{
template <HasReactionTemperatureModes G, typename Source, std::size_t N>
struct ModelOtherTemperature {
    static constexpr std::size_t mode_count = G::mode_count;
    const G& gas;
    TemperatureSpecies<Source, N> species;

    OptionalModeIndex checked(OptionalModeIndex mode, std::string_view component) const
    {
        if (!mode) {
            throw ParameterError(
              std::format("temperature component is not assigned to a mode: {}", component));
        }
        if ((*mode).value() >= mode_count) {
            throw ParameterError(
              std::format("temperature mode is outside the model mode range: {}", component));
        }
        return mode;
    }
    OptionalModeIndex operator()(TranslationTemperatureInput) const { return {}; }
    OptionalModeIndex operator()(RotationTemperatureInput) const { return {}; }
    OptionalModeIndex operator()(VibrationTemperatureInput input) const
    {
        const auto i = species.vibration(input.species);
        return checked(species_vibration_mode(gas, i),
                       std::format("vibration of {}", input.species));
    }
    OptionalModeIndex operator()(BoundElectronicTemperatureInput input) const
    {
        const auto i = species.electronic(input.species);
        return checked(species_electronic_mode(gas, i),
                       std::format("bound electronic of {}", input.species));
    }
    OptionalModeIndex operator()(FreeElectronTemperatureInput) const
    {
        (void)species.electron();
        return checked(electron_energy_mode(gas), "free electron");
    }
};
}  // namespace detail

// Resolve each physical requirement against gas, in the host species order.
// Translation and rotation use the primary temperature. Vibration, bound
// electronic excitation and free electrons require their configured components.
// Even a one-mode model requires explicit active Arrhenius and equilibrium
// sources. Unused supplied requirements are validated. Invalid species metadata,
// missing components and out-of-range model queries throw ParameterError with
// reaction, rate and component context. No reactant inference is performed.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          std::size_t K,
          typename Source,
          HasReactionTemperatureModes G>
    requires(N > 0 && R > 0 && K > 0) && std::same_as<typename G::units, U>
ReactionMechanism<P, U, N, R> make_reaction_mechanism(const MechanismInput<R, K>& input,
                                                      const ReactionSpeciesView<Source, N>& species,
                                                      const UnitScales<P, U>& scales,
                                                      Pressure<P, U> standard_pressure,
                                                      const G& gas)
{
    return detail::build_mechanism(
      input,
      species,
      scales,
      standard_pressure,
      detail::ModelOtherTemperature<G, Source, N> {.gas = gas, .species = {species}});
}

// Owned host species definitions are borrowed only during construction.
template <typename... Explicit, typename Input, typename Source, typename... Rest>
auto make_reaction_mechanism(const Input& input,
                             const detail::OwnedReactionSpecies<Source>& species,
                             const Rest&... rest)
{
    return detail::with_inputs(
      [&](const auto& borrowed, const auto& definitions) {
          return make_reaction_mechanism<Explicit...>(borrowed,
                                                      reaction_species(definitions),
                                                      rest...);
      },
      input,
      species.input);
}

// A thermo model of N species whose species_properties(gas, temp) gives the
// SpeciesProperties of every species at a temperature with value type T.
template <typename G, typename T, std::size_t N>
concept HasSpeciesProperties =
  GasModelType<G> && Scalar<T> && requires(const G& gas, Temperature<T, typename G::units> temp) {
      {
          species_properties(gas, temp)
      } -> std::convertible_to<SpeciesProperties<T, typename G::units, N>>;
  };

// Forward and backward rates of progress of R reactions, in particles per
// unit volume and time.
template <Scalar T, UnitSystem U, std::size_t R>
    requires(R > 0)
struct RatesOfProgress {
    static constexpr std::size_t reaction_count = R;

  private:
    std::array<NumberRateDensity<T, U>, R> forward_;
    std::array<NumberRateDensity<T, U>, R> backward_;

  public:
    // Views borrow this output. Reacquire them after moving or replacing the owner.
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto forward() & noexcept
    {
        return as_indexed<ReactionIndex>(
          Kokkos::mdspan<NumberRateDensity<T, U>, Kokkos::extents<std::size_t, R>>(
            forward_.data()));
    }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto forward() const& noexcept
    {
        return as_indexed<ReactionIndex>(
          Kokkos::mdspan<const NumberRateDensity<T, U>, Kokkos::extents<std::size_t, R>>(
            forward_.data()));
    }
    auto forward() && = delete;
    auto forward() const&& = delete;
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto backward() & noexcept
    {
        return as_indexed<ReactionIndex>(
          Kokkos::mdspan<NumberRateDensity<T, U>, Kokkos::extents<std::size_t, R>>(
            backward_.data()));
    }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr auto backward() const& noexcept
    {
        return as_indexed<ReactionIndex>(
          Kokkos::mdspan<const NumberRateDensity<T, U>, Kokkos::extents<std::size_t, R>>(
            backward_.data()));
    }
    auto backward() && = delete;
    auto backward() const&& = delete;
};

// A thermo model and value type for which rates_of_progress and
// production_rates of a mechanism of N species with parameter type P can be
// evaluated.
template <typename G, typename T, typename P, std::size_t N>
concept ReactingModel = HasSpeciesProperties<G, T, N> && std::floating_point<P> &&
                        ScalarWith<T, P> && HasExp<T> && HasLog<T>;

namespace detail
{

// temp clamped to [temperature_lower, temperature_upper] of mech, on its
// value part; a clamped temperature has derivative zero.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
CANARY_HOST_DEVICE T clamped_temperature(const ReactionMechanism<P, U, N, R>& mech,
                                         Temperature<T, U> temp)
{
    using V = ScalarTraits<T>::value_type;
    const V t_value = ScalarTraits<T>::value(temp.value());
    const auto lower = static_cast<V>(mech.temperature_lower.value());
    const auto upper = static_cast<V>(mech.temperature_upper.value());
    if (t_value < lower) { return T(lower); }
    if (upper < t_value) { return T(upper); }
    return temp.value();
}

// ln(k / A) of an Arrhenius rate at T_l, from ln T_l and 1 / T_l:
// n ln T_l - theta / T_l.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE T arrhenius_exponent(const ArrheniusRate<P, U>& rate, T ln_t, T inv_t)
{ return T(rate.n * ln_t) - T(rate.activation_temperature.value() * inv_t); }

// ln(k / A) of an Arrhenius rate at T_c = T_l^s T_ve^(1 - s), from ln T_l,
// 1 / T_l and ln_ratio = ln(T_ve / T_l): with w = (1 - s) ln_ratio,
// ln T_c = ln T_l + w and 1 / T_c = exp(-w) / T_l. w = 0 gives the value of
// the expression at T_l.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P> && HasExp<T>
CANARY_HOST_DEVICE T
  arrhenius_exponent(const ArrheniusRate<P, U>& rate, T ln_t, T inv_t, T ln_ratio)
{
    const T w = T((P(1) - rate.translational_exponent) * ln_ratio);
    return T(rate.n * T(ln_t + w)) -
           T(rate.activation_temperature.value() * T(inv_t * math::exp(T(-w))));
}

// ln(k / A) of an Arrhenius rate in rates_of_progress_at: at T_l with no
// other temperature (M = 0); otherwise at T_c for s != 1, from ln_ratio =
// ln(T_o / T_l) of the rate's other temperature, and at T_l for s = 1, so
// that a rate that depends on T only computes no exponential of w.
template <std::size_t M, Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P> && HasExp<T>
CANARY_HOST_DEVICE T
  rate_exponent(const ArrheniusRate<P, U>& rate, T ln_t, T inv_t, const std::array<T, M>& ln_ratio)
{
    if constexpr (M > 0) {
        if (rate.other_mode && rate.translational_exponent != P(1)) {
            return arrhenius_exponent(rate, ln_t, inv_t, ln_ratio[(*rate.other_mode).value()]);
        }
    }
    return arrhenius_exponent(rate, ln_t, inv_t);
}

// ln(p0 / (k_B T)) and the standard Gibbs energies per particle over k_B T,
// g_i m_i / (k_B T), at a temperature T, for the equilibrium constants.
template <Scalar T, std::size_t N> struct EquilibriumTerms {
    T ln_n0;
    std::array<T, N> g;
};

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE EquilibriumTerms<T, N>
  equilibrium_terms(const ReactionMechanism<P, U, N, R>& mech, const G& gas, T t)
{
    using V = ScalarTraits<T>::value_type;
    const SpeciesProperties<T, U, N> props = species_properties(gas, Temperature<T, U> {t});
    const T inv_kt = T(V(1)) / T(mech.boltzmann.value() * t);
    EquilibriumTerms<T, N> out {};
    out.ln_n0 = math::log(T(mech.standard_pressure.value() * inv_kt));
    for (const SpeciesIndex i : species_indices(0, N)) {
        out.g[i.value()] =
          T(props.gibbs_energy[i.value()].value() * mech.mass()[i].value() * inv_kt);
    }
    return out;
}

// ln(k_b / A_f) of an equilibrium backward rate, from ln(k_f / A_f) at the
// equilibrium temperature and the terms of equilibrium_terms there:
// ln(k_f / A_f) + sum dnu_i g_i m_i / (k_B T) - sum(dnu) ln(p0 / (k_B T)).
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE T equilibrium_exponent(const Reaction<P, U, N>& r,
                                          T forward_exponent,
                                          const EquilibriumTerms<T, N>& terms)
{
    using V = ScalarTraits<T>::value_type;
    T exponent = forward_exponent - T(V(r.coefficient_change) * terms.ln_n0);
    for (const SpeciesIndex i : species_indices(0, N)) {
        const int dnu = r.products()[i] - r.reactants()[i];
        if (dnu != 0) { exponent = exponent + T(V(dnu) * terms.g[i.value()]); }
    }
    return exponent;
}

// Explicit reverse blocks carry no equilibrium storage. An equilibrium source
// borrows the single terms object owned by its source traversal.
struct ArrheniusReverseSource {};
template <Scalar T, std::size_t N> struct EquilibriumRateSource {
    const EquilibriumTerms<T, N>& terms;
    T ln_temperature;
    T inverse_temperature;
};

template <bool ForwardOther,
          bool BackwardOther,
          bool ThirdBody,
          bool EquilibriumOther,
          std::size_t M,
          Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          typename Source>
CANARY_HOST_DEVICE void evaluate_rate_block(const ReactionMechanism<P, U, N, R>& mech,
                                            const RateSchedule<R>& schedule,
                                            std::size_t begin,
                                            std::size_t end,
                                            T ln_t,
                                            T inv_t,
                                            const std::array<T, M>& ln_ratio,
                                            const std::array<T, N>& n,
                                            const Source& source,
                                            RatesOfProgress<T, U, R>& out)
{
    constexpr bool equilibrium = !std::same_as<Source, ArrheniusReverseSource>;
    const auto records = mech.reactions();
    for (std::size_t position = begin; position < end; ++position) {
        const ReactionIndex k = schedule.order[position];
        const auto& r = records[k];
        const T forward_exponent = arrhenius_exponent(r.forward, ln_t, inv_t);
        T controlled_exponent = forward_exponent;
        if constexpr (ForwardOther && M > 0) {
            controlled_exponent =
              arrhenius_exponent(r.forward, ln_t, inv_t, ln_ratio[(*r.forward.other_mode).value()]);
        }
        const T k_f = T(r.forward.A * math::exp(controlled_exponent));
        T k_b = k_f;
        if constexpr (equilibrium) {
            T exponent = forward_exponent;
            if constexpr (EquilibriumOther && M > 0) {
                exponent = equilibrium_exponent(
                  r,
                  arrhenius_exponent(r.forward, source.ln_temperature, source.inverse_temperature),
                  source.terms);
            } else {
                exponent = equilibrium_exponent(r, forward_exponent, source.terms);
            }
            k_b = T(r.forward.A * math::exp(exponent));
        } else {
            const T exponent = [&] {
                if constexpr (BackwardOther && M > 0) {
                    return arrhenius_exponent(r.backward,
                                              ln_t,
                                              inv_t,
                                              ln_ratio[(*r.backward.other_mode).value()]);
                } else {
                    return arrhenius_exponent(r.backward, ln_t, inv_t);
                }
            }();
            k_b = T(r.backward.A * math::exp(exponent));
        }
        T q_f = k_f;
        T q_b = k_b;
        if constexpr (ThirdBody) {
            T m = T(r.efficiencies()[SpeciesIndex {0}] * n[0]);
            for (const SpeciesIndex i : species_indices(1, N)) {
                m = m + T(r.efficiencies()[i] * n[i.value()]);
            }
            q_f = q_f * m;
            q_b = q_b * m;
        }
        for (const SpeciesIndex i : species_indices(0, N)) {
            for (int c = 0; c < r.reactants()[i]; ++c) { q_f = q_f * n[i.value()]; }
            for (int c = 0; c < r.products()[i]; ++c) { q_b = q_b * n[i.value()]; }
        }
        out.forward()[k] = NumberRateDensity<T, U> {q_f};
        out.backward()[k] = NumberRateDensity<T, U> {q_b};
    }
}

template <bool EquilibriumOther,
          std::size_t M,
          Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          typename Source>
CANARY_HOST_DEVICE void dispatch_rate_block(const ReactionMechanism<P, U, N, R>& mech,
                                            const RateSchedule<R>& schedule,
                                            std::size_t block,
                                            T ln_t,
                                            T inv_t,
                                            const std::array<T, M>& ln_ratio,
                                            const std::array<T, N>& n,
                                            const Source& source,
                                            RatesOfProgress<T, U, R>& out)
{
    const std::size_t begin = block == 0 ? 0 : schedule.blocks[block - 1].end;
    const std::size_t end = schedule.blocks[block].end;
    const auto evaluate = [&]<bool F, bool B, bool Third>() {
        evaluate_rate_block<F && (M > 0), B && (M > 0), Third, EquilibriumOther>(mech,
                                                                                 schedule,
                                                                                 begin,
                                                                                 end,
                                                                                 ln_t,
                                                                                 inv_t,
                                                                                 ln_ratio,
                                                                                 n,
                                                                                 source,
                                                                                 out);
    };
    if constexpr (std::same_as<Source, ArrheniusReverseSource>) {
        switch (schedule.blocks[block].behavior) {
        case RateBehavior::ArrheniusTranslationTranslation:
            evaluate.template operator()<false, false, false>();
            break;
        case RateBehavior::ArrheniusTranslationTranslationThirdBody:
            evaluate.template operator()<false, false, true>();
            break;
        case RateBehavior::ArrheniusTranslationOther:
            evaluate.template operator()<false, true, false>();
            break;
        case RateBehavior::ArrheniusTranslationOtherThirdBody:
            evaluate.template operator()<false, true, true>();
            break;
        case RateBehavior::ArrheniusOtherTranslation:
            evaluate.template operator()<true, false, false>();
            break;
        case RateBehavior::ArrheniusOtherTranslationThirdBody:
            evaluate.template operator()<true, false, true>();
            break;
        case RateBehavior::ArrheniusOtherOther:
            evaluate.template operator()<true, true, false>();
            break;
        case RateBehavior::ArrheniusOtherOtherThirdBody:
            evaluate.template operator()<true, true, true>();
            break;
        default: CANARY_ASSERT(false); break;
        }
    } else {
        switch (schedule.blocks[block].behavior) {
        case RateBehavior::EquilibriumForwardTranslation:
            evaluate.template operator()<false, false, false>();
            break;
        case RateBehavior::EquilibriumForwardTranslationThirdBody:
            evaluate.template operator()<false, false, true>();
            break;
        case RateBehavior::EquilibriumForwardOther:
            evaluate.template operator()<true, false, false>();
            break;
        case RateBehavior::EquilibriumForwardOtherThirdBody:
            evaluate.template operator()<true, false, true>();
            break;
        default: CANARY_ASSERT(false); break;
        }
    }
}

// Rates at clamped primary temperature t and sampled, unclamped internal
// temperatures. Only prepared dependencies are evaluated. Each equilibrium
// source owns one Gibbs array while all of its behavior blocks are consumed.
template <std::size_t M,
          Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE RatesOfProgress<T, U, R>
  rates_of_progress_at(const ReactionMechanism<P, U, N, R>& mech,
                       const G& gas,
                       const NormalisedComposition<T, U, N>& mix,
                       Density<T, U> rho,
                       T t,
                       const std::array<T, M>& other)
{
    using V = ScalarTraits<T>::value_type;
    const auto& schedule = reaction_schedule<M>(mech);
    const T ln_t = math::log(t);
    const T inv_t = T(V(1)) / t;
    std::array<T, N> n {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        n[i.value()] = T(rho.value() * mix.fractions[i.value()].value() / mech.mass()[i].value());
    }
    std::array<T, M> ln_ratio {};
    if constexpr (M > 0) {
        const auto& references = reaction_mode_references(mech);
        if constexpr (M == 1) {
            if (references.forward_count != 0 || references.backward_only_count != 0) {
                ln_ratio[0] = math::log(T(other[0] / t));
            }
        } else {
            for (std::size_t j = 0; j < references.forward_count; ++j) {
                const ModeIndex m = references.forward[j];
                ln_ratio[m.value()] = math::log(T(other[m.value()] / t));
            }
            for (std::size_t j = 0; j < references.backward_only_count; ++j) {
                const ModeIndex m = references.backward_only[j];
                ln_ratio[m.value()] = math::log(T(other[m.value()] / t));
            }
        }
    }
    RatesOfProgress<T, U, R> out {};
    for (std::size_t b = 0; b < schedule.explicit_end; ++b) {
        dispatch_rate_block<false>(mech,
                                   schedule,
                                   b,
                                   ln_t,
                                   inv_t,
                                   ln_ratio,
                                   n,
                                   ArrheniusReverseSource {},
                                   out);
    }
    // Explicit single-temperature evaluation shares primary properties across
    // all equilibrium blocks, retaining their original prepared permutation.
    const std::size_t primary_end = M == 0 ? schedule.block_count : schedule.translation_end;
    if (schedule.explicit_end < primary_end) {
        const auto terms = equilibrium_terms(mech, gas, t);
        const EquilibriumRateSource<T, N> source {terms, ln_t, inv_t};
        for (std::size_t b = schedule.explicit_end; b < primary_end; ++b) {
            dispatch_rate_block<false>(mech, schedule, b, ln_t, inv_t, ln_ratio, n, source, out);
        }
    }
    if constexpr (M > 0) {
        std::size_t b = schedule.translation_end;
        while (b < schedule.block_count) {
            const auto mode = schedule.blocks[b].equilibrium_mode;
            const T temperature = other[(*mode).value()];
            const T ln_other = math::log(temperature);
            const T inv_other = T(V(1)) / temperature;
            const auto terms = equilibrium_terms(mech, gas, temperature);
            const EquilibriumRateSource<T, N> source {terms, ln_other, inv_other};
            do {
                dispatch_rate_block<true>(mech, schedule, b, ln_t, inv_t, ln_ratio, n, source, out);
                ++b;
            } while (b < schedule.block_count && schedule.blocks[b].equilibrium_mode == mode);
        }
    }
    return out;
}

// omega_i = m_i sum_r (nu''_ir - nu'_ir) (q_f,r - q_b,r) from the rates of
// progress q of mech.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE CANARY_FORCE_INLINE std::array<MassProductionRate<T, U>, N>
  mass_production_rates(const ReactionMechanism<P, U, N, R>& mech,
                        const RatesOfProgress<T, U, R>& q)
{
    using V = ScalarTraits<T>::value_type;
    const auto records = mech.reactions();
    std::array<T, N> net {};
    net.fill(T(V(0)));
    for (const ReactionIndex k : reaction_indices(0, R)) {
        const Reaction<P, U, N>& r = records[k];
        const T rate = q.forward()[k].value() - q.backward()[k].value();
        for (const SpeciesIndex i : species_indices(0, N)) {
            const int dnu = r.products()[i] - r.reactants()[i];
            if (dnu != 0) { net[i.value()] = net[i.value()] + T(V(dnu) * rate); }
        }
    }
    std::array<MassProductionRate<T, U>, N> out {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        out[i.value()] = MassProductionRate<T, U> {T(mech.mass()[i].value() * net[i.value()])};
    }
    return out;
}

}  // namespace detail

// Rates of progress of every reaction of mech in a cell of density rho,
// translational temperature temp and vibrational-electronic temperature
// temp_ve, with the normalised mass fractions of mix. Requires a one-mode gas
// and a mechanism constructed with a one-mode gas. Every rate constant
// and the Gibbs energies of the equilibrium constants use T_l, temp clamped
// to [temperature_lower, temperature_upper] of mech (on its value part; a
// clamped temperature has derivative zero), except that temp_ve is T_o of
// rates resolved to the internal mode: it is not clamped and enters through the rate-controlling
// temperatures of the Arrhenius rates and the equilibrium backward rates
// at T_o (see the top of this file). The Arrhenius expressions compute on
// the values of the quantities, since A has no fixed kind. The gas model is
// called through species_properties at T_l, and at T_ve when a reaction
// has its equilibrium backward rate at T_o. The status of mix is not
// checked.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U> &&
             detail::HasModeCount<G> && (G::mode_count == 1)
CANARY_HOST_DEVICE RatesOfProgress<T, U, R>
  rates_of_progress(const ReactionMechanism<P, U, N, R>& mech,
                    const G& gas,
                    const NormalisedComposition<T, U, N>& mix,
                    Density<T, U> rho,
                    Temperature<T, U> temp,
                    Temperature<T, U> temp_ve)
{
    CANARY_ASSERT(detail::reaction_mode_count(mech) == 1);
    return detail::rates_of_progress_at<1>(mech,
                                           gas,
                                           mix,
                                           rho,
                                           detail::clamped_temperature(mech, temp),
                                           std::array<T, 1> {temp_ve.value()});
}

// Mode temperatures borrow the caller's storage. The extent must equal
// G::mode_count and the mechanism's configured mode count must agree before any
// input read. The model determines the size of the numerical scratch arrays.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          typename G,
          typename V>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U> &&
             detail::HasModeCount<G> && ModeRangeOf<V, kind::mode_temperature, T, U>
CANARY_HOST_DEVICE RatesOfProgress<T, U, R>
  rates_of_progress(const ReactionMechanism<P, U, N, R>& mech,
                    const G& gas,
                    const NormalisedComposition<T, U, N>& mix,
                    Density<T, U> rho,
                    Temperature<T, U> temp,
                    const V& mode_temp)
{
    constexpr std::size_t modes = G::mode_count;
    CANARY_ASSERT(static_cast<std::size_t>(mode_temp.extent(0)) == modes);
    CANARY_ASSERT(detail::reaction_mode_count(mech) == modes);
    std::array<T, modes> other {};
    for (const ModeIndex m : mode_indices(0, modes)) {
        other[m.value()] = detail::mode_at(mode_temp, m).value();
    }
    return detail::rates_of_progress_at<modes>(mech,
                                               gas,
                                               mix,
                                               rho,
                                               detail::clamped_temperature(mech, temp),
                                               other);
}

// Rates of progress with the M mode temperatures mode_temp of the cell
// (for example MultiTemperatureState::mode_temp): as above, with T_o of
// each rate the temperature of its other mode (ArrheniusRate::other_mode;
// resolved during construction). The gas model is called through
// species_properties at T_l and at each mode temperature at which a
// reaction has its equilibrium backward rate.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          typename G,
          std::size_t M>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U> &&
             detail::HasModeCount<G> && (M == G::mode_count)
CANARY_HOST_DEVICE RatesOfProgress<T, U, R>
  rates_of_progress(const ReactionMechanism<P, U, N, R>& mech,
                    const G& gas,
                    const NormalisedComposition<T, U, N>& mix,
                    Density<T, U> rho,
                    Temperature<T, U> temp,
                    const std::array<ModeTemperature<T, U>, M>& mode_temp)
{
    const auto view = as_indexed<ModeIndex>(
      Kokkos::mdspan<const ModeTemperature<T, U>, Kokkos::extents<std::size_t, M>>(
        mode_temp.data()));
    return rates_of_progress(mech, gas, mix, rho, temp, view);
}

// Rates of progress at one temperature: every T_o is T_l, so that every
// rate-controlling temperature and every equilibrium temperature is T_l.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE RatesOfProgress<T, U, R>
  rates_of_progress(const ReactionMechanism<P, U, N, R>& mech,
                    const G& gas,
                    const NormalisedComposition<T, U, N>& mix,
                    Density<T, U> rho,
                    Temperature<T, U> temp)
{
    return detail::rates_of_progress_at<0>(mech,
                                           gas,
                                           mix,
                                           rho,
                                           detail::clamped_temperature(mech, temp),
                                           std::array<T, 0> {});
}

// Mass production rate of every species, omega_i = m_i sum_r (nu''_ir -
// nu'_ir) (q_f,r - q_b,r), from the rates of progress of rates_of_progress
// with the same arguments.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U> &&
             detail::HasModeCount<G> && (G::mode_count == 1)
CANARY_HOST_DEVICE std::array<MassProductionRate<T, U>, N>
  production_rates(const ReactionMechanism<P, U, N, R>& mech,
                   const G& gas,
                   const NormalisedComposition<T, U, N>& mix,
                   Density<T, U> rho,
                   Temperature<T, U> temp,
                   Temperature<T, U> temp_ve)
{
    return detail::mass_production_rates(mech,
                                         rates_of_progress(mech, gas, mix, rho, temp, temp_ve));
}

// Production rates using a borrowed mode-temperature view of G::mode_count entries.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          typename G,
          typename V>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U> &&
             detail::HasModeCount<G> && ModeRangeOf<V, kind::mode_temperature, T, U>
CANARY_HOST_DEVICE std::array<MassProductionRate<T, U>, N>
  production_rates(const ReactionMechanism<P, U, N, R>& mech,
                   const G& gas,
                   const NormalisedComposition<T, U, N>& mix,
                   Density<T, U> rho,
                   Temperature<T, U> temp,
                   const V& mode_temp)
{
    return detail::mass_production_rates(mech,
                                         rates_of_progress(mech, gas, mix, rho, temp, mode_temp));
}

// As above with the M mode temperatures of the cell.
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t R,
          typename G,
          std::size_t M>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U> &&
             detail::HasModeCount<G> && (M == G::mode_count)
CANARY_HOST_DEVICE std::array<MassProductionRate<T, U>, N>
  production_rates(const ReactionMechanism<P, U, N, R>& mech,
                   const G& gas,
                   const NormalisedComposition<T, U, N>& mix,
                   Density<T, U> rho,
                   Temperature<T, U> temp,
                   const std::array<ModeTemperature<T, U>, M>& mode_temp)
{
    return detail::mass_production_rates(mech,
                                         rates_of_progress(mech, gas, mix, rho, temp, mode_temp));
}

// As above at one temperature.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t R, typename G>
    requires ReactingModel<G, T, P, N> && std::same_as<typename G::units, U>
CANARY_HOST_DEVICE std::array<MassProductionRate<T, U>, N>
  production_rates(const ReactionMechanism<P, U, N, R>& mech,
                   const G& gas,
                   const NormalisedComposition<T, U, N>& mix,
                   Density<T, U> rho,
                   Temperature<T, U> temp)
{ return detail::mass_production_rates(mech, rates_of_progress(mech, gas, mix, rho, temp)); }

namespace detail
{

struct OwnedVibrationTemperature {
    OwnedSpeciesName species;
};
struct OwnedBoundElectronicTemperature {
    OwnedSpeciesName species;
};
struct StorePhysicalTemperature {
    using Storage = std::variant<TranslationTemperatureInput,
                                 RotationTemperatureInput,
                                 OwnedVibrationTemperature,
                                 OwnedBoundElectronicTemperature,
                                 FreeElectronTemperatureInput>;
    Storage operator()(TranslationTemperatureInput v) const { return v; }
    Storage operator()(RotationTemperatureInput v) const { return v; }
    Storage operator()(VibrationTemperatureInput v) const
    { return OwnedVibrationTemperature {OwnedSpeciesName {v.species}}; }
    Storage operator()(BoundElectronicTemperatureInput v) const
    { return OwnedBoundElectronicTemperature {OwnedSpeciesName {v.species}}; }
    Storage operator()(FreeElectronTemperatureInput v) const { return v; }
};
struct BorrowPhysicalTemperature {
    PhysicalTemperatureInput operator()(TranslationTemperatureInput v) const { return v; }
    PhysicalTemperatureInput operator()(RotationTemperatureInput v) const { return v; }
    PhysicalTemperatureInput operator()(const OwnedVibrationTemperature& v) const
    { return VibrationTemperatureInput {v.species.c_str()}; }
    PhysicalTemperatureInput operator()(const OwnedBoundElectronicTemperature& v) const
    { return BoundElectronicTemperatureInput {v.species.c_str()}; }
    PhysicalTemperatureInput operator()(FreeElectronTemperatureInput v) const { return v; }
};
struct OwnedPhysicalTemperature {
    std::optional<StorePhysicalTemperature::Storage> value;
    OwnedPhysicalTemperature() = default;
    explicit OwnedPhysicalTemperature(const std::optional<PhysicalTemperatureInput>& input)
    {
        if (input) { value = std::visit(StorePhysicalTemperature {}, *input); }
    }
    std::optional<PhysicalTemperatureInput> borrow() const
    {
        if (value) { return std::visit(BorrowPhysicalTemperature {}, *value); }
        return std::nullopt;
    }
};

template <> struct InputStorage<ArrheniusInput> {
    double A {};
    double n {};
    double activation_temperature {};
    double translational_exponent {};
    OwnedPhysicalTemperature other_temperature {};
    using Scratch = std::tuple<>;
    InputStorage() = default;
    explicit InputStorage(const ArrheniusInput& input) :
        A(input.A), n(input.n), activation_temperature(input.activation_temperature),
        translational_exponent(input.translational_exponent),
        other_temperature(input.other_temperature)
    {
    }
    ArrheniusInput borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .A = A,
          .n = n,
          .activation_temperature = activation_temperature,
          .translational_exponent = translational_exponent,
          .other_temperature = other_temperature.borrow(),
        };
    }
};

template <> struct InputStorage<SpeciesCoefficient> {
    OwnedSpeciesName name {};
    double value {};
    using Scratch = std::tuple<>;
    InputStorage() = default;
    explicit InputStorage(const SpeciesCoefficient& input) : name(input.name), value(input.value) {}
    SpeciesCoefficient borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .name = name.c_str(),
          .value = value,
        };
    }
};

template <std::size_t K> struct InputStorage<SpeciesCoefficients<K>> {
    int count {};
    InputStorage<std::array<SpeciesCoefficient, K>> entries {};
    using Scratch = std::tuple<typename InputStorage<std::array<SpeciesCoefficient, K>>::Scratch>;
    InputStorage() = default;
    explicit InputStorage(const SpeciesCoefficients<K>& input) :
        count(input.count), entries(input.entries)
    {
    }
    SpeciesCoefficients<K> borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .count = count,
          .entries = entries.borrow(std::get<0>(scratch)),
        };
    }
};

template <std::size_t K> struct InputStorage<ReactionInput<K>> {
    OwnedReactionEquation equation {};
    InputStorage<SpeciesCoefficients<K>> reactants {};
    InputStorage<SpeciesCoefficients<K>> products {};
    InputStorage<ArrheniusInput> forward {};
    BackwardRateModel backward_model {};
    InputStorage<ArrheniusInput> backward {};
    OwnedPhysicalTemperature equilibrium_temperature {};
    bool third_body {};
    InputStorage<SpeciesCoefficients<K>> efficiencies {};
    using Scratch = std::tuple<typename InputStorage<SpeciesCoefficients<K>>::Scratch,
                               typename InputStorage<SpeciesCoefficients<K>>::Scratch,
                               typename InputStorage<ArrheniusInput>::Scratch,
                               typename InputStorage<ArrheniusInput>::Scratch,
                               typename InputStorage<SpeciesCoefficients<K>>::Scratch>;
    InputStorage() = default;
    explicit InputStorage(const ReactionInput<K>& input) :
        equation(input.equation), reactants(input.reactants), products(input.products),
        forward(input.forward), backward_model(input.backward_model), backward(input.backward),
        equilibrium_temperature(input.equilibrium_temperature), third_body(input.third_body),
        efficiencies(input.efficiencies)
    {
    }
    ReactionInput<K> borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .equation = equation.c_str(),
          .reactants = reactants.borrow(std::get<0>(scratch)),
          .products = products.borrow(std::get<1>(scratch)),
          .forward = forward.borrow(std::get<2>(scratch)),
          .backward_model = backward_model,
          .backward = backward.borrow(std::get<3>(scratch)),
          .equilibrium_temperature = equilibrium_temperature.borrow(),
          .third_body = third_body,
          .efficiencies = efficiencies.borrow(std::get<4>(scratch)),
        };
    }
};

template <std::size_t R, std::size_t K> struct InputStorage<MechanismInput<R, K>> {
    double temperature_lower {};
    double temperature_upper {};
    InputStorage<std::array<ReactionInput<K>, R>> reactions {};
    using Scratch = std::tuple<typename InputStorage<std::array<ReactionInput<K>, R>>::Scratch>;
    InputStorage() = default;
    explicit InputStorage(const MechanismInput<R, K>& input) :
        temperature_lower(input.temperature_lower), temperature_upper(input.temperature_upper),
        reactions(input.reactions)
    {
    }
    MechanismInput<R, K> borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .temperature_lower = temperature_lower,
          .temperature_upper = temperature_upper,
          .reactions = reactions.borrow(std::get<0>(scratch)),
        };
    }
};

}  // namespace detail

// Builds from owned and direct inputs, keeping borrowed views inside this call.
template <typename... Explicit, typename... Inputs>
    requires detail::HasOwnedInput<Inputs...>
auto make_reaction_mechanism(const Inputs&... inputs)
{
    return detail::with_inputs(
      [](const auto&... borrowed) { return make_reaction_mechanism<Explicit...>(borrowed...); },
      inputs...);
}

}  // namespace canary

#endif  // CANARY_REACTIONS_HPP
