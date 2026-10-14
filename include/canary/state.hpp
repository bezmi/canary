#ifndef CANARY_STATE_HPP
#define CANARY_STATE_HPP

#include <array>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/status.hpp>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <limits>
#include <type_traits>

namespace canary
{

// An internal energy mode identity, separate from species identities and
// array extents. The containing model determines the mode-count bound.
// Default construction denotes mode 0; value() is the storage-access boundary.
class ModeIndex
{
  public:
    CANARY_HOST_DEVICE constexpr ModeIndex() noexcept {}
    CANARY_HOST_DEVICE constexpr explicit ModeIndex(std::size_t value) noexcept : value_(value) {}
    [[nodiscard]] CANARY_HOST_DEVICE constexpr std::size_t value() const noexcept { return value_; }
    friend CANARY_HOST_DEVICE constexpr bool operator==(ModeIndex, ModeIndex) noexcept = default;

  private:
    std::size_t value_ {0};
};

// An optional mode identity in one size_t. Default construction denotes
// absence. A present index must be less than SIZE_MAX; dereference requires
// presence. The containing model determines the mode-count bound.
class OptionalModeIndex
{
  public:
    CANARY_HOST_DEVICE constexpr OptionalModeIndex() noexcept {}
    CANARY_HOST_DEVICE constexpr explicit OptionalModeIndex(ModeIndex index) noexcept :
        value_(index.value())
    { CANARY_ASSERT(value_ < std::numeric_limits<std::size_t>::max()); }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr bool has_value() const noexcept
    { return value_ != std::numeric_limits<std::size_t>::max(); }
    CANARY_HOST_DEVICE constexpr explicit operator bool() const noexcept { return has_value(); }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr ModeIndex operator*() const noexcept
    {
        CANARY_ASSERT(has_value());
        return ModeIndex {value_};
    }
    friend CANARY_HOST_DEVICE constexpr bool operator==(OptionalModeIndex,
                                                        OptionalModeIndex) noexcept = default;

  private:
    std::size_t value_ {std::numeric_limits<std::size_t>::max()};
};

// Consecutive mode identities in [first, last). The endpoint may be SIZE_MAX;
// it is never yielded. Iteration retains a native-size position.
class ModeIndexRange
{
  public:
    struct sentinel {
        std::size_t last {0};
    };
    class iterator
    {
      public:
        using value_type = ModeIndex;
        using difference_type = std::ptrdiff_t;
        using iterator_concept = std::forward_iterator_tag;
        CANARY_HOST_DEVICE constexpr iterator() noexcept {}
        [[nodiscard]] CANARY_HOST_DEVICE constexpr ModeIndex operator*() const noexcept
        { return ModeIndex {position_}; }
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
        friend class ModeIndexRange;
        std::size_t position_ {0};
    };
    [[nodiscard]] CANARY_HOST_DEVICE constexpr iterator begin() const noexcept
    { return iterator {first_}; }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr sentinel end() const noexcept { return {last_}; }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr std::size_t size() const noexcept
    { return last_ - first_; }

  private:
    CANARY_HOST_DEVICE constexpr ModeIndexRange(std::size_t first, std::size_t last) noexcept :
        first_(first), last_(last)
    {
    }
    friend CANARY_HOST_DEVICE constexpr ModeIndexRange mode_indices(std::size_t first,
                                                                    std::size_t last) noexcept;
    std::size_t first_;
    std::size_t last_;
};

// Mode identities from first up to, excluding, last. Requires first <= last.
// Equal bounds produce an empty range. No allocation or throwing access.
[[nodiscard]] CANARY_HOST_DEVICE constexpr ModeIndexRange mode_indices(std::size_t first,
                                                                       std::size_t last) noexcept
{
    CANARY_ASSERT(first <= last);
    return {first, last};
}

// A thermodynamic state with value type T in unit system U: the members
// rho, e (internal energy per unit mass), p and temp, and the member types
// value_type and units. Every model's state type satisfies it; models with
// more state variables add members.
template <typename S, typename T, typename U>
concept ThermoStateOf = Scalar<T> && UnitSystem<U> && std::same_as<typename S::value_type, T> &&
                        std::same_as<typename S::units, U> && requires(const S& s) {
                            { s.rho } -> std::convertible_to<Density<T, U>>;
                            { s.e } -> std::convertible_to<SpecificEnergy<T, U>>;
                            { s.p } -> std::convertible_to<Pressure<T, U>>;
                            { s.temp } -> std::convertible_to<Temperature<T, U>>;
                        };

// A state type of any value type and unit system.
template <typename S>
concept ThermoStateType = ThermoStateOf<S, typename S::value_type, typename S::units>;

// State of a single-temperature model.
template <Scalar T, UnitSystem U> struct ThermoState {
    using value_type = T;
    using units = U;

    Density<T, U> rho;
    SpecificEnergy<T, U> e;
    Pressure<T, U> p;
    Temperature<T, U> temp;
};

// State of a model with M internal energy modes. temp is the temperature of
// the translational and rotational energy; mode_e[k] is the energy per unit
// mass of mode k, which follows the temperature mode_temp[k]. e is the total
// internal energy per unit mass, the mode energies included.
template <Scalar T, UnitSystem U, std::size_t M>
    requires(M > 0)
struct MultiTemperatureState {
    using value_type = T;
    using units = U;

    Density<T, U> rho;
    SpecificEnergy<T, U> e;
    Pressure<T, U> p;
    Temperature<T, U> temp;
    std::array<SpecificEnergy<T, U>, M> mode_e;
    std::array<ModeTemperature<T, U>, M> mode_temp;
};

// Acoustic response under the model's documented pressure-work and energy
// constraints. sound_speed is authoritative; acoustic_exponent = rho*a*a/p.
// An acoustic exponent alone does not specify a pressure/internal-energy law.
template <Scalar T, UnitSystem U> struct AcousticProperties {
    SoundSpeed<T, U> sound_speed;
    AcousticExponent<T> acoustic_exponent;
};

// Single-temperature heat capacities and gas constant are per unit mass.
// capacity_ratio is cp/cv; acoustics describes the model's pressure response.
template <Scalar T, UnitSystem U> struct SingleTemperatureProperties {
    SpecificHeat<T, U> cv;
    SpecificHeat<T, U> cp;
    HeatCapacityRatio<T> capacity_ratio;
    SpecificGasConstant<T, U> gas_constant;
    AcousticProperties<T, U> acoustics;
};

// Common-velocity mixture properties with M internal energy modes. All heat
// capacities and gas constants are per mixture mass. mode_cv follows the same
// mode ordering as MultiTemperatureState. The heavy and electron gas constants
// determine their respective partial pressures at their own temperatures.
template <Scalar T, UnitSystem U, std::size_t M>
    requires(M > 0)
struct MultiTemperatureProperties {
    SpecificHeat<T, U> cv_heavy;
    std::array<SpecificHeat<T, U>, M> mode_cv;
    SpecificGasConstant<T, U> gas_constant_heavy;
    SpecificGasConstant<T, U> gas_constant_electron;
    AcousticProperties<T, U> acoustics;
};

// Structural property contracts permit caller-defined records with quantities
// of the required kind, scalar and unit system.
template <typename P, typename T, typename U>
concept AcousticPropertiesOf = Scalar<T> && UnitSystem<U> && requires(const P& p) {
    { p.sound_speed } -> std::convertible_to<SoundSpeed<T, U>>;
    { p.acoustic_exponent } -> std::convertible_to<AcousticExponent<T>>;
};

template <typename P, typename T, typename U>
concept SingleTemperaturePropertiesOf = Scalar<T> && UnitSystem<U> && requires(const P& p) {
    { p.cv } -> std::convertible_to<SpecificHeat<T, U>>;
    { p.cp } -> std::convertible_to<SpecificHeat<T, U>>;
    { p.capacity_ratio } -> std::convertible_to<HeatCapacityRatio<T>>;
    { p.gas_constant } -> std::convertible_to<SpecificGasConstant<T, U>>;
    requires AcousticPropertiesOf<std::remove_cvref_t<decltype(p.acoustics)>, T, U>;
};

// mode_cv must contain one capacity for every mode in the corresponding model.
// The concept checks indexed access and quantity types; the producer is
// responsible for the mode count and ordering.
template <typename P, typename T, typename U>
concept MultiTemperaturePropertiesOf = Scalar<T> && UnitSystem<U> && requires(const P& p) {
    { p.cv_heavy } -> std::convertible_to<SpecificHeat<T, U>>;
    { p.mode_cv.size() } -> std::convertible_to<std::size_t>;
    { p.mode_cv[std::size_t {}] } -> std::convertible_to<SpecificHeat<T, U>>;
    { p.gas_constant_heavy } -> std::convertible_to<SpecificGasConstant<T, U>>;
    { p.gas_constant_electron } -> std::convertible_to<SpecificGasConstant<T, U>>;
    requires AcousticPropertiesOf<std::remove_cvref_t<decltype(p.acoustics)>, T, U>;
};

// Initial temperature for a closure that solves for the temperature
// iteratively. Passed to thermo_state after the state variables; models
// that accept it say so through AcceptsInput. Only its value part is used.
template <Scalar T, UnitSystem U> struct TemperatureGuess {
    Temperature<T, U> temp;
};

// Result of a closure: the state of the model and the status of the
// calculation. The state must not be used when the status is Invalid.
template <ThermoStateType S> struct ClosureOf {
    S state;
    Status status;
};

// Result of a closure of a single-temperature model.
template <Scalar T, UnitSystem U> using Closure = ClosureOf<ThermoState<T, U>>;

}  // namespace canary

#endif  // CANARY_STATE_HPP
