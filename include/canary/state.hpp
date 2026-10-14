#ifndef CANARY_STATE_HPP
#define CANARY_STATE_HPP

#include <array>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/status.hpp>
#include <concepts>
#include <cstddef>
#include <type_traits>

namespace canary
{

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

// Derived properties at a state.
template <Scalar T, UnitSystem U> struct Properties {
    SpecificHeat<T, U> cv;
    SpecificHeat<T, U> cp;
    HeatCapacityRatio<T> gamma;
    SoundSpeed<T, U> sound_speed;
    SpecificGasConstant<T, U> gas_constant;
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
