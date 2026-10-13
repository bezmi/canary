#ifndef CANARY_STATE_HPP
#define CANARY_STATE_HPP

#include <array>
#include <concepts>
#include <cstddef>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/status.hpp>

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
