#ifndef GASLIB_STATE_HPP
#define GASLIB_STATE_HPP

#include <concepts>
#include <gaslib/quantity.hpp>
#include <gaslib/scalar.hpp>
#include <gaslib/status.hpp>

namespace gaslib
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

// Derived properties at a state.
template <Scalar T, UnitSystem U> struct Properties {
    SpecificHeat<T, U> cv;
    SpecificHeat<T, U> cp;
    HeatCapacityRatio<T> gamma;
    SoundSpeed<T, U> sound_speed;
    SpecificGasConstant<T, U> gas_constant;
};

// Result of a closure: the state of the model and the status of the
// calculation. The state must not be used when the status is Invalid.
template <ThermoStateType S> struct ClosureOf {
    S state;
    Status status;
};

// Result of a closure of a single-temperature model.
template <Scalar T, UnitSystem U> using Closure = ClosureOf<ThermoState<T, U>>;

}  // namespace gaslib

#endif  // GASLIB_STATE_HPP
