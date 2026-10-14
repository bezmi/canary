#ifndef CANARY_TRANSPORT_TRANSPORT_HPP
#define CANARY_TRANSPORT_TRANSPORT_HPP

#include <array>
#include <concepts>
#include <cstddef>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/state.hpp>
#include <canary/thermo.hpp>
#include <type_traits>

namespace canary
{

// A viscosity law: viscosity(law, temp) for state value type T in unit
// system U.
template <typename V, typename T, typename U>
concept ViscosityModel = std::is_trivially_copyable_v<V> && Scalar<T> && UnitSystem<U> &&
                         requires(const V& law, Temperature<T, U> temp) {
                             { viscosity(law, temp) } -> std::same_as<DynamicViscosity<T, U>>;
                         };

// A conductivity law: conductivity(law, mu, properties) for state value
// type T in unit system U. The law reads what it needs from the viscosity
// and the gas properties at the state.
template <typename K, typename T, typename U>
concept ConductivityModel =
  std::is_trivially_copyable_v<K> && Scalar<T> && UnitSystem<U> &&
  requires(
    const K& law, DynamicViscosity<T, U> mu, const SingleTemperatureProperties<T, U>& properties) {
      { conductivity(law, mu, properties) } -> std::same_as<ThermalConductivity<T, U>>;
  };

// A law used in a Transport pair.
template <typename L>
concept TransportLaw = std::is_trivially_copyable_v<L>;

// A viscosity law and a conductivity law used together.
template <TransportLaw V, TransportLaw K> struct Transport {
    V viscosity;
    K conductivity;
};

// Transport coefficients at a state.
template <Scalar T, UnitSystem U> struct TransportCoefficients {
    DynamicViscosity<T, U> viscosity;
    ThermalConductivity<T, U> conductivity;
};

// Transport coefficients at a state with M internal energy modes:
// conductivity carries the energy that follows the temperature temp of the
// state, mode_conductivity[k] the energy of mode k, which follows
// mode_temp[k].
template <Scalar T, UnitSystem U, std::size_t M>
    requires(M > 0)
struct MultiTemperatureTransportCoefficients {
    DynamicViscosity<T, U> viscosity;
    ThermalConductivity<T, U> conductivity;
    std::array<ThermalConductivity<T, U>, M> mode_conductivity;
};

namespace detail
{
template <TransportLaw V, TransportLaw K, Scalar T, UnitSystem U, typename Props>
    requires ViscosityModel<V, T, U> && PropertiesOf<Props, T, U> &&
             requires(const K& law, DynamicViscosity<T, U> mu, const Props& properties) {
                 { conductivity(law, mu, properties) } -> std::same_as<ThermalConductivity<T, U>>;
             }
CANARY_HOST_DEVICE TransportCoefficients<T, U>
  transport_from(const Transport<V, K>& model, Temperature<T, U> temp, const Props& properties)
{
    const DynamicViscosity<T, U> mu = viscosity(model.viscosity, temp);
    return {.viscosity = mu, .conductivity = conductivity(model.conductivity, mu, properties)};
}
}  // namespace detail

// Viscosity and conductivity at a state of a single-species gas model.
template <TransportLaw V, TransportLaw K, GasModelType G, ThermoStateType S>
    requires ViscosityModel<V, typename S::value_type, typename S::units> &&
             ConductivityModel<K, typename S::value_type, typename S::units> &&
             detail::PropertiesAlone<G, typename S::value_type> &&
             requires(const Transport<V, K>& model, const G& gas, const S& state) {
                 detail::transport_from(model, state.temp, properties(gas, state));
             }
CANARY_HOST_DEVICE auto transport(const Transport<V, K>& model, const G& gas, const S& state)
{ return detail::transport_from(model, state.temp, properties(gas, state)); }

// Viscosity and conductivity at a state of a mixture gas model.
template <TransportLaw V, TransportLaw K, GasModelType G, ThermoStateType S>
    requires HasMixture<G, typename S::value_type> &&
             ViscosityModel<V, typename S::value_type, typename S::units> &&
             ConductivityModel<K, typename S::value_type, typename S::units> &&
             detail::PropertiesWithMixture<G, typename S::value_type> &&
             requires(const Transport<V, K>& model,
                      const G& gas,
                      const mixture_t<G, typename S::value_type>& mix,
                      const S& state) {
                 detail::transport_from(model, state.temp, properties(gas, mix, state));
             }
CANARY_HOST_DEVICE auto transport(const Transport<V, K>& model,
                                  const G& gas,
                                  const mixture_t<G, typename S::value_type>& mix,
                                  const S& state)
{ return detail::transport_from(model, state.temp, properties(gas, mix, state)); }

namespace detail
{
// The coefficients with one conductivity per mode of a state type S;
// void for a state without modes.
template <typename S> struct mode_coefficients {
    using type = void;
};

template <Scalar T, UnitSystem U, std::size_t M>
struct mode_coefficients<MultiTemperatureState<T, U, M>> {
    using type = MultiTemperatureTransportCoefficients<T, U, M>;
};

// C is what transport may return at a state of type S: TransportCoefficients,
// or for a state with modes also MultiTemperatureTransportCoefficients.
template <typename C, typename S>
concept TransportCoefficientsFor =
  ThermoStateType<S> &&
  (std::same_as<C, TransportCoefficients<typename S::value_type, typename S::units>> ||
   std::same_as<C, typename mode_coefficients<S>::type>);

template <typename M, typename G, typename T>
concept TransportAlone = requires(const M& model, const G& gas, const state_t<G, T>& state) {
    { transport(model, gas, state) } -> TransportCoefficientsFor<state_t<G, T>>;
};

template <typename M, typename G, typename T>
concept TransportWithMixture =
  requires(const M& model, const G& gas, const mixture_t<G, T>& mix, const state_t<G, T>& state) {
      { transport(model, gas, mix, state) } -> TransportCoefficientsFor<state_t<G, T>>;
  };
}  // namespace detail

// Transport model M gives viscosity and conductivity at a state of gas model
// G with value type T (with the mixture value for a mixture model): as
// TransportCoefficients, or, at a state with internal energy modes, as
// TransportCoefficients or MultiTemperatureTransportCoefficients.
template <typename M, typename G, typename T>
concept TransportModel = GasModel<G, T> && std::is_trivially_copyable_v<M> &&
                         ((HasMixture<G, T> && detail::TransportWithMixture<M, G, T>) ||
                          (!HasMixture<G, T> && detail::TransportAlone<M, G, T>));

}  // namespace canary

#endif  // CANARY_TRANSPORT_TRANSPORT_HPP
