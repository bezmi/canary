#ifndef GASLIB_TRANSPORT_TRANSPORT_HPP
#define GASLIB_TRANSPORT_TRANSPORT_HPP

#include <concepts>
#include <gaslib/config.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/scalar.hpp>
#include <gaslib/state.hpp>
#include <gaslib/thermo.hpp>
#include <type_traits>

namespace gaslib
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
  requires(const K& law, DynamicViscosity<T, U> mu, const Properties<T, U>& properties) {
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

namespace detail
{
template <TransportLaw V, TransportLaw K, Scalar T, UnitSystem U>
    requires ViscosityModel<V, T, U> && ConductivityModel<K, T, U>
GASLIB_HOST_DEVICE TransportCoefficients<T, U> transport_from(const Transport<V, K>& model,
                                                              Temperature<T, U> temp,
                                                              const Properties<T, U>& properties)
{
    const DynamicViscosity<T, U> mu = viscosity(model.viscosity, temp);
    return {.viscosity = mu, .conductivity = conductivity(model.conductivity, mu, properties)};
}
}  // namespace detail

// Viscosity and conductivity at a state of a single-species gas model.
template <TransportLaw V, TransportLaw K, GasModelType G, ThermoStateType S>
    requires ViscosityModel<V, typename S::value_type, typename S::units> &&
             ConductivityModel<K, typename S::value_type, typename S::units> &&
             detail::PropertiesAlone<G, typename S::value_type>
GASLIB_HOST_DEVICE auto transport(const Transport<V, K>& model, const G& gas, const S& state)
{ return detail::transport_from(model, state.temp, properties(gas, state)); }

// Viscosity and conductivity at a state of a mixture gas model.
template <TransportLaw V, TransportLaw K, GasModelType G, ThermoStateType S>
    requires HasMixture<G, typename S::value_type> &&
             ViscosityModel<V, typename S::value_type, typename S::units> &&
             ConductivityModel<K, typename S::value_type, typename S::units> &&
             detail::PropertiesWithMixture<G, typename S::value_type>
GASLIB_HOST_DEVICE auto transport(const Transport<V, K>& model,
                                  const G& gas,
                                  const mixture_t<G, typename S::value_type>& mix,
                                  const S& state)
{ return detail::transport_from(model, state.temp, properties(gas, mix, state)); }

namespace detail
{
template <typename M, typename G, typename T>
concept TransportAlone = requires(const M& model, const G& gas, const state_t<G, T>& state) {
    { transport(model, gas, state) } -> std::same_as<TransportCoefficients<T, typename G::units>>;
};

template <typename M, typename G, typename T>
concept TransportWithMixture =
  requires(const M& model, const G& gas, const mixture_t<G, T>& mix, const state_t<G, T>& state) {
      {
          transport(model, gas, mix, state)
      } -> std::same_as<TransportCoefficients<T, typename G::units>>;
  };
}  // namespace detail

// Transport model M gives viscosity and conductivity at a state of gas model
// G with value type T (with the mixture value for a mixture model).
template <typename M, typename G, typename T>
concept TransportModel = GasModel<G, T> && std::is_trivially_copyable_v<M> &&
                         ((HasMixture<G, T> && detail::TransportWithMixture<M, G, T>) ||
                          (!HasMixture<G, T> && detail::TransportAlone<M, G, T>));

}  // namespace gaslib

#endif  // GASLIB_TRANSPORT_TRANSPORT_HPP
