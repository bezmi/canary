#ifndef GASLIB_MODELS_TWO_TEMPERATURE_HPP
#define GASLIB_MODELS_TWO_TEMPERATURE_HPP

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>
#include <gaslib/build.hpp>
#include <gaslib/config.hpp>
#include <gaslib/detail/solve.hpp>
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/nasa9.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/scalar.hpp>
#include <gaslib/species.hpp>
#include <gaslib/state.hpp>
#include <gaslib/status.hpp>
#include <gaslib/thermo.hpp>
#include <gaslib/units.hpp>
#include <string_view>

namespace gaslib
{

// Ideal gas mixture of N species with two temperatures (Park's
// two-temperature model). The translational and rotational energy of every
// species follows the temperature T (State::temp); the vibrational and
// electronic energy follows the temperature T_ve (State::mode_temp[0]). For
// species i, with gas constant R_i, cp_i(T) and h_i(T) per unit mass from
// its NASA 9-coefficient curve, and T_ref the reference temperature:
// - cp_tr,i = (5/2 + d_i/2) R_i and cv_tr,i = cp_tr,i - R_i, where d_i is
//   the number of rotational degrees of freedom: 0 for an atom, 2 for a
//   linear molecule, 3 for a nonlinear molecule;
// - dh_f,i = h_i(T_ref);
// - e_ve,i(T_ve) = h_i(T_ve) - cp_tr,i (T_ve - T_ref) - dh_f,i, which is 0
//   at T_ref, and cv_ve,i(T_ve) = cp_i(T_ve) - cp_tr,i;
// - h_i(T, T_ve) = cp_tr,i (T - T_ref) + dh_f,i + e_ve,i(T_ve), which is
//   h_i(T) when T_ve = T, and e_i = h_i - R_i T.
// For mass fractions y_i and R = sum(y_i R_i): the internal energy per unit
// mass e = sum(y_i e_i) includes e_ve = sum(y_i e_ve,i) (State::mode_e[0]),
// and p = rho R T. Free electrons are not part of the model. All N species
// are used. The fields are in unit system U and are derived by
// make_two_temperature_mixture:
// - thermally_perfect: the single-temperature model of the same species. Its
//   gas constants, curves and reference temperature define the species
//   above; its temperature limits, bracket width, tolerance and iteration
//   limit are those of the T_ve solve; it gives the closures without mode
//   input;
// - translational_rotational_cp: cp_tr,i;
// - formation_enthalpy: dh_f,i.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S = 4>
    requires(N > 0 && S > 0)
struct TwoTemperatureMixture {
    using parameter_type = P;
    using units = U;
    template <Scalar T> using State = MultiTemperatureState<T, U, 1>;
    template <Scalar T> using Mixture = NormalisedComposition<T, U, N>;

    ThermallyPerfectMixture<P, U, N, S> thermally_perfect;
    std::array<SpecificHeat<P, U>, N> translational_rotational_cp;
    std::array<SpecificEnthalpy<P, U>, N> formation_enthalpy;
};

// Two-temperature model of the given species, with d_i in
// rotational_degrees, in the unit system of the scales, with its species
// data. The single-temperature model and the species data are those of
// make_thermally_perfect_mixture, with its temperature limits and solve
// settings. Throws ParameterError, with the species name in the message,
// when a species is named e-, a number of rotational degrees of freedom is
// not 0, 2 or 3, or make_thermally_perfect_mixture rejects a species.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires(N > 0 && S > 0)
Built<TwoTemperatureMixture<P, U, N, S>, SpeciesData<P, U, N>>
  make_two_temperature_mixture(const std::array<Nasa9Species<S>, N>& species,
                               const std::array<int, N>& rotational_degrees,
                               const UnitScales<P, U>& scales)
{
    for (std::size_t i = 0; i < N; ++i) {
        if (std::string_view(species[i].name) == "e-") {
            throw ParameterError(
              std::format("species {}: free electrons are not supported", species[i].name));
        }
        const int d = rotational_degrees[i];
        if (d != 0 && d != 2 && d != 3) {
            throw ParameterError(
              std::format("species {}: rotational degrees of freedom must be 0, 2 or 3, got {}",
                          species[i].name,
                          d));
        }
    }
    const auto single = make_thermally_perfect_mixture(species, scales);
    Built<TwoTemperatureMixture<P, U, N, S>, SpeciesData<P, U, N>> built {};
    auto& model = built.model;
    model.thermally_perfect = single.model;
    built.species = single.species;
    const P t_ref = single.model.reference_temperature.value();
    const P ln_t_ref = std::log(t_ref);
    for (std::size_t i = 0; i < N; ++i) {
        const P r = single.model.gas_constant[i].value();
        const auto d = static_cast<P>(rotational_degrees[i]);
        model.translational_rotational_cp[i] = SpecificHeat<P, U> {((P(5) + d) / P(2)) * r};
        model.formation_enthalpy[i] =
          SpecificEnthalpy<P, U> {r * nasa9_cp_h(single.model.curve[i], t_ref, ln_t_ref).h_over_r};
    }
    return built;
}

namespace detail
{

// sum(y_i cp_tr,i) and sum(y_i dh_f,i), with plain values of type W.
template <Scalar W> struct TranslationalRotationalSums {
    W cp;
    W formation_enthalpy;
};

template <Scalar W, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<W, P>
GASLIB_HOST_DEVICE TranslationalRotationalSums<W>
  translational_rotational_sums(const TwoTemperatureMixture<P, U, N, S>& gas,
                                const std::array<MassFraction<W>, N>& y)
{
    using V = ScalarTraits<W>::value_type;
    W cp = W(V(0));
    W h = W(V(0));
    for (std::size_t i = 0; i < N; ++i) {
        cp = cp + W(y[i].value() * gas.translational_rotational_cp[i].value());
        h = h + W(y[i].value() * gas.formation_enthalpy[i].value());
    }
    return {.cp = cp, .formation_enthalpy = h};
}

// e_ve = sum(y_i e_ve,i(t)) and cv_ve = sum(y_i cv_ve,i(t)) at the
// vibrational-electronic temperature t (model units), with plain values of
// type W. The curves are empirical fits and are evaluated on plain values,
// with one logarithm.
template <Scalar W, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<W, P> && HasLog<W>
GASLIB_HOST_DEVICE EnergyAndCv<W> vibrational_electronic_energy(
  const TwoTemperatureMixture<P, U, N, S>& gas, const std::array<MassFraction<W>, N>& y, W t)
{
    using V = ScalarTraits<W>::value_type;
    const auto& single = gas.thermally_perfect;
    const P t_ref = single.reference_temperature.value();
    const W ln_t = math::log(t);
    W e = W(V(0));
    W cv = W(V(0));
    for (std::size_t i = 0; i < N; ++i) {
        const Nasa9CpH<W> c = nasa9_cp_h(single.curve[i], t, ln_t);
        const P r = single.gas_constant[i].value();
        const P cp_tr = gas.translational_rotational_cp[i].value();
        const P h_f = gas.formation_enthalpy[i].value();
        e = e + W(y[i].value() * ((r * c.h_over_r) - (cp_tr * (t - t_ref)) - h_f));
        cv = cv + W(y[i].value() * ((r * c.cp_over_r) - cp_tr));
    }
    return {.e = e, .cv = cv};
}

// e = sum(y_i cp_tr,i) (T - T_ref) + sum(y_i dh_f,i) - R T + e_ve.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE SpecificEnergy<T, U>
  two_temperature_energy(const TwoTemperatureMixture<P, U, N, S>& gas,
                         const NormalisedComposition<T, U, N>& mix,
                         const TranslationalRotationalSums<T>& sums,
                         Temperature<T, U> temp,
                         SpecificEnergy<T, U> e_ve)
{
    const P t_ref = gas.thermally_perfect.reference_temperature.value();
    const T t = temp.value();
    const T e = T((sums.cp * T(t - t_ref)) + sums.formation_enthalpy -
                  (mix.gas_constant.value() * t) + e_ve.value());
    return SpecificEnergy<T, U> {e};
}

// The inverse of two_temperature_energy for T:
// T = (e - e_ve - sum(y_i dh_f,i) + T_ref sum(y_i cp_tr,i)) /
//     (sum(y_i cp_tr,i) - R).
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE Temperature<T, U>
  two_temperature_temperature(const TwoTemperatureMixture<P, U, N, S>& gas,
                              const NormalisedComposition<T, U, N>& mix,
                              const TranslationalRotationalSums<T>& sums,
                              SpecificEnergy<T, U> e,
                              SpecificEnergy<T, U> e_ve)
{
    const P t_ref = gas.thermally_perfect.reference_temperature.value();
    const T numerator = T(e.value() - e_ve.value() - sums.formation_enthalpy + T(t_ref * sums.cp));
    return Temperature<T, U> {T(numerator / (sums.cp - mix.gas_constant.value()))};
}

// Valid when the density and both temperatures have finite, positive value
// parts and both energies finite value parts (they may be negative);
// Invalid otherwise.
template <Scalar T, UnitSystem U>
GASLIB_HOST_DEVICE Status two_temperature_status(Density<T, U> rho,
                                                 SpecificEnergy<T, U> e,
                                                 Temperature<T, U> temp,
                                                 SpecificEnergy<T, U> e_ve,
                                                 ModeTemperature<T, U> t_ve)
{
    return combine(thermally_perfect_status(rho, e, temp),
                   thermally_perfect_status(rho, e_ve, Temperature<T, U> {t_ve.value()}));
}

// Closure of the state with the given members. The status combines
// two_temperature_status with status.
template <Scalar T, UnitSystem U>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  two_temperature_closure(Density<T, U> rho,
                          SpecificEnergy<T, U> e,
                          Pressure<T, U> p,
                          Temperature<T, U> temp,
                          SpecificEnergy<T, U> e_ve,
                          ModeTemperature<T, U> t_ve,
                          Status status)
{
    return {
      .state =
        {
          .rho = rho,
          .e = e,
          .p = p,
          .temp = temp,
          .mode_e = {e_ve},
          .mode_temp = {t_ve},
        },
      .status = combine(two_temperature_status(rho, e, temp, e_ve, t_ve), status),
    };
}

// The element of a mode range of one mode.
template <SpeciesRange V> GASLIB_HOST_DEVICE range_quantity_t<V> single_mode(const V& modes)
{
    // NOLINTNEXTLINE(misc-static-assert): run-time extent
    GASLIB_ASSERT(static_cast<std::size_t>(modes.extent(0)) == 1);
    return modes[0];
}

// T_ve from e_ve by solve_temperature of the single-temperature model on
// e_ve(T_ve) - e_ve = 0, with the derivative cv_ve(T_ve), starting from
// guess.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
GASLIB_HOST_DEVICE SolveResult<T>
  solve_mode_temperature(const TwoTemperatureMixture<P, U, N, S>& gas,
                         const NormalisedComposition<T, U, N>& mix,
                         SpecificEnergy<T, U> e_ve,
                         typename ScalarTraits<T>::value_type guess)
{
    using V = ScalarTraits<T>::value_type;
    const std::array<MassFraction<V>, N> y = fraction_values(mix);
    const V target = ScalarTraits<T>::value(e_ve.value());
    const auto f_value = [&gas, &y, target](V t) {
        const EnergyAndCv<V> r = vibrational_electronic_energy(gas, y, t);
        return Residual<V> {.value = r.e - target, .derivative = r.cv};
    };
    const auto f_full = [&gas, &mix, e_ve](T t) {
        const EnergyAndCv<T> r = vibrational_electronic_energy(gas, mix.fractions, t);
        return Residual<T> {.value = r.e - e_ve.value(), .derivative = r.cv};
    };
    return solve_temperature<T>(gas.thermally_perfect, f_value, f_full, guess);
}

// State from rho, T and T_ve: e_ve(T_ve), e from two_temperature_energy,
// p = rho R T.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  two_temperature_from_temperatures(const TwoTemperatureMixture<P, U, N, S>& gas,
                                    const NormalisedComposition<T, U, N>& mix,
                                    Density<T, U> rho,
                                    Temperature<T, U> temp,
                                    ModeTemperature<T, U> t_ve)
{
    const T e_ve_value = vibrational_electronic_energy(gas, mix.fractions, t_ve.value()).e;
    const SpecificEnergy<T, U> e_ve {e_ve_value};
    const TranslationalRotationalSums<T> sums = translational_rotational_sums(gas, mix.fractions);
    return two_temperature_closure(rho,
                                   two_temperature_energy(gas, mix, sums, temp, e_ve),
                                   Pressure<T, U> {pressure(rho, mix.gas_constant, temp)},
                                   temp,
                                   e_ve,
                                   t_ve,
                                   mix.status);
}

// State from rho, T and e_ve: T_ve by solve_mode_temperature starting from
// T, e from two_temperature_energy, p = rho R T.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  two_temperature_from_mode_energy(const TwoTemperatureMixture<P, U, N, S>& gas,
                                   const NormalisedComposition<T, U, N>& mix,
                                   Density<T, U> rho,
                                   Temperature<T, U> temp,
                                   SpecificEnergy<T, U> e_ve)
{
    const SolveResult<T> r =
      solve_mode_temperature(gas, mix, e_ve, ScalarTraits<T>::value(temp.value()));
    const TranslationalRotationalSums<T> sums = translational_rotational_sums(gas, mix.fractions);
    return two_temperature_closure(rho,
                                   two_temperature_energy(gas, mix, sums, temp, e_ve),
                                   Pressure<T, U> {pressure(rho, mix.gas_constant, temp)},
                                   temp,
                                   e_ve,
                                   ModeTemperature<T, U> {r.root},
                                   combine(r.status, mix.status));
}

// e_ve = (rho e_ve) / rho from the mode energy density.
template <Scalar T, UnitSystem U>
GASLIB_HOST_DEVICE SpecificEnergy<T, U> mode_specific_energy(ModeEnergyDensity<T, U> rho_e_ve,
                                                             Density<T, U> rho)
{ return SpecificEnergy<T, U> {T(rho_e_ve.value() / rho.value())}; }

// The closure c of the single-temperature model at thermal equilibrium:
// T_ve = T and e_ve = e_ve(T), with the state and status of c.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  at_thermal_equilibrium(const TwoTemperatureMixture<P, U, N, S>& gas,
                         const NormalisedComposition<T, U, N>& mix,
                         const Closure<T, U>& c)
{
    const T e_ve_value = vibrational_electronic_energy(gas, mix.fractions, c.state.temp.value()).e;
    const SpecificEnergy<T, U> e_ve {e_ve_value};
    return {
      .state =
        {
          .rho = c.state.rho,
          .e = c.state.e,
          .p = c.state.p,
          .temp = c.state.temp,
          .mode_e = {e_ve},
          .mode_temp = {ModeTemperature<T, U> {c.state.temp.value()}},
        },
      .status = c.status,
    };
}

}  // namespace detail

// Per-cell composition from mass fractions, as for the single-temperature
// model.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V>
    requires SpeciesWeights<V, kind::mass_fraction, unitless, P> &&
             ScalarWith<typename range_quantity_t<V>::value_type, P>
GASLIB_HOST_DEVICE auto mixture(const TwoTemperatureMixture<P, U, N, S>& gas, const V& fractions)
{ return mixture(gas.thermally_perfect, fractions); }

// Per-cell composition from partial densities and the density, as for the
// single-temperature model.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V, Scalar T>
    requires SpeciesWeights<V, kind::density, U, P> &&
             std::same_as<typename range_quantity_t<V>::value_type, T> && ScalarWith<T, P>
GASLIB_HOST_DEVICE auto mixture(const TwoTemperatureMixture<P, U, N, S>& gas,
                                const V& partial_densities,
                                Density<T, U> rho)
{ return mixture(gas.thermally_perfect, partial_densities, rho); }

// Closures at thermal equilibrium (no mode input): T, p, e and the status
// are those of the single-temperature model; T_ve = T and e_ve = e_ve(T).

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        SpecificEnergy<T, U> e)
{ return detail::at_thermal_equilibrium(gas, mix, close(gas.thermally_perfect, mix, rho, e)); }

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Pressure<T, U> p)
{ return detail::at_thermal_equilibrium(gas, mix, close(gas.thermally_perfect, mix, rho, p)); }

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Temperature<T, U> temp)
{ return detail::at_thermal_equilibrium(gas, mix, close(gas.thermally_perfect, mix, rho, temp)); }

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Pressure<T, U> p,
        Temperature<T, U> temp)
{ return detail::at_thermal_equilibrium(gas, mix, close(gas.thermally_perfect, mix, p, temp)); }

// Closures with the vibrational-electronic temperature T_ve given as a mode
// temperature range of one element. e_ve = e_ve(T_ve). The status combines
// detail::two_temperature_status and the composition status.

// T in closed form from e and e_ve (detail::two_temperature_temperature),
// p = rho R T; e is kept.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V>
    requires ScalarWith<T, P> && HasLog<T> && ModeRangeOf<V, kind::mode_temperature, T, U>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        SpecificEnergy<T, U> e,
        const V& mode_temperatures)
{
    const ModeTemperature<T, U> t_ve = detail::single_mode(mode_temperatures);
    const T e_ve_value = detail::vibrational_electronic_energy(gas, mix.fractions, t_ve.value()).e;
    const SpecificEnergy<T, U> e_ve {e_ve_value};
    const auto sums = detail::translational_rotational_sums(gas, mix.fractions);
    const Temperature<T, U> temp = detail::two_temperature_temperature(gas, mix, sums, e, e_ve);
    return detail::two_temperature_closure(rho,
                                           e,
                                           Pressure<T, U> {pressure(rho, mix.gas_constant, temp)},
                                           temp,
                                           e_ve,
                                           t_ve,
                                           mix.status);
}

// T = p / (rho R); p is kept.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V>
    requires ScalarWith<T, P> && HasLog<T> && ModeRangeOf<V, kind::mode_temperature, T, U>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Pressure<T, U> p,
        const V& mode_temperatures)
{
    const Temperature<T, U> temp {temperature(p, rho, mix.gas_constant)};
    auto c = detail::two_temperature_from_temperatures(gas,
                                                       mix,
                                                       rho,
                                                       temp,
                                                       detail::single_mode(mode_temperatures));
    c.state.p = p;
    return c;
}

// p = rho R T.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V>
    requires ScalarWith<T, P> && HasLog<T> && ModeRangeOf<V, kind::mode_temperature, T, U>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Temperature<T, U> temp,
        const V& mode_temperatures)
{
    return detail::two_temperature_from_temperatures(gas,
                                                     mix,
                                                     rho,
                                                     temp,
                                                     detail::single_mode(mode_temperatures));
}

// rho = p / (R T); p is kept.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V>
    requires ScalarWith<T, P> && HasLog<T> && ModeRangeOf<V, kind::mode_temperature, T, U>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Pressure<T, U> p,
        Temperature<T, U> temp,
        const V& mode_temperatures)
{
    const Density<T, U> rho {density(p, mix.gas_constant, temp)};
    auto c = detail::two_temperature_from_temperatures(gas,
                                                       mix,
                                                       rho,
                                                       temp,
                                                       detail::single_mode(mode_temperatures));
    c.state.p = p;
    return c;
}

// Closures with the vibrational-electronic energy given as a mode energy
// density range of one element, rho e_ve; e_ve = (rho e_ve) / rho. T_ve is
// found by detail::solve_mode_temperature, starting from T. The status
// combines detail::two_temperature_status, the solve status (NotConverged,
// or Invalid when no bracket has a sign change) and the composition status.

// T in closed form from e and e_ve (detail::two_temperature_temperature),
// p = rho R T; e is kept.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V>
    requires ScalarWith<T, P> && HasLog<T> && ModeRangeOf<V, kind::mode_energy_density, T, U>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        SpecificEnergy<T, U> e,
        const V& mode_energy_densities)
{
    const SpecificEnergy<T, U> e_ve =
      detail::mode_specific_energy(detail::single_mode(mode_energy_densities), rho);
    const auto sums = detail::translational_rotational_sums(gas, mix.fractions);
    const Temperature<T, U> temp = detail::two_temperature_temperature(gas, mix, sums, e, e_ve);
    const detail::SolveResult<T> r =
      detail::solve_mode_temperature(gas, mix, e_ve, ScalarTraits<T>::value(temp.value()));
    return detail::two_temperature_closure(rho,
                                           e,
                                           Pressure<T, U> {pressure(rho, mix.gas_constant, temp)},
                                           temp,
                                           e_ve,
                                           ModeTemperature<T, U> {r.root},
                                           combine(r.status, mix.status));
}

// T = p / (rho R); p is kept.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V>
    requires ScalarWith<T, P> && HasLog<T> && ModeRangeOf<V, kind::mode_energy_density, T, U>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Pressure<T, U> p,
        const V& mode_energy_densities)
{
    const Temperature<T, U> temp {temperature(p, rho, mix.gas_constant)};
    auto c = detail::two_temperature_from_mode_energy(
      gas,
      mix,
      rho,
      temp,
      detail::mode_specific_energy(detail::single_mode(mode_energy_densities), rho));
    c.state.p = p;
    return c;
}

// p = rho R T.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V>
    requires ScalarWith<T, P> && HasLog<T> && ModeRangeOf<V, kind::mode_energy_density, T, U>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Density<T, U> rho,
        Temperature<T, U> temp,
        const V& mode_energy_densities)
{
    return detail::two_temperature_from_mode_energy(
      gas,
      mix,
      rho,
      temp,
      detail::mode_specific_energy(detail::single_mode(mode_energy_densities), rho));
}

// rho = p / (R T), then e_ve = (rho e_ve) / rho; p is kept.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V>
    requires ScalarWith<T, P> && HasLog<T> && ModeRangeOf<V, kind::mode_energy_density, T, U>
GASLIB_HOST_DEVICE ClosureOf<MultiTemperatureState<T, U, 1>>
  close(const TwoTemperatureMixture<P, U, N, S>& gas,
        const NormalisedComposition<T, U, N>& mix,
        Pressure<T, U> p,
        Temperature<T, U> temp,
        const V& mode_energy_densities)
{
    const Density<T, U> rho {density(p, mix.gas_constant, temp)};
    auto c = detail::two_temperature_from_mode_energy(
      gas,
      mix,
      rho,
      temp,
      detail::mode_specific_energy(detail::single_mode(mode_energy_densities), rho));
    c.state.p = p;
    return c;
}

// cv = sum(y_i (cv_tr,i + cv_ve,i(T_ve))), cp = cv + R, gamma = cp / cv,
// the frozen sound speed a = sqrt(gamma p / rho), and R.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T> && HasSqrt<T>
GASLIB_HOST_DEVICE Properties<T, U> properties(const TwoTemperatureMixture<P, U, N, S>& gas,
                                               const NormalisedComposition<T, U, N>& mix,
                                               const MultiTemperatureState<T, U, 1>& state)
{
    const auto sums = detail::translational_rotational_sums(gas, mix.fractions);
    const T cv_ve =
      detail::vibrational_electronic_energy(gas, mix.fractions, state.mode_temp[0].value()).cv;
    const SpecificHeat<T, U> cv {T((sums.cp - mix.gas_constant.value()) + cv_ve)};
    const SpecificHeat<T, U> cp {cv.value() + mix.gas_constant.value()};
    const HeatCapacityRatio<T> gamma {heat_capacity_ratio(cp, cv)};
    return {
      .cv = cv,
      .cp = cp,
      .gamma = gamma,
      .sound_speed = SoundSpeed<T, U> {sound_speed(gamma, state.p, state.rho)},
      .gas_constant = mix.gas_constant,
    };
}

}  // namespace gaslib

#endif  // GASLIB_MODELS_TWO_TEMPERATURE_HPP
