#ifndef CANARY_MODELS_THERMALLY_PERFECT_HPP
#define CANARY_MODELS_THERMALLY_PERFECT_HPP

#include <algorithm>
#include <array>
#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/detail/solve.hpp>
#include <canary/nasa9.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/status.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>

namespace canary
{

// Host input for one species of a thermally perfect mixture: its name (for
// error messages), molar mass [kg/mol], charge [elementary charges] and fit.
template <std::size_t S>
    requires(S > 0)
struct Nasa9Species {
    const char* name;
    double molar_mass;
    double charge;
    Nasa9Data<S> fit;
};

// Per-cell composition of a thermally perfect mixture: the mass fractions,
// each weight clamped to [0, total] and divided by the sum of the clamped
// weights, the mixture gas constant R = sum(y_i R_i), and the status of the
// composition they came from.
template <Scalar T, UnitSystem U, std::size_t N>
    requires(N > 0)
struct NormalisedComposition {
    int n;
    std::array<MassFraction<T>, N> fractions;
    SpecificGasConstant<T, U> gas_constant;
    Status status;
};

// Ideal gas mixture of N species at one temperature, with the heat capacity
// of each species a function of temperature given by a NASA 9-coefficient
// curve of up to S segments. All N species are used; n equals N. The
// fields are in unit system U and are derived by
// make_thermally_perfect_mixture:
// - gas_constant: R_i = R_u / M_i;
// - curve: the fits folded into the temperature unit of U;
// - reference_temperature (298.15 K), and the species energy and cv there,
//   from which the closure from (rho, e) forms its default initial guess;
// - temperature_min and temperature_max: the full bracket of the
//   temperature solve; bracket_width: the width of the first bracket
//   around the guess;
// - relative_tolerance and max_iterations: settings of the solve.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S = 4>
    requires(N > 0 && S > 0)
struct ThermallyPerfectMixture {
    using parameter_type = P;
    using units = U;
    template <Scalar T> using State = ThermoState<T, U>;
    template <Scalar T> using Mixture = NormalisedComposition<T, U, N>;

    int n;
    std::array<SpecificGasConstant<P, U>, N> gas_constant;
    std::array<Nasa9Curve<P, S>, N> curve;
    Temperature<P, U> reference_temperature;
    std::array<SpecificEnergy<P, U>, N> reference_energy;
    std::array<SpecificHeat<P, U>, N> reference_cv;
    Temperature<P, U> temperature_min;
    Temperature<P, U> temperature_max;
    Temperature<P, U> bracket_width;
    P relative_tolerance;
    int max_iterations;
};

namespace detail
{

// Fills entry i of the model and of the species data from one species, after
// the model's reference temperature is set. Throws ParameterError.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
void fill_thermally_perfect_species(ThermallyPerfectMixture<P, U, N, S>& model,
                                    SpeciesData<P, U, N>& species,
                                    std::size_t i,
                                    const Nasa9Species<S>& input,
                                    const UnitScales<P, U>& scales)
{
    if (!std::isfinite(input.molar_mass) || !(input.molar_mass > 0.0)) {
        throw ParameterError(
          std::format("molar mass must be finite and positive, got {}", input.molar_mass));
    }
    if (!std::isfinite(input.charge)) {
        throw ParameterError(std::format("charge must be finite, got {}", input.charge));
    }
    const P tau = scale_factor<kind::temperature>(scales);
    const auto gas_constant =
      static_cast<P>(constants::universal_gas_constant<double> / input.molar_mass);
    const auto mass = static_cast<P>(input.molar_mass / constants::avogadro<double>);
    const auto charge = static_cast<P>(input.charge * constants::elementary_charge<double>);
    model.gas_constant[i] = convert<U>(SpecificGasConstant<P, si> {gas_constant}, scales);
    model.curve[i] = make_nasa9_curve(input.fit, tau);
    species.mass[i] = convert<U>(ParticleMass<P, si> {mass}, scales);
    species.charge[i] = convert<U>(Charge<P, si> {charge}, scales);

    const P t_ref = model.reference_temperature.value();
    const Nasa9CpH<P> ref = nasa9_cp_h(model.curve[i], t_ref, std::log(t_ref));
    const P r = model.gas_constant[i].value();
    model.reference_energy[i] = SpecificEnergy<P, U> {r * (ref.h_over_r - t_ref)};
    model.reference_cv[i] = SpecificHeat<P, U> {r * (ref.cp_over_r - P(1))};
}

}  // namespace detail

// Mixture model of the given species in the unit system of the scales, with
// its species data (particle mass M / N_A and charge q e, in U). The
// temperature limits are 10 K and 100000 K, the first bracket is 1000 K
// wide, the relative tolerance is 1e-10 and the iteration limit 30. Throws
// ParameterError, with the species name in the message, when a molar mass
// is not finite and positive, a charge is not finite, or a fit is rejected
// by make_nasa9_curve.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires(N > 0 && S > 0)
Built<ThermallyPerfectMixture<P, U, N, S>, SpeciesData<P, U, N>>
  make_thermally_perfect_mixture(const std::array<Nasa9Species<S>, N>& species,
                                 const UnitScales<P, U>& scales)
{
    Built<ThermallyPerfectMixture<P, U, N, S>, SpeciesData<P, U, N>> built {};
    auto& model = built.model;
    model.n = static_cast<int>(N);
    built.species.n = static_cast<int>(N);
    model.reference_temperature = convert<U>(Temperature<P, si> {P(298.15)}, scales);
    model.temperature_min = convert<U>(Temperature<P, si> {P(10)}, scales);
    model.temperature_max = convert<U>(Temperature<P, si> {P(100000)}, scales);
    model.bracket_width = convert<U>(Temperature<P, si> {P(1000)}, scales);
    model.relative_tolerance = P(1.0e-10);
    model.max_iterations = 30;
    for (std::size_t i = 0; i < N; ++i) {
        try {
            detail::fill_thermally_perfect_species(model, built.species, i, species[i], scales);
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("species {}: {}", species[i].name, error.what()));
        }
    }
    return built;
}

namespace detail
{

// Normalised fractions w_i / sum(w) of the clamped weights and the
// mixture gas constant.
template <std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t S,
          SpeciesRange V,
          Scalar T>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE NormalisedComposition<T, U, N> normalise(
  const ThermallyPerfectMixture<P, U, N, S>& gas, const V& weights, T total, Status status)
{
    // NOLINTNEXTLINE(misc-static-assert): run-time extent
    CANARY_ASSERT(static_cast<std::size_t>(weights.extent(0)) == N);
    std::array<T, N> w {};
    T sum = T(typename ScalarTraits<T>::value_type(0));
    for (std::size_t i = 0; i < N; ++i) {
        w[i] = clamped_weight(weights, i, total);
        sum = sum + w[i];
    }
    NormalisedComposition<T, U, N> mix {};
    mix.n = gas.n;
    mix.status = status;
    for (std::size_t i = 0; i < N; ++i) { mix.fractions[i] = MassFraction<T> {w[i] / sum}; }
    auto gas_constant = gas.gas_constant[0] * mix.fractions[0].value();
    for (std::size_t i = 1; i < N; ++i) {
        gas_constant = gas_constant + (gas.gas_constant[i] * mix.fractions[i].value());
    }
    mix.gas_constant = SpecificGasConstant<T, U> {gas_constant};
    return mix;
}

// Energy per unit mass and cv of a mixture at one temperature.
template <Scalar W> struct EnergyAndCv {
    W e;
    W cv;
};

// e = sum(y_i R_i (h_i/R_i - t)) and cv = sum(y_i R_i (cp_i/R_i - 1)) at
// temperature t (model units), with plain values of type W. Every species'
// curve is evaluated at t, with one logarithm.

template <Scalar W, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<W, P> && HasLog<W>
CANARY_HOST_DEVICE EnergyAndCv<W> mixture_energy(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                                 const std::array<MassFraction<W>, N>& y,
                                                 W t)
{
    using V = ScalarTraits<W>::value_type;
    const W ln_t = math::log(t);
    W e = W(V(0));
    W cv = W(V(0));
    for (std::size_t i = 0; i < N; ++i) {
        const Nasa9CpH<W> c = nasa9_cp_h(gas.curve[i], t, ln_t);
        const P r = gas.gas_constant[i].value();
        e = e + W(y[i].value() * (r * (c.h_over_r - t)));
        cv = cv + W(y[i].value() * (r * (c.cp_over_r - P(1))));
    }
    return {.e = e, .cv = cv};
}

// Valid when density and temperature have finite, positive value parts and
// the energy a finite value part (it may be negative); Invalid otherwise.
template <Scalar T, UnitSystem U>
CANARY_HOST_DEVICE Status thermally_perfect_status(Density<T, U> rho,
                                                   SpecificEnergy<T, U> e,
                                                   Temperature<T, U> temp)
{
    using V = ScalarTraits<T>::value_type;
    const auto positive = [](T x) {
        const V v = ScalarTraits<T>::value(x);
        return std::isfinite(v) && v > V(0);
    };
    return positive(rho.value()) && positive(temp.value()) && math::isfinite(e.value())
             ? Status::Valid
             : Status::Invalid;
}

// State validation includes the actual returned pressure. The three-field
// overload also serves the independently evaluated energy modes.
template <Scalar T, UnitSystem U>
CANARY_HOST_DEVICE Status thermally_perfect_status(Density<T, U> rho,
                                                   SpecificEnergy<T, U> e,
                                                   Temperature<T, U> temp,
                                                   Pressure<T, U> p)
{
    using V = ScalarTraits<T>::value_type;
    const V value = ScalarTraits<T>::value(p.value());
    const Status pressure_status =
      std::isfinite(value) && value > V(0) ? Status::Valid : Status::Invalid;
    return combine(thermally_perfect_status(rho, e, temp), pressure_status);
}

template <Scalar T, UnitSystem U, std::size_t N>
CANARY_HOST_DEVICE std::array<MassFraction<typename ScalarTraits<T>::value_type>, N>
  fraction_values(const NormalisedComposition<T, U, N>& mix)
{
    using V = ScalarTraits<T>::value_type;
    std::array<MassFraction<V>, N> y {};
    for (std::size_t i = 0; i < N; ++i) {
        y[i] = MassFraction<V> {ScalarTraits<T>::value(mix.fractions[i].value())};
    }
    return y;
}

// Default initial temperature, from the value parts:
// T_ref + (e - sum(y_i e_i(T_ref))) / sum(y_i cv_i(T_ref)).
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, std::floating_point V>
CANARY_HOST_DEVICE V default_temperature_guess(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                               const std::array<MassFraction<V>, N>& y,
                                               V e)
{
    V e_ref = V(0);
    V cv_ref = V(0);
    for (std::size_t i = 0; i < N; ++i) {
        e_ref += y[i].value() * static_cast<V>(gas.reference_energy[i].value());
        cv_ref += y[i].value() * static_cast<V>(gas.reference_cv[i].value());
    }
    return static_cast<V>(gas.reference_temperature.value()) + ((e - e_ref) / cv_ref);
}

// Temperature root of a residual in temperature (model units) by
// detail::solve, with the model's limits and solve settings. The guess
// (value part) is clamped to [temperature_min, temperature_max]; a guess
// that is not a number becomes temperature_min. The first bracket is
// [lower, upper], where lower = guess - bracket_width / 2 and
// upper = lower + bracket_width, each clamped to the configured limits.
// When it has no sign change, the solve is repeated on
// [temperature_min, temperature_max].
template <Scalar T,
          std::floating_point P,
          UnitSystem U,
          std::size_t N,
          std::size_t S,
          typename FValue,
          typename FFull>
    requires ResidualFunction<FValue, typename ScalarTraits<T>::value_type> &&
             ResidualFunction<FFull, T>
CANARY_HOST_DEVICE SolveResult<T> solve_temperature(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                                    const FValue& f_value,
                                                    const FFull& f_full,
                                                    typename ScalarTraits<T>::value_type guess)
{
    using V = ScalarTraits<T>::value_type;
    const auto t_min = static_cast<V>(gas.temperature_min.value());
    const auto t_max = static_cast<V>(gas.temperature_max.value());
    const auto width = static_cast<V>(gas.bracket_width.value());
    V start = guess > t_max ? t_max : guess;
    if (!(start >= t_min)) { start = t_min; }
    const V lower = std::clamp(start - (width / V(2)), t_min, t_max);
    const V upper = std::clamp(lower + width, t_min, t_max);
    const SolveSettings<V> settings {
      .relative_tolerance = static_cast<V>(gas.relative_tolerance),
      .max_iterations = gas.max_iterations,
    };
    SolveResult<T> r =
      solve<T>(f_value, f_full, Bracket<V> {.lower = lower, .upper = upper}, start, settings);
    if (r.status == Status::Invalid) {
        r = solve<T>(f_value, f_full, Bracket<V> {.lower = t_min, .upper = t_max}, start, settings);
    }
    return r;
}

// Temperature from (rho, e) by solve_temperature on e(T) - e = 0, with the
// derivative cv(T). Then p = rho R T, and e is kept.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE Closure<T, U>
  thermally_perfect_from_rho_e(const ThermallyPerfectMixture<P, U, N, S>& gas,
                               const NormalisedComposition<T, U, N>& mix,
                               Density<T, U> rho,
                               SpecificEnergy<T, U> e,
                               typename ScalarTraits<T>::value_type guess)
{
    using V = ScalarTraits<T>::value_type;
    const std::array<MassFraction<V>, N> y = fraction_values(mix);
    const V target = ScalarTraits<T>::value(e.value());
    const auto f_value = [&gas, &y, target](V t) {
        const EnergyAndCv<V> r = mixture_energy(gas, y, t);
        return Residual<V> {.value = r.e - target, .derivative = r.cv};
    };
    const auto f_full = [&gas, &mix, e](T t) {
        const EnergyAndCv<T> r = mixture_energy(gas, mix.fractions, t);
        return Residual<T> {.value = r.e - e.value(), .derivative = r.cv};
    };
    const SolveResult<T> r = solve_temperature<T>(gas, f_value, f_full, guess);
    const Temperature<T, U> temp {r.root};
    const Pressure<T, U> p {pressure(rho, mix.gas_constant, temp)};
    return {
      .state =
        {
          .rho = rho,
          .e = e,
          .p = p,
          .temp = temp,
        },
      .status = combine(combine(thermally_perfect_status(rho, e, temp, p), r.status), mix.status),
    };
}

// State from rho and T, retaining the supplied pressure; e comes from the curves.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE Closure<T, U>
  thermally_perfect_from_rho_temp(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                  const NormalisedComposition<T, U, N>& mix,
                                  Density<T, U> rho,
                                  Temperature<T, U> temp,
                                  Pressure<T, U> p)
{
    const SpecificEnergy<T, U> e {mixture_energy(gas, mix.fractions, temp.value()).e};
    return {
      .state =
        {
          .rho = rho,
          .e = e,
          .p = p,
          .temp = temp,
        },
      .status = combine(thermally_perfect_status(rho, e, temp, p), mix.status),
    };
}

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE Closure<T, U>
  thermally_perfect_from_rho_temp(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                  const NormalisedComposition<T, U, N>& mix,
                                  Density<T, U> rho,
                                  Temperature<T, U> temp)
{
    return thermally_perfect_from_rho_temp(gas,
                                           mix,
                                           rho,
                                           temp,
                                           Pressure<T, U> {pressure(rho, mix.gas_constant, temp)});
}

}  // namespace detail

// Per-cell composition from mass fractions.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V>
    requires SpeciesWeights<V, kind::mass_fraction, unitless, P> &&
             ScalarWith<typename range_quantity_t<V>::value_type, P>
CANARY_HOST_DEVICE auto mixture(const ThermallyPerfectMixture<P, U, N, S>& gas, const V& fractions)
{
    using T = range_quantity_t<V>::value_type;
    return detail::normalise(gas,
                             fractions,
                             T(typename ScalarTraits<T>::value_type(1)),
                             composition_status(fractions));
}

// Per-cell composition from partial densities and the density.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V, Scalar T>
    requires SpeciesWeights<V, kind::density, U, P> &&
             std::same_as<typename range_quantity_t<V>::value_type, T> && ScalarWith<T, P>
CANARY_HOST_DEVICE auto mixture(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                const V& partial_densities,
                                Density<T, U> rho)
{
    return detail::normalise(gas,
                             partial_densities,
                             rho.value(),
                             composition_status(partial_densities, rho));
}

// Closure from density and energy per unit mass, with the default initial
// guess (see detail::default_temperature_guess and
// detail::thermally_perfect_from_rho_e). The status combines the state
// status, the solve status (NotConverged, or Invalid when no bracket has a
// sign change) and the composition status.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE Closure<T, U> close(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                       const NormalisedComposition<T, U, N>& mix,
                                       Density<T, U> rho,
                                       SpecificEnergy<T, U> e)
{
    using V = ScalarTraits<T>::value_type;
    const V guess = detail::default_temperature_guess(gas,
                                                      detail::fraction_values(mix),
                                                      ScalarTraits<T>::value(e.value()));
    return detail::thermally_perfect_from_rho_e(gas, mix, rho, e, guess);
}

// As above, starting from the value part of the given temperature.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE Closure<T, U> close(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                       const NormalisedComposition<T, U, N>& mix,
                                       Density<T, U> rho,
                                       SpecificEnergy<T, U> e,
                                       const TemperatureGuess<T, U>& guess)
{
    return detail::thermally_perfect_from_rho_e(gas,
                                                mix,
                                                rho,
                                                e,
                                                ScalarTraits<T>::value(guess.temp.value()));
}

// T = p / (rho R), e from the curves.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE Closure<T, U> close(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                       const NormalisedComposition<T, U, N>& mix,
                                       Density<T, U> rho,
                                       Pressure<T, U> p)
{
    const Temperature<T, U> temp {temperature(p, rho, mix.gas_constant)};
    auto c = detail::thermally_perfect_from_rho_temp(gas, mix, rho, temp, p);
    return c;
}

// p = rho R T, e from the curves.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE Closure<T, U> close(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                       const NormalisedComposition<T, U, N>& mix,
                                       Density<T, U> rho,
                                       Temperature<T, U> temp)
{ return detail::thermally_perfect_from_rho_temp(gas, mix, rho, temp); }

// rho = p / (R T), e from the curves.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE Closure<T, U> close(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                       const NormalisedComposition<T, U, N>& mix,
                                       Pressure<T, U> p,
                                       Temperature<T, U> temp)
{
    const Density<T, U> rho {density(p, mix.gas_constant, temp)};
    auto c = detail::thermally_perfect_from_rho_temp(gas, mix, rho, temp, p);
    return c;
}

// cv = sum(y_i cv_i(T)), cp = cv + R, gamma = cp / cv, the frozen sound
// speed a = sqrt(gamma p / rho), and R.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T> && HasSqrt<T>
CANARY_HOST_DEVICE Properties<T, U> properties(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                               const NormalisedComposition<T, U, N>& mix,
                                               const ThermoState<T, U>& state)
{
    const SpecificHeat<T, U> cv {detail::mixture_energy(gas, mix.fractions, state.temp.value()).cv};
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

// Per-species properties per unit mass at one temperature, in model units:
// energy e_i, cv_i, enthalpy h_i, entropy s_i and Gibbs energy
// g_i = h_i - T s_i, with s_i and g_i at the standard pressure of the fit.
template <Scalar T, UnitSystem U, std::size_t N>
    requires(N > 0)
struct SpeciesProperties {
    int n;
    std::array<SpecificEnergy<T, U>, N> energy;
    std::array<SpecificHeat<T, U>, N> cv;
    std::array<SpecificEnthalpy<T, U>, N> enthalpy;
    std::array<SpecificEntropy<T, U>, N> entropy;
    std::array<SpecificGibbsEnergy<T, U>, N> gibbs_energy;
};

// Properties of every species at temperature temp, with one logarithm.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE SpeciesProperties<T, U, N>
  species_properties(const ThermallyPerfectMixture<P, U, N, S>& gas, Temperature<T, U> temp)
{
    const T t = temp.value();
    const T ln_t = math::log(t);
    SpeciesProperties<T, U, N> out {};
    out.n = gas.n;
    for (std::size_t i = 0; i < N; ++i) {
        const Nasa9CpH<T> c = nasa9_cp_h(gas.curve[i], t, ln_t);
        const T s_over_r = nasa9_s(gas.curve[i], t, ln_t);
        const P r = gas.gas_constant[i].value();
        const T h = T(r * c.h_over_r);
        const T s = T(r * s_over_r);
        out.energy[i] = SpecificEnergy<T, U> {T(r * (c.h_over_r - t))};
        out.cv[i] = SpecificHeat<T, U> {T(r * (c.cp_over_r - P(1)))};
        out.enthalpy[i] = SpecificEnthalpy<T, U> {h};
        out.entropy[i] = SpecificEntropy<T, U> {s};
        out.gibbs_energy[i] = SpecificGibbsEnergy<T, U> {h - (t * s)};
    }
    return out;
}

// Isobaric specific heat cp_i = R_i (cp/R)_i(T) of every species at
// temperature temp, in model units, with one logarithm.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires ScalarWith<T, P> && HasLog<T>
CANARY_HOST_DEVICE std::array<SpecificHeat<T, U>, N>
  species_isobaric_specific_heat(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                 Temperature<T, U> temp)
{
    const T t = temp.value();
    const T ln_t = math::log(t);
    std::array<SpecificHeat<T, U>, N> out {};
    for (std::size_t i = 0; i < N; ++i) {
        const T cp_over_r = nasa9_cp_h(gas.curve[i], t, ln_t).cp_over_r;
        out[i] = SpecificHeat<T, U> {T(gas.gas_constant[i].value() * cp_over_r)};
    }
    return out;
}

}  // namespace canary

#endif  // CANARY_MODELS_THERMALLY_PERFECT_HPP
