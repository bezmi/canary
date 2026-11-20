// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_MODELS_THERMALLY_PERFECT_HPP
#define CANARY_MODELS_THERMALLY_PERFECT_HPP

#include <algorithm>
#include <array>
#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/detail/closure.hpp>
#include <canary/detail/solve.hpp>
#include <canary/layout.hpp>
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
#include <limits>
#include <utility>

namespace canary
{

// Host input for one species of a thermally perfect mixture: its name (for
// error messages), molar mass [kg/mol], charge [elementary charges], metadata
// and fit. The name and borrowed metadata must outlive construction.
template <std::size_t S>
    requires(S > 0)
struct Nasa9Species {
    const char* name;
    double molar_mass;
    double charge;
    SpeciesMetadataView metadata;
    Nasa9Data<S> fit;
};

// Per-cell composition of a thermally perfect mixture: the mass fractions,
// each weight clamped to [0, total] and divided by the sum of the clamped
// weights, the mixture gas constant R = sum(y_i R_i), and the status of the
// composition they came from.
template <Scalar T, UnitSystem U, std::size_t N>
    requires(N > 0)
struct NormalisedComposition {
    std::array<MassFraction<T>, N> fractions;
    SpecificGasConstant<T, U> gas_constant;
    Status status;
    CompositionReason reason;
};

// Ideal gas mixture of N species at one temperature, with the heat capacity
// of each species a function of temperature given by a NASA 9-coefficient
// curve of up to S segments. All N species are used. The
// fields are in unit system U and are derived by
// make_thermally_perfect_mixture:
// - gas_constant: R_i = R_u / M_i;
// - curve: the fits folded into the temperature unit of U;
// - reference_temperature (298.15 K), and the species energy and cv there,
//   from which the closure from (rho, e) forms its default initial guess;
// - temperature_min and temperature_max: the full bracket of the
//   temperature solve; bracket_width: the width of the first bracket
//   around the guess;
// - relative_tolerance, absolute_tolerance and max_iterations: settings of
//   the solve; absolute_tolerance is in the model temperature units.
// State scalars with a single-precision value type are experimental: a
// closure that solves for a temperature can return Status::NotConverged or a
// diagnostic with ClosureReason::PrecisionLimit, so callers check the status of
// every closure.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S = 4>
    requires(N > 0 && S > 0)
struct ThermallyPerfectMixture {
    using parameter_type = P;
    using units = U;
    template <Scalar T> using State = ThermoState<T, U>;
    template <Scalar T> using Mixture = NormalisedComposition<T, U, N>;

    std::array<SpecificGasConstant<P, U>, N> gas_constant;
    std::array<Nasa9Curve<P, S>, N> curve;
    Temperature<P, U> reference_temperature;
    std::array<SpecificEnergy<P, U>, N> reference_energy;
    std::array<SpecificHeat<P, U>, N> reference_cv;
    Temperature<P, U> temperature_min;
    Temperature<P, U> temperature_max;
    Temperature<P, U> bracket_width;
    P relative_tolerance;
    std::size_t max_iterations;
    Temperature<P, U> absolute_tolerance {P(0)};

    static consteval bool padding_free()
    {
        return layout::members_fill<ThermallyPerfectMixture,
                                    decltype(gas_constant),
                                    decltype(curve),
                                    decltype(reference_temperature),
                                    decltype(reference_energy),
                                    decltype(reference_cv),
                                    decltype(temperature_min),
                                    decltype(temperature_max),
                                    decltype(bracket_width),
                                    decltype(relative_tolerance),
                                    decltype(max_iterations),
                                    decltype(absolute_tolerance)>();
    }
};

namespace detail
{

// Fills entry i of the model and of the species data from one species, after
// the model's reference temperature is set. Throws ParameterError.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
void fill_thermally_perfect_species(ThermallyPerfectMixture<P, U, N, S>& model,
                                    SpeciesData<P, U, N>& species,
                                    SpeciesIndex index,
                                    const Nasa9Species<S>& input,
                                    const UnitScales<P, U>& scales)
{
    const auto i = index.value();
    if (!std::isfinite(input.molar_mass) || !(input.molar_mass > 0.0)) {
        throw ParameterError(
          std::format("molar mass must be finite and positive, got {}", input.molar_mass));
    }
    if (!std::isfinite(input.charge)) {
        throw ParameterError(std::format("charge must be finite, got {}", input.charge));
    }
    check_species_metadata(input.metadata, input.charge);
    const auto molar_mass = static_cast<long double>(input.molar_mass);
    const SpecificGasConstant<long double, si> gas_constant {
      constants::universal_gas_constant<long double> / molar_mass};
    const Charge<long double, si> charge {static_cast<long double>(input.charge) *
                                          constants::elementary_charge<long double>};
    if (input.charge != 0.0 && charge.value() == 0.0L) {
        throw ParameterError("nonzero charge rounds to zero in the model's units");
    }
    model.gas_constant[i] =
      positive_thermodynamic_parameter<P>(convert<U>(gas_constant, scales), "stored gas constant");
    species.mass[i] = particle_mass_parameter(MolarMass<double> {input.molar_mass}, scales);
    species.charge[i] = species_charge_parameter(charge, scales);

    const P t_ref = model.reference_temperature.value();
    const Nasa9CpH<P> ref = nasa9_cp_h(model.curve[i], t_ref, std::log(t_ref));
    const P r = model.gas_constant[i].value();
    static_cast<void>(positive_thermodynamic_parameter<P>(SpecificHeat<P, U> {r * ref.cp_over_r},
                                                          "reference isobaric specific heat"));
    model.reference_energy[i] =
      thermodynamic_parameter<P>(SpecificEnergy<P, U> {r * (ref.h_over_r - t_ref)},
                                 "reference specific energy");
    model.reference_cv[i] =
      positive_thermodynamic_parameter<P>(SpecificHeat<P, U> {r * (ref.cp_over_r - P(1))},
                                          "reference isochoric specific heat");
}

}  // namespace detail

// Mixture model of the given species in the unit system of the scales, with
// its species data (particle mass M / N_A and charge q e, in U). The
// temperature limits are 10 K and 100000 K, the first bracket is 1000 K
// wide, the relative tolerance is 1e-10, the absolute temperature tolerance
// is zero and the iteration limit 30. Throws
// ParameterError, with the species name in the message, when a molar mass
// is not finite and positive, a charge is not finite, or a fit is rejected
// by make_nasa9_curve. Metadata requires a specified, known role, nonempty
// distinct element identifiers and positive element counts. A free electron
// requires charge -1 and empty composition when supplied. At most one species
// may have the electron role. Names must be non-null and nonempty; a name
// error gives the species index. Stored masses, gas constants and reference
// cp/cv must be finite and positive; reference energies must be finite. Nonzero
// charges must remain nonzero and finite in the model's units. Default
// temperatures and bracket width must remain finite and positive after conversion.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S>
    requires(N > 0 && S > 0)
Built<ThermallyPerfectMixture<P, U, N, S>, SpeciesData<P, U, N>>
  make_thermally_perfect_mixture(const std::array<Nasa9Species<S>, N>& species,
                                 const UnitScales<P, U>& scales)
{
    const auto reference_temperature = detail::positive_thermodynamic_parameter<P>(
      convert<U>(Temperature<long double, si> {298.15L}, scales),
      "reference temperature");
    const auto temperature_min = detail::positive_thermodynamic_parameter<P>(
      convert<U>(Temperature<long double, si> {10}, scales),
      "minimum temperature");
    const auto temperature_max = detail::positive_thermodynamic_parameter<P>(
      convert<U>(Temperature<long double, si> {100000}, scales),
      "maximum temperature");
    const auto bracket_width = detail::positive_thermodynamic_parameter<P>(
      convert<U>(Temperature<long double, si> {1000}, scales),
      "temperature bracket width");
    // Each array entry is constructed from a complete validated curve.
    const auto make_curve = [&](std::size_t i) {
        if (species[i].name == nullptr || species[i].name[0] == '\0') {
            throw ParameterError(std::format("species {}: name must not be null or empty", i));
        }
        try {
            return make_nasa9_curve(species[i].fit,
                                    static_cast<P>(scale_factor<kind::temperature>(scales)));
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("species {}: {}", species[i].name, error.what()));
        }
    };
    const auto curves = [&]<std::size_t... I>(std::index_sequence<I...>) {
        return std::array<Nasa9Curve<P, S>, N> {make_curve(I)...};
    }(std::make_index_sequence<N> {});
    Built<ThermallyPerfectMixture<P, U, N, S>, SpeciesData<P, U, N>> built {
      .model = {.gas_constant = {},
                .curve = curves,
                .reference_temperature = reference_temperature,
                .reference_energy = {},
                .reference_cv = {},
                .temperature_min = temperature_min,
                .temperature_max = temperature_max,
                .bracket_width = bracket_width,
                .relative_tolerance = P(1.0e-10),
                .max_iterations = 30},
      .species = {}};
    auto& model = built.model;
    built.species.n = N;
    bool has_electron = false;
    for (const SpeciesIndex i : species_indices(0, N)) {
        if (species[i.value()].name == nullptr || species[i.value()].name[0] == '\0') {
            throw ParameterError(
              std::format("species {}: name must not be null or empty", i.value()));
        }
        try {
            if (species[i.value()].metadata.role == SpeciesRole::electron) {
                if (has_electron) { throw ParameterError("more than one electron-role species"); }
                has_electron = true;
            }
            detail::fill_thermally_perfect_species(model,
                                                   built.species,
                                                   i,
                                                   species[i.value()],
                                                   scales);
        } catch (const ParameterError& error) {
            throw ParameterError(
              std::format("species {}: {}", species[i.value()].name, error.what()));
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
  const ThermallyPerfectMixture<P, U, N, S>& gas, const V& weights, T total, CompositionCheck check)
{
    // NOLINTNEXTLINE(misc-static-assert): run-time extent
    CANARY_ASSERT(static_cast<std::size_t>(weights.extent(0)) == N);
    std::array<T, N> w {};
    T sum = T(typename ScalarTraits<T>::value_type(0));
    for (const SpeciesIndex i : species_indices(0, N)) {
        w[i.value()] = clamped_weight(weights, i, total);
        sum = sum + w[i.value()];
    }
    NormalisedComposition<T, U, N> mix {};
    mix.status = check.status;
    mix.reason = check.reason;
    for (const SpeciesIndex i : species_indices(0, N)) {
        mix.fractions[i.value()] = MassFraction<T> {w[i.value()] / sum};
    }
    auto gas_constant = gas.gas_constant[0] * mix.fractions[0].value();
    for (const SpeciesIndex i : species_indices(1, N)) {
        gas_constant =
          gas_constant + (gas.gas_constant[i.value()] * mix.fractions[i.value()].value());
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
    for (const SpeciesIndex i : species_indices(0, N)) {
        const Nasa9CpH<W> c = nasa9_cp_h(gas.curve[i.value()], t, ln_t);
        const P r = gas.gas_constant[i.value()].value();
        e = e + W(y[i.value()].value() * (r * (c.h_over_r - t)));
        cv = cv + W(y[i.value()].value() * (r * (c.cp_over_r - P(1))));
    }
    return {.e = e, .cv = cv};
}

// Valid when density and temperature have finite, positive value parts and
// the energy a finite value part (it may be negative); Invalid otherwise.
template <Scalar T, UnitSystem U>
CANARY_HOST_DEVICE ClosureCheck thermally_perfect_check(Density<T, U> rho,
                                                        SpecificEnergy<T, U> e,
                                                        Temperature<T, U> temp)
{
    auto check = positive_check(rho.value(),
                                ClosureReason::NonfiniteDensity,
                                ClosureReason::NonpositiveDensity);
    check = combine(check, finite_check(e.value(), ClosureReason::NonfiniteEnergy));
    return combine(check,
                   positive_check(temp.value(),
                                  ClosureReason::NonfiniteTemperature,
                                  ClosureReason::NonpositiveTemperature));
}

// State validation includes the actual returned pressure. The three-field
// overload also serves the independently evaluated energy modes.
template <Scalar T, UnitSystem U>
CANARY_HOST_DEVICE ClosureCheck thermally_perfect_check(Density<T, U> rho,
                                                        SpecificEnergy<T, U> e,
                                                        Temperature<T, U> temp,
                                                        Pressure<T, U> p)
{
    return combine(thermally_perfect_check(rho, e, temp),
                   positive_check(p.value(),
                                  ClosureReason::NonfinitePressure,
                                  ClosureReason::NonpositivePressure));
}

template <Scalar T, UnitSystem U, std::size_t N>
CANARY_HOST_DEVICE std::array<MassFraction<typename ScalarTraits<T>::value_type>, N>
  fraction_values(const NormalisedComposition<T, U, N>& mix)
{
    using V = ScalarTraits<T>::value_type;
    std::array<MassFraction<V>, N> y {};
    for (const SpeciesIndex i : species_indices(0, N)) {
        y[i.value()] = MassFraction<V> {ScalarTraits<T>::value(mix.fractions[i.value()].value())};
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
    for (const SpeciesIndex i : species_indices(0, N)) {
        e_ref += y[i.value()].value() * static_cast<V>(gas.reference_energy[i.value()].value());
        cv_ref += y[i.value()].value() * static_cast<V>(gas.reference_cv[i.value()].value());
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
      .absolute_tolerance = static_cast<V>(gas.absolute_tolerance.value()),
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
    const auto check = combine(combine(composition_check(mix), solve_check(r)),
                               thermally_perfect_check(rho, e, temp, p));
    return {
      .state =
        {
          .rho = rho,
          .e = e,
          .p = p,
          .temp = temp,
        },
      .status = check.status,
      .diagnostic = check.diagnostic,
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
    const auto check = combine(composition_check(mix), thermally_perfect_check(rho, e, temp, p));
    return {
      .state =
        {
          .rho = rho,
          .e = e,
          .p = p,
          .temp = temp,
        },
      .status = check.status,
      .diagnostic = check.diagnostic,
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
                             check_composition(fractions));
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
                             check_composition(partial_densities, rho));
}

// Per-cell composition from partial densities and a previously computed composition check.
// check must be from check_composition(partial_densities, rho) for the same,
// unchanged inputs, with status other than Invalid. Normalisation divides by
// the sum of clamped weights.
template <std::floating_point P, UnitSystem U, std::size_t N, std::size_t S, typename V, Scalar T>
    requires SpeciesWeights<V, kind::density, U, P> &&
             std::same_as<typename range_quantity_t<V>::value_type, T> && ScalarWith<T, P>
CANARY_HOST_DEVICE auto mixture(const ThermallyPerfectMixture<P, U, N, S>& gas,
                                const V& partial_densities,
                                Density<T, U> rho,
                                CompositionCheck check) noexcept
{
    return detail::normalise(gas, partial_densities, rho.value(), check);
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
CANARY_HOST_DEVICE SingleTemperatureProperties<T, U>
  properties(const ThermallyPerfectMixture<P, U, N, S>& gas,
             const NormalisedComposition<T, U, N>& mix,
             const ThermoState<T, U>& state) noexcept
{
    const SpecificHeat<T, U> cv {detail::mixture_energy(gas, mix.fractions, state.temp.value()).cv};
    const SpecificHeat<T, U> cp {cv.value() + mix.gas_constant.value()};
    const HeatCapacityRatio<T> gamma {heat_capacity_ratio(cp, cv)};
    return {
      .cv = cv,
      .cp = cp,
      .capacity_ratio = gamma,
      .gas_constant = mix.gas_constant,
      .acoustics = {.sound_speed =
                      SoundSpeed<T, U> {
                        sound_speed(AcousticExponent<T> {gamma.value()}, state.p, state.rho)},
                    .acoustic_exponent = AcousticExponent<T> {gamma.value()}},
    };
}

// Per-species properties per unit mass at one temperature, in model units:
// energy e_i, cv_i, enthalpy h_i, entropy s_i and Gibbs energy
// g_i = h_i - T s_i, with s_i and g_i at the standard pressure of the fit.
template <Scalar T, UnitSystem U, std::size_t N>
    requires(N > 0)
struct SpeciesProperties {
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
    for (const SpeciesIndex i : species_indices(0, N)) {
        const Nasa9CpH<T> c = nasa9_cp_h(gas.curve[i.value()], t, ln_t);
        const T s_over_r = nasa9_s(gas.curve[i.value()], t, ln_t);
        const P r = gas.gas_constant[i.value()].value();
        const T h = T(r * c.h_over_r);
        const T s = T(r * s_over_r);
        out.energy[i.value()] = SpecificEnergy<T, U> {T(r * (c.h_over_r - t))};
        out.cv[i.value()] = SpecificHeat<T, U> {T(r * (c.cp_over_r - P(1)))};
        out.enthalpy[i.value()] = SpecificEnthalpy<T, U> {h};
        out.entropy[i.value()] = SpecificEntropy<T, U> {s};
        out.gibbs_energy[i.value()] = SpecificGibbsEnergy<T, U> {h - (t * s)};
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
    for (const SpeciesIndex i : species_indices(0, N)) {
        const T cp_over_r = nasa9_cp_h(gas.curve[i.value()], t, ln_t).cp_over_r;
        out[i.value()] = SpecificHeat<T, U> {T(gas.gas_constant[i.value()].value() * cp_over_r)};
    }
    return out;
}

namespace detail
{

template <std::size_t S> struct InputStorage<Nasa9Species<S>> {
    OwnedSpeciesName name {};
    double molar_mass {};
    double charge {};
    InputStorage<SpeciesMetadataView> metadata {};
    Nasa9Data<S> fit {};
    using Scratch = std::tuple<typename InputStorage<SpeciesMetadataView>::Scratch>;
    InputStorage() = default;
    explicit InputStorage(const Nasa9Species<S>& input) :
        name(input.name), molar_mass(input.molar_mass), charge(input.charge),
        metadata(input.metadata), fit(input.fit)
    {
    }
    Nasa9Species<S> borrow([[maybe_unused]] Scratch& scratch) const
    {
        return {
          .name = name.c_str(),
          .molar_mass = molar_mass,
          .charge = charge,
          .metadata = metadata.borrow(std::get<0>(scratch)),
          .fit = fit,
        };
    }
};

}  // namespace detail

// Builds from owned and direct inputs, keeping borrowed views inside this call.
template <typename... Explicit, typename... Inputs>
    requires detail::HasOwnedInput<Inputs...>
auto make_thermally_perfect_mixture(const Inputs&... inputs)
{
    return detail::with_inputs(
      [](const auto&... borrowed) {
          return make_thermally_perfect_mixture<Explicit...>(borrowed...);
      },
      inputs...);
}

}  // namespace canary

#endif  // CANARY_MODELS_THERMALLY_PERFECT_HPP
