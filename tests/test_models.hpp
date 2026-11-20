// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_TESTS_TEST_MODELS_HPP
#define CANARY_TESTS_TEST_MODELS_HPP

// Minimal gas models for tests of the closure front end and the model
// concepts. Each is an ideal gas with constant cv, in SI units, with values
// of type double for parameters.

#include <concepts>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/status.hpp>

namespace canary::test
{

// Parameters shared by the test models.
struct IdealParameters {
    SpecificGasConstant<double, si> gas_constant;
    SpecificHeat<double, si> cv;
};

template <Scalar T>
Closure<T, si>
  ideal_from_rho_e(const IdealParameters& g, Density<T, si> rho, SpecificEnergy<T, si> e)
{
    const Temperature<T, si> temp {temperature(e, g.cv)};
    return {
      .state =
        {
          .rho = rho,
          .e = e,
          .p = Pressure<T, si> {pressure(rho, g.gas_constant, temp)},
          .temp = temp,
        },
      .status = Status::Valid,
      .diagnostic = {},
    };
}

template <Scalar T>
Closure<T, si>
  ideal_from_rho_temp(const IdealParameters& g, Density<T, si> rho, Temperature<T, si> temp)
{
    return {
      .state =
        {
          .rho = rho,
          .e = SpecificEnergy<T, si> {specific_energy(g.cv, temp)},
          .p = Pressure<T, si> {pressure(rho, g.gas_constant, temp)},
          .temp = temp,
        },
      .status = Status::Valid,
      .diagnostic = {},
    };
}

// A caller-owned record exercises the structural property contract.
template <Scalar T> struct IdealProperties {
    SpecificHeat<T, si> cv;
    SpecificHeat<T, si> cp;
    HeatCapacityRatio<T> capacity_ratio;
    SpecificGasConstant<T, si> gas_constant;
    AcousticProperties<T, si> acoustics;
};

template <Scalar T>
IdealProperties<T> ideal_properties(const IdealParameters& g, const ThermoState<T, si>& s)
{
    const double cp = g.cv.value() + g.gas_constant.value();
    const HeatCapacityRatio<double> gamma {cp / g.cv.value()};
    const auto a = sound_speed(AcousticExponent<double> {gamma.value()}, s.p, s.rho);
    return {
      .cv = SpecificHeat<T, si> {T(g.cv.value())},
      .cp = SpecificHeat<T, si> {T(cp)},
      .capacity_ratio = HeatCapacityRatio<T> {T(gamma.value())},
      .gas_constant = SpecificGasConstant<T, si> {T(g.gas_constant.value())},
      .acoustics = {.sound_speed = SoundSpeed<T, si> {a},
                    .acoustic_exponent = AcousticExponent<T> {T(gamma.value())}},
    };
}

// A complete single-species model.
struct IdealGas {
    using parameter_type = double;
    using units = si;
    template <Scalar T> using State = ThermoState<T, si>;

    IdealParameters parameters;
};

template <Scalar T>
Closure<T, si> close(const IdealGas& g, Density<T, si> rho, SpecificEnergy<T, si> e)
{ return ideal_from_rho_e(g.parameters, rho, e); }

template <Scalar T> Closure<T, si> close(const IdealGas& g, Density<T, si> rho, Pressure<T, si> p)
{
    return ideal_from_rho_temp(g.parameters,
                               rho,
                               Temperature<T, si> {temperature(p, rho, g.parameters.gas_constant)});
}

template <Scalar T>
Closure<T, si> close(const IdealGas& g, Density<T, si> rho, Temperature<T, si> temp)
{ return ideal_from_rho_temp(g.parameters, rho, temp); }

template <Scalar T>
Closure<T, si> close(const IdealGas& g, Pressure<T, si> p, Temperature<T, si> temp)
{
    return ideal_from_rho_temp(g.parameters,
                               Density<T, si> {density(p, g.parameters.gas_constant, temp)},
                               temp);
}

template <Scalar T> IdealProperties<T> properties(const IdealGas& g, const ThermoState<T, si>& s)
{ return ideal_properties(g.parameters, s); }

// An initial temperature for an iterative closure.
struct TemperatureGuess {
    Temperature<double, si> temp;
};

// A complete model that also has a (rho, e) closure taking a guess. The
// closure here ignores the guess and marks its result NotConverged, so tests
// can tell that the overload was selected.
struct GuessGas {
    using parameter_type = double;
    using units = si;
    template <Scalar T> using State = ThermoState<T, si>;

    IdealParameters parameters;
};

template <Scalar T>
Closure<T, si> close(const GuessGas& g, Density<T, si> rho, SpecificEnergy<T, si> e)
{ return ideal_from_rho_e(g.parameters, rho, e); }

template <Scalar T>
Closure<T, si> close(const GuessGas& g,
                     Density<T, si> rho,
                     SpecificEnergy<T, si> e,
                     const TemperatureGuess& /*guess*/)
{
    auto c = ideal_from_rho_e(g.parameters, rho, e);
    c.status = Status::NotConverged;
    c.diagnostic = {ClosureReason::IterationLimit, {}};
    return c;
}

template <Scalar T> Closure<T, si> close(const GuessGas& g, Density<T, si> rho, Pressure<T, si> p)
{
    return ideal_from_rho_temp(g.parameters,
                               rho,
                               Temperature<T, si> {temperature(p, rho, g.parameters.gas_constant)});
}

template <Scalar T>
Closure<T, si> close(const GuessGas& g, Density<T, si> rho, Temperature<T, si> temp)
{ return ideal_from_rho_temp(g.parameters, rho, temp); }

template <Scalar T>
Closure<T, si> close(const GuessGas& g, Pressure<T, si> p, Temperature<T, si> temp)
{
    return ideal_from_rho_temp(g.parameters,
                               Density<T, si> {density(p, g.parameters.gas_constant, temp)},
                               temp);
}

template <Scalar T> IdealProperties<T> properties(const GuessGas& g, const ThermoState<T, si>& s)
{ return ideal_properties(g.parameters, s); }

// A model with its own (rho, rho e) closure, which marks its result
// CompositionViolation so tests can tell that it was selected over the
// generic one.
struct OwnEnergyDensityGas {
    using parameter_type = double;
    using units = si;
    template <Scalar T> using State = ThermoState<T, si>;

    IdealParameters parameters;
};

template <Scalar T>
Closure<T, si> close(const OwnEnergyDensityGas& g, Density<T, si> rho, SpecificEnergy<T, si> e)
{ return ideal_from_rho_e(g.parameters, rho, e); }

template <Scalar T>
Closure<T, si> close(const OwnEnergyDensityGas& g, Density<T, si> rho, EnergyDensity<T, si> rho_e)
{
    auto c = ideal_from_rho_e(g.parameters, rho, specific_energy(rho_e, rho));
    c.status = Status::CompositionViolation;
    c.diagnostic = {ClosureReason::CompositionSumMismatch, {}};
    return c;
}

// A model with closures from (rho, e), (rho, p) and (rho, T) only.
struct NoPressureTemperatureGas {
    using parameter_type = double;
    using units = si;
    template <Scalar T> using State = ThermoState<T, si>;

    IdealParameters parameters;
};

template <Scalar T>
Closure<T, si> close(const NoPressureTemperatureGas& g, Density<T, si> rho, SpecificEnergy<T, si> e)
{ return ideal_from_rho_e(g.parameters, rho, e); }

template <Scalar T>
Closure<T, si> close(const NoPressureTemperatureGas& g, Density<T, si> rho, Pressure<T, si> p)
{
    return ideal_from_rho_temp(g.parameters,
                               rho,
                               Temperature<T, si> {temperature(p, rho, g.parameters.gas_constant)});
}

template <Scalar T>
Closure<T, si> close(const NoPressureTemperatureGas& g, Density<T, si> rho, Temperature<T, si> temp)
{ return ideal_from_rho_temp(g.parameters, rho, temp); }

template <Scalar T>
IdealProperties<T> properties(const NoPressureTemperatureGas& g, const ThermoState<T, si>& s)
{ return ideal_properties(g.parameters, s); }

// A model whose (rho, e) closure requires a guess.
struct GuessRequiredGas {
    using parameter_type = double;
    using units = si;
    template <Scalar T> using State = ThermoState<T, si>;

    IdealParameters parameters;
};

template <Scalar T>
Closure<T, si> close(const GuessRequiredGas& g,
                     Density<T, si> rho,
                     SpecificEnergy<T, si> e,
                     const TemperatureGuess& /*guess*/)
{ return ideal_from_rho_e(g.parameters, rho, e); }

// A two-species mixture whose per-cell mixture value is the mass-weighted
// gas constant and cv.
struct TwoSpeciesGas {
    using parameter_type = double;
    using units = si;
    template <Scalar T> using State = ThermoState<T, si>;
    template <Scalar T> struct Mixture {
        SpecificGasConstant<T, si> gas_constant;
        SpecificHeat<T, si> cv;
    };

    IdealParameters species[2];
};

template <typename V>
    requires SpeciesRangeOf<V, kind::mass_fraction>
auto mixture(const TwoSpeciesGas& g, const V& y)
{
    using T = range_quantity_t<V>::value_type;
    const T y0 = y[0].value();
    const T y1 = y[1].value();
    const T r = y0 * g.species[0].gas_constant.value() + y1 * g.species[1].gas_constant.value();
    const T cv = y0 * g.species[0].cv.value() + y1 * g.species[1].cv.value();
    return TwoSpeciesGas::Mixture<T> {SpecificGasConstant<T, si> {r}, SpecificHeat<T, si> {cv}};
}

template <Scalar T>
IdealParameters as_parameters(const TwoSpeciesGas::Mixture<T>& m)
    requires std::same_as<T, double>
{ return {m.gas_constant, m.cv}; }

template <Scalar T>
Closure<T, si> close(const TwoSpeciesGas& /*g*/,
                     const TwoSpeciesGas::Mixture<T>& m,
                     Density<T, si> rho,
                     SpecificEnergy<T, si> e)
{ return ideal_from_rho_e(as_parameters(m), rho, e); }

template <Scalar T>
Closure<T, si> close(const TwoSpeciesGas& /*g*/,
                     const TwoSpeciesGas::Mixture<T>& m,
                     Density<T, si> rho,
                     Pressure<T, si> p)
{ return ideal_from_rho_temp(as_parameters(m), rho, temperature(p, rho, m.gas_constant)); }

template <Scalar T>
Closure<T, si> close(const TwoSpeciesGas& /*g*/,
                     const TwoSpeciesGas::Mixture<T>& m,
                     Density<T, si> rho,
                     Temperature<T, si> temp)
{ return ideal_from_rho_temp(as_parameters(m), rho, temp); }

template <Scalar T>
Closure<T, si> close(const TwoSpeciesGas& /*g*/,
                     const TwoSpeciesGas::Mixture<T>& m,
                     Pressure<T, si> p,
                     Temperature<T, si> temp)
{ return ideal_from_rho_temp(as_parameters(m), density(p, m.gas_constant, temp), temp); }

template <Scalar T>
IdealProperties<T> properties(const TwoSpeciesGas& /*g*/,
                              const TwoSpeciesGas::Mixture<T>& m,
                              const ThermoState<T, si>& s)
{ return ideal_properties(as_parameters(m), s); }

}  // namespace canary::test

#endif  // CANARY_TESTS_TEST_MODELS_HPP
