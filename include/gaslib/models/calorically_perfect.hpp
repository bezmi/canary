#ifndef GASLIB_MODELS_CALORICALLY_PERFECT_HPP
#define GASLIB_MODELS_CALORICALLY_PERFECT_HPP

#include <cmath>
#include <concepts>
#include <format>
#include <gaslib/build.hpp>
#include <gaslib/config.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/scalar.hpp>
#include <gaslib/species.hpp>
#include <gaslib/state.hpp>
#include <gaslib/status.hpp>

namespace gaslib
{

// Ideal gas of one species with constant heat capacities. The fields are
// derived once by make_calorically_perfect, in unit system U.
template <std::floating_point P, UnitSystem U> struct CaloricallyPerfect {
    using parameter_type = P;
    using units = U;
    template <Scalar T> using State = ThermoState<T, U>;

    SpecificGasConstant<P, U> gas_constant;
    SpecificHeat<P, U> cv;
    SpecificHeat<P, U> cp;
    HeatCapacityRatio<P> gamma;
};

namespace detail
{

template <std::floating_point P> void require_positive(const char* name, P value)
{
    if (!std::isfinite(value) || !(value > P(0))) {
        throw ParameterError(std::format("{} must be finite and positive, got {}", name, value));
    }
}

template <std::floating_point P, UnitSystem U>
CaloricallyPerfect<P, U> calorically_perfect_from(SpecificGasConstant<P, U> gas_constant,
                                                  HeatCapacityRatio<P> gamma)
{
    require_positive("gas constant", gas_constant.value());
    if (!std::isfinite(gamma.value()) || !(gamma.value() > P(1))) {
        throw ParameterError(
          std::format("heat capacity ratio must be finite and greater than 1, got {}",
                      gamma.value()));
    }
    const SpecificHeat<P, U> cv = isochoric_specific_heat(gas_constant, gamma);
    return {
      .gas_constant = gas_constant,
      .cv = cv,
      .cp = isobaric_specific_heat(cv, gamma),
      .gamma = gamma,
    };
}

// Valid when density, energy per unit volume and temperature have finite,
// positive value parts; Invalid otherwise.
template <Scalar T, UnitSystem U>
GASLIB_HOST_DEVICE constexpr Status
  state_status(Density<T, U> rho, EnergyDensity<T, U> rho_e, Temperature<T, U> temp)
{
    using V = ScalarTraits<T>::value_type;
    const auto positive = [](T x) {
        const V v = ScalarTraits<T>::value(x);
        return std::isfinite(v) && v > V(0);
    };
    return positive(rho.value()) && positive(rho_e.value()) && positive(temp.value())
             ? Status::Valid
             : Status::Invalid;
}

}  // namespace detail

// Model from particle mass, charge, Boltzmann constant (exactly 1 in scaled
// units, see boltzmann_constant) and heat capacity ratio, with its species
// data. Throws ParameterError when the mass or the Boltzmann constant is not
// positive, the charge is not finite, or gamma is not greater than 1.
template <std::floating_point P, UnitSystem U>
Built<CaloricallyPerfect<P, U>, SpeciesData<P, U, 1>>
  make_calorically_perfect(ParticleMass<P, U> mass,
                           Charge<P, U> charge,
                           BoltzmannConstant<P, U> boltzmann,
                           HeatCapacityRatio<P> gamma)
{
    detail::require_positive("particle mass", mass.value());
    detail::require_positive("Boltzmann constant", boltzmann.value());
    if (!std::isfinite(charge.value())) {
        throw ParameterError(std::format("charge must be finite, got {}", charge.value()));
    }
    return {
      .model = detail::calorically_perfect_from(specific_gas_constant(boltzmann, mass), gamma),
      .species = {.n = 1, .mass = {mass}, .charge = {charge}},
    };
}

// Model from the specific gas constant and the heat capacity ratio, for
// codes that do not need particle masses. Throws ParameterError as above.
template <std::floating_point P, UnitSystem U>
CaloricallyPerfect<P, U> make_calorically_perfect(SpecificGasConstant<P, U> gas_constant,
                                                  HeatCapacityRatio<P> gamma)
{ return detail::calorically_perfect_from(gas_constant, gamma); }

// Closure from density and energy per unit volume:
// p = (gamma - 1) rho_e, T = p / (rho R), e = rho_e / rho.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE Closure<T, U>
  close(const CaloricallyPerfect<P, U>& gas, Density<T, U> rho, EnergyDensity<T, U> rho_e)
{
    const Pressure<T, U> p {pressure(gas.gamma, rho_e)};
    const Temperature<T, U> temp {temperature(p, rho, gas.gas_constant)};
    return {
      .state = {.rho = rho, .e = specific_energy(rho_e, rho), .p = p, .temp = temp},
      .status = detail::state_status(rho, rho_e, temp),
    };
}

// Closure from density and energy per unit mass, through rho_e = rho e.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE Closure<T, U>
  close(const CaloricallyPerfect<P, U>& gas, Density<T, U> rho, SpecificEnergy<T, U> e)
{
    auto c = close(gas, rho, energy_density(rho, e));
    c.state.e = e;
    return c;
}

// Closure from density and pressure: rho_e = p / (gamma - 1), T = p / (rho R).
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE Closure<T, U>
  close(const CaloricallyPerfect<P, U>& gas, Density<T, U> rho, Pressure<T, U> p)
{
    const EnergyDensity<T, U> rho_e {energy_density(gas.gamma, p)};
    const Temperature<T, U> temp {temperature(p, rho, gas.gas_constant)};
    return {
      .state = {.rho = rho, .e = specific_energy(rho_e, rho), .p = p, .temp = temp},
      .status = detail::state_status(rho, rho_e, temp),
    };
}

// Closure from density and temperature: p = rho R T, e = cv T.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE Closure<T, U>
  close(const CaloricallyPerfect<P, U>& gas, Density<T, U> rho, Temperature<T, U> temp)
{
    const SpecificEnergy<T, U> e {specific_energy(gas.cv, temp)};
    return {
      .state =
        {
          .rho = rho,
          .e = e,
          .p = Pressure<T, U> {pressure(rho, gas.gas_constant, temp)},
          .temp = temp,
        },
      .status = detail::state_status(rho, energy_density(rho, e), temp),
    };
}

// Closure from pressure and temperature: rho = p / (R T), e = cv T.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE Closure<T, U>
  close(const CaloricallyPerfect<P, U>& gas, Pressure<T, U> p, Temperature<T, U> temp)
{
    const Density<T, U> rho {density(p, gas.gas_constant, temp)};
    const SpecificEnergy<T, U> e {specific_energy(gas.cv, temp)};
    return {
      .state = {.rho = rho, .e = e, .p = p, .temp = temp},
      .status = detail::state_status(rho, energy_density(rho, e), temp),
    };
}

// Heat capacities, gamma and gas constant of the model, and the sound speed
// a = sqrt(gamma p / rho) at the state.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P> && HasSqrt<T>
GASLIB_HOST_DEVICE Properties<T, U> properties(const CaloricallyPerfect<P, U>& gas,
                                               const ThermoState<T, U>& state)
{
    return {
      .cv = SpecificHeat<T, U> {gas.cv},
      .cp = SpecificHeat<T, U> {gas.cp},
      .gamma = HeatCapacityRatio<T> {gas.gamma},
      .sound_speed = SoundSpeed<T, U> {sound_speed(gas.gamma, state.p, state.rho)},
      .gas_constant = SpecificGasConstant<T, U> {gas.gas_constant},
    };
}

}  // namespace gaslib

#endif  // GASLIB_MODELS_CALORICALLY_PERFECT_HPP
