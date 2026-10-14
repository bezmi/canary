#ifndef CANARY_MODELS_CALORICALLY_PERFECT_HPP
#define CANARY_MODELS_CALORICALLY_PERFECT_HPP

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>
#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/quantity.hpp>
#include <canary/scalar.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/status.hpp>
#include <ranges>

namespace canary
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
    const SpecificHeat<P, U> cp = isobaric_specific_heat(cv, gamma);
    require_positive("isochoric specific heat", cv.value());
    require_positive("isobaric specific heat", cp.value());
    return {
      .gas_constant = gas_constant,
      .cv = cv,
      .cp = cp,
      .gamma = gamma,
    };
}

}  // namespace detail

// Model from particle mass, charge, Boltzmann constant (exactly 1 in scaled
// units, see boltzmann_constant) and heat capacity ratio, with its species
// data. Throws ParameterError when the mass or the Boltzmann constant is not
// positive, the charge is not finite, gamma is not greater than 1, or either
// derived specific heat is not finite and positive.
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

namespace detail
{

// Coefficients a calorically perfect closure reads: those of the
// one-species model (value type P) or of a per-cell mixture (value type T).
template <Scalar C, UnitSystem U> struct PerfectCoefficients {
    HeatCapacityRatio<C> gamma;
    SpecificGasConstant<C, U> gas_constant;
    SpecificHeat<C, U> cv;
    SpecificHeat<C, U> cp;
};

// Every state value, energy density and coefficient must have a finite,
// positive value part; gamma must exceed one. Derivatives do not set validity.
template <Scalar T, Scalar C, UnitSystem U>
CANARY_HOST_DEVICE constexpr Status state_status(const PerfectCoefficients<C, U>& k,
                                                 const ThermoState<T, U>& state,
                                                 EnergyDensity<T, U> rho_e)
{
    const auto positive = []<Scalar S>(S x) {
        using V = ScalarTraits<S>::value_type;
        const V v = ScalarTraits<S>::value(x);
        return std::isfinite(v) && v > V(0);
    };
    return positive(state.rho.value()) && positive(state.e.value()) && positive(state.p.value()) &&
               positive(state.temp.value()) && positive(rho_e.value()) &&
               positive(k.gas_constant.value()) && positive(k.cv.value()) &&
               positive(k.cp.value()) && positive(k.gamma.value() - C(1))
             ? Status::Valid
             : Status::Invalid;
}

// p = (gamma - 1) rho_e, T = p / (rho R), e = rho_e / rho.
template <Scalar T, Scalar C, UnitSystem U>
    requires ScalarPair<T, C>
CANARY_HOST_DEVICE Closure<T, U> perfect_from_rho_rho_e(const PerfectCoefficients<C, U>& k,
                                                        Density<T, U> rho,
                                                        EnergyDensity<T, U> rho_e)
{
    const Pressure<T, U> p {pressure(k.gamma, rho_e)};
    const Temperature<T, U> temp {temperature(p, rho, k.gas_constant)};
    const ThermoState<T, U> state {.rho = rho,
                                   .e = specific_energy(rho_e, rho),
                                   .p = p,
                                   .temp = temp};
    return {.state = state, .status = state_status(k, state, rho_e)};
}

// rho_e = rho e, then as above; the given e is kept.
template <Scalar T, Scalar C, UnitSystem U>
    requires ScalarPair<T, C>
CANARY_HOST_DEVICE Closure<T, U>
  perfect_from_rho_e(const PerfectCoefficients<C, U>& k, Density<T, U> rho, SpecificEnergy<T, U> e)
{
    const EnergyDensity<T, U> rho_e {energy_density(rho, e)};
    const Pressure<T, U> p {pressure(k.gamma, rho_e)};
    const Temperature<T, U> temp {temperature(p, rho, k.gas_constant)};
    const ThermoState<T, U> state {.rho = rho, .e = e, .p = p, .temp = temp};
    return {.state = state, .status = state_status(k, state, rho_e)};
}

// rho_e = p / (gamma - 1), T = p / (rho R).
template <Scalar T, Scalar C, UnitSystem U>
    requires ScalarPair<T, C>
CANARY_HOST_DEVICE Closure<T, U>
  perfect_from_rho_p(const PerfectCoefficients<C, U>& k, Density<T, U> rho, Pressure<T, U> p)
{
    const EnergyDensity<T, U> rho_e {energy_density(k.gamma, p)};
    const Temperature<T, U> temp {temperature(p, rho, k.gas_constant)};
    const ThermoState<T, U> state {.rho = rho,
                                   .e = specific_energy(rho_e, rho),
                                   .p = p,
                                   .temp = temp};
    return {.state = state, .status = state_status(k, state, rho_e)};
}

// p = rho R T, e = cv T.
template <Scalar T, Scalar C, UnitSystem U>
    requires ScalarPair<T, C>
CANARY_HOST_DEVICE Closure<T, U> perfect_from_rho_temp(const PerfectCoefficients<C, U>& k,
                                                       Density<T, U> rho,
                                                       Temperature<T, U> temp)
{
    const SpecificEnergy<T, U> e {specific_energy(k.cv, temp)};
    const ThermoState<T, U> state {
      .rho = rho,
      .e = e,
      .p = Pressure<T, U> {pressure(rho, k.gas_constant, temp)},
      .temp = temp,
    };
    return {.state = state, .status = state_status(k, state, energy_density(rho, e))};
}

// rho = p / (R T), e = cv T.
template <Scalar T, Scalar C, UnitSystem U>
    requires ScalarPair<T, C>
CANARY_HOST_DEVICE Closure<T, U>
  perfect_from_p_temp(const PerfectCoefficients<C, U>& k, Pressure<T, U> p, Temperature<T, U> temp)
{
    const Density<T, U> rho {density(p, k.gas_constant, temp)};
    const SpecificEnergy<T, U> e {specific_energy(k.cv, temp)};
    const ThermoState<T, U> state {.rho = rho, .e = e, .p = p, .temp = temp};
    return {.state = state, .status = state_status(k, state, energy_density(rho, e))};
}

// cv, cp, gamma and R, and the sound speed a = sqrt(gamma p / rho).
template <Scalar T, Scalar C, UnitSystem U>
    requires ScalarPair<T, C> && HasSqrt<T>
CANARY_HOST_DEVICE SingleTemperatureProperties<T, U>
  perfect_properties(const PerfectCoefficients<C, U>& k, const ThermoState<T, U>& state) noexcept
{
    return {
      .cv = SpecificHeat<T, U> {k.cv},
      .cp = SpecificHeat<T, U> {k.cp},
      .capacity_ratio = HeatCapacityRatio<T> {k.gamma},
      .gas_constant = SpecificGasConstant<T, U> {k.gas_constant},
      .acoustics = {.sound_speed =
                      SoundSpeed<T, U> {
                        sound_speed(AcousticExponent<C> {k.gamma.value()}, state.p, state.rho)},
                    .acoustic_exponent = AcousticExponent<T> {k.gamma.value()}},
    };
}

template <std::floating_point P, UnitSystem U>
CANARY_HOST_DEVICE constexpr PerfectCoefficients<P, U>
  coefficients(const CaloricallyPerfect<P, U>& gas)
{ return {.gamma = gas.gamma, .gas_constant = gas.gas_constant, .cv = gas.cv, .cp = gas.cp}; }

}  // namespace detail

// Closures of the one-species model. See the detail functions above for the
// formulas.
template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Closure<T, U>
  close(const CaloricallyPerfect<P, U>& gas, Density<T, U> rho, EnergyDensity<T, U> rho_e)
{ return detail::perfect_from_rho_rho_e(detail::coefficients(gas), rho, rho_e); }

template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Closure<T, U>
  close(const CaloricallyPerfect<P, U>& gas, Density<T, U> rho, SpecificEnergy<T, U> e)
{ return detail::perfect_from_rho_e(detail::coefficients(gas), rho, e); }

template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Closure<T, U>
  close(const CaloricallyPerfect<P, U>& gas, Density<T, U> rho, Pressure<T, U> p)
{ return detail::perfect_from_rho_p(detail::coefficients(gas), rho, p); }

template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Closure<T, U>
  close(const CaloricallyPerfect<P, U>& gas, Density<T, U> rho, Temperature<T, U> temp)
{ return detail::perfect_from_rho_temp(detail::coefficients(gas), rho, temp); }

template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Closure<T, U>
  close(const CaloricallyPerfect<P, U>& gas, Pressure<T, U> p, Temperature<T, U> temp)
{ return detail::perfect_from_p_temp(detail::coefficients(gas), p, temp); }

template <Scalar T, std::floating_point P, UnitSystem U>
    requires ScalarWith<T, P> && HasSqrt<T>
CANARY_HOST_DEVICE SingleTemperatureProperties<T, U>
  properties(const CaloricallyPerfect<P, U>& gas, const ThermoState<T, U>& state) noexcept
{ return detail::perfect_properties(detail::coefficients(gas), state); }

// Per-cell coefficients of a calorically perfect mixture: the weighted means
// of the species gas constants and heat capacities, and the status of the
// composition they came from.
template <Scalar T, UnitSystem U> struct Mixture {
    SpecificGasConstant<T, U> gas_constant;
    SpecificHeat<T, U> cv;
    SpecificHeat<T, U> cp;
    Status status;
};

// Ideal gas mixture of n species, each with constant heat capacities. Only
// the first n entries of each array are used.
template <std::floating_point P, UnitSystem U, std::size_t Capacity = CANARY_MAX_SPECIES>
    requires(Capacity > 0)
struct CaloricallyPerfectMixture {
    using parameter_type = P;
    using units = U;
    template <Scalar T> using State = ThermoState<T, U>;
    template <Scalar T> using Mixture = canary::Mixture<T, U>;

    int n;
    std::array<SpecificGasConstant<P, U>, Capacity> gas_constant;
    std::array<SpecificHeat<P, U>, Capacity> cv;
    std::array<SpecificHeat<P, U>, Capacity> cp;
    std::array<HeatCapacityRatio<P>, Capacity> gamma;
};

// A sized range whose elements have type Q, for building models.
template <typename R, typename Q>
concept SizedRangeOf =
  std::ranges::sized_range<R> && std::same_as<std::ranges::range_value_t<R>, Q>;

// Mixture model from the particle mass, charge and heat capacity ratio of
// each species and the Boltzmann constant, with its species data. Throws
// ParameterError when the ranges differ in length, are empty or longer than
// Capacity, or when a species parameter is out of range (see
// make_calorically_perfect).
template <std::size_t Capacity = CANARY_MAX_SPECIES,
          typename Masses,
          typename Charges,
          std::floating_point P,
          UnitSystem U,
          typename Gammas>
    requires SizedRangeOf<Masses, ParticleMass<P, U>> && SizedRangeOf<Charges, Charge<P, U>> &&
             SizedRangeOf<Gammas, HeatCapacityRatio<P>>
Built<CaloricallyPerfectMixture<P, U, Capacity>, SpeciesData<P, U, Capacity>>
  make_calorically_perfect_mixture(const Masses& masses,
                                   const Charges& charges,
                                   BoltzmannConstant<P, U> boltzmann,
                                   const Gammas& gammas)
{
    const std::size_t n = std::ranges::size(masses);
    if (std::ranges::size(charges) != n || std::ranges::size(gammas) != n) {
        throw ParameterError(
          std::format("species counts differ: {} masses, {} charges, {} heat capacity ratios",
                      n,
                      std::ranges::size(charges),
                      std::ranges::size(gammas)));
    }
    if (n == 0 || n > Capacity) {
        throw ParameterError(
          std::format("species count must be between 1 and {}, got {}", Capacity, n));
    }
    Built<CaloricallyPerfectMixture<P, U, Capacity>, SpeciesData<P, U, Capacity>> built {};
    built.model.n = static_cast<int>(n);
    built.species.n = static_cast<int>(n);
    auto mass = std::ranges::begin(masses);
    auto charge = std::ranges::begin(charges);
    auto gamma = std::ranges::begin(gammas);
    for (const SpeciesIndex i : species_indices(0, n)) {
        try {
            const auto one = make_calorically_perfect(*mass, *charge, boltzmann, *gamma);
            built.model.gas_constant[i.value()] = one.model.gas_constant;
            built.model.cv[i.value()] = one.model.cv;
            built.model.cp[i.value()] = one.model.cp;
            built.model.gamma[i.value()] = one.model.gamma;
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("species {}: {}", i.value(), error.what()));
        }
        built.species.mass[i.value()] = *mass;
        built.species.charge[i.value()] = *charge;
        ++mass;
        ++charge;
        ++gamma;
    }
    return built;
}

namespace detail
{

// Weighted means of R, cv and cp over the clamped weights, normalised by
// their sum.
template <std::floating_point P, UnitSystem U, std::size_t Capacity, SpeciesRange V, Scalar T>
    requires ScalarPair<T, P>
CANARY_HOST_DEVICE Mixture<T, U> mix(const CaloricallyPerfectMixture<P, U, Capacity>& gas,
                                     const V& weights,
                                     T total,
                                     Status status)
{
    // NOLINTNEXTLINE(misc-static-assert): run-time species count and extent
    CANARY_ASSERT(gas.n >= 1 &&
                  static_cast<std::size_t>(gas.n) == static_cast<std::size_t>(weights.extent(0)));
    const T w0 = clamped_weight(weights, SpeciesIndex {0}, total);
    T sum = w0;
    auto gas_constant = gas.gas_constant[0] * w0;
    auto cv = gas.cv[0] * w0;
    auto cp = gas.cp[0] * w0;
    for (const SpeciesIndex i : species_indices(1, static_cast<std::size_t>(gas.n))) {
        const T w = clamped_weight(weights, i, total);
        sum = sum + w;
        gas_constant = gas_constant + gas.gas_constant[i.value()] * w;
        cv = cv + gas.cv[i.value()] * w;
        cp = cp + gas.cp[i.value()] * w;
    }
    return {
      .gas_constant = SpecificGasConstant<T, U> {gas_constant / sum},
      .cv = SpecificHeat<T, U> {cv / sum},
      .cp = SpecificHeat<T, U> {cp / sum},
      .status = status,
    };
}

template <Scalar T, UnitSystem U>
CANARY_HOST_DEVICE constexpr PerfectCoefficients<T, U> coefficients(const Mixture<T, U>& mix)
{
    return {
      .gamma = heat_capacity_ratio(mix.cp, mix.cv),
      .gas_constant = mix.gas_constant,
      .cv = mix.cv,
      .cp = mix.cp,
    };
}

template <Scalar T, UnitSystem U>
CANARY_HOST_DEVICE constexpr Closure<T, U> with_composition_status(Closure<T, U> c,
                                                                   const Mixture<T, U>& mix)
{
    c.status = combine(c.status, mix.status);
    return c;
}

}  // namespace detail

// Per-cell mixture coefficients from mass fractions.
template <std::floating_point P, UnitSystem U, std::size_t Capacity, typename V>
    requires SpeciesWeights<V, kind::mass_fraction, unitless, P>
CANARY_HOST_DEVICE auto mixture(const CaloricallyPerfectMixture<P, U, Capacity>& gas,
                                const V& fractions)
{
    using T = range_quantity_t<V>::value_type;
    return detail::mix(gas,
                       fractions,
                       T(typename ScalarTraits<T>::value_type(1)),
                       composition_status(fractions));
}

// Per-cell mixture coefficients from partial densities and the density.
template <std::floating_point P, UnitSystem U, std::size_t Capacity, typename V, Scalar T>
    requires SpeciesWeights<V, kind::density, U, P> &&
             std::same_as<typename range_quantity_t<V>::value_type, T>
CANARY_HOST_DEVICE auto mixture(const CaloricallyPerfectMixture<P, U, Capacity>& gas,
                                const V& partial_densities,
                                Density<T, U> rho)
{
    return detail::mix(gas,
                       partial_densities,
                       rho.value(),
                       composition_status(partial_densities, rho));
}

// Per-cell mixture coefficients from partial densities and a previously computed status.
// status must equal composition_status(partial_densities, rho) for these inputs
// and must not be Invalid. Normalisation divides by the sum of clamped weights.
template <std::floating_point P, UnitSystem U, std::size_t Capacity, typename V, Scalar T>
    requires SpeciesWeights<V, kind::density, U, P> &&
             std::same_as<typename range_quantity_t<V>::value_type, T>
CANARY_HOST_DEVICE auto mixture(const CaloricallyPerfectMixture<P, U, Capacity>& gas,
                                const V& partial_densities,
                                Density<T, U> rho,
                                Status status) noexcept
{
    return detail::mix(gas, partial_densities, rho.value(), status);
}

// Closures of the mixture model, with the coefficients of a per-cell mixture
// value and the same formulas as the one-species model. The status combines
// the composition status with the state status.
template <Scalar T, std::floating_point P, UnitSystem U, std::size_t Capacity>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Closure<T, U> close(const CaloricallyPerfectMixture<P, U, Capacity>& /*gas*/,
                                       const Mixture<T, U>& mix,
                                       Density<T, U> rho,
                                       EnergyDensity<T, U> rho_e)
{
    return detail::with_composition_status(
      detail::perfect_from_rho_rho_e(detail::coefficients(mix), rho, rho_e),
      mix);
}

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t Capacity>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Closure<T, U> close(const CaloricallyPerfectMixture<P, U, Capacity>& /*gas*/,
                                       const Mixture<T, U>& mix,
                                       Density<T, U> rho,
                                       SpecificEnergy<T, U> e)
{
    return detail::with_composition_status(
      detail::perfect_from_rho_e(detail::coefficients(mix), rho, e),
      mix);
}

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t Capacity>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Closure<T, U> close(const CaloricallyPerfectMixture<P, U, Capacity>& /*gas*/,
                                       const Mixture<T, U>& mix,
                                       Density<T, U> rho,
                                       Pressure<T, U> p)
{
    return detail::with_composition_status(
      detail::perfect_from_rho_p(detail::coefficients(mix), rho, p),
      mix);
}

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t Capacity>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Closure<T, U> close(const CaloricallyPerfectMixture<P, U, Capacity>& /*gas*/,
                                       const Mixture<T, U>& mix,
                                       Density<T, U> rho,
                                       Temperature<T, U> temp)
{
    return detail::with_composition_status(
      detail::perfect_from_rho_temp(detail::coefficients(mix), rho, temp),
      mix);
}

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t Capacity>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE Closure<T, U> close(const CaloricallyPerfectMixture<P, U, Capacity>& /*gas*/,
                                       const Mixture<T, U>& mix,
                                       Pressure<T, U> p,
                                       Temperature<T, U> temp)
{
    return detail::with_composition_status(
      detail::perfect_from_p_temp(detail::coefficients(mix), p, temp),
      mix);
}

template <Scalar T, std::floating_point P, UnitSystem U, std::size_t Capacity>
    requires ScalarWith<T, P> && HasSqrt<T>
CANARY_HOST_DEVICE SingleTemperatureProperties<T, U>
  properties(const CaloricallyPerfectMixture<P, U, Capacity>& /*gas*/,
             const Mixture<T, U>& mix,
             const ThermoState<T, U>& state) noexcept
{ return detail::perfect_properties(detail::coefficients(mix), state); }

}  // namespace canary

#endif  // CANARY_MODELS_CALORICALLY_PERFECT_HPP
