#ifndef CANARY_TESTS_DEVICE_THERMODYNAMICS_INPUTS_HPP
#define CANARY_TESTS_DEVICE_THERMODYNAMICS_INPUTS_HPP

// Host fixtures of the thermodynamic checks: models and the declared states
// they are evaluated at. States are input data only, in SI units; inputs
// derived from them by a forward closure (energies, pressures, energy
// densities) are generated offline into derived_inputs.hpp.

#include "host_inputs.hpp"

#include <array>
#include <canary/data/air11.hpp>
#include <canary/data/air11_species.hpp>
#include <canary/data/air5.hpp>
#include <canary/data/air5_electronic.hpp>
#include <canary/models/calorically_perfect.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/nasa9.hpp>
#include <canary/quantity.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <vector>

namespace canary::test
{

// Perfect gas: R = 287 J/(kg K), gamma = 1.4, at rho = 1 kg/m^3 and T = 300
// and 4000 K.
inline constexpr double perfect_gas_constant = 287.0;
inline constexpr double perfect_gamma = 1.4;
inline constexpr std::array<FixtureState<1, 0>, 2> perfect_states {{
  {.y = {1.0}, .rho = 1.0, .temp = 300.0, .modes = {}},
  {.y = {1.0}, .rho = 1.0, .temp = 4000.0, .modes = {}},
}};

template <std::floating_point P, UnitSystem U>
CaloricallyPerfect<P, U> make_perfect_gas(const UnitScales<P, U>& scales)
{
    return make_calorically_perfect(
      convert<U>(SpecificGasConstant<P, si> {P(perfect_gas_constant)}, scales),
      HeatCapacityRatio<P> {P(perfect_gamma)});
}

// Perfect mixture of n species (the benchmark mixture): particle masses
// 4.65e-26 (1 + 0.1 s) kg, no charge, gamma 1.4 for even and 5/3 for odd s.
template <std::size_t Capacity, std::floating_point P, UnitSystem U>
CaloricallyPerfectMixture<P, U, Capacity> make_perfect_mixture(std::size_t n,
                                                               const UnitScales<P, U>& scales)
{
    std::vector<ParticleMass<P, U>> masses;
    std::vector<Charge<P, U>> charges;
    std::vector<HeatCapacityRatio<P>> gammas;
    for (std::size_t s = 0; s < n; ++s) {
        masses.push_back(
          convert<U>(ParticleMass<P, si> {P(4.65e-26 * (1.0 + (0.1 * static_cast<double>(s))))},
                     scales));
        charges.push_back(convert<U>(Charge<P, si> {P(0)}, scales));
        gammas.emplace_back(P(s % 2 == 0 ? 1.4 : 5.0 / 3.0));
    }
    return make_calorically_perfect_mixture<Capacity>(masses,
                                                      charges,
                                                      boltzmann_constant(scales),
                                                      gammas)
      .model;
}

// Raw weights 1 + ((3 s) mod 7) of species s, normalised in binary64.
template <std::size_t N> std::array<double, N> perfect_mixture_fractions()
{
    std::array<double, N> y {};
    double total = 0.0;
    for (std::size_t s = 0; s < N; ++s) {
        y[s] = 1.0 + static_cast<double>((3 * s) % 7);
        total += y[s];
    }
    for (double& weight : y) { weight /= total; }
    return y;
}

template <std::size_t N> std::array<FixtureState<N, 0>, 2> perfect_mixture_states()
{
    const auto y = perfect_mixture_fractions<N>();
    return {{{.y = y, .rho = 1.0, .temp = 300.0, .modes = {}},
             {.y = y, .rho = 1.0, .temp = 4000.0, .modes = {}}}};
}

// The NASA9 fit of one segment on [50, 100000] K with constant cp/R = a2 and
// the integration constants a7 = 100 K and a8 = 2.
inline constexpr Nasa9Data<1> constant_fit(double a2)
{
    return {.segments = 1,
            .breaks = {50.0, 100000.0},
            .blends = {},
            .a = {{{0.0, 0.0, a2, 0.0, 0.0, 0.0, 0.0, 100.0, 2.0}}}};
}

// One-species thermally perfect gas with the N2 record of air5 and the
// constant fit cp/R = 4.
inline constexpr std::array<Nasa9Species<1>, 1> controlled_nitrogen {{
  {.name = air5_species[0].name,
   .molar_mass = air5_species[0].molar_mass,
   .charge = air5_species[0].charge,
   .metadata = air5_species[0].metadata,
   .fit = constant_fit(4.0)},
}};
inline constexpr std::array<FixtureState<1, 0>, 3> controlled_nitrogen_states {{
  {.y = {1.0}, .rho = 1.0, .temp = 300.0, .modes = {}},
  {.y = {1.0}, .rho = 1.0, .temp = 1000.0, .modes = {}},
  {.y = {1.0}, .rho = 1.0, .temp = 8000.0, .modes = {}},
}};

// One-species thermally perfect N2 of air5 with its named fit.
inline constexpr std::array<Nasa9Species<4>, 1> named_nitrogen {{air5_species[0]}};

// Named thermally perfect states: densities 0.01 and 1 kg/m^3 at 300, 1000,
// 6000, 10000 and 20000 K.
template <std::size_t N>
constexpr std::array<FixtureState<N, 0>, 10> named_states(const std::array<double, N>& y)
{
    std::array<FixtureState<N, 0>, 10> out {};
    std::size_t k = 0;
    for (const double rho : {0.01, 1.0}) {
        for (const double temp : {300.0, 1000.0, 6000.0, 10000.0, 20000.0}) {
            out[k++] = {.y = y, .rho = rho, .temp = temp, .modes = {}};
        }
    }
    return out;
}

// The relative temperature guesses of the thermally perfect energy closures,
// as factors of the generating temperature.
inline constexpr std::array<double, 2> guess_factors {0.9, 1.1};

// Two-temperature air5: air and mixed at rho = 0.01 kg/m^3.
inline constexpr std::array<FixtureState<5, 1>, 6> two_temperature_states {{
  {.y = air5_air, .rho = 0.01, .temp = 4000.0, .modes = {3000.0}},
  {.y = air5_air, .rho = 0.01, .temp = 8000.0, .modes = {8000.0}},
  {.y = air5_air, .rho = 0.01, .temp = 12000.0, .modes = {5000.0}},
  {.y = air5_mixed, .rho = 0.01, .temp = 4000.0, .modes = {3000.0}},
  {.y = air5_mixed, .rho = 0.01, .temp = 8000.0, .modes = {8000.0}},
  {.y = air5_mixed, .rho = 0.01, .temp = 12000.0, .modes = {5000.0}},
}};

// Four-temperature air5: mixed and air (absent NO mode) at rho = 0.01 kg/m^3.
inline constexpr std::array<FixtureState<5, 3>, 4> four_temperature_states {{
  {.y = air5_mixed, .rho = 0.01, .temp = 8000.0, .modes = {3000.0, 4000.0, 5000.0}},
  {.y = air5_mixed, .rho = 0.01, .temp = 8000.0, .modes = {8000.0, 8000.0, 8000.0}},
  {.y = air5_air, .rho = 0.01, .temp = 8000.0, .modes = {3000.0, 4000.0, 5000.0}},
  {.y = air5_air, .rho = 0.01, .temp = 8000.0, .modes = {8000.0, 8000.0, 8000.0}},
}};

// Three-temperature air11: weakly ionised, neutral and electron-dominated at
// rho = 0.01 kg/m^3.
inline constexpr std::array<double, 11> air11_neutral {
  0.3, 0.1, 0.3, 0.25, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
inline constexpr std::array<double, 11> air11_electron_dominated {
  0.999, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001};
inline constexpr std::array<FixtureState<11, 2>, 6> three_temperature_states {{
  {.y = air11_weakly_ionised, .rho = 0.01, .temp = 4000.0, .modes = {3000.0, 12000.0}},
  {.y = air11_weakly_ionised, .rho = 0.01, .temp = 8000.0, .modes = {8000.0, 8000.0}},
  {.y = air11_neutral, .rho = 0.01, .temp = 4000.0, .modes = {3000.0, 12000.0}},
  {.y = air11_neutral, .rho = 0.01, .temp = 8000.0, .modes = {8000.0, 8000.0}},
  {.y = air11_electron_dominated, .rho = 0.01, .temp = 4000.0, .modes = {3000.0, 12000.0}},
  {.y = air11_electron_dominated, .rho = 0.01, .temp = 8000.0, .modes = {8000.0, 8000.0}},
}};

// Controlled two-species modal model: N2 and N of air5 with constant fits
// (cp/R 4.5 and 2.5), N electronic levels (g, theta) = (2, 0) and (4, 4000 K)
// or the ground level alone, and modes 0 = N2 vibration, 1 = N bound
// electronic.
inline constexpr std::array<Nasa9Species<1>, 2> controlled_modal_species {{
  {.name = air5_species[0].name,
   .molar_mass = air5_species[0].molar_mass,
   .charge = air5_species[0].charge,
   .metadata = air5_species[0].metadata,
   .fit = constant_fit(4.5)},
  {.name = air5_species[2].name,
   .molar_mass = air5_species[2].molar_mass,
   .charge = air5_species[2].charge,
   .metadata = air5_species[2].metadata,
   .fit = constant_fit(2.5)},
}};
inline constexpr std::array<int, 2> controlled_modal_rotation {air5_rotational_degrees[0],
                                                               air5_rotational_degrees[2]};
inline constexpr std::array<ElectronicLevelsInput<2>, 1> controlled_two_levels {{
  {.name = "N", .levels = 2, .degeneracies = {2, 4}, .temperatures = {0.0, 4000.0}},
}};
inline constexpr std::array<ElectronicLevelsInput<2>, 1> controlled_ground_level {{
  {.name = "N", .levels = 1, .degeneracies = {2, 0}, .temperatures = {0.0, 0.0}},
}};
inline constexpr std::array<EnergyModeInput<1>, 2> controlled_modes {{
  {.count = 1, .components = {{{.species = "N2", .component = ModeComponent::vibration}}}},
  {.count = 1, .components = {{{.species = "N", .component = ModeComponent::electronic}}}},
}};

template <std::floating_point P, UnitSystem U>
auto make_controlled_modal(const std::array<ElectronicLevelsInput<2>, 1>& levels,
                           const UnitScales<P, U>& scales)
{
    return make_multi_temperature_mixture(controlled_modal_species,
                                          controlled_modal_rotation,
                                          levels,
                                          controlled_modes,
                                          scales)
      .model;
}

inline constexpr std::array<FixtureState<2, 2>, 4> controlled_modal_states {{
  {.y = {0.75, 0.25}, .rho = 0.01, .temp = 4000.0, .modes = {3000.0, 5000.0}},
  {.y = {0.75, 0.25}, .rho = 0.01, .temp = 8000.0, .modes = {8000.0, 8000.0}},
  {.y = {1.0, 0.0}, .rho = 0.01, .temp = 4000.0, .modes = {3000.0, 5000.0}},
  {.y = {1.0, 0.0}, .rho = 0.01, .temp = 8000.0, .modes = {8000.0, 8000.0}},
}};

// The ground-level-only model at the mixed composition.
inline constexpr std::array<FixtureState<2, 2>, 2> controlled_ground_level_states {
  {controlled_modal_states[0], controlled_modal_states[1]}};

// Direct NASA9 temperatures [K] of the constant one-segment curve.
inline constexpr std::array<double, 6> constant_curve_temperatures {
  25.0, 50.0, 300.0, 1000.0, 100000.0, 120000.0};

// The seven term-isolated two-segment curves: breaks {50, 1000, 3000} K,
// transition width 1800 K, a2 = 3 on both segments, upper a[k] increased by
// 1000^(2 - k), lower and upper a7 = {100, 200}, a8 = {2, 5}.
inline Nasa9Data<2> term_isolated_fit(std::size_t k)
{
    Nasa9Data<2> fit {.segments = 2,
                      .breaks = {50.0, 1000.0, 3000.0},
                      .blends = {1800.0, 0.0},
                      .a = {}};
    fit.a[0][2] = 3.0;
    fit.a[1][2] = 3.0;
    fit.a[1][k] += std::pow(1000.0, 2.0 - static_cast<double>(k));
    fit.a[0][7] = 100.0;
    fit.a[1][7] = 200.0;
    fit.a[0][8] = 2.0;
    fit.a[1][8] = 5.0;
    return fit;
}
inline constexpr std::array<double, 16> term_isolated_temperatures {25.0,
                                                                    50.0,
                                                                    99.0,
                                                                    100.0,
                                                                    101.0,
                                                                    199.0,
                                                                    201.0,
                                                                    298.15,
                                                                    800.0,
                                                                    1000.0,
                                                                    1800.0,
                                                                    1899.0,
                                                                    1900.0,
                                                                    1901.0,
                                                                    3000.0,
                                                                    3100.0};

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_THERMODYNAMICS_INPUTS_HPP
