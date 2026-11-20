// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_TESTS_DEVICE_EXCHANGE_INPUTS_HPP
#define CANARY_TESTS_DEVICE_EXCHANGE_INPUTS_HPP

// Host fixtures of the energy-exchange checks: exchange inputs with every
// process configured, none, or exactly one, and the declared states they are
// evaluated at, in SI units. The full states (energies, pressures) and the
// rates of progress supplied to the standalone chemistry-coupled terms at the
// named states are derived offline (derived_inputs.hpp).

#include "host_inputs.hpp"
#include "reactions_inputs.hpp"
#include "thermodynamics_inputs.hpp"

#include <array>
#include <canary/data/air11_kim_jo_exchange.hpp>
#include <canary/data/air11_species.hpp>
#include <canary/data/air5.hpp>
#include <canary/data/air5_park_exchange.hpp>
#include <canary/energy_exchange.hpp>
#include <cstddef>
#include <string_view>

namespace canary::test
{

// The names of species in their model order.
template <std::size_t S, std::size_t N>
constexpr std::array<std::string_view, N>
  species_names(const std::array<Nasa9Species<S>, N>& species)
{
    std::array<std::string_view, N> names {};
    for (std::size_t i = 0; i < N; ++i) { names[i] = species[i].name; }
    return names;
}

// Exchange inputs with no process configured.
inline constexpr EnergyExchangeInput<0, 0> no_exchange {};
inline constexpr ElectronExchangeInput<0, 0, 0> no_electron_exchange {};

// The first vibration-translation pair of a named input with the
// high-temperature correction model and sigma_0 = 1e-21 m^2,
// T_ref = 50000 K and omega = 0.5 (each used only by its model).
template <std::size_t V, std::size_t C>
constexpr EnergyExchangeInput<1, 0>
  first_vibration_translation(const EnergyExchangeInput<V, C>& input,
                              RelaxationCorrection correction)
{
    VibrationTranslationInput pair = input.vibration_translation[0];
    pair.correction = {.model = correction,
                       .cross_section = 1.0e-21,
                       .reference_temperature = 50000.0,
                       .exponent = 0.5};
    return {.vibration_translation = {pair}, .chemistry_vibration = {}};
}

inline constexpr std::array<RelaxationCorrection, 4> corrections {RelaxationCorrection::none,
                                                                  RelaxationCorrection::park_1985,
                                                                  RelaxationCorrection::park_1993,
                                                                  RelaxationCorrection::kim_2021};

// The first chemistry-vibration coupling of a named input alone.
template <std::size_t V, std::size_t C>
constexpr EnergyExchangeInput<0, 1>
  first_chemistry_vibration(const EnergyExchangeInput<V, C>& input)
{ return {.vibration_translation = {}, .chemistry_vibration = {input.chemistry_vibration[0]}}; }

// One electron process of air11 alone: the elastic pair `pair`, the first
// vibration-electron exchange or the first electron-impact ionisation.
inline constexpr ElectronExchangeInput<1, 0, 0> air11_electron_translation(std::size_t pair)
{
    return {.electron_translation = {air11_kim_jo_exchange_electrons.electron_translation[pair]},
            .vibration_electron = {},
            .electron_impact_ionisation = {}};
}
inline constexpr ElectronExchangeInput<0, 1, 0> air11_vibration_electron {
  .electron_translation = {},
  .vibration_electron = {air11_kim_jo_exchange_electrons.vibration_electron[0]},
  .electron_impact_ionisation = {}};
inline constexpr ElectronExchangeInput<0, 0, 1> air11_electron_impact_ionisation {
  .electron_translation = {},
  .vibration_electron = {},
  .electron_impact_ionisation = {air11_kim_jo_exchange_electrons.electron_impact_ionisation[0]}};

// The neutral (first) and charged (sixth) elastic pairs of air11.
inline constexpr std::size_t neutral_pair = 0;
inline constexpr std::size_t charged_pair = 5;

// Two-temperature process states: air5 mixed at rho = 0.01 kg/m^3 with
// (T, T_ve) = (10000, 3000) K and (8000, 8000) K.
inline constexpr std::array<FixtureState<5, 1>, 2> two_temperature_process_states {{
  {.y = air5_mixed, .rho = 0.01, .temp = 10000.0, .modes = {3000.0}},
  {.y = air5_mixed, .rho = 0.01, .temp = 8000.0, .modes = {8000.0}},
}};

// Multi-temperature process states of air11 at rho = 1e-3 kg/m^3: weakly
// ionised with (T, T_v, T_e) = (12000, 5000, 8000) K and all at 8000 K; the
// same composition with the electron weight moved to N2; and the
// electron-temperature series of the elastic and vibration-electron terms.
inline constexpr std::array<double, 11> air11_no_electrons = [] {
    std::array<double, 11> y = air11_weakly_ionised;
    y[0] += y[10];
    y[10] = 0.0;
    return y;
}();

template <std::size_t C>
constexpr std::array<FixtureState<11, 2>, C>
  electron_temperature_series(const std::array<double, C>& t_e)
{
    std::array<FixtureState<11, 2>, C> out {};
    for (std::size_t k = 0; k < C; ++k) {
        out[k] = {.y = air11_weakly_ionised,
                  .rho = 1.0e-3,
                  .temp = 12000.0,
                  .modes = {5000.0, t_e[k]}};
    }
    return out;
}

inline constexpr std::array<FixtureState<11, 2>, 3> multi_temperature_process_states {{
  {.y = air11_weakly_ionised, .rho = 1.0e-3, .temp = 12000.0, .modes = {5000.0, 8000.0}},
  {.y = air11_weakly_ionised, .rho = 1.0e-3, .temp = 8000.0, .modes = {8000.0, 8000.0}},
  {.y = air11_no_electrons, .rho = 1.0e-3, .temp = 12000.0, .modes = {5000.0, 8000.0}},
}};
inline constexpr auto neutral_electron_states =
  electron_temperature_series(std::array<double, 3> {8000.0, 30000.0, 45000.0});
inline constexpr auto charged_electron_states =
  electron_temperature_series(std::array<double, 2> {8000.0, 45000.0});
inline constexpr auto vibration_electron_states =
  electron_temperature_series(std::array<double, 5> {6000.0, 6950.0, 6980.0, 7050.0, 8000.0});

// Independent rates of progress of the isolated chemistry-coupled terms
// [particles/(m^3 s)]: every reaction zero except the coupled one, with
// (forward, backward) = (2, 1), (1, 2) and, for ionisation, (1, 1).
inline constexpr std::array<std::array<double, 2>, 2> coupling_rates {{{2.0, 1.0}, {1.0, 2.0}}};
inline constexpr std::array<std::array<double, 2>, 3> ionisation_rates {
  {{2.0, 1.0}, {1.0, 2.0}, {1.0, 1.0}}};

// The rates of a mechanism of R reactions that are zero except reaction k.
template <std::size_t R>
constexpr DerivedRates<R> isolated_rates(std::size_t k, const std::array<double, 2>& rates)
{
    DerivedRates<R> out {};
    out.forward[k] = rates[0];
    out.backward[k] = rates[1];
    return out;
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_EXCHANGE_INPUTS_HPP
