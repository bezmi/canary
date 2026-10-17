#ifndef CANARY_TESTS_DEVICE_REACTIONS_INPUTS_HPP
#define CANARY_TESTS_DEVICE_REACTIONS_INPUTS_HPP

// Host fixtures of the reaction checks: input states and the models and
// mechanisms they are evaluated with. States are input data only, in SI
// units; no expected result is recorded here. Every mechanism uses the
// standard pressure 101325 Pa.

#include "host_inputs.hpp"
#include "thermodynamics_inputs.hpp"

#include <array>
#include <canary/data/air11.hpp>
#include <canary/data/air11_kim_jo.hpp>
#include <canary/data/air11_species.hpp>
#include <canary/data/air5.hpp>
#include <canary/data/air5_kim_jo.hpp>
#include <canary/data/air5_park.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/units.hpp>
#include <concepts>
#include <cstddef>

namespace canary::test
{

// The input records of the single-temperature air5 Kim-Jo reference states
// (tests/test_reactions.cpp).
inline constexpr std::array<FixtureState<5, 0>, 5> air5_kim_jo_states {{
  {.y = {0.767, 0.233, 0.0, 0.0, 0.0}, .rho = 0.006939829108798786, .temp = 5000.0, .modes = {}},
  {.y = {0.6, 0.1, 0.1, 0.15, 0.05}, .rho = 0.007037619261282653, .temp = 8000.0, .modes = {}},
  {.y = {0.3, 0.1, 0.3, 0.25, 0.05}, .rho = 0.0009537284289470741, .temp = 12000.0, .modes = {}},
  {.y = {0.6, 0.1, 0.1, 0.15, 0.05}, .rho = 1.1260190818052245, .temp = 250.0, .modes = {}},
  {.y = {0.3, 0.1, 0.3, 0.25, 0.05}, .rho = 0.00028611852868412217, .temp = 40000.0, .modes = {}},
}};

// The input records of the two-temperature air5 Park reference states
// (tests/test_energy_exchange.cpp); the mode temperature is T_ve.
inline constexpr std::array<FixtureState<5, 1>, 8> air5_park_states {{
  {.y = air5_mixed, .rho = 0.01, .temp = 10000.0, .modes = {3000.0}},
  {.y = air5_mixed, .rho = 0.001, .temp = 20000.0, .modes = {5000.0}},
  {.y = air5_dissociated, .rho = 0.01, .temp = 1500.0, .modes = {4000.0}},
  {.y = air5_dissociated, .rho = 0.01, .temp = 6000.0, .modes = {6000.000006}},
  {.y = air5_mixed, .rho = 0.01, .temp = 8000.0, .modes = {8000.0}},
  {.y = air5_air, .rho = 0.1, .temp = 5000.0, .modes = {1000.0}},
  {.y = air5_mixed, .rho = 0.001, .temp = 60000.0, .modes = {10000.0}},
  {.y = air5_dissociated, .rho = 0.1, .temp = 250.0, .modes = {2000.0}},
}};

// The input records of the three-temperature air11 Kim-Jo reference states
// (tests/test_multi_temperature_kinetics.cpp); the modes are T_v and T_e.
inline constexpr std::array<FixtureState<11, 2>, 4> air11_kim_jo_states {{
  {.y = air11_weakly_ionised, .rho = 0.001, .temp = 12000.0, .modes = {5000.0, 8000.0}},
  {.y = air11_strongly_ionised, .rho = 0.0001, .temp = 10000.0, .modes = {9000.0, 14000.0}},
  {.y = air11_weakly_ionised, .rho = 0.01, .temp = 7000.0, .modes = {7000.0, 7000.0}},
  {.y = air11_strongly_ionised, .rho = 0.0001, .temp = 3000.0, .modes = {6000.0, 5000.0}},
}};

// The states of the interleaved four-mode mechanism.
inline constexpr std::array<FixtureState<11, 4>, 2> interleaved_states {{
  {.y = air11_weakly_ionised,
   .rho = 1.0e-3,
   .temp = 12000.0,
   .modes = {3000.0, 5000.0, 7000.0, 11000.0}},
  {.y = air11_weakly_ionised,
   .rho = 1.0e-3,
   .temp = 8000.0,
   .modes = {8000.0, 8000.0, 8000.0, 8000.0}},
}};

// Four energy modes of eleven-species air: the vibration of N2; the vibration
// of O2, NO, N2+, O2+ and NO+; the bound electronic energy of the ten heavy
// species; and the energy of the free electrons.
inline constexpr std::array<EnergyModeInput<11>, 4> interleaved_modes {{
  {.count = 1, .components = {{{.species = "N2", .component = ModeComponent::vibration}}}},
  {.count = 5,
   .components = {{
     {.species = "O2", .component = ModeComponent::vibration},
     {.species = "NO", .component = ModeComponent::vibration},
     {.species = "N2+", .component = ModeComponent::vibration},
     {.species = "O2+", .component = ModeComponent::vibration},
     {.species = "NO+", .component = ModeComponent::vibration},
   }}},
  {.count = 10,
   .components = {{
     {.species = "N2", .component = ModeComponent::electronic},
     {.species = "O2", .component = ModeComponent::electronic},
     {.species = "N", .component = ModeComponent::electronic},
     {.species = "O", .component = ModeComponent::electronic},
     {.species = "NO", .component = ModeComponent::electronic},
     {.species = "N+", .component = ModeComponent::electronic},
     {.species = "O+", .component = ModeComponent::electronic},
     {.species = "N2+", .component = ModeComponent::electronic},
     {.species = "O2+", .component = ModeComponent::electronic},
     {.species = "NO+", .component = ModeComponent::electronic},
   }}},
  {.count = 1, .components = {{{.species = "e-", .component = ModeComponent::electron}}}},
}};

// The interleaved twelve-reaction mechanism: every forward rate (A, n, theta)
// = (1e10, 0, 4000 K) and every explicit backward rate (2e9, 0, 2000 K), SI
// per mole, with the rate-controlling and equilibrium temperatures of each
// reaction drawn from every physical source. With all_explicit, every
// equilibrium backward rate is replaced by the explicit rate with s = 0.5 and
// the N2 vibration temperature. Controlled numerical inputs, not a physical
// kinetics dataset.
inline MechanismInput<12, 3> interleaved_mechanism_input(bool all_explicit)
{
    using T = PhysicalTemperatureInput;
    const auto coefficients = [](std::initializer_list<SpeciesCoefficient> entries) {
        SpeciesCoefficients<3> out {.count = static_cast<int>(entries.size()), .entries = {}};
        std::size_t i = 0;
        for (const SpeciesCoefficient& entry : entries) { out.entries[i++] = entry; }
        return out;
    };
    const auto forward = [](double s, std::optional<T> other) {
        return ArrheniusInput {.A = 1.0e10,
                               .n = 0.0,
                               .activation_temperature = 4000.0,
                               .translational_exponent = s,
                               .other_temperature = other};
    };
    const auto explicit_backward = [](double s, std::optional<T> other) {
        return ArrheniusInput {.A = 2.0e9,
                               .n = 0.0,
                               .activation_temperature = 2000.0,
                               .translational_exponent = s,
                               .other_temperature = other};
    };
    const auto exchange_n2_o = [&] {
        return std::pair {coefficients({{"N2", 1.0}, {"O", 1.0}}),
                          coefficients({{"NO", 1.0}, {"N", 1.0}})};
    };
    const auto exchange_no_o = [&] {
        return std::pair {coefficients({{"NO", 1.0}, {"O", 1.0}}),
                          coefficients({{"O2", 1.0}, {"N", 1.0}})};
    };
    const auto exchange_o2_n = [&] {
        return std::pair {coefficients({{"O2", 1.0}, {"N", 1.0}}),
                          coefficients({{"NO", 1.0}, {"O", 1.0}})};
    };
    const auto ionisation = [&](const char* atom, const char* ion) {
        return std::pair {coefficients({{atom, 1.0}, {"e-", 1.0}}),
                          coefficients({{ion, 1.0}, {"e-", 2.0}})};
    };
    const auto dissociation = [&](const char* molecule, const char* atom) {
        return std::pair {coefficients({{molecule, 1.0}}), coefficients({{atom, 2.0}})};
    };
    const SpeciesCoefficients<3> none {.count = 0, .entries = {}};
    const SpeciesCoefficients<3> efficiencies =
      coefficients({{"N2", 2.0}, {"N", 3.0}, {"e-", 0.0}});
    const T translation = TranslationTemperatureInput {};
    const T rotation = RotationTemperatureInput {};
    const T free_electron = FreeElectronTemperatureInput {};
    const T n2_vibration = VibrationTemperatureInput {"N2"};
    const T o2_vibration = VibrationTemperatureInput {"O2"};
    const T n_electronic = BoundElectronicTemperatureInput {"N"};
    const T o_electronic = BoundElectronicTemperatureInput {"O"};

    // Each reaction: equation, species, forward rate, backward rate (an
    // explicit rate or the equilibrium temperature) and third body.
    const auto reaction = [&](const char* equation,
                              std::pair<SpeciesCoefficients<3>, SpeciesCoefficients<3>> species,
                              ArrheniusInput forward_rate,
                              std::optional<ArrheniusInput> backward_rate,
                              std::optional<T> equilibrium_temperature,
                              bool third_body) {
        ReactionInput<3> r {.equation = equation,
                            .reactants = species.first,
                            .products = species.second,
                            .forward = forward_rate,
                            .backward_model = BackwardRateModel::arrhenius,
                            .backward = {},
                            .equilibrium_temperature = {},
                            .third_body = third_body,
                            .efficiencies = third_body ? efficiencies : none};
        if (backward_rate) {
            r.backward = *backward_rate;
        } else if (all_explicit) {
            r.backward = explicit_backward(0.5, n2_vibration);
        } else {
            r.backward_model = BackwardRateModel::equilibrium;
            r.equilibrium_temperature = equilibrium_temperature;
        }
        return r;
    };
    const std::optional<ArrheniusInput> equilibrium {};
    return {
      .temperature_lower = 300.0,
      .temperature_upper = 30000.0,
      .reactions = {{
        reaction("N2 + O <=> NO + N",
                 exchange_n2_o(),
                 forward(0.5, n2_vibration),
                 equilibrium,
                 translation,
                 false),
        reaction("O + e- <=> O+ + 2 e-",
                 ionisation("O", "O+"),
                 forward(0.0, free_electron),
                 equilibrium,
                 free_electron,
                 false),
        reaction("NO + O <=> O2 + N",
                 exchange_no_o(),
                 forward(1.0, std::nullopt),
                 explicit_backward(0.5, o2_vibration),
                 std::nullopt,
                 false),
        reaction("O2 + N <=> NO + O",
                 exchange_o2_n(),
                 forward(0.5, o2_vibration),
                 equilibrium,
                 o2_vibration,
                 false),
        reaction("N2 + M <=> 2 N + M",
                 dissociation("N2", "N"),
                 forward(0.5, n2_vibration),
                 equilibrium,
                 translation,
                 true),
        reaction("N + e- <=> N+ + 2 e-",
                 ionisation("N", "N+"),
                 forward(0.25, n_electronic),
                 equilibrium,
                 n_electronic,
                 false),
        reaction("NO + O <=> O2 + N",
                 exchange_no_o(),
                 forward(0.0, rotation),
                 explicit_backward(0.0, free_electron),
                 std::nullopt,
                 false),
        reaction("N2 + O <=> NO + N",
                 exchange_n2_o(),
                 forward(0.0, n2_vibration),
                 equilibrium,
                 n2_vibration,
                 false),
        reaction("O2 + M <=> 2 O + M",
                 dissociation("O2", "O"),
                 forward(0.5, o2_vibration),
                 equilibrium,
                 o2_vibration,
                 true),
        reaction("N + e- <=> N+ + 2 e-",
                 ionisation("N", "N+"),
                 forward(0.0, free_electron),
                 equilibrium,
                 free_electron,
                 false),
        reaction("NO + O <=> O2 + N",
                 exchange_no_o(),
                 forward(0.5, o_electronic),
                 explicit_backward(1.0, std::nullopt),
                 std::nullopt,
                 false),
        reaction("N2 + O <=> NO + N",
                 exchange_n2_o(),
                 forward(1.0, translation),
                 equilibrium,
                 translation,
                 false),
      }},
    };
}

// The standard pressure of every fixture mechanism, 101325 Pa, in U.
template <std::floating_point P, UnitSystem U>
Pressure<P, U> standard_pressure(const UnitScales<P, U>& scales)
{ return convert<U>(Pressure<P, si> {P(101325)}, scales); }

// The interleaved four-mode model of air11.
template <std::floating_point P, UnitSystem U>
auto make_interleaved_model(const UnitScales<P, U>& scales)
{
    return make_multi_temperature_mixture(air11_species,
                                          air11_rotational_degrees,
                                          air11_electronic_levels,
                                          interleaved_modes,
                                          scales);
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_REACTIONS_INPUTS_HPP
