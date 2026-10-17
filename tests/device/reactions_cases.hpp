#ifndef CANARY_TESTS_DEVICE_REACTIONS_CASES_HPP
#define CANARY_TESTS_DEVICE_REACTIONS_CASES_HPP

// Reaction checks for the host and native compilations (see cell_inputs.hpp).
// rates_of_progress and production_rates are separate checks: production uses
// its own net-rate evaluation.

#include "cell_inputs.hpp"
#include "check_groups.hpp"
#include "host_inputs.hpp"
#include "native_runtime.hpp"
#include "reactions_inputs.hpp"

#include <array>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/state.hpp>
#include <concepts>
#include <cstddef>
#include <string>

namespace canary::test
{

// The including file defines CANARY_CHECK_ENTRY (see thermodynamics_cases.hpp).

// Rates of progress of the selected cell at one temperature, or with the scalar
// T_ve of a one-mode model as the extra argument.
template <Storage S, std::size_t N, typename M, typename G, typename T, typename... Extra>
CANARY_CHECK_ENTRY void
  rates_check(Output<RatesOfProgress<T, typename G::units, M::reaction_count>> out,
              const M* mechanism,
              const G* gas,
              const T* composition,
              Density<T, typename G::units> rho,
              Temperature<T, typename G::units> temp,
              Extra... extra)
{
    out[0] = rates_of_progress(*mechanism,
                               *gas,
                               cell_composition<Composition::Mixture, S, N>(*gas, composition),
                               rho,
                               temp,
                               extra...);
}

template <Storage S, std::size_t N, typename M, typename G, typename T, typename... Extra>
CANARY_CHECK_ENTRY void
  production_check(Output<std::array<MassProductionRate<T, typename G::units>, N>> out,
                   const M* mechanism,
                   const G* gas,
                   const T* composition,
                   Density<T, typename G::units> rho,
                   Temperature<T, typename G::units> temp,
                   Extra... extra)
{
    out[0] = production_rates(*mechanism,
                              *gas,
                              cell_composition<Composition::Mixture, S, N>(*gas, composition),
                              rho,
                              temp,
                              extra...);
}

// As above with the mode temperatures of the selected cell in form F.
template <ModeForm F, Storage S, std::size_t N, typename M, typename G, typename T>
CANARY_CHECK_ENTRY void
  modal_rates_check(Output<RatesOfProgress<T, typename G::units, M::reaction_count>> out,
                    const M* mechanism,
                    const G* gas,
                    const T* composition,
                    Density<T, typename G::units> rho,
                    Temperature<T, typename G::units> temp,
                    const T* modes)
{
    out[0] = rates_of_progress(
      *mechanism,
      *gas,
      cell_composition<Composition::Mixture, S, N>(*gas, composition),
      rho,
      temp,
      cell_modes<F, S, kind::mode_temperature, typename G::units, G::mode_count>(modes));
}

template <ModeForm F, Storage S, std::size_t N, typename M, typename G, typename T>
CANARY_CHECK_ENTRY void
  modal_production_check(Output<std::array<MassProductionRate<T, typename G::units>, N>> out,
                         const M* mechanism,
                         const G* gas,
                         const T* composition,
                         Density<T, typename G::units> rho,
                         Temperature<T, typename G::units> temp,
                         const T* modes)
{
    out[0] = production_rates(
      *mechanism,
      *gas,
      cell_composition<Composition::Mixture, S, N>(*gas, composition),
      rho,
      temp,
      cell_modes<F, S, kind::mode_temperature, typename G::units, G::mode_count>(modes));
}

// Runners. Each evaluates the checks of one model and mechanism at every
// state of a fixture and every seed (see host_inputs.hpp).

// The inputs of one reaction state: the raw fractions, density, temperature
// and, when Modes, the mode temperatures; and the fixture description.
template <typename T, UnitSystem U, std::size_t N, std::size_t K> struct ReactionInputs {
    std::array<T, N> y;
    Density<T, U> rho;
    Temperature<T, U> temp;
    std::array<T, K> modes;
    std::string fixture;
};

// Calls evaluate(inputs) for every state of states and every seed. Mode
// temperatures are inputs only when Modes.
template <typename T,
          bool Modes,
          std::size_t N,
          std::size_t K,
          std::size_t C,
          typename P,
          typename U,
          typename F>
void for_each_reaction_input(const std::string& name,
                             const std::array<FixtureState<N, K>, C>& states,
                             const UnitScales<P, U>& scales,
                             F&& evaluate)
{
    for (std::size_t s = 0; s < C; ++s) {
        for_each_seed<T>([&](Seeds<T>& seeds) {
            const FixtureState<N, K>& state = states[s];
            ReactionInputs<T, U, N, K> inputs {
              .y = seeds.fractions("y", state.y),
              .rho = seeds.template quantity<kind::density>("rho", state.rho, scales),
              .temp = seeds.template quantity<kind::temperature>("temp", state.temp, scales),
              .modes = {},
              .fixture = {}};
            if constexpr (Modes) {
                inputs.modes =
                  seeds.template values<kind::mode_temperature>("modes", state.modes, scales);
            }
            inputs.fixture = name + " state " + std::to_string(s) + ", " + seeds.label();
            evaluate(inputs);
        });
    }
}

// Rates and production at one temperature (no mode argument).
template <typename T,
          typename M,
          typename G,
          std::size_t N,
          std::size_t K,
          std::size_t C,
          typename P,
          typename U>
void single_rates(Report& report,
                  const std::string& name,
                  const M& mechanism,
                  const G& gas,
                  const std::array<FixtureState<N, K>, C>& states,
                  const UnitScales<P, U>& scales)
{
    for_each_reaction_input<T, false>(name, states, scales, [&](const auto& in) {
        const auto y = cell_storage<Storage::Contiguous>(in.y);
        evaluate<&rates_check<Storage::Contiguous, N, M, G, T>>(report,
                                                                in.fixture,
                                                                mechanism,
                                                                gas,
                                                                y,
                                                                in.rho,
                                                                in.temp);
        evaluate<&production_check<Storage::Contiguous, N, M, G, T>>(report,
                                                                     in.fixture,
                                                                     mechanism,
                                                                     gas,
                                                                     y,
                                                                     in.rho,
                                                                     in.temp);
    });
}

// Rates and production with mode temperatures in form F, and the composition
// and mode temperatures in storage S.
template <ModeForm F,
          Storage S,
          typename T,
          typename M,
          typename G,
          std::size_t N,
          std::size_t K,
          std::size_t C,
          typename P,
          typename U>
void modal_rates(Report& report,
                 const std::string& name,
                 const M& mechanism,
                 const G& gas,
                 const std::array<FixtureState<N, K>, C>& states,
                 const UnitScales<P, U>& scales)
{
    for_each_reaction_input<T, true>(name, states, scales, [&](const auto& in) {
        const auto y = cell_storage<S>(in.y);
        const auto modes = cell_storage<S>(in.modes);
        evaluate<&modal_rates_check<F, S, N, M, G, T>>(report,
                                                       in.fixture,
                                                       mechanism,
                                                       gas,
                                                       y,
                                                       in.rho,
                                                       in.temp,
                                                       modes);
        evaluate<&modal_production_check<F, S, N, M, G, T>>(report,
                                                            in.fixture,
                                                            mechanism,
                                                            gas,
                                                            y,
                                                            in.rho,
                                                            in.temp,
                                                            modes);
    });
}

// Rates and production with each borrowed mode form and grid storage once.
template <typename T,
          typename M,
          typename G,
          std::size_t N,
          std::size_t K,
          std::size_t C,
          typename P,
          typename U>
void borrowed_rates(Report& report,
                    const std::string& name,
                    const M& mechanism,
                    const G& gas,
                    const std::array<FixtureState<N, K>, C>& states,
                    const UnitScales<P, U>& scales)
{
    modal_rates<ModeForm::Quantities, Storage::LayoutRight, T>(report,
                                                               name,
                                                               mechanism,
                                                               gas,
                                                               states,
                                                               scales);
    modal_rates<ModeForm::Indexed, Storage::LayoutLeft, T>(report,
                                                           name,
                                                           mechanism,
                                                           gas,
                                                           states,
                                                           scales);
}

template <std::floating_point P, typename T, UnitSystem U> Report reaction_checks()
{
    Report report;
    const auto scales = fixture_scales<P, U>();
    const auto p0 = standard_pressure(scales);
    // Single temperature: air5 with the 13-reaction mechanism.
    {
        const auto gas = make_air5(scales).model;
        single_rates<T>(report,
                        "air5 Kim-Jo",
                        make_air5_kim_jo(scales, p0),
                        gas,
                        air5_kim_jo_states,
                        scales);
    }
    // Two temperature: Park's 5 reactions with scalar T_ve and owning modes,
    // and the explicit single-temperature mechanism.
    {
        const auto gas = make_air5_two_temperature(scales).model;
        const auto mechanism = make_air5_park(scales, p0, gas);
        using M = std::remove_const_t<decltype(mechanism)>;
        using G = std::remove_const_t<decltype(gas)>;
        for_each_reaction_input<T, true>(
          "air5 Park two-temperature",
          air5_park_states,
          scales,
          [&](const auto& in) {
              const auto y = cell_storage<Storage::Contiguous>(in.y);
              const Temperature<T, U> temp_ve {in.modes[0]};
              evaluate<&rates_check<Storage::Contiguous, 5, M, G, T, Temperature<T, U>>>(report,
                                                                                         in.fixture,
                                                                                         mechanism,
                                                                                         gas,
                                                                                         y,
                                                                                         in.rho,
                                                                                         in.temp,
                                                                                         temp_ve);
              evaluate<&production_check<Storage::Contiguous, 5, M, G, T, Temperature<T, U>>>(
                report,
                in.fixture,
                mechanism,
                gas,
                y,
                in.rho,
                in.temp,
                temp_ve);
          });
        modal_rates<ModeForm::Owning, Storage::Contiguous, T>(report,
                                                              "air5 Park two-temperature",
                                                              mechanism,
                                                              gas,
                                                              air5_park_states,
                                                              scales);
        single_rates<T>(report,
                        "air5 Park single-temperature",
                        make_air5_park(scales, p0),
                        gas,
                        air5_park_states,
                        scales);
    }
    // Three temperature: air11 with 40 reactions, every mode form.
    {
        const auto gas = make_air11_three_temperature(scales).model;
        const auto mechanism = make_air11_kim_jo(scales, p0, gas);
        const std::string name = "air11 Kim-Jo three-temperature";
        modal_rates<ModeForm::Owning, Storage::Contiguous, T>(report,
                                                              name,
                                                              mechanism,
                                                              gas,
                                                              air11_kim_jo_states,
                                                              scales);
        borrowed_rates<T>(report, name, mechanism, gas, air11_kim_jo_states, scales);
        single_rates<T>(report,
                        "air11 Kim-Jo single-temperature",
                        make_air11_kim_jo(scales, p0),
                        gas,
                        air11_kim_jo_states,
                        scales);
    }
    // Four internal modes: the 12-reaction interleaved mechanism, with
    // equilibrium and all-explicit backward rates.
    {
        const auto gas = make_interleaved_model(scales).model;
        for (const bool all_explicit : {false, true}) {
            const auto input = interleaved_mechanism_input(all_explicit);
            const auto species = reaction_species(air11_species);
            const std::string name =
              all_explicit ? "interleaved all-explicit" : "interleaved equilibrium";
            modal_rates<ModeForm::Owning, Storage::Contiguous, T>(
              report,
              name,
              make_reaction_mechanism(input, species, scales, p0, gas),
              gas,
              interleaved_states,
              scales);
            single_rates<T>(report,
                            name + " single-temperature",
                            make_reaction_mechanism(input, species, scales, p0),
                            gas,
                            interleaved_states,
                            scales);
        }
    }
    // Park's reactions against four-temperature air5 (used without electrons).
    {
        const auto gas = make_air5_four_temperature(scales).model;
        modal_rates<ModeForm::Owning, Storage::Contiguous, T>(report,
                                                              "air5 Park four-temperature",
                                                              make_air5_park(scales, p0, gas),
                                                              gas,
                                                              four_temperature_states,
                                                              scales);
    }
    return report;
}

// Multi-temperature rates and production with P and T of different precision.
template <std::floating_point P, typename T, UnitSystem U> Report reaction_mixed_checks()
{
    Report report;
    const auto scales = fixture_scales<P, U>();
    const auto p0 = standard_pressure(scales);
    {
        const auto gas = make_air11_three_temperature(scales).model;
        modal_rates<ModeForm::Owning, Storage::Contiguous, T>(report,
                                                              "air11 Kim-Jo three-temperature",
                                                              make_air11_kim_jo(scales, p0, gas),
                                                              gas,
                                                              air11_kim_jo_states,
                                                              scales);
    }
    {
        const auto gas = make_interleaved_model(scales).model;
        const auto species = reaction_species(air11_species);
        for (const bool all_explicit : {false, true}) {
            modal_rates<ModeForm::Owning, Storage::Contiguous, T>(
              report,
              all_explicit ? "interleaved all-explicit" : "interleaved equilibrium",
              make_reaction_mechanism(interleaved_mechanism_input(all_explicit),
                                      species,
                                      scales,
                                      p0,
                                      gas),
              gas,
              interleaved_states,
              scales);
        }
    }
    return report;
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_REACTIONS_CASES_HPP
