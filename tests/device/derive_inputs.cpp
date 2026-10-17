// Generator of tests/device/derived_inputs.hpp: the inputs of the device-check
// fixtures that a forward closure derives from a declared state (density,
// specific energy, pressure and mode specific energies). Each fixture state of
// the {subsystem}_inputs.hpp headers is closed in binary64 and SI units with
// Canary's host closure of the fixture model, and each derived value is
// written as a hexadecimal binary64 literal. The values are inputs only; no
// expected result is derived here.
//
// Usage: canary_device_derive_inputs <source description> > derived_inputs.hpp
// The source description (for example the Canary commit) is recorded in the
// generated header. Exit status 1 if a closure is not valid.

#define CANARY_HOST_DEVICE
#include "exchange_inputs.hpp"
#include "host_inputs.hpp"
#include "reactions_inputs.hpp"
#include "thermodynamics_inputs.hpp"
#include "transport_inputs.hpp"

#include <array>
#include <canary/data/air11.hpp>
#include <canary/data/air5.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/thermo.hpp>
#include <canary/units.hpp>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <mdspan/mdspan.hpp>
#include <string>
#include <string_view>
#include <vector>

namespace
{

using namespace canary;
using namespace canary::test;

const UnitScales<double, si> scales {};

template <typename Closure> void require_valid(const Closure& closure, std::string_view fixture)
{
    if (closure.status != Status::Valid) {
        std::fprintf(stderr,
                     "%.*s: closure status %d, reason %d\n",
                     static_cast<int>(fixture.size()),
                     fixture.data(),
                     static_cast<int>(closure.status),
                     static_cast<int>(closure.diagnostic.reason));
        std::exit(1);
    }
}

// The derived inputs of a closed state of a model with M modes.
template <std::size_t M, typename State> DerivedState<M> derived(const State& state)
{
    DerivedState<M> out {.rho = state.rho.value(),
                         .e = state.e.value(),
                         .p = state.p.value(),
                         .mode_e = {}};
    if constexpr (M > 0) {
        for (std::size_t m = 0; m < M; ++m) { out.mode_e[m] = state.mode_e[m].value(); }
    }
    return out;
}

// Closes every state of states from (rho, T) with the mass fractions of the
// state and, when Modal, its mode temperatures; otherwise the modes of a modal
// model follow T.
template <std::size_t M, bool Modal, typename G, std::size_t N, std::size_t K, std::size_t C>
std::array<DerivedState<M>, C>
  derive(const G& gas, const std::array<FixtureState<N, K>, C>& states, std::string_view fixture)
{
    std::array<DerivedState<M>, C> out {};
    for (std::size_t s = 0; s < C; ++s) {
        const FixtureState<N, K>& state = states[s];
        const auto y = as_kind<kind::mass_fraction>(
          Kokkos::mdspan<const double, Kokkos::extents<std::size_t, N>>(state.y.data()));
        const Density<double, si> rho {state.rho};
        const Temperature<double, si> temp {state.temp};
        if constexpr (Modal) {
            const auto modes = as_kind<kind::mode_temperature, si>(
              Kokkos::mdspan<const double, Kokkos::extents<std::size_t, K>>(state.modes.data()));
            const auto closure = thermo_state(gas, y, rho, temp, modes);
            require_valid(closure, fixture);
            out[s] = derived<M>(closure.state);
        } else {
            const auto closure = thermo_state(gas, y, rho, temp);
            require_valid(closure, fixture);
            out[s] = derived<M>(closure.state);
        }
    }
    return out;
}

// As derive, for states declared by pressure: closes (p, T) with the mass
// fractions and, when Modal, the mode temperatures of the state.
template <std::size_t M, bool Modal, typename G, std::size_t N, std::size_t K, std::size_t C>
std::array<DerivedState<M>, C> derive_from_pressure(
  const G& gas, const std::array<PressureState<N, K>, C>& states, std::string_view fixture)
{
    std::array<DerivedState<M>, C> out {};
    for (std::size_t s = 0; s < C; ++s) {
        const PressureState<N, K>& state = states[s];
        const auto y = as_kind<kind::mass_fraction>(
          Kokkos::mdspan<const double, Kokkos::extents<std::size_t, N>>(state.y.data()));
        const Pressure<double, si> p {state.p};
        const Temperature<double, si> temp {state.temp};
        if constexpr (Modal) {
            const auto modes = as_kind<kind::mode_temperature, si>(
              Kokkos::mdspan<const double, Kokkos::extents<std::size_t, K>>(state.modes.data()));
            const auto closure = thermo_state(gas, y, p, temp, modes);
            require_valid(closure, fixture);
            out[s] = derived<M>(closure.state);
        } else {
            const auto closure = thermo_state(gas, y, p, temp);
            require_valid(closure, fixture);
            out[s] = derived<M>(closure.state);
        }
    }
    return out;
}

// The perfect gas closed from (rho, T) at every state.
template <std::size_t N, std::size_t C>
std::array<DerivedState<0>, C> derive_perfect(const std::array<FixtureState<N, 0>, C>& states)
{
    const auto gas = make_perfect_gas(scales);
    std::array<DerivedState<0>, C> out {};
    for (std::size_t s = 0; s < C; ++s) {
        const auto closure = thermo_state(gas,
                                          Density<double, si> {states[s].rho},
                                          Temperature<double, si> {states[s].temp});
        require_valid(closure, "perfect gas");
        out[s] = derived<0>(closure.state);
    }
    return out;
}

void print_value(double value) { std::printf("%a", value); }

// Writes derived as an array named name.
template <std::size_t M, std::size_t C>
void print(std::string_view name, const std::array<DerivedState<M>, C>& derived)
{
    std::printf("inline constexpr std::array<DerivedState<%zu>, %zu> %.*s {{\n",
                M,
                C,
                static_cast<int>(name.size()),
                name.data());
    for (const DerivedState<M>& d : derived) {
        std::printf("  {.rho = ");
        print_value(d.rho);
        std::printf(", .e = ");
        print_value(d.e);
        std::printf(", .p = ");
        print_value(d.p);
        std::printf(", .mode_e = {");
        for (std::size_t m = 0; m < M; ++m) {
            if (m > 0) { std::printf(", "); }
            print_value(d.mode_e[m]);
        }
        std::printf("}},\n");
    }
    std::printf("}};\n\n");
}

void thermodynamics()
{
    print("perfect_derived", derive_perfect(perfect_states));
    print("perfect_mixture5_derived",
          derive<0, false>(make_perfect_mixture<5>(5, scales),
                           perfect_mixture_states<5>(),
                           "perfect mixture 5"));
    print("perfect_mixture11_derived",
          derive<0, false>(make_perfect_mixture<11>(11, scales),
                           perfect_mixture_states<11>(),
                           "perfect mixture 11"));
    print("controlled_nitrogen_derived",
          derive<0, false>(make_thermally_perfect_mixture(controlled_nitrogen, scales).model,
                           controlled_nitrogen_states,
                           "controlled nitrogen"));
    print("named_nitrogen_derived",
          derive<0, false>(make_thermally_perfect_mixture(named_nitrogen, scales).model,
                           named_states<1>({1.0}),
                           "named nitrogen"));
    {
        const auto gas = make_air5(scales).model;
        print("air5_air_derived", derive<0, false>(gas, named_states(air5_air), "air5 air"));
        print("air5_mixed_derived", derive<0, false>(gas, named_states(air5_mixed), "air5 mixed"));
    }
    {
        const auto gas = make_air5_two_temperature(scales).model;
        print("two_temperature_equilibrium_derived",
              derive<1, false>(gas, two_temperature_states, "two-temperature equilibrium"));
        print("two_temperature_modal_derived",
              derive<1, true>(gas, two_temperature_states, "two-temperature"));
    }
    {
        const auto gas = make_air5_four_temperature(scales).model;
        print("four_temperature_equilibrium_derived",
              derive<3, false>(gas, four_temperature_states, "four-temperature equilibrium"));
        print("four_temperature_modal_derived",
              derive<3, true>(gas, four_temperature_states, "four-temperature"));
    }
    {
        const auto gas = make_air11_three_temperature(scales).model;
        print("three_temperature_equilibrium_derived",
              derive<2, false>(gas, three_temperature_states, "three-temperature equilibrium"));
        print("three_temperature_modal_derived",
              derive<2, true>(gas, three_temperature_states, "three-temperature"));
    }
    {
        const auto gas = make_controlled_modal(controlled_two_levels, scales);
        print("controlled_modal_equilibrium_derived",
              derive<2, false>(gas, controlled_modal_states, "controlled modal equilibrium"));
        print("controlled_modal_derived",
              derive<2, true>(gas, controlled_modal_states, "controlled modal"));
        const auto ground = make_controlled_modal(controlled_ground_level, scales);
        print("controlled_ground_level_derived",
              derive<2, true>(ground, controlled_ground_level_states, "controlled ground level"));
    }
}

// The rates of progress of mechanism with gas at every state, with the mode
// temperatures of the state as the scalar T_ve (Scalar) or as mode inputs.
template <bool Scalar, typename M, typename G, std::size_t N, std::size_t K, std::size_t C>
std::array<DerivedRates<M::reaction_count>, C>
  derive_rates(const M& mechanism, const G& gas, const std::array<FixtureState<N, K>, C>& states)
{
    std::array<DerivedRates<M::reaction_count>, C> out {};
    for (std::size_t s = 0; s < C; ++s) {
        const FixtureState<N, K>& state = states[s];
        const auto mix =
          mixture(gas,
                  as_kind<kind::mass_fraction>(
                    Kokkos::mdspan<const double, Kokkos::extents<std::size_t, N>>(state.y.data())));
        const Density<double, si> rho {state.rho};
        const Temperature<double, si> temp {state.temp};
        const auto rates = [&] {
            if constexpr (Scalar) {
                return rates_of_progress(mechanism,
                                         gas,
                                         mix,
                                         rho,
                                         temp,
                                         Temperature<double, si> {state.modes[0]});
            } else {
                return rates_of_progress(
                  mechanism,
                  gas,
                  mix,
                  rho,
                  temp,
                  as_kind<kind::mode_temperature, si>(
                    Kokkos::mdspan<const double, Kokkos::extents<std::size_t, K>>(
                      state.modes.data())));
            }
        }();
        for (const ReactionIndex k : reaction_indices(0, M::reaction_count)) {
            out[s].forward[k.value()] = rates.forward()[k].value();
            out[s].backward[k.value()] = rates.backward()[k].value();
        }
    }
    return out;
}

// Writes rates as an array named name.
template <std::size_t R, std::size_t C>
void print_rates(std::string_view name, const std::array<DerivedRates<R>, C>& rates)
{
    std::printf("inline constexpr std::array<DerivedRates<%zu>, %zu> %.*s {{\n",
                R,
                C,
                static_cast<int>(name.size()),
                name.data());
    for (const DerivedRates<R>& r : rates) {
        for (const auto* values : {&r.forward, &r.backward}) {
            std::printf(values == &r.forward ? "  {.forward = {" : ", .backward = {");
            for (std::size_t k = 0; k < R; ++k) {
                if (k > 0) { std::printf(", "); }
                print_value((*values)[k]);
            }
            std::printf("}");
        }
        std::printf("},\n");
    }
    std::printf("}};\n\n");
}

void exchange()
{
    const Pressure<double, si> p0 = standard_pressure(scales);
    {
        const auto gas = make_air5_two_temperature(scales).model;
        print_rates("air5_park_rates_derived",
                    derive_rates<true>(make_air5_park(scales, p0, gas), gas, air5_park_states));
    }
    const auto gas = make_air11_three_temperature(scales).model;
    print("air11_kim_jo_derived", derive<2, true>(gas, air11_kim_jo_states, "air11 Kim-Jo"));
    print_rates("air11_kim_jo_rates_derived",
                derive_rates<false>(make_air11_kim_jo(scales, p0, gas), gas, air11_kim_jo_states));
    print("multi_temperature_process_derived",
          derive<2, true>(gas, multi_temperature_process_states, "multi-temperature process"));
    print("neutral_electron_derived",
          derive<2, true>(gas, neutral_electron_states, "neutral electron-translation"));
    print("charged_electron_derived",
          derive<2, true>(gas, charged_electron_states, "charged electron-translation"));
    print("vibration_electron_derived",
          derive<2, true>(gas, vibration_electron_states, "vibration-electron"));
}

// The structural states of every order of N species, order by order.
template <std::size_t N> void structural(std::string_view name)
{
    const auto orders = structural_orders<N>();
    std::vector<DerivedState<0>> all;
    for (const auto& order : orders) {
        const auto models = make_structural_models<N>(order, scales);
        for (const DerivedState<0>& d :
             derive_from_pressure<0, false>(models.gas, structural_states(order), name)) {
            all.push_back(d);
        }
    }
    std::printf("inline constexpr std::array<DerivedState<0>, %zu> %.*s {{\n",
                all.size(),
                static_cast<int>(name.size()),
                name.data());
    for (const DerivedState<0>& d : all) {
        std::printf("  {.rho = ");
        print_value(d.rho);
        std::printf(", .e = ");
        print_value(d.e);
        std::printf(", .p = ");
        print_value(d.p);
        std::printf(", .mode_e = {}},\n");
    }
    std::printf("}};\n\n");
}

void transport()
{
    print("simple_law_perfect_derived", derive_perfect(simple_law_states));
    print("simple_law_mixture_derived",
          derive<0, false>(make_perfect_mixture<5>(5, scales),
                           simple_law_mixture_states(),
                           "simple-law mixture"));
    print("named_transport_derived",
          derive_from_pressure<0, false>(make_air5(scales).model,
                                         named_transport_states,
                                         "named transport"));
    print("named_transport_two_temperature_derived",
          derive_from_pressure<1, true>(make_air5_two_temperature(scales).model,
                                        named_transport_states,
                                        "named transport two-temperature"));
    structural<1>("structural1_derived");
    structural<2>("structural2_derived");
    structural<16>("structural16_derived");
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::fprintf(stderr, "usage: %s <source description>\n", argv[0]);
        return 2;
    }
    std::printf("#ifndef CANARY_TESTS_DEVICE_DERIVED_INPUTS_HPP\n"
                "#define CANARY_TESTS_DEVICE_DERIVED_INPUTS_HPP\n\n"
                "// Generated by tests/device/derive_inputs.cpp from %s, compiled by\n"
                "// %s; do not edit.\n"
                "// Inputs derived from the declared fixture states by binary64 SI forward\n"
                "// closures, in the order of the fixture arrays (see host_inputs.hpp).\n\n"
                "#include \"host_inputs.hpp\"\n\n"
                "#include <array>\n\n"
                "namespace canary::test\n{\n\n",
                argv[1],
                __VERSION__);
    thermodynamics();
    transport();
    exchange();
    std::printf(
      "}  // namespace canary::test\n\n#endif  // CANARY_TESTS_DEVICE_DERIVED_INPUTS_HPP\n");
    return 0;
}
