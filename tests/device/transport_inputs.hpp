// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_TESTS_DEVICE_TRANSPORT_INPUTS_HPP
#define CANARY_TESTS_DEVICE_TRANSPORT_INPUTS_HPP

// Host fixtures of the transport checks: transport laws and models, and the
// declared states they are evaluated at, in SI units. The densities and
// energies of the owning states are derived offline (derived_inputs.hpp).

#include "host_inputs.hpp"
#include "thermodynamics_inputs.hpp"

#include <algorithm>
#include <array>
#include <canary/data/air5.hpp>
#include <canary/data/air5_transport.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/transport/constant.hpp>
#include <canary/transport/mixture.hpp>
#include <canary/transport/power_law.hpp>
#include <canary/transport/prandtl.hpp>
#include <canary/transport/sutherland.hpp>
#include <canary/transport/transport.hpp>
#include <canary/units.hpp>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

namespace canary::test
{

// Simple laws: constant viscosity 1.8e-5 Pa s; Sutherland and power law with
// mu_0 = 1.716e-5 Pa s at T_0 = 273.15 K, S = 110.4 K and exponent 0.76; a
// constant Prandtl number 0.71.
template <std::floating_point P, UnitSystem U>
ConstantViscosity<P, U> make_constant_law(const UnitScales<P, U>& scales)
{ return make_constant_viscosity(convert<U>(DynamicViscosity<P, si> {P(1.8e-5)}, scales)); }

template <std::floating_point P, UnitSystem U>
Sutherland<P, U> make_sutherland_law(const UnitScales<P, U>& scales)
{
    return make_sutherland(convert<U>(DynamicViscosity<P, si> {P(1.716e-5)}, scales),
                           convert<U>(Temperature<P, si> {P(273.15)}, scales),
                           convert<U>(Temperature<P, si> {P(110.4)}, scales));
}

template <std::floating_point P, UnitSystem U>
PowerLaw<P, U> make_power_law_law(const UnitScales<P, U>& scales)
{
    return make_power_law(convert<U>(DynamicViscosity<P, si> {P(1.716e-5)}, scales),
                          convert<U>(Temperature<P, si> {P(273.15)}, scales),
                          P(0.76));
}

template <std::floating_point P> ConstantPrandtl<P> make_prandtl()
{ return make_constant_prandtl(PrandtlNumber<P> {P(0.71)}); }

// The simple-law states: rho = 1 kg/m^3 at 273.15, 300 and 3000 K, for the
// perfect gas and the five-species perfect mixture.
inline constexpr std::array<FixtureState<1, 0>, 3> simple_law_states {{
  {.y = {1.0}, .rho = 1.0, .temp = 273.15, .modes = {}},
  {.y = {1.0}, .rho = 1.0, .temp = 300.0, .modes = {}},
  {.y = {1.0}, .rho = 1.0, .temp = 3000.0, .modes = {}},
}};

inline std::array<FixtureState<5, 0>, 3> simple_law_mixture_states()
{
    const auto y = perfect_mixture_fractions<5>();
    return {{{.y = y, .rho = 1.0, .temp = 273.15, .modes = {}},
             {.y = y, .rho = 1.0, .temp = 300.0, .modes = {}},
             {.y = y, .rho = 1.0, .temp = 3000.0, .modes = {}}}};
}

// A state declared by pressure [Pa] instead of density.
template <std::size_t N, std::size_t M> struct PressureState {
    std::array<double, N> y;
    double p;
    double temp;
    std::array<double, M> modes;
};

// Named air5 states: air, mixed, pure N2 and a trace of O2, at p = 10000 Pa and
// T = 300, 3000 and 10000 K, with T_ve = 0.5 T + 150 K for the two-temperature
// model.
inline constexpr std::array<double, 5> air5_pure_nitrogen {1.0, 0.0, 0.0, 0.0, 0.0};
inline constexpr std::array<double, 5> air5_trace_oxygen {1.0 - 0x1p-20, 0x1p-20, 0.0, 0.0, 0.0};
inline constexpr std::array<PressureState<5, 1>, 12> named_transport_states = [] {
    std::array<PressureState<5, 1>, 12> out {};
    std::size_t k = 0;
    for (const auto& y : {air5_air, air5_mixed, air5_pure_nitrogen, air5_trace_oxygen}) {
        for (const double temp : {300.0, 3000.0, 10000.0}) {
            out[k++] = {.y = y, .p = 10000.0, .temp = temp, .modes = {(0.5 * temp) + 150.0}};
        }
    }
    return out;
}();

// Structural transport models of N species: species i repeats the air5
// species, transport and rotation records of species order[i] mod 5 under the
// name "species<order[i]>", with the collision pairs of the matching air5
// species renamed alike. The orders are every permutation for N <= 5, and for
// larger N the N cyclic shifts and their reversals.
template <std::size_t N> std::vector<std::array<std::size_t, N>> structural_orders()
{
    std::array<std::size_t, N> order {};
    std::iota(order.begin(), order.end(), std::size_t {0});
    std::vector<std::array<std::size_t, N>> out;
    if constexpr (N <= 5) {
        do {
            out.push_back(order);
        } while (std::next_permutation(order.begin(), order.end()));
    } else {
        for (std::size_t shift = 0; shift < N; ++shift) {
            for (std::size_t i = 0; i < N; ++i) { order[i] = (i + shift) % N; }
            out.push_back(order);
            std::ranges::reverse(order);
            out.push_back(order);
        }
    }
    return out;
}

// The owned names and the renamed inputs of one structural order.
template <std::size_t N> struct StructuralInputs {
    std::array<std::string, N> owned_names;
    std::array<std::string_view, N> names;
    std::array<Nasa9Species<4>, N> species;
    std::array<SpeciesTransportInput<3>, N> transport;
    std::array<int, N> rotation;
    std::array<CollisionPairInput, (N * (N + 1)) / 2> pairs;
};

// Fills inputs for order; names in the inputs point into inputs.owned_names,
// which must not move while the inputs are used.
template <std::size_t N>
void fill_structural_inputs(StructuralInputs<N>& inputs, const std::array<std::size_t, N>& order)
{
    for (std::size_t i = 0; i < N; ++i) {
        inputs.owned_names[i] = "species" + std::to_string(order[i]);
        inputs.names[i] = inputs.owned_names[i];
        inputs.species[i] = air5_species[order[i] % 5];
        inputs.species[i].name = inputs.owned_names[i].c_str();
        inputs.transport[i] = air5_species_transport[order[i] % 5];
        inputs.transport[i].name = inputs.owned_names[i].c_str();
        inputs.rotation[i] = air5_rotational_degrees[order[i] % 5];
    }
    std::size_t k = 0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = i; j < N; ++j) {
            const std::string_view a = air5_species[order[i] % 5].name;
            const std::string_view b = air5_species[order[j] % 5].name;
            const auto pair = std::ranges::find_if(air5_collision_integrals, [&](const auto& p) {
                return (p.first == a && p.second == b) || (p.first == b && p.second == a);
            });
            if (pair == air5_collision_integrals.end()) { std::abort(); }
            inputs.pairs[k] = *pair;
            inputs.pairs[k].first = inputs.owned_names[i].c_str();
            inputs.pairs[k].second = inputs.owned_names[j].c_str();
            ++k;
        }
    }
}

// The thermally perfect gas and the Wilke and Gupta-Yos models of one order.
template <std::size_t N, std::floating_point P, UnitSystem U> struct StructuralModels {
    ThermallyPerfectMixture<P, U, N, 4> gas;
    WilkeTransport<P, U, N, 3> wilke;
    GuptaYosTransport<P, U, N> gupta_yos;
};

template <std::size_t N, std::floating_point P, UnitSystem U>
StructuralModels<N, P, U> make_structural_models(const std::array<std::size_t, N>& order,
                                                 const UnitScales<P, U>& scales)
{
    StructuralInputs<N> inputs {};
    fill_structural_inputs(inputs, order);
    const auto built = make_thermally_perfect_mixture(inputs.species, scales);
    return {.gas = built.model,
            .wilke = make_wilke_transport(inputs.transport, inputs.names, built.species, scales),
            .gupta_yos = make_gupta_yos_transport(inputs.pairs,
                                                  inputs.names,
                                                  built.species,
                                                  inputs.rotation,
                                                  scales)};
}

// The structural states of one order: normalised weights 1 + id, pure id 0, and
// id 0 with a trace 2^-20 of id 1 (id = order[i]), each at p = 10000 Pa and
// T = 300 and 10000 K; weights normalised in binary64.
template <std::size_t N>
std::array<PressureState<N, 0>, 6> structural_states(const std::array<std::size_t, N>& order)
{
    std::array<PressureState<N, 0>, 6> out {};
    std::size_t k = 0;
    for (int composition = 0; composition < 3; ++composition) {
        std::array<double, N> y {};
        double total = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            const std::size_t id = order[i];
            y[i] = composition == 0   ? 1.0 + static_cast<double>(id)
                   : composition == 1 ? (id == 0 ? 1.0 : 0.0)
                                      : (id == 0 ? 1.0 : (id == 1 ? 0x1p-20 : 0.0));
            total += y[i];
        }
        for (double& weight : y) { weight /= total; }
        for (const double temp : {300.0, 10000.0}) {
            out[k++] = {.y = y, .p = 10000.0, .temp = temp, .modes = {}};
        }
    }
    return out;
}

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_TRANSPORT_INPUTS_HPP
