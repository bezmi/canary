// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// Values that cerberus's calorically perfect gas produces in its contract
// tests, reproduced with canary in scaled units (k_B = 1). The reference
// values and inputs are copied from cerberus
// Source/states/Eulerian/hydro/tests/test_gas_contract.cpp at commit
// 9d45fc5; each test names the line of the cerberus test it follows.

#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <cstddef>
#include <canary/models/calorically_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/status.hpp>
#include <canary/thermo.hpp>
#include <limits>
#include <mdspan/mdspan.hpp>

namespace
{

using canary::scaled;
using canary::Status;
using Rho = canary::Density<double, scaled>;
using P = canary::Pressure<double, scaled>;
using Temp = canary::Temperature<double, scaled>;
using Mass = canary::ParticleMass<double, scaled>;
using Q = canary::Charge<double, scaled>;
using Gamma = canary::HeatCapacityRatio<double>;
using View = Kokkos::mdspan<const double, Kokkos::dextents<std::size_t, 1>>;

constexpr double eps = std::numeric_limits<double>::epsilon();
constexpr canary::BoltzmannConstant<double, scaled> k_b {1.0};

bool within_ulp(double a, double b, double ulps)
{ return std::abs(a - b) <= ulps * eps * std::abs(b); }

template <std::size_t N> auto fractions(const std::array<double, N>& y)
{ return canary::as_kind<canary::kind::mass_fraction>(View(y.data(), N)); }

// mass={2,4}, charge={1,-1}, gamma={1.4,5/3}
auto two_species()
{
    return canary::make_calorically_perfect_mixture<2>(std::array {Mass {2.0}, Mass {4.0}},
                                                       std::array {Q {1.0}, Q {-1.0}},
                                                       k_b,
                                                       std::array {Gamma {1.4}, Gamma {5.0 / 3.0}});
}

}  // namespace

void test_cerberus_parity()
{
    using namespace boost::ut;
    using canary::thermo_state;

    // Line 244: Y = {0.25, 0.75}, rho = 8, p = 5.
    "two species oracle: mass, charge, cv, cp, gamma, T, rho_e, sound speed"_test = [] {
        const auto gas = two_species();
        const std::array y {0.25, 0.75};
        const auto mix = canary::mixture(gas.model, fractions(y));
        expect(within_ulp(canary::mean_particle_mass(gas.species, fractions(y)).value(), 3.2, 4.0));
        expect(within_ulp(canary::mean_charge(gas.species, fractions(y)).value(), -0.2, 4.0));
        expect(within_ulp(mix.cp.value(), 0.90625, 4.0));
        expect(within_ulp(mix.cv.value(), 0.59375, 4.0));
        const auto c = thermo_state(gas.model, mix, Rho {8.0}, P {5.0});
        expect(c.status == Status::Valid);
        expect(within_ulp(c.state.temp.value(), 2.0, 4.0));
        expect(within_ulp(c.state.e.value() * 8.0, 9.5, 4.0));
        const auto props = canary::properties(gas.model, mix, c.state);
        expect(within_ulp(props.capacity_ratio.value(), 29.0 / 19, 4.0));
        expect(
          within_ulp(props.acoustics.sound_speed.value(), std::sqrt((29.0 / 19) * 5 / 8), 4.0));
        // Partial densities rho Y = {2, 6} (cerberus U[5], U[6]).
        const std::array rho_i {2.0, 6.0};
        const auto from_rho =
          canary::mixture(gas.model,
                          canary::as_kind<canary::kind::density, scaled>(View(rho_i.data(), 2)),
                          Rho {8.0});
        expect(from_rho.cp == mix.cp && from_rho.cv == mix.cv &&
               from_rho.gas_constant == mix.gas_constant);
    };

    // Line 274: rho in {1e-16, 1e-14, 1e-12}, p = rho * 5 / 8.
    "states at small density give T = 2 and the same sound speed"_test = [] {
        const auto gas = two_species();
        const std::array y {0.25, 0.75};
        for (const double rho : {1.0e-16, 1.0e-14, 1.0e-12}) {
            const auto c = thermo_state(gas.model, fractions(y), Rho {rho}, P {rho * 5 / 8});
            expect(c.status == Status::Valid);
            expect(within_ulp(c.state.temp.value(), 2.0, 4.0)) << rho;
            const auto props =
              canary::properties(gas.model, canary::mixture(gas.model, fractions(y)), c.state);
            expect(
              within_ulp(props.acoustics.sound_speed.value(), std::sqrt((29.0 / 19) * 5 / 8), 4.0))
              << rho;
        }
    };

    // Line 298: mass={2,4,8}, charge={1,-1,2}, gamma={1.4,5/3,1.5}.
    "pure endpoints and an unequal three-species mixture match the analytic weights"_test = [] {
        const auto gas = canary::make_calorically_perfect_mixture<3>(
          std::array {Mass {2.0}, Mass {4.0}, Mass {8.0}},
          std::array {Q {1.0}, Q {-1.0}, Q {2.0}},
          k_b,
          std::array {Gamma {1.4}, Gamma {5.0 / 3.0}, Gamma {1.5}});
        const std::array<std::array<double, 3>, 4> mixtures {
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.125, 0.375, 0.5}},
        };
        const std::array expected_mass {2.0, 4.0, 8.0, 32.0 / 7};
        const std::array expected_charge {1.0, -1.0, 2.0, 3.0 / 7};
        const std::array expected_cv {1.25, 0.375, 0.25, 27.0 / 64};
        const std::array expected_cp {1.75, 0.625, 0.375, 41.0 / 64};
        for (std::size_t i = 0; i < mixtures.size(); ++i) {
            const auto y = fractions(mixtures[i]);
            const auto mix = canary::mixture(gas.model, y);
            expect(
              within_ulp(canary::mean_particle_mass(gas.species, y).value(), expected_mass[i], 4.0))
              << i;
            expect(within_ulp(canary::mean_charge(gas.species, y).value(), expected_charge[i], 4.0))
              << i;
            expect(within_ulp(mix.cv.value(), expected_cv[i], 4.0)) << i;
            expect(within_ulp(mix.cp.value(), expected_cp[i], 4.0)) << i;
            expect(within_ulp(canary::heat_capacity_ratio(mix.cp, mix.cv).value(),
                              expected_cp[i] / expected_cv[i],
                              4.0))
              << i;
        }
    };

    // Line 354: equal or unequal masses with equal charge and gamma.
    "equal charges and gammas give the constant values"_test = [] {
        for (const double second_mass : {2.0, 4.0}) {
            const auto gas = canary::make_calorically_perfect_mixture<2>(
              std::array {Mass {2.0}, Mass {second_mass}},
              std::array {Q {-1.0}, Q {-1.0}},
              k_b,
              std::array {Gamma {1.4}, Gamma {1.4}});
            const std::array y {0.25, 0.75};
            const double mass = second_mass == 2.0 ? 2.0 : 3.2;
            const auto mix = canary::mixture(gas.model, fractions(y));
            expect(
              within_ulp(canary::mean_particle_mass(gas.species, fractions(y)).value(), mass, 4.0));
            expect(within_ulp(canary::mean_charge(gas.species, fractions(y)).value(), -1.0, 4.0));
            expect(within_ulp(canary::heat_capacity_ratio(mix.cp, mix.cv).value(), 1.4, 4.0));
            expect(within_ulp(mix.cp.value(), 3.5 / mass, 4.0));
        }
    };

    // Line 771: each pair of density, pressure and temperature closes to
    // rho = 8, p = 5, T = 2, rho_e = 9.5.
    "each pair of density, pressure and temperature gives the same state"_test = [] {
        const auto gas = two_species();
        const std::array y {0.25, 0.75};
        const auto mix = canary::mixture(gas.model, fractions(y));
        for (const auto& c : {
               thermo_state(gas.model, mix, Rho {8.0}, P {5.0}),
               thermo_state(gas.model, mix, Rho {8.0}, Temp {2.0}),
               thermo_state(gas.model, mix, P {5.0}, Temp {2.0}),
             }) {
            expect(within_ulp(c.state.rho.value(), 8.0, 4.0));
            expect(within_ulp(c.state.p.value(), 5.0, 4.0));
            expect(within_ulp(c.state.temp.value(), 2.0, 4.0));
            expect(within_ulp(c.state.rho.value() * c.state.e.value(), 9.5, 4.0));
        }
    };

    // Line 487: mass=2, charge=-1, gamma=1.4, one species.
    "the one-species model has R = 1/m, cv = 1.25 and cp = 1.75"_test = [] {
        const auto gas = canary::make_calorically_perfect(Mass {2.0}, Q {-1.0}, k_b, Gamma {1.4});
        expect(gas.model.gas_constant.value() == 0.5);
        expect(within_ulp(gas.model.cv.value(), 1.25, 4.0));
        expect(within_ulp(gas.model.cp.value(), 1.75, 4.0));
    };
}
