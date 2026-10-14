#include "dual.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/models/calorically_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/status.hpp>
#include <canary/thermo.hpp>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <string>
#include <type_traits>
#include <vector>

namespace
{

using canary::scaled;
using canary::Status;
using D = canary::test::Dual<double>;
using Rho = canary::Density<double, scaled>;
using E = canary::SpecificEnergy<double, scaled>;
using RhoE = canary::EnergyDensity<double, scaled>;
using P = canary::Pressure<double, scaled>;
using Temp = canary::Temperature<double, scaled>;
using Mass = canary::ParticleMass<double, scaled>;
using Q = canary::Charge<double, scaled>;
using Gamma = canary::HeatCapacityRatio<double>;
using Extents = Kokkos::dextents<std::size_t, 1>;
using View = Kokkos::mdspan<const double, Extents>;
using StridedView = Kokkos::mdspan<const double, Extents, Kokkos::layout_stride>;

constexpr double eps = std::numeric_limits<double>::epsilon();
constexpr canary::BoltzmannConstant<double, scaled> k_b {1.0};

bool within_ulp(double a, double b, double ulps)
{ return std::abs(a - b) <= ulps * eps * std::abs(b); }

template <std::size_t N> auto fractions(const std::array<double, N>& y)
{ return canary::as_kind<canary::kind::mass_fraction>(View(y.data(), N)); }

template <std::size_t N> auto partial_densities(const std::array<double, N>& rho_i)
{ return canary::as_kind<canary::kind::density, scaled>(View(rho_i.data(), N)); }

bool same_state(const canary::Closure<double, scaled>& a,
                const canary::Closure<double, scaled>& b,
                double ulps)
{
    return within_ulp(a.state.rho.value(), b.state.rho.value(), ulps) &&
           within_ulp(a.state.e.value(), b.state.e.value(), ulps) &&
           within_ulp(a.state.p.value(), b.state.p.value(), ulps) &&
           within_ulp(a.state.temp.value(), b.state.temp.value(), ulps);
}

bool identical(const canary::Closure<double, scaled>& a, const canary::Closure<double, scaled>& b)
{
    return a.state.rho == b.state.rho && a.state.e == b.state.e && a.state.p == b.state.p &&
           a.state.temp == b.state.temp && a.status == b.status;
}

using Mixture3 = canary::CaloricallyPerfectMixture<double, scaled, 3>;

// Electrons, protons and a heavier ion with different gammas.
auto three_species()
{
    return canary::make_calorically_perfect_mixture<3>(
      std::array {Mass {1.0 / 1836.0}, Mass {1.0}, Mass {4.0}},
      std::array {Q {-1.0}, Q {1.0}, Q {2.0}},
      k_b,
      std::array {Gamma {5.0 / 3.0}, Gamma {5.0 / 3.0}, Gamma {1.4}});
}

}  // namespace

static_assert(canary::MixtureGasModel<Mixture3, double>);
static_assert(canary::MixtureGasModel<Mixture3, D>);
static_assert(canary::MixtureGasModel<Mixture3, float>);
static_assert(canary::MixtureGasModel<canary::CaloricallyPerfectMixture<float, scaled>, double>);
static_assert(std::is_trivially_copyable_v<Mixture3>);
static_assert(std::is_trivially_copyable_v<canary::Mixture<double, scaled>>);
static_assert(std::same_as<canary::CaloricallyPerfectMixture<double, scaled>,
                           canary::CaloricallyPerfectMixture<double, scaled, CANARY_MAX_SPECIES>>);

void test_calorically_perfect_mixture()
{
    using namespace boost::ut;
    using canary::mixture;
    using canary::thermo_state;

    "N identical species give the one-species results"_test = [] {
        const auto one = canary::make_calorically_perfect(Mass {2.0}, Q {1.0}, k_b, Gamma {1.4});
        const std::array<Mass, 5> masses {
          Mass {2.0},
          Mass {2.0},
          Mass {2.0},
          Mass {2.0},
          Mass {2.0},
        };
        const std::array<Q, 5> charges {Q {1.0}, Q {1.0}, Q {1.0}, Q {1.0}, Q {1.0}};
        const std::array<Gamma, 5> gammas {
          Gamma {1.4},
          Gamma {1.4},
          Gamma {1.4},
          Gamma {1.4},
          Gamma {1.4},
        };
        const auto five = canary::make_calorically_perfect_mixture<5>(masses, charges, k_b, gammas);
        const std::array y {0.1, 0.3, 0.05, 0.25, 0.3};
        const auto mix = mixture(five.model, fractions(y));
        expect(mix.status == Status::Valid);
        const Rho rho {3.0};
        const E e {1.7};
        const auto ref = thermo_state(one.model, rho, e);
        const auto c = thermo_state(five.model, mix, rho, e);
        expect(same_state(c, ref, 4.0));
        expect(same_state(thermo_state(five.model, mix, rho, ref.state.p),
                          thermo_state(one.model, rho, ref.state.p),
                          4.0));
        expect(same_state(thermo_state(five.model, mix, rho, ref.state.temp),
                          thermo_state(one.model, rho, ref.state.temp),
                          4.0));
        expect(same_state(thermo_state(five.model, mix, ref.state.p, ref.state.temp),
                          thermo_state(one.model, ref.state.p, ref.state.temp),
                          4.0));
        expect(same_state(thermo_state(five.model, mix, rho, RhoE {3.0 * 1.7}),
                          thermo_state(one.model, rho, RhoE {3.0 * 1.7}),
                          4.0));
        const auto props = canary::properties(five.model, mix, c.state);
        const auto ref_props = canary::properties(one.model, ref.state);
        expect(within_ulp(props.acoustics.sound_speed.value(),
                          ref_props.acoustics.sound_speed.value(),
                          4.0));
        expect(within_ulp(props.capacity_ratio.value(), 1.4, 4.0));
        expect(
          within_ulp(canary::mean_particle_mass(five.species, fractions(y)).value(), 2.0, 4.0));
        expect(within_ulp(canary::mean_charge(five.species, fractions(y)).value(), 1.0, 4.0));
    };

    "mixture coefficients are the weighted means of the species values"_test = [] {
        const auto built = three_species();
        const std::array y {0.2, 0.5, 0.3};
        const auto mix = mixture(built.model, fractions(y));
        double r = 0;
        double cv = 0;
        double cp = 0;
        for (std::size_t i = 0; i < 3; ++i) {
            r += y[i] * built.model.gas_constant[i].value();
            cv += y[i] * built.model.cv[i].value();
            cp += y[i] * built.model.cp[i].value();
        }
        expect(within_ulp(mix.gas_constant.value(), r, 4.0));
        expect(within_ulp(mix.cv.value(), cv, 4.0));
        expect(within_ulp(mix.cp.value(), cp, 4.0));
        // R of the mixture is k_B over the mean particle mass.
        expect(within_ulp(mix.gas_constant.value(),
                          1.0 / canary::mean_particle_mass(built.species, fractions(y)).value(),
                          8.0));
        // cp - cv = R for the mixture as for each species.
        expect(within_ulp(mix.cp.value() - mix.cv.value(), mix.gas_constant.value(), 16.0));
    };

    "mixture results do not change when the species are permuted"_test = [] {
        const std::array<Mass, 3> masses {Mass {1.0 / 1836.0}, Mass {1.0}, Mass {4.0}};
        const std::array<Q, 3> charges {Q {-1.0}, Q {1.0}, Q {2.0}};
        const std::array<Gamma, 3> gammas {Gamma {5.0 / 3.0}, Gamma {5.0 / 3.0}, Gamma {1.4}};
        const std::array<double, 3> y {0.2, 0.5, 0.3};
        const std::array<std::array<std::size_t, 3>, 6> orders {
          {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}},
        };
        const auto reference =
          thermo_state(three_species().model, fractions(y), Rho {3.0}, E {1.7});
        for (const auto& order : orders) {
            std::array<Mass, 3> m {};
            std::array<Q, 3> q {};
            std::array<Gamma, 3> g {};
            std::array<double, 3> yp {};
            for (std::size_t i = 0; i < 3; ++i) {
                m[i] = masses[order[i]];
                q[i] = charges[order[i]];
                g[i] = gammas[order[i]];
                yp[i] = y[order[i]];
            }
            const auto built = canary::make_calorically_perfect_mixture<3>(m, q, k_b, g);
            expect(same_state(thermo_state(built.model, fractions(yp), Rho {3.0}, E {1.7}),
                              reference,
                              4.0));
        }
    };

    "strided and contiguous views give identical results"_test = [] {
        const auto built = three_species();
        const std::array y {0.2, 0.5, 0.3};
        std::vector<double> storage(21, -99.0);
        for (std::size_t i = 0; i < 3; ++i) { storage[7 * i] = y[i]; }
        const StridedView strided(
          storage.data(),
          Kokkos::layout_stride::mapping<Extents>(Extents(3), std::array<std::size_t, 1> {7}));
        const auto ys = canary::as_kind<canary::kind::mass_fraction>(strided);
        expect(identical(thermo_state(built.model, ys, Rho {3.0}, E {1.7}),
                         thermo_state(built.model, fractions(y), Rho {3.0}, E {1.7})));
    };

    "partial densities give the same mixture as the equivalent fractions"_test = [] {
        const auto built = three_species();
        const std::array y {0.25, 0.5, 0.25};
        const std::array rho_i {0.5, 1.0, 0.5};
        const Rho rho {2.0};
        const auto from_y = mixture(built.model, fractions(y));
        const auto from_rho = mixture(built.model, partial_densities(rho_i), rho);
        expect(from_y.gas_constant == from_rho.gas_constant && from_y.cv == from_rho.cv &&
               from_y.cp == from_rho.cp);
        expect(from_rho.status == Status::Valid);
        // The front end finds the density among the state variables.
        expect(identical(thermo_state(built.model, partial_densities(rho_i), E {1.7}, rho),
                         thermo_state(built.model, fractions(y), rho, E {1.7})));
    };

    "a composition view passed to thermo_state equals the mixture value"_test = [] {
        const auto built = three_species();
        const std::array y {0.2, 0.5, 0.3};
        const auto mix = mixture(built.model, fractions(y));
        expect(identical(thermo_state(built.model, fractions(y), E {1.7}, Rho {3.0}),
                         thermo_state(built.model, mix, Rho {3.0}, E {1.7})));
    };

    "the closure status combines the composition and state statuses"_test = [] {
        const auto built = three_species();
        const std::array violated {0.2, 0.5, 0.3 + 1.0e-10};
        const std::array not_finite {std::numeric_limits<double>::quiet_NaN(), 0.5, 0.5};
        expect(thermo_state(built.model, fractions(violated), Rho {3.0}, E {1.7}).status ==
               Status::CompositionViolation);
        expect(thermo_state(built.model, fractions(violated), Rho {3.0}, E {-1.7}).status ==
               Status::Invalid);
        expect(mixture(built.model, fractions(not_finite)).status == Status::Invalid);
    };

    "dual-number states carry the derivative through the mixture closure"_test = [] {
        const auto built = three_species();
        const std::array<D, 3> y {D(0.2), D(0.5), D(0.3)};
        const auto yd = canary::as_kind<canary::kind::mass_fraction>(
          Kokkos::mdspan<const D, Extents>(y.data(), 3));
        // dp/drho at constant e = R T.
        const auto c = thermo_state(built.model,
                                    yd,
                                    canary::Density<D, scaled> {D(3.0, 1.0)},
                                    canary::SpecificEnergy<D, scaled> {D(1.7)});
        const auto mix = mixture(built.model, fractions(std::array {0.2, 0.5, 0.3}));
        expect(within_ulp(c.state.p.value().dual(),
                          mix.gas_constant.value() * c.state.temp.value().real(),
                          4.0));
    };

    "make_calorically_perfect_mixture rejects inconsistent or out-of-range species"_test = [] {
        const std::array<Mass, 2> m {Mass {1.0}, Mass {2.0}};
        const std::array<Q, 2> q {Q {0.0}, Q {0.0}};
        const std::array<Gamma, 2> g {Gamma {1.4}, Gamma {1.4}};
        const std::array<Gamma, 1> g1 {Gamma {1.4}};
        const std::array<Gamma, 2> bad_g {Gamma {1.4}, Gamma {0.9}};
        const std::vector<Mass> none;
        const std::vector<Q> no_q;
        const std::vector<Gamma> no_g;
        expect(throws<canary::ParameterError>(
          [&] { canary::make_calorically_perfect_mixture<2>(m, q, k_b, g1); }));
        expect(throws<canary::ParameterError>(
          [&] { canary::make_calorically_perfect_mixture<2>(none, no_q, k_b, no_g); }));
        expect(throws<canary::ParameterError>(
          [&] { canary::make_calorically_perfect_mixture<1>(m, q, k_b, g); }));
        try {
            canary::make_calorically_perfect_mixture<2>(m, q, k_b, bad_g);
            expect(false) << "no exception";
        } catch (const canary::ParameterError& error) {
            expect(std::string(error.what()).starts_with("species 1: heat capacity ratio"))
              << error.what();
        }
        const std::vector<Mass> mv(m.begin(), m.end());
        const std::vector<Q> qv(q.begin(), q.end());
        const std::vector<Gamma> gv(g.begin(), g.end());
        expect(canary::make_calorically_perfect_mixture(mv, qv, k_b, gv).model.n == 2_i);
    };
}
