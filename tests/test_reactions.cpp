#include "dual.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/data/air5.hpp>
#include <canary/data/air5_kim_jo.hpp>
#include <canary/models/calorically_perfect.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/species.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <random>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace
{

using canary::scaled;
using canary::si;
using D = canary::test::Dual<double>;
using Air = canary::Air5<double, si>;
using AirScaled = canary::Air5<double, scaled>;
using Mechanism = canary::Air5KimJo<double, si>;
using Extents = Kokkos::dextents<std::size_t, 1>;

constexpr std::size_t species_count = 5;
constexpr std::size_t reaction_count = 13;
using Fractions = std::array<double, species_count>;
using PerReaction = std::array<double, reaction_count>;

constexpr double atmosphere = 101325.0;
constexpr double avogadro = canary::constants::avogadro<double>;
// Universal gas constant of the reference implementation [J/(mol K)].
constexpr double reference_gas_constant = 8.31451;

const canary::UnitScales<double, scaled> reference_scales {
  .mass = 4.65e-26,
  .density = 1.0e-3,
  .temperature = 1000.0,
  .length = 1.0,
};

template <typename T> auto fractions(const std::array<T, species_count>& y)
{
    return canary::as_kind<canary::kind::mass_fraction>(
      Kokkos::mdspan<const T, Extents>(y.data(), species_count));
}

bool relative(double a, double b, double bound) { return std::abs(a - b) <= bound * std::abs(b); }

// A state of the reference implementation of the mechanism: its rates of
// progress in mol/(m^3 s) and its mass production rates in kg/(m^3 s),
// computed with R_u = 8.31451 J/(mol K) and the standard pressure 1 atm.
struct ReferenceState {
    Fractions y;
    double temp;
    double rho;
    PerReaction forward;
    PerReaction backward;
    Fractions source;
};

// clang-format off
const std::array<ReferenceState, 5> reference_states {
 {
  {
    .y = {0.767, 0.233, 0.0, 0.0, 0.0}, .temp = 5000.0, .rho = 0.006939829108798786,
    .forward = {0.02086957258136873, 0.01193076483970515, 0.0, 0.0, 0.0, 21.452051125326943, 53.11942775274426, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    .backward = {},
    .source = {-0.0009188489723115111, -2.386197838323625, 0.0009188489723115111, 2.386197838323625, 0.0},
  },
  {
    .y = {0.6, 0.1, 0.1, 0.15, 0.05}, .temp = 8000.0, .rho = 0.007037619261282653,
    .forward = {36.11917338782714, 9.449099298657492, 5.0383061883730536, 31.919883547238125, 121.48841955416773, 564.9614172603348, 258.8716541637941, 43.95353546523503, 9092.468950324308, 3582.0073365831445, 362.6249273227603, 571.4791326314848, 6547.884340817446},
    .backward = {0.0465333278431809, 0.012173535389804535, 0.006490988902777724, 0.04112326685528125, 0.15651688358320112, 0.014252179409990634, 0.006530508361429085, 0.001108807883185507, 0.22937406838895197, 0.09036265070356325, 0.34483357978179385, 21542.18651615003, 8871.057891464243},
    .source = {59.372202435036435, -1104.3627016168407, -315.4882659828913, 811.8096129430239, 548.6691522216718},
  },
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .temp = 12000.0, .rho = 0.0009537284289470741,
    .forward = {11.3324353641388, 5.070311373421975, 2.703514944892124, 59.79775804381441, 108.64952943047089, 55.42078097247583, 19.88725849679765, 8.623382721593023, 3253.5127831374225, 712.0732261767776, 165.10652060828053, 39.337479941326585, 330.90062474954647},
    .backward = {0.00028215932058449665, 0.00012624255654735965, 6.731315163329544e-05, 0.0014888675064065548, 0.0027051976403005633, 5.7203579795290745e-05, 2.0527000131884224e-05, 8.90078329760392e-06, 0.0033581731408229513, 0.0007349794950367644, 0.004831933515966059, 1115.3539068743562, 281.5217883242919},
    .source = {-6.637150873842631, -164.01079908131675, -6.813393138677023, 148.6466780625252, 28.81466503131116},
  },
  // Below the lower temperature limit: rates at 300 K.
  {
    .y = {0.6, 0.1, 0.1, 0.15, 0.05}, .temp = 250.0, .rho = 1.1260190818052245,
    .forward = {9.302211786680948e-151, 8.642965192338661e-151, 4.608471520734223e-151, 8.551087362336059e-151, 1.1112383818721133e-149, 4.405252399611495e-76, 4.0035624018478834e-73, 3.4272502805339773e-77, 3.98903507851722e-73, 1.5714931769578992e-73, 1.5378558164223767e-98, 1.1597380101060866e-20, 1.397106543456571e-44},
    .backward = {24152302.333758928, 22440631.666159894, 11965455.100114929, 22202079.67671047, 288522407.1363414, 186304.83039025188, 169316748.86525002, 14494.363189635562, 168702366.24567926, 66460838.842832245, 4497691.490702478, 8020.115369279894, 1240303008.3030038},
    .source = {45089973.208899, 12949041.84974347, -27780531.21325408, 6822973.419029424, -37081457.26441781},
  },
  // Above the upper temperature limit: rates at 30000 K.
  {
    .y = {0.3, 0.1, 0.3, 0.25, 0.05}, .temp = 40000.0, .rho = 0.00028611852868412217,
    .forward = {96.28826408351549, 30.24683197649958, 16.12775947306739, 502.52737869697563, 648.1463994964196, 76.11369194089052, 3.2822268167112036, 11.843165766388294, 1451.145961431145, 317.60200598112453, 638.195866106026, 9.362685332361984, 81.2541390796447},
    .backward = {2.371158872684948e-07, 7.44847201207008e-08, 3.9715618860802984e-08, 1.2375051768831707e-06, 1.5961011454435885e-06, 2.1313739978806003e-07, 9.191057106673343e-09, 3.316383021687169e-08, 4.063572125401226e-06, 8.893651588318524e-07, 1.6100815331684422e-05, 84.35239825157483, 7.931200106984255},
    .source = {-38.28478118462953, -61.91693433659379, 45.14642817564361, 69.75477157648831, -14.699484230908613},
  },
 },
};
// clang-format on

// Change of the coefficient of species i in reaction k.
int coefficient_change(const Mechanism& mech, std::size_t k, std::size_t i)
{ return mech.reactions[k].products[i] - mech.reactions[k].reactants[i]; }

struct Rates {
    PerReaction forward;
    PerReaction backward;
    Fractions production;
};

template <typename M, typename G>
Rates evaluate(const M& mech, const G& gas, const Fractions& y, double rho, double temp)
{
    using U = M::units;
    const auto mix = canary::mixture(gas, fractions(y));
    const canary::Density<double, U> rho_q {rho};
    const canary::Temperature<double, U> temp_q {temp};
    const auto q = canary::rates_of_progress(mech, gas, mix, rho_q, temp_q);
    const auto omega = canary::production_rates(mech, gas, mix, rho_q, temp_q);
    Rates out {};
    for (std::size_t k = 0; k < reaction_count; ++k) {
        out.forward[k] = q.forward[k].value();
        out.backward[k] = q.backward[k].value();
    }
    for (std::size_t i = 0; i < species_count; ++i) { out.production[i] = omega[i].value(); }
    return out;
}

// Equilibrium constant of reaction k per particle in SI, from the
// definition: exp(-sum dnu_i g_i m_i / (k_B T)) (p0 / (k_B T))^sum(dnu).
double
  equilibrium_constant(const Mechanism& mech, const Air& gas, std::size_t k, double temp, double p0)
{
    const double kb = canary::constants::boltzmann<double>;
    const auto sp = canary::species_properties(gas, canary::Temperature<double, si> {temp});
    double exponent = 0;
    int change = 0;
    for (std::size_t i = 0; i < species_count; ++i) {
        const int dnu = coefficient_change(mech, k, i);
        exponent -= dnu * sp.gibbs_energy[i].value() * mech.mass[i].value() / (kb * temp);
        change += dnu;
    }
    return std::exp(exponent) * std::pow(p0 / (kb * temp), change);
}

template <typename F> std::string error_message(const F& f)
{
    try {
        f();
    } catch (const canary::ParameterError& error) {
        return error.what();
    }
    return "no exception";
}

}  // namespace

static_assert(canary::HasSpeciesProperties<Air, double, 5>);
static_assert(canary::HasSpeciesProperties<Air, float, 5>);
static_assert(canary::HasSpeciesProperties<Air, D, 5>);
static_assert(canary::HasSpeciesProperties<AirScaled, D, 5>);
static_assert(!canary::HasSpeciesProperties<Air, double, 4>);
static_assert(!canary::HasSpeciesProperties<canary::CaloricallyPerfect<double, si>, double, 1>);
static_assert(canary::ReactingModel<Air, double, double, 5>);
static_assert(canary::ReactingModel<Air, float, double, 5>);
static_assert(canary::ReactingModel<Air, D, double, 5>);
static_assert(std::is_trivially_copyable_v<Mechanism>);
static_assert(std::is_trivially_copyable_v<canary::Air5KimJo<float, scaled>>);
static_assert(std::is_trivially_copyable_v<canary::RatesOfProgress<D, si, 13>>);
static_assert(canary::air5_kim_jo.reactions.size() == 13);

void test_reactions()
{
    using namespace boost::ut;
    const canary::UnitScales<double, si> si_scales {};
    const auto built = canary::make_air5<double>(si_scales);
    const Air& gas = built.model;
    const Mechanism mech =
      canary::make_air5_kim_jo(built.species, si_scales, canary::Pressure<double, si> {atmosphere});

    "the builder converts rate constants to per-particle units"_test = [&] {
        // N2 + N2 <=> N + N + N2: order 2. NO + M <=> N + O + M: order 2.
        expect(relative(mech.reactions[0].forward.A, 1.216e14 / avogadro, 1.0e-15));
        expect(mech.reactions[0].reactants[0] == 2_i && mech.reactions[0].products[2] == 2_i);
        expect(mech.reactions[0].coefficient_change == 1_i);
        expect(mech.reactions[10].third_body);
        expect(relative(mech.reactions[10].forward.A, 1.45e9 / avogadro, 1.0e-15));
        expect(mech.reactions[10].efficiencies[0] == 1.0_d);
        expect(mech.reactions[10].efficiencies[4] == 0.664827586_d);
        expect(mech.reactions[11].coefficient_change == 0_i);
        // The translational exponent defaults to 1.
        expect(mech.reactions[0].forward.translational_exponent == 1.0_d);
        expect(mech.temperature_lower.value() == 300.0_d);
        expect(mech.temperature_upper.value() == 30000.0_d);
        expect(mech.standard_pressure.value() == atmosphere);
    };

    "rates of progress and production rates equal the reference at five states"_test = [&] {
        // The reference values come from GDTk (gdtk-uq/gdtk commit bf559c3f,
        // libgas.so through its C interface), whose mechanism preparation
        // writes the third-body efficiency of NO, N and O in
        // NO + M <=> N + O + M to 7 significant digits (0.6648276); the file
        // keeps the source value 0.664827586.
        auto input = canary::air5_kim_jo;
        for (std::size_t e = 2; e < 5; ++e) {
            input.reactions[10].efficiencies.entries[e].value = 0.6648276;
        }
        const Mechanism reference_mech = canary::make_reaction_mechanism(
          input,
          std::array<std::string_view, 5> {"N2", "O2", "N", "O", "NO"},
          built.species,
          si_scales,
          canary::Pressure<double, si> {atmosphere});
        const double r = canary::constants::universal_gas_constant<double> / reference_gas_constant;
        for (const ReferenceState& s : reference_states) {
            const Rates got = evaluate(reference_mech, gas, s.y, s.rho, s.temp);
            PerReaction forward {};
            PerReaction backward {};
            for (std::size_t k = 0; k < reaction_count; ++k) {
                forward[k] = avogadro * s.forward[k];
                backward[k] =
                  avogadro * s.backward[k] * std::pow(r, mech.reactions[k].coefficient_change);
                expect(relative(got.forward[k], forward[k], 1.0e-12))
                  << s.temp << k << got.forward[k] << forward[k];
                expect(relative(got.backward[k], backward[k], 1.0e-10))
                  << s.temp << k << got.backward[k] << backward[k];
            }
            for (std::size_t i = 0; i < species_count; ++i) {
                double expected = 0;
                double scale = 0;
                for (std::size_t k = 0; k < reaction_count; ++k) {
                    const double term = built.species.mass[i].value() *
                                        coefficient_change(mech, k, i) * (forward[k] - backward[k]);
                    expected += term;
                    scale += std::abs(term);
                }
                expect(std::abs(got.production[i] - expected) <= 1.0e-10 * scale)
                  << s.temp << i << got.production[i] << expected;
            }
        }
    };

    "production rates conserve mass"_test = [&] {
        const auto check = [&](const Fractions& y, double rho, double temp) {
            const Rates got = evaluate(mech, gas, y, rho, temp);
            double sum = 0;
            double scale = 0;
            for (const double w : got.production) {
                sum += w;
                scale += std::abs(w);
            }
            expect(std::abs(sum) <= 1.0e-13 * scale) << temp << sum << scale;
        };
        for (const ReferenceState& s : reference_states) { check(s.y, s.rho, s.temp); }
        // A fixed seed gives the same compositions in every run.
        std::mt19937 engine(12345);  // NOLINT(bugprone-random-generator-seed)
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        for (int n = 0; n < 200; ++n) {
            Fractions y {};
            double total = 0;
            for (double& v : y) {
                v = unit(engine);
                total += v;
            }
            for (double& v : y) { v /= total; }
            const double rho = std::pow(10.0, -4.0 + (4.0 * unit(engine)));
            const double temp = 300.0 + (29700.0 * unit(engine));
            check(y, rho, temp);
        }
    };

    "forward and backward rates balance at the equilibrium composition of each reaction"_test =
      [&] {
          for (const double temp : {3000.0, 8000.0, 15000.0}) {
              for (std::size_t k = 0; k < reaction_count; ++k) {
                  const double kc = equilibrium_constant(mech, gas, k, temp, atmosphere);
                  // Every species at n_base, except the first one whose
                  // coefficient increases, which is set so that
                  // prod n_i^dnu_i = K_c.
                  const double n_base = 1.0e22;
                  std::size_t j = species_count;
                  double others = 1.0;
                  for (std::size_t i = 0; i < species_count; ++i) {
                      const int dnu = coefficient_change(mech, k, i);
                      if (dnu > 0 && j == species_count) {
                          j = i;
                      } else {
                          others *= std::pow(n_base, dnu);
                      }
                  }
                  std::array<double, species_count> n {};
                  n.fill(n_base);
                  n[j] = std::pow(kc / others, 1.0 / coefficient_change(mech, k, j));
                  double rho = 0;
                  for (std::size_t i = 0; i < species_count; ++i) {
                      rho += n[i] * built.species.mass[i].value();
                  }
                  Fractions y {};
                  for (std::size_t i = 0; i < species_count; ++i) {
                      y[i] = n[i] * built.species.mass[i].value() / rho;
                  }
                  const Rates got = evaluate(mech, gas, y, rho, temp);
                  expect(relative(got.backward[k], got.forward[k], 1.0e-12))
                    << temp << k << got.forward[k] << got.backward[k];
              }
          }
      };

    "the default standard pressure is 1 bar"_test = [&] {
        const Mechanism bar = canary::make_air5_kim_jo(built.species, si_scales);
        expect(bar.standard_pressure.value() == 1.0e5);
        const ReferenceState& s = reference_states[1];
        const Rates a = evaluate(bar, gas, s.y, s.rho, s.temp);
        const Rates b = evaluate(mech, gas, s.y, s.rho, s.temp);
        for (std::size_t k = 0; k < reaction_count; ++k) {
            expect(a.forward[k] == b.forward[k]);
            // K_c(1 bar) / K_c(1 atm) = (q_f / q_b)(1 bar) / (q_f / q_b)(1 atm).
            const double ratio = b.backward[k] / a.backward[k];
            const double expected =
              std::pow(1.0e5 / atmosphere, bar.reactions[k].coefficient_change);
            expect(relative(ratio, expected, 1.0e-13)) << k << ratio << expected;
        }
    };

    "rates outside the temperature limits are those at the limits"_test = [&] {
        const canary::Density<double, si> rho {0.01};
        for (const auto& [outside, limit] :
             {std::array {250.0, 300.0}, std::array {4.0e4, 3.0e4}}) {
            const Rates a = evaluate(mech, gas, reference_states[1].y, rho.value(), outside);
            const Rates b = evaluate(mech, gas, reference_states[1].y, rho.value(), limit);
            expect(a.forward == b.forward && a.backward == b.backward &&
                   a.production == b.production)
              << outside;
            const std::array<D, species_count> y_d {0.6, 0.1, 0.1, 0.15, 0.05};
            const auto mix_d = canary::mixture(gas, fractions(y_d));
            const canary::Density<D, si> rho_d {D(rho.value())};
            const auto q = canary::rates_of_progress(mech,
                                                     gas,
                                                     mix_d,
                                                     rho_d,
                                                     canary::Temperature<D, si> {D(outside, 1.0)});
            const auto omega =
              canary::production_rates(mech,
                                       gas,
                                       mix_d,
                                       rho_d,
                                       canary::Temperature<D, si> {D(outside, 1.0)});
            for (std::size_t k = 0; k < reaction_count; ++k) {
                expect(q.forward[k].value().dual() == 0.0 && q.backward[k].value().dual() == 0.0);
                expect(q.forward[k].value().real() == a.forward[k]);
            }
            for (std::size_t i = 0; i < species_count; ++i) {
                expect(omega[i].value().dual() == 0.0);
            }
        }
    };

    "dual-number derivatives equal central differences"_test = [&] {
        for (const std::size_t state : {std::size_t {1}, std::size_t {2}}) {
            const ReferenceState& s = reference_states[state];
            Fractions rho_i {};
            for (std::size_t j = 0; j < species_count; ++j) { rho_i[j] = s.rho * s.y[j]; }
            // omega from partial densities and T, with value type T.
            const auto omega = [&]<typename T>(const std::array<T, species_count>& partial,
                                               T temp) {
                T rho = partial[0];
                for (std::size_t j = 1; j < species_count; ++j) { rho = rho + partial[j]; }
                const auto view = canary::as_kind<canary::kind::density, si>(
                  Kokkos::mdspan<const T, Extents>(partial.data(), species_count));
                const canary::Density<T, si> rho_q {rho};
                const auto mix = canary::mixture(gas, view, rho_q);
                return canary::production_rates(mech,
                                                gas,
                                                mix,
                                                rho_q,
                                                canary::Temperature<T, si> {temp});
            };
            // Columns 0 to 4: d/d(rho_j); column 5: d/dT.
            for (std::size_t c = 0; c <= species_count; ++c) {
                std::array<D, species_count> partial_d {};
                for (std::size_t j = 0; j < species_count; ++j) {
                    partial_d[j] = D(rho_i[j], j == c ? 1.0 : 0.0);
                }
                const auto ad = omega(partial_d, D(s.temp, c == species_count ? 1.0 : 0.0));
                const double h = 1.0e-5 * (c == species_count ? s.temp : rho_i[c]);
                Fractions up = rho_i;
                Fractions down = rho_i;
                double t_up = s.temp;
                double t_down = s.temp;
                if (c == species_count) {
                    t_up += h;
                    t_down -= h;
                } else {
                    up[c] += h;
                    down[c] -= h;
                }
                const auto plus = omega(up, t_up);
                const auto minus = omega(down, t_down);
                for (std::size_t i = 0; i < species_count; ++i) {
                    const double fd = (plus[i].value() - minus[i].value()) / (2.0 * h);
                    const double derivative = ad[i].value().dual();
                    expect(relative(fd, derivative, 1.0e-6)) << state << c << i << derivative << fd;
                }
            }
        }
    };

    "a scaled mechanism gives the SI rates after conversion"_test = [&] {
        const auto built_s = canary::make_air5<double>(reference_scales);
        const auto mech_s = canary::make_air5_kim_jo(
          built_s.species,
          reference_scales,
          canary::convert<scaled>(canary::Pressure<double, si> {atmosphere}, reference_scales));
        for (const ReferenceState& s : reference_states) {
            const auto mix = canary::mixture(gas, fractions(s.y));
            const auto mix_s = canary::mixture(built_s.model, fractions(s.y));
            const canary::Density<double, si> rho {s.rho};
            const canary::Temperature<double, si> temp {s.temp};
            const auto rho_s = canary::convert<scaled>(rho, reference_scales);
            const auto temp_s = canary::convert<scaled>(temp, reference_scales);
            const auto q = canary::rates_of_progress(mech, gas, mix, rho, temp);
            const auto q_s = canary::rates_of_progress(mech_s, built_s.model, mix_s, rho_s, temp_s);
            // The Gibbs energies of the scaled thermo model differ from the
            // SI ones by up to 6e-13 relative at 30000 K; in K_c they are
            // multiplied by m / (k_B T), about 40.
            for (std::size_t k = 0; k < reaction_count; ++k) {
                const double f = canary::convert<si>(q_s.forward[k], reference_scales).value();
                const double b = canary::convert<si>(q_s.backward[k], reference_scales).value();
                expect(relative(f, q.forward[k].value(), 1.0e-12)) << s.temp << k;
                expect(relative(b, q.backward[k].value(), 1.0e-10)) << s.temp << k;
            }
            const auto omega = canary::production_rates(mech, gas, mix, rho, temp);
            const auto omega_s =
              canary::production_rates(mech_s, built_s.model, mix_s, rho_s, temp_s);
            double scale = 0;
            for (const auto w : omega) { scale = std::max(scale, std::abs(w.value())); }
            for (std::size_t i = 0; i < species_count; ++i) {
                const double w = canary::convert<si>(omega_s[i], reference_scales).value();
                expect(std::abs(w - omega[i].value()) <= 1.0e-12 * scale)
                  << s.temp << i << w << omega[i].value();
            }
        }
    };

    "float states give the double rates within float precision"_test = [&] {
        for (const std::size_t state : {std::size_t {1}, std::size_t {2}}) {
            const ReferenceState& s = reference_states[state];
            std::array<float, species_count> y {};
            for (std::size_t i = 0; i < species_count; ++i) { y[i] = static_cast<float>(s.y[i]); }
            const auto mix = canary::mixture(gas, fractions(y));
            const auto omega =
              canary::production_rates(mech,
                                       gas,
                                       mix,
                                       canary::Density<float, si> {static_cast<float>(s.rho)},
                                       canary::Temperature<float, si> {static_cast<float>(s.temp)});
            const Rates ref = evaluate(mech, gas, s.y, s.rho, s.temp);
            double scale = 0;
            for (const double w : ref.production) { scale = std::max(scale, std::abs(w)); }
            for (std::size_t i = 0; i < species_count; ++i) {
                expect(std::abs(omega[i].value() - ref.production[i]) <= 1.0e-4 * scale)
                  << state << i << omega[i].value() << ref.production[i];
            }
        }
    };

    "a float-parameter scaled mechanism gives the SI rates within float precision"_test = [&] {
        const canary::UnitScales<float, scaled> scales_f {
          .mass = 4.65e-26F,
          .density = 1.0e-3F,
          .temperature = 1000.0F,
          .length = 1.0F,
        };
        const auto built_f = canary::make_air5<float>(scales_f);
        const auto mech_f = canary::make_air5_kim_jo(
          built_f.species,
          scales_f,
          canary::convert<scaled>(canary::Pressure<float, si> {101325.0F}, scales_f));
        for (const std::size_t state : {std::size_t {1}, std::size_t {2}}) {
            const ReferenceState& s = reference_states[state];
            std::array<float, species_count> y {};
            for (std::size_t i = 0; i < species_count; ++i) { y[i] = static_cast<float>(s.y[i]); }
            const auto omega = canary::production_rates(
              mech_f,
              built_f.model,
              canary::mixture(built_f.model, fractions(y)),
              canary::convert<scaled>(canary::Density<float, si> {static_cast<float>(s.rho)},
                                      scales_f),
              canary::convert<scaled>(canary::Temperature<float, si> {static_cast<float>(s.temp)},
                                      scales_f));
            const Rates ref = evaluate(mech, gas, s.y, s.rho, s.temp);
            double scale = 0;
            for (const double w : ref.production) { scale = std::max(scale, std::abs(w)); }
            for (std::size_t i = 0; i < species_count; ++i) {
                const double w = canary::convert<si>(omega[i], scales_f).value();
                expect(std::abs(w - ref.production[i]) <= 1.0e-3 * scale)
                  << state << i << w << ref.production[i];
            }
        }
    };

    "reaction orders and dimension exponents are checked before narrowing"_test = [&] {
        const std::intmax_t maximum = std::numeric_limits<int>::max();
        expect(canary::detail::reaction_order(std::array {std::numeric_limits<int>::max(), 1},
                                              true,
                                              "reactant") == maximum + 2);
        const auto limit = maximum / 3 + 1;
        const auto dimension = canary::detail::rate_dimension(limit, "forward");
        expect(dimension.length == 3 * (limit - 1));
        expect(dimension.count == -(limit - 1));
        for (const auto order : {std::intmax_t(0), limit + 1, maximum}) {
            expect(error_message([&] {
                       (void)canary::detail::rate_dimension(order, "forward");
                   }).starts_with("forward rate: reaction order"));
        }
        const auto build = [&](const auto& input) {
            (void)canary::make_reaction_mechanism(
              input,
              std::array<std::string_view, 5> {"N2", "O2", "N", "O", "NO"},
              built.species,
              si_scales);
        };
        auto input = canary::air5_kim_jo;
        auto& reaction = input.reactions[0];
        reaction.reactants = {.count = 2, .entries = {{{"N2", double(maximum)}, {"O2", 1.0}}}};
        reaction.products = reaction.reactants;
        expect(error_message([&] { build(input); }).find("forward rate: reaction order") !=
               std::string::npos);
        reaction.third_body = true;
        expect(error_message([&] { build(input); }).find("forward rate: reaction order") !=
               std::string::npos);
        reaction.third_body = false;
        reaction.reactants = {.count = 1, .entries = {{{"N2", 1.0}}}};
        reaction.products = {.count = 2,
                             .entries = {{{"N", double(maximum)}, {"O", double(maximum)}}}};
        expect(error_message([&] { build(input); }).find("net stoichiometric coefficient") !=
               std::string::npos);
        std::swap(reaction.reactants, reaction.products);
        expect(error_message([&] { build(input); }).find("net stoichiometric coefficient") !=
               std::string::npos);
    };

    "the builder rejects invalid input with the reaction in the message"_test = [&] {
        const auto build = [&](const canary::MechanismInput<13, 5>& input) {
            (void)canary::make_reaction_mechanism(
              input,
              std::array<std::string_view, 5> {"N2", "O2", "N", "O", "NO"},
              built.species,
              si_scales);
        };
        const auto starts = [](const std::string& message, std::string_view prefix) {
            return message.starts_with(prefix);
        };
        auto input = canary::air5_kim_jo;
        input.reactions[3].reactants.entries[1].name = "Ar";
        std::string message = error_message([&] { build(input); });
        expect(starts(message, "reaction 3 (N2 + N <=> N + N + N): unknown species \"Ar\""))
          << message;

        input = canary::air5_kim_jo;
        input.reactions[5].forward.A = -1.0;
        message = error_message([&] { build(input); });
        expect(starts(message, "reaction 5 (O2 + N2 <=> O + O + N2): forward rate: A must be"))
          << message;
        input.reactions[5].forward.A = 0.0;
        expect(starts(error_message([&] { build(input); }), "reaction 5"));
        input.reactions[5].forward.A = std::numeric_limits<double>::infinity();
        expect(starts(error_message([&] { build(input); }), "reaction 5"));

        input = canary::air5_kim_jo;
        input.reactions[1].forward.activation_temperature = -1.0;
        expect(starts(error_message([&] { build(input); }), "reaction 1"));
        input = canary::air5_kim_jo;
        input.reactions[1].forward.n = std::numeric_limits<double>::quiet_NaN();
        expect(starts(error_message([&] { build(input); }), "reaction 1"));
        input = canary::air5_kim_jo;
        input.reactions[0].reactants.entries[0].value = 1.5;
        expect(starts(error_message([&] { build(input); }), "reaction 0"));
        input = canary::air5_kim_jo;
        input.reactions[0].products.count = 0;
        expect(starts(error_message([&] { build(input); }), "reaction 0"));
        input = canary::air5_kim_jo;
        input.reactions[1].reactants.entries[1].name = "N2";
        expect(starts(error_message([&] { build(input); }), "reaction 1"));
        input = canary::air5_kim_jo;
        input.reactions[10].efficiencies.entries[0].value = -1.0;
        expect(starts(error_message([&] { build(input); }), "reaction 10"));
        input = canary::air5_kim_jo;
        input.reactions[10].third_body = false;
        expect(starts(error_message([&] { build(input); }), "reaction 10"));
        input = canary::air5_kim_jo;
        input.reactions[11].backward_model = canary::BackwardRateModel::arrhenius;
        expect(starts(error_message([&] { build(input); }),
                      "reaction 11 (NO + O <=> O2 + N): backward rate"));

        for (const double s : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN()}) {
            input = canary::air5_kim_jo;
            input.reactions[2].forward.translational_exponent = s;
            message = error_message([&] { build(input); });
            expect(starts(message,
                          "reaction 2 (N2 + NO <=> N + N + NO): forward rate: translational "
                          "exponent must be finite and between 0 and 1"))
              << message;
            input = canary::air5_kim_jo;
            input.reactions[11].backward_model = canary::BackwardRateModel::arrhenius;
            input.reactions[11].backward = input.reactions[11].forward;
            input.reactions[11].backward.translational_exponent = s;
            message = error_message([&] { build(input); });
            expect(starts(message,
                          "reaction 11 (NO + O <=> O2 + N): backward rate: translational exponent"))
              << message;
        }
        // An equilibrium backward rate ignores its Arrhenius input.
        input = canary::air5_kim_jo;
        input.reactions[11].backward.translational_exponent = 2.0;
        build(input);

        input = canary::air5_kim_jo;
        input.temperature_lower = 40000.0;
        expect(starts(error_message([&] { build(input); }), "temperature limits"));
        expect(starts(error_message([&] {
                          (void)canary::make_air5_kim_jo(built.species,
                                                         si_scales,
                                                         canary::Pressure<double, si> {0.0});
                      }),
                      "standard pressure"));
        auto species = built.species;
        species.n = 4;
        expect(starts(error_message([&] { (void)canary::make_air5_kim_jo(species, si_scales); }),
                      "species data"));
        expect(starts(error_message([&] {
                          (void)canary::make_reaction_mechanism(
                            canary::air5_kim_jo,
                            std::array<std::string_view, 5> {"N2", "O2", "N", "N2", "NO"},
                            built.species,
                            si_scales);
                      }),
                      "species \"N2\" is named twice"));
    };
}
