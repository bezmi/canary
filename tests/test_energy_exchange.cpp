#include "dual.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <cstddef>
#include <gaslib/build.hpp>
#include <gaslib/data/air5.hpp>
#include <gaslib/data/air5_park.hpp>
#include <gaslib/data/air5_park_exchange.hpp>
#include <gaslib/energy_exchange.hpp>
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/models/two_temperature.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/reactions.hpp>
#include <gaslib/species.hpp>
#include <gaslib/status.hpp>
#include <gaslib/thermo.hpp>
#include <gaslib/units.hpp>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <numbers>
#include <random>
#include <string>
#include <string_view>
#include <type_traits>

namespace
{

using gaslib::scaled;
using gaslib::si;
using D = gaslib::test::Dual<double>;
using Air = gaslib::Air5TwoTemperature<double, si>;
using AirScaled = gaslib::Air5TwoTemperature<double, scaled>;
using AirSingle = gaslib::Air5<double, si>;
using Mechanism = gaslib::Air5Park<double, si>;
using Exchange = gaslib::Air5ParkExchange<double, si>;
using Rho = gaslib::Density<double, si>;
using Temp = gaslib::Temperature<double, si>;
using Extents = Kokkos::dextents<std::size_t, 1>;

constexpr std::size_t species_count = 5;
constexpr std::size_t reaction_count = 5;
using Fractions = std::array<double, species_count>;
using PerReaction = std::array<double, reaction_count>;

constexpr double atmosphere = 101325.0;
constexpr double avogadro = gaslib::constants::avogadro<double>;
constexpr double boltzmann = gaslib::constants::boltzmann<double>;
// Universal gas constant of the reference implementation [J/(mol K)].
constexpr double reference_gas_constant = 8.31451;

const gaslib::UnitScales<double, scaled> reference_scales {
  .mass = 4.65e-26,
  .density = 1.0e-3,
  .temperature = 1000.0,
  .length = 1.0,
};

template <typename T> auto fractions(const std::array<T, species_count>& y)
{
    return gaslib::as_kind<gaslib::kind::mass_fraction>(
      Kokkos::mdspan<const T, Extents>(y.data(), species_count));
}

bool relative(double a, double b, double bound) { return std::abs(a - b) <= bound * std::abs(b); }

// A state of a reference implementation of the two-temperature air model,
// the mechanism air5_park and the exchange air5_park_exchange, computed with
// R_u = 8.31451 J/(mol K), k_B = 1.380658e-23 J/K and the standard pressure
// 1 atm: rates of progress in mol/(m^3 s), mass production rates in
// kg/(m^3 s), and the sources of vibrational-electronic energy Q = Q_VT +
// Q_CV in J/(m^3 s). The values come from GDTk (gdtk-uq/gdtk commit
// bf559c3f, libgas.so built with -version=multi_T_gas, through its C
// interface). rates_listed is false for a state whose listed rates of
// progress were evaluated with other temperature limits than its sources.
struct ReferenceState {
    Fractions y;
    double rho;
    double temp;
    double t_ve;
    PerReaction forward;
    PerReaction backward;
    Fractions source;
    double total;
    double vt;
    double cv;
    bool rates_listed = true;
};

constexpr Fractions mixed {0.6, 0.1, 0.1, 0.15, 0.05};
constexpr Fractions dissociated {0.3, 0.1, 0.3, 0.25, 0.05};
constexpr Fractions cold_air {0.767, 0.233, 0.0, 0.0, 0.0};

// clang-format off
const std::array<ReferenceState, 8> reference_states {
 {
  {
    .y = mixed, .rho = 0.01, .temp = 10000.0, .t_ve = 3000.0,
    .forward = {1.6038512832126857, 3212.233354819928, 306.4285531904488, 1876.3961515915175, 27622.447324252233},
    .backward = {0.715679576844746, 1.6643762994328262, 32.96268696228066, 42133.758577092616, 13348.95451299561},
    .source = {-399.8739442281355, -1390.9216434110062, -360.09303146806093, 522.8359679161556, 1628.052651191047},
    .total = 60227278346.57714, .vt = 60311142698.43336, .cv = -83864351.8562156,
  },
  {
    .y = mixed, .rho = 0.001, .temp = 20000.0, .t_ve = 5000.0,
    .forward = {70.20764414259725, 1771.7729507439935, 1562.6923332997578, 49.62218884836747, 942.0578300944791},
    .backward = {0.0001482166476693414, 0.0006167115915839476, 0.027006752087634416, 391.3994131968709, 52.965633125358245},
    .source = {-26.873246016946766, -67.63104960871893, 31.520611722732454, 72.93958501485652, -9.955901111923275},
    .total = 4794024843.332647, .vt = 4949330548.142285, .cv = -155305704.8096371,
  },
  {
    .y = dissociated, .rho = 0.01, .temp = 1500.0, .t_ve = 4000.0,
    .forward = {4.19924939241379e-11, 0.029209560525757138, 2.590128440894305e-05, 0.051112163981894596, 5.4417893045302135e-05},
    .backward = {260.10349575368025, 75.49258341184914, 264.1670437130537, 26293.81998235093, 256709.92677607664},
    .source = {7198.604244492461, -838.9543139161444, -3974.9343544826347, 4521.248076060205, -6905.963652153886},
    .total = -443868639.6751198, .vt = -661509592.4220978, .cv = 217640952.74697796,
  },
  // T_ve = 6000 (1 + 1e-9).
  {
    .y = dissociated, .rho = 0.01, .temp = 6000.0, .t_ve = 6000.000006,
    .forward = {7.539094013086201, 13320.078136616828, 2160.7153006598637, 855.141739678003, 2965.7585011250517},
    .backward = {28.198111937994238, 13.235054575841255, 321.0484793376328, 127876.93303414612, 74655.8056053904},
    .source = {2008.8606948833647, -4490.347905687039, -2758.1081749660634, 3634.5067630313724, 1605.0886227383662},
    .total = -3739566326.40284, .vt = -35.436764796991206, .cv = -3739566290.9660754,
  },
  {
    .y = mixed, .rho = 0.01, .temp = 8000.0, .t_ve = 8000.0,
    .forward = {592.0890608030162, 55944.841118592776, 23666.969270122936, 1153.8486379368017, 13220.54821709806},
    .backward = {1.083895914765931, 2.0053804105740167, 31.97928689432072, 43494.89097772273, 17911.16679622102},
    .source = {114.8441104183226, -3144.9661578412865, -311.1540867087124, 2920.7276266626036, 420.5485074690731},
    .total = -20864549082.012775, .vt = 0.0, .cv = -20864549082.012775,
  },
  {
    .y = cold_air, .rho = 0.1, .temp = 5000.0, .t_ve = 1000.0,
    .forward = {3.0006878820365215e-11, 0.13262844379992447, 0.0, 0.0, 0.0},
    .backward = {0.0, 0.0, 0.0, 0.0, 0.0},
    .source = {-8.405946991464189e-13, -0.004243951047465023, 8.405946991464189e-13, 0.004243951047465023, 0.0},
    .total = 105511841694.30983, .vt = 105511842179.75139, .cv = -485.4415485772803,
  },
  // T above the upper limit 50000 K of the mechanism. The reference's
  // stand-alone rates of progress clamp T to 30000 K instead of the
  // mechanism's limits, so they give m_N2 (-(q_f,0 - q_b,0) - (q_f,4 -
  // q_b,4)) = -130 kg/(m^3 s) for N2, not the listed -316 kg/(m^3 s), which
  // uses 50000 K; the species and energy sources are compared, the rates are
  // not.
  {
    .y = mixed, .rho = 0.001, .temp = 60000.0, .t_ve = 10000.0,
    .forward = {3487.576348254868, 9609.83827637619, 37993.52141008766, 68.62142173516973, 1191.0631079824614},
    .backward = {1.0965662903628135e-05, 6.344260726703739e-05, 0.004228635665655557, 343.4663649334981, 32.29428595713462},
    .source = {-315.78569603087834, -466.8418568147578, 1715.7398959191787, 2065.9642215182835, -2999.076564591826},
    .total = -4717349647.927643, .vt = 6292928353.119854, .cv = -11010278001.047497,
    .rates_listed = false,
  },
  // T below the lower limit 300 K of the mechanism.
  {
    .y = dissociated, .rho = 0.1, .temp = 250.0, .t_ve = 2000.0,
    .forward = {1.0604386645370377e-51, 2.5404098502789777e-22, 2.9427373301523874e-32, 1.5244638855451262e-22, 9.182412110315756e-47},
    .backward = {7888073.289522587, 1845602.6278664698, 412594.52436420403, 189.76248978040806, 29346583.696194474},
    .source = {1043069.3400036863, 59050.997196615615, -637802.2918168769, 403872.4330731866, -868190.4784566117},
    .total = 4184945841820.315, .vt = -130923021.33504651, .cv = 4185076764841.6504,
  },
 },
};
// clang-format on

// Change of the coefficient of species i in reaction k.
int coefficient_change(const Mechanism& mech, std::size_t k, std::size_t i)
{ return mech.reactions[k].products[i] - mech.reactions[k].reactants[i]; }

// Rates of progress, mass production rates and energy sources of gaslib in
// one cell, as plain values.
struct Values {
    PerReaction forward;
    PerReaction backward;
    Fractions production;
    double total;
    double vt;
    double cv;
};

template <typename M, typename X, typename G>
Values evaluate(const M& mech,
                const X& exchange,
                const G& gas,
                const Fractions& y,
                double rho,
                double temp,
                double t_ve)
{
    using U = M::units;
    const auto mix = gaslib::mixture(gas, fractions(y));
    const gaslib::Density<double, U> rho_q {rho};
    const gaslib::Temperature<double, U> temp_q {temp};
    const gaslib::Temperature<double, U> t_ve_q {t_ve};
    const auto q = gaslib::rates_of_progress(mech, gas, mix, rho_q, temp_q, t_ve_q);
    const auto sources =
      gaslib::two_temperature_sources(mech, exchange, gas, mix, rho_q, temp_q, t_ve_q);
    Values out {};
    for (std::size_t k = 0; k < reaction_count; ++k) {
        out.forward[k] = q.forward[k].value();
        out.backward[k] = q.backward[k].value();
    }
    for (std::size_t i = 0; i < species_count; ++i) {
        out.production[i] = sources.species[i].value();
    }
    out.total = sources.vibrational_electronic.value();
    out.vt =
      gaslib::vibration_translation_source(exchange, gas, mix, rho_q, temp_q, t_ve_q).value();
    out.cv = gaslib::chemistry_vibration_source(exchange, mech, q, temp_q, t_ve_q).value();
    return out;
}

// G(x) / k_B in long double (64-bit significand): the definition with expm1
// for |theta_d x| >= 1e-3, where the cancellation between the two terms
// costs at most 4 of the 19 digits; below, the series
// (theta_d - theta_v) / 2 + sum_k B_2k / (2k)! (theta_v^2k - theta_d^2k) x^(2k - 1)
// to k = 4, whose first omitted term is below 1e-30 of the result.
long double precise_g(long double theta_v, long double theta_d, long double x)
{
    if (std::abs(theta_d * x) >= 1.0e-3L) {
        return (theta_v / std::expm1(theta_v * x)) - (theta_d / std::expm1(theta_d * x));
    }
    const std::array<long double, 4> b {
      1.0L / 12.0L,
      -1.0L / 720.0L,
      1.0L / 30240.0L,
      -1.0L / 1209600.0L,
    };
    long double sum = (theta_d - theta_v) / 2.0L;
    long double power_v = theta_v;
    long double power_d = theta_d;
    long double power_x = 1.0L;
    for (const long double c : b) {
        power_v *= theta_v;
        power_d *= theta_d;
        sum += c * (power_v - power_d) * power_x * x;
        power_v *= theta_v;
        power_d *= theta_d;
        power_x *= x * x;
    }
    return sum;
}

template <typename F> std::string error_message(const F& f)
{
    try {
        f();
    } catch (const gaslib::ParameterError& error) {
        return error.what();
    }
    return "no exception";
}

}  // namespace

static_assert(std::numeric_limits<long double>::digits >= 64);
static_assert(gaslib::HasSpeciesProperties<Air, double, 5>);
static_assert(gaslib::HasSpeciesProperties<Air, D, 5>);
static_assert(gaslib::ReactingModel<Air, double, double, 5>);
static_assert(gaslib::ReactingModel<AirScaled, D, double, 5>);
static_assert(gaslib::HasVibrationalElectronicEnergy<Air, double, 5>);
static_assert(gaslib::HasVibrationalElectronicEnergy<Air, float, 5>);
static_assert(gaslib::HasVibrationalElectronicEnergy<AirScaled, D, 5>);
static_assert(!gaslib::HasVibrationalElectronicEnergy<Air, double, 4>);
static_assert(!gaslib::HasVibrationalElectronicEnergy<AirSingle, double, 5>);
static_assert(gaslib::RelaxingModel<Air, double, double, 5>);
static_assert(gaslib::RelaxingModel<Air, float, double, 5>);
static_assert(gaslib::RelaxingModel<Air, D, double, 5>);
static_assert(!gaslib::RelaxingModel<AirSingle, double, double, 5>);
static_assert(std::is_trivially_copyable_v<Exchange>);
static_assert(std::is_trivially_copyable_v<gaslib::Air5ParkExchange<float, scaled>>);
static_assert(std::is_trivially_copyable_v<gaslib::EnergyExchange<double, si, 5, 0, 0>>);
static_assert(std::is_trivially_copyable_v<gaslib::TwoTemperatureSources<D, si, 5>>);
static_assert(gaslib::air5_park.reactions.size() == 5);
static_assert(gaslib::air5_park_exchange.vibration_translation.size() == 15);

void test_energy_exchange()
{
    using namespace boost::ut;
    const gaslib::UnitScales<double, si> si_scales {};
    const auto built = gaslib::make_air5_two_temperature<double>(si_scales);
    const Air& gas = built.model;
    const Mechanism mech =
      gaslib::make_air5_park(built.species, si_scales, gaslib::Pressure<double, si> {atmosphere});
    const Exchange exchange = gaslib::make_air5_park_exchange(built.species, mech, si_scales);
    const double ru = gaslib::constants::universal_gas_constant<double>;

    "energy exchange: the generated mechanism and exchange hold their constants"_test = [&] {
        expect(mech.reactions[0].forward.translational_exponent == 0.5_d);
        expect(mech.reactions[2].forward.translational_exponent == 0.5_d);
        expect(mech.reactions[3].forward.translational_exponent == 1.0_d);
        expect(mech.reactions[4].forward.translational_exponent == 1.0_d);
        // N2 + M <=> N + N + M: order 2.
        expect(relative(mech.reactions[0].forward.A, 7.0e15 / avogadro, 1.0e-15));
        expect(mech.temperature_upper.value() == 50000.0_d);

        // Pair 2 is N2-N2 and pair 12 NO-N.
        const auto& n2_n2 = exchange.vibration_translation[2];
        expect(n2_n2.species == 0_ul && n2_n2.partner == 0_ul);
        expect(exchange.vibration_translation[12].species == 4_ul);
        expect(exchange.vibration_translation[12].partner == 2_ul);
        expect(n2_n2.millikan_white_a == 221.0_d);
        expect(relative(n2_n2.millikan_white_offset,
                        std::log(atmosphere) - (221.0 * 0.029) - 18.42,
                        1.0e-15));
        // 1 / K = 1 / (sigma_ref T_ref^2 sqrt(8 k_B / (pi mu))), mu = m_N2 / 2.
        const double mu = built.species.mass[0].value() / 2.0;
        const double k =
          1.0e-21 * 50000.0 * 50000.0 * std::sqrt(8.0 * boltzmann / (std::numbers::pi * mu));
        expect(relative(n2_n2.correction_coefficient, 1.0 / k, 1.0e-15));
        expect(n2_n2.correction_exponent == 1.5_d);

        const auto& n2 = exchange.chemistry_vibration[0];
        expect(n2.reaction == 0_ul && n2.species == 0_ul);
        expect(relative(n2.dissociation_temperature.value(), 945360.0 / ru, 1.0e-15));
        expect(n2.vibrational_temperature.value() == 3393.44_d);
        expect(exchange.chemistry_vibration[2].reaction == 2_ul);
        expect(exchange.chemistry_vibration[2].species == 4_ul);
        expect(exchange.boltzmann.value() == boltzmann);
    };

    "energy exchange: rates, species sources and energy sources equal the reference"_test = [&] {
        const double r = ru / reference_gas_constant;
        for (const ReferenceState& s : reference_states) {
            const Values got = evaluate(mech, exchange, gas, s.y, s.rho, s.temp, s.t_ve);
            PerReaction forward {};
            PerReaction backward {};
            for (std::size_t k = 0; k < reaction_count; ++k) {
                forward[k] = avogadro * s.forward[k];
                backward[k] =
                  avogadro * s.backward[k] * std::pow(r, mech.reactions[k].coefficient_change);
                if (!s.rates_listed) { continue; }
                // Forward rates contain no constant that differs between
                // the implementations; backward rates contain R_u through
                // the Gibbs energies, corrected above, and differ by
                // rounding amplified by the exponent of K_c (up to about
                // 50).
                expect(relative(got.forward[k], forward[k], 1.0e-12))
                  << s.temp << k << got.forward[k] << forward[k];
                expect(relative(got.backward[k], backward[k], 1.0e-10))
                  << s.temp << k << got.backward[k] << backward[k];
            }
            for (std::size_t i = 0; i < species_count; ++i) {
                // The bound is relative to the sum of the magnitudes of the
                // terms, the rounding scale of the backward rates.
                double expected = s.source[i];
                double scale = 0;
                for (std::size_t k = 0; k < reaction_count; ++k) {
                    const double m = built.species.mass[i].value();
                    const int dnu = coefficient_change(mech, k, i);
                    scale += m * std::abs(dnu) * (got.forward[k] + got.backward[k]);
                }
                if (s.rates_listed) {
                    expected = 0;
                    for (std::size_t k = 0; k < reaction_count; ++k) {
                        expected += built.species.mass[i].value() * coefficient_change(mech, k, i) *
                                    (forward[k] - backward[k]);
                    }
                }
                expect(std::abs(got.production[i] - expected) <= 1.0e-10 * scale)
                  << s.temp << i << got.production[i] << expected;
            }
            // k_B, R_u and N_A of the reference differ from the exact values
            // by up to 7e-6 relative; they enter e_ve, tau and G linearly or
            // through exponents of order 1.
            const double scale = std::abs(s.vt) + std::abs(s.cv);
            expect(std::abs(got.total - s.total) <= 3.0e-5 * scale)
              << s.temp << got.total << s.total;
            expect(std::abs(got.vt - s.vt) <= 3.0e-5 * scale) << s.temp << got.vt << s.vt;
            expect(std::abs(got.cv - s.cv) <= 3.0e-5 * scale) << s.temp << got.cv << s.cv;
        }
    };

    "energy exchange: the forward rate with s = 1/2 is the Arrhenius rate at sqrt(T T_ve)"_test =
      [&] {
          for (const auto& [temp, t_ve] : {
                 std::array {10000.0, 3000.0},
                 std::array {1500.0, 4000.0},
                 std::array {250.0, 250.0},
               }) {
              const auto mix = gaslib::mixture(gas, fractions(mixed));
              const double rho = 0.01;
              const auto q =
                gaslib::rates_of_progress(mech, gas, mix, Rho {rho}, Temp {temp}, Temp {t_ve});
              // Below the lower limit T_l = 300 K; T_ve is not clamped.
              const double t_l = std::max(temp, 300.0);
              const double t_c = std::sqrt(t_l * t_ve);
              std::array<double, species_count> n {};
              for (std::size_t i = 0; i < species_count; ++i) {
                  n[i] = rho * mixed[i] / built.species.mass[i].value();
              }
              for (std::size_t k = 0; k < 3; ++k) {
                  const auto& rate = mech.reactions[k].forward;
                  double third_body = 0;
                  for (std::size_t i = 0; i < species_count; ++i) {
                      third_body += mech.reactions[k].efficiencies[i] * n[i];
                  }
                  // Reactions 0 to 2 dissociate species 0, 1 and 4.
                  const std::size_t molecule = k == 2 ? 4 : k;
                  const double expected = rate.A * std::pow(t_c, rate.n) *
                                          std::exp(-rate.activation_temperature.value() / t_c) *
                                          third_body * n[molecule];
                  // The exponent theta / T_c is at most about 410 here, so a
                  // few ulp in T_c give up to about 2e-13 in k_f.
                  expect(relative(q.forward[k].value(), expected, 1.0e-12))
                    << temp << t_ve << k << q.forward[k].value() << expected;
              }
          }
      };

    "energy exchange: rates at T_ve = T equal the single-temperature rates"_test = [&] {
        for (const ReferenceState& s : reference_states) {
            if (s.temp < 300.0 || s.temp > 50000.0) { continue; }
            const auto mix = gaslib::mixture(gas, fractions(s.y));
            const auto one = gaslib::rates_of_progress(mech, gas, mix, Rho {s.rho}, Temp {s.temp});
            const auto two =
              gaslib::rates_of_progress(mech, gas, mix, Rho {s.rho}, Temp {s.temp}, Temp {s.temp});
            const auto omega_one =
              gaslib::production_rates(mech, gas, mix, Rho {s.rho}, Temp {s.temp});
            const auto omega_two =
              gaslib::production_rates(mech, gas, mix, Rho {s.rho}, Temp {s.temp}, Temp {s.temp});
            for (std::size_t k = 0; k < reaction_count; ++k) {
                expect(one.forward[k] == two.forward[k] && one.backward[k] == two.backward[k])
                  << s.temp << k;
            }
            for (std::size_t i = 0; i < species_count; ++i) {
                expect(omega_one[i] == omega_two[i]) << s.temp << i;
            }
        }
    };

    "energy exchange: rates with s = 1 and equilibrium backward rates do not depend on T_ve"_test =
      [&] {
          const auto mix = gaslib::mixture(gas, fractions(dissociated));
          const Rho rho {0.01};
          for (const double temp : {250.0, 1500.0, 8000.0, 60000.0}) {
              const auto a =
                gaslib::rates_of_progress(mech, gas, mix, rho, Temp {temp}, Temp {1000.0});
              const auto b =
                gaslib::rates_of_progress(mech, gas, mix, rho, Temp {temp}, Temp {9000.0});
              for (std::size_t k = 0; k < reaction_count; ++k) {
                  expect(a.backward[k] == b.backward[k]) << temp << k;
              }
              expect(a.forward[3] == b.forward[3] && a.forward[4] == b.forward[4]) << temp;
              expect(!(a.forward[0] == b.forward[0])) << temp;
              // Derivatives with respect to T_ve.
              const std::array<D, species_count> y_d {0.3, 0.1, 0.3, 0.25, 0.05};
              const auto q = gaslib::rates_of_progress(mech,
                                                       gas,
                                                       gaslib::mixture(gas, fractions(y_d)),
                                                       gaslib::Density<D, si> {D(0.01)},
                                                       gaslib::Temperature<D, si> {D(temp)},
                                                       gaslib::Temperature<D, si> {D(3000.0, 1.0)});
              for (std::size_t k = 0; k < reaction_count; ++k) {
                  expect(q.backward[k].value().dual() == 0.0) << temp << k;
              }
              expect(q.forward[3].value().dual() == 0.0 && q.forward[4].value().dual() == 0.0);
              expect(q.forward[0].value().dual() > 0.0);
          }
      };

    "energy exchange: below the lower limit T is the limit and T_ve is not clamped"_test = [&] {
        const auto mix = gaslib::mixture(gas, fractions(dissociated));
        const Rho rho {0.1};
        for (const double t_ve : {250.0, 2000.0}) {
            const auto below =
              gaslib::rates_of_progress(mech, gas, mix, rho, Temp {250.0}, Temp {t_ve});
            const auto limit =
              gaslib::rates_of_progress(mech, gas, mix, rho, Temp {300.0}, Temp {t_ve});
            for (std::size_t k = 0; k < reaction_count; ++k) {
                expect(below.forward[k] == limit.forward[k] &&
                       below.backward[k] == limit.backward[k])
                  << t_ve << k;
            }
        }
        const auto clamped =
          gaslib::rates_of_progress(mech, gas, mix, rho, Temp {250.0}, Temp {300.0});
        const auto unclamped =
          gaslib::rates_of_progress(mech, gas, mix, rho, Temp {250.0}, Temp {250.0});
        expect(unclamped.forward[0].value() < clamped.forward[0].value());
    };

    "energy exchange: two-temperature production rates conserve mass"_test = [&] {
        const auto check = [&](const Fractions& y, double rho, double temp, double t_ve) {
            const Values got = evaluate(mech, exchange, gas, y, rho, temp, t_ve);
            double sum = 0;
            double scale = 0;
            for (const double w : got.production) {
                sum += w;
                scale += std::abs(w);
            }
            // The sum of terms of opposite signs of size scale, in double.
            expect(std::abs(sum) <= 1.0e-13 * scale) << temp << t_ve << sum << scale;
        };
        for (const ReferenceState& s : reference_states) { check(s.y, s.rho, s.temp, s.t_ve); }
        // A fixed seed gives the same states in every run.
        std::mt19937 engine(2468);  // NOLINT(bugprone-random-generator-seed)
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        for (int n = 0; n < 100; ++n) {
            Fractions y {};
            double total = 0;
            for (double& v : y) {
                v = unit(engine);
                total += v;
            }
            for (double& v : y) { v /= total; }
            const double rho = std::pow(10.0, -4.0 + (4.0 * unit(engine)));
            const double temp = 300.0 + (29700.0 * unit(engine));
            const double t_ve = 300.0 + (29700.0 * unit(engine));
            check(y, rho, temp, t_ve);
        }
    };

    "energy exchange: Q_VT is zero at T = T_ve and has the sign of T - T_ve"_test = [&] {
        for (const auto& y : {mixed, dissociated, cold_air}) {
            const auto mix = gaslib::mixture(gas, fractions(y));
            const Rho rho {0.01};
            for (const double temp : {500.0, 3000.0, 9000.0, 25000.0}) {
                const auto same = gaslib::vibration_translation_source(exchange,
                                                                       gas,
                                                                       mix,
                                                                       rho,
                                                                       Temp {temp},
                                                                       Temp {temp});
                expect(same.value() == 0.0) << temp;
                for (const double factor : {0.5, 0.99, 1.01, 2.0}) {
                    const double t_ve = factor * temp;
                    const double q = gaslib::vibration_translation_source(exchange,
                                                                          gas,
                                                                          mix,
                                                                          rho,
                                                                          Temp {temp},
                                                                          Temp {t_ve})
                                       .value();
                    expect(std::isfinite(q));
                    expect(temp > t_ve ? q > 0.0 : q < 0.0) << temp << t_ve << q;
                }
            }
        }
    };

    "energy exchange: an absent partner contributes zero and an empty exchange gives zero"_test =
      [&] {
          // N2 relaxing with O only, in air without atoms.
          gaslib::EnergyExchangeInput<1, 0> input {};
          input.vibration_translation[0] = gaslib::air5_park_exchange.vibration_translation[1];
          const auto absent = gaslib::make_energy_exchange(input,
                                                           gaslib::detail::air5_park_species_names,
                                                           built.species,
                                                           mech,
                                                           si_scales);
          const auto mix = gaslib::mixture(gas, fractions(cold_air));
          const double q = gaslib::vibration_translation_source(absent,
                                                                gas,
                                                                mix,
                                                                Rho {0.1},
                                                                Temp {5000.0},
                                                                Temp {1000.0})
                             .value();
          expect(q == 0.0) << q;

          const gaslib::EnergyExchangeInput<0, 0> nothing {};
          const auto empty = gaslib::make_energy_exchange(nothing,
                                                          gaslib::detail::air5_park_species_names,
                                                          built.species,
                                                          mech,
                                                          si_scales);
          const auto sources = gaslib::two_temperature_sources(mech,
                                                               empty,
                                                               gas,
                                                               mix,
                                                               Rho {0.1},
                                                               Temp {5000.0},
                                                               Temp {1000.0});
          expect(sources.vibrational_electronic.value() == 0.0);
          expect(sources.species[1].value() < 0.0);
      };

    "energy exchange: Q_VT of pure N2 is the Landau-Teller rate with tau of the N2-N2 pair"_test =
      [&] {
          const Fractions y {1.0, 0.0, 0.0, 0.0, 0.0};
          const auto mix = gaslib::mixture(gas, fractions(y));
          const double m = built.species.mass[0].value();
          for (const double rho : {1.0e-3, 0.1}) {
              for (const double temp : {2000.0, 6000.0, 15000.0}) {
                  const double t_ve = temp * (1.0 - 1.0e-4);
                  const double q = gaslib::vibration_translation_source(exchange,
                                                                        gas,
                                                                        mix,
                                                                        Rho {rho},
                                                                        Temp {temp},
                                                                        Temp {t_ve})
                                     .value();
                  const auto e_t = gaslib::species_vibrational_electronic_energy(gas, Temp {temp});
                  const auto e_ve = gaslib::species_vibrational_electronic_energy(gas, Temp {t_ve});
                  const double rate = q / (rho * (e_t[0].value() - e_ve[0].value()));
                  // tau from the formulas, with x_q = 1 and p = n k_B T.
                  const double n = rho / m;
                  const double p = n * boltzmann * temp;
                  const double tau_mw =
                    (atmosphere / p) * std::exp((221.0 * (std::cbrt(1.0 / temp) - 0.029)) - 18.42);
                  const double sigma = 1.0e-21 * (50000.0 / temp) * (50000.0 / temp);
                  const double speed =
                    std::sqrt(8.0 * boltzmann * temp / (std::numbers::pi * (m / 2.0)));
                  const double tau_p = 1.0 / (sigma * speed * n);
                  // The exponent of tau_MW is about 10, so rounding in
                  // T^(-1/3) gives a few 1e-15.
                  expect(relative(rate, 1.0 / (tau_mw + tau_p), 1.0e-12))
                    << rho << temp << rate << 1.0 / (tau_mw + tau_p);
              }
          }
      };

    "energy exchange: each relaxation-time correction adds 1 / (sigma c n_q) to tau"_test = [&] {
        using gaslib::RelaxationCorrection;
        const Fractions y {1.0, 0.0, 0.0, 0.0, 0.0};
        const auto mix = gaslib::mixture(gas, fractions(y));
        const double m = built.species.mass[0].value();
        // sigma_0, T_ref and omega of each model; a negative Millikan-White a
        // with the last.
        const std::array<gaslib::RelaxationCorrectionInput, 4> corrections {
          {
            {
              .model = RelaxationCorrection::none,
              .cross_section = 0.0,
              .reference_temperature = 0.0,
              .exponent = 0.0,
            },
            {
              .model = RelaxationCorrection::park_1985,
              .cross_section = 5.8e-20,
              .reference_temperature = 0.0,
              .exponent = 0.0,
            },
            {
              .model = RelaxationCorrection::park_1993,
              .cross_section = 1.0e-21,
              .reference_temperature = 50000.0,
              .exponent = 0.0,
            },
            {
              .model = RelaxationCorrection::kim_2021,
              .cross_section = 1.7e-23,
              .reference_temperature = 0.0,
              .exponent = 0.55,
            },
          },
        };
        for (std::size_t c = 0; c < corrections.size(); ++c) {
            const double a = c == 3 ? -27.32 : 221.0;
            const double b = c == 3 ? 0.1465 : 0.029;
            gaslib::EnergyExchangeInput<1, 0> input {};
            input.vibration_translation[0] = {
              .species = "N2",
              .partner = "N2",
              .millikan_white = {.a = a, .b = b},
              .correction = corrections[c],
            };
            const auto single =
              gaslib::make_energy_exchange(input,
                                           gaslib::detail::air5_park_species_names,
                                           built.species,
                                           mech,
                                           si_scales);
            for (const double temp : {2000.0, 6000.0, 15000.0}) {
                const double rho = 0.01;
                const double t_ve = temp * (1.0 - 1.0e-4);
                const double q = gaslib::vibration_translation_source(single,
                                                                      gas,
                                                                      mix,
                                                                      Rho {rho},
                                                                      Temp {temp},
                                                                      Temp {t_ve})
                                   .value();
                const auto e_t = gaslib::species_vibrational_electronic_energy(gas, Temp {temp});
                const auto e_ve = gaslib::species_vibrational_electronic_energy(gas, Temp {t_ve});
                const double rate = q / (rho * (e_t[0].value() - e_ve[0].value()));
                const double n = rho / m;
                const double tau_mw = (atmosphere / (n * boltzmann * temp)) *
                                      std::exp((a * (std::cbrt(1.0 / temp) - b)) - 18.42);
                const gaslib::RelaxationCorrectionInput& k = corrections[c];
                double sigma = 0.0;
                if (k.model == RelaxationCorrection::park_1985) { sigma = k.cross_section; }
                if (k.model == RelaxationCorrection::park_1993) {
                    sigma = k.cross_section * std::pow(k.reference_temperature / temp, 2);
                }
                if (k.model == RelaxationCorrection::kim_2021) {
                    sigma = k.cross_section * std::pow(temp, k.exponent);
                }
                const double speed =
                  std::sqrt(8.0 * boltzmann * temp / (std::numbers::pi * (m / 2.0)));
                const double tau_corr = sigma > 0.0 ? 1.0 / (sigma * speed * n) : 0.0;
                // As for the N2-N2 pair above: a few 1e-15 from T^(-1/3) in an
                // exponent of about 10, and T^beta through exp(beta ln T).
                expect(relative(rate, 1.0 / (tau_mw + tau_corr), 1.0e-12))
                  << c << temp << rate << 1.0 / (tau_mw + tau_corr);
            }
        }
    };

    "energy exchange: relaxation at fixed energy reaches the equilibrium temperature"_test = [&] {
        const Fractions y {1.0, 0.0, 0.0, 0.0, 0.0};
        const auto mix = gaslib::mixture(gas, fractions(y));
        const double rho = 0.01;
        const std::array t_ve_start {1000.0};
        const auto start =
          gaslib::thermo_state(gas,
                               mix,
                               Rho {rho},
                               Temp {8000.0},
                               gaslib::as_kind<gaslib::kind::mode_temperature, si>(
                                 Kokkos::mdspan<const double, Extents>(t_ve_start.data(), 1)));
        const double rho_e = rho * start.state.e.value();
        const auto equilibrium =
          gaslib::thermo_state(gas, mix, Rho {rho}, gaslib::EnergyDensity<double, si> {rho_e});
        const double t_eq = equilibrium.state.temp.value();
        std::array rho_e_ve {rho * start.state.mode_e[0].value()};
        double temp = 8000.0;
        double t_ve = 1000.0;
        bool monotonic = true;
        bool converged = true;
        int steps = 0;
        // Explicit steps of 1 % of the current relaxation time
        // rho (e_ve(T) - e_ve(T_ve)) / Q_VT, until T - T_ve < 1e-6 T; T and
        // T_ve from the closure of (rho, rho e, rho e_ve).
        while (temp - t_ve > 1.0e-6 * temp && steps < 100000) {
            const double q = gaslib::vibration_translation_source(exchange,
                                                                  gas,
                                                                  mix,
                                                                  Rho {rho},
                                                                  Temp {temp},
                                                                  Temp {t_ve})
                               .value();
            const auto e_t = gaslib::species_vibrational_electronic_energy(gas, Temp {temp});
            const auto e_ve = gaslib::species_vibrational_electronic_energy(gas, Temp {t_ve});
            const double tau = rho * (e_t[0].value() - e_ve[0].value()) / q;
            rho_e_ve[0] += 0.01 * tau * q;
            const auto next =
              gaslib::thermo_state(gas,
                                   mix,
                                   Rho {rho},
                                   gaslib::EnergyDensity<double, si> {rho_e},
                                   gaslib::as_kind<gaslib::kind::mode_energy_density, si>(
                                     Kokkos::mdspan<const double, Extents>(rho_e_ve.data(), 1)));
            converged = converged && next.status == gaslib::Status::Valid;
            const double next_temp = next.state.temp.value();
            const double next_t_ve = next.state.mode_temp[0].value();
            monotonic =
              monotonic && next_temp <= temp && next_t_ve >= t_ve && next_t_ve <= next_temp;
            temp = next_temp;
            t_ve = next_t_ve;
            ++steps;
        }
        expect(converged);
        expect(monotonic);
        expect(steps < 100000);
        // e = e_tr(T) + e_ve(T_ve) with both increasing, so T_eq lies
        // between T_ve and T, which differ by less than 1e-6 T; the
        // closures solve to 1e-10.
        expect(t_ve <= t_eq * (1.0 + 1.0e-10) && t_eq <= temp * (1.0 + 1.0e-10))
          << t_ve << t_eq << temp;
        expect(relative(temp, t_eq, 1.0e-6)) << temp << t_eq;
    };

    "energy exchange: G is accurate for every x"_test = [&] {
        std::array<double, 60> xs {};
        std::size_t count = 0;
        xs[count++] = 0.0;
        for (const double sign : {1.0, -1.0}) {
            xs[count++] = sign * 1.0e-300;
            for (int e = -12; e <= -2; ++e) { xs[count++] = sign * std::pow(10.0, e); }
        }
        for (const auto& coupling : exchange.chemistry_vibration) {
            const double theta_v = coupling.vibrational_temperature.value();
            const double theta_d = coupling.dissociation_temperature.value();
            std::size_t total = count;
            for (const double sign : {1.0, -1.0}) {
                for (const double w :
                     {0.1, 0.4999, 0.5, 0.5001, 1.0, 3.0, 10.0, 100.0, 1000.0, 1.0e4}) {
                    xs[total++] = sign * w / theta_d;
                }
            }
            for (std::size_t j = 0; j < total; ++j) {
                const double x = xs[j];
                const double g =
                  gaslib::detail::truncated_oscillator_energy(theta_v, theta_d, x).value;
                const auto precise = static_cast<double>(precise_g(theta_v, theta_d, x));
                expect(relative(g, precise, 1.0e-12)) << theta_d << x << g << precise;
            }
            const auto g0 = gaslib::detail::truncated_oscillator_energy(theta_v, theta_d, 0.0);
            expect(g0.value == (theta_d - theta_v) / 2.0);
        }
    };

    "energy exchange: G_vanish equals G_appear at T = T_ve and Q_CV < 0 without products"_test =
      [&] {
          for (const ReferenceState& s : reference_states) {
              const auto mix = gaslib::mixture(gas, fractions(s.y));
              const auto q = gaslib::rates_of_progress(mech,
                                                       gas,
                                                       mix,
                                                       Rho {s.rho},
                                                       Temp {s.temp},
                                                       Temp {s.temp});
              const double cv =
                gaslib::chemistry_vibration_source(exchange, mech, q, Temp {s.temp}, Temp {s.temp})
                  .value();
              // With G_vanish = G_appear = k_B (theta_D - theta_v) / 2,
              // Q_CV = sum_c dnu (q_f - q_b) G_appear.
              double expected = 0;
              double scale = 0;
              for (const auto& c : exchange.chemistry_vibration) {
                  const double g =
                    boltzmann *
                    (c.dissociation_temperature.value() - c.vibrational_temperature.value()) / 2.0;
                  const double dnu = coefficient_change(mech, c.reaction, c.species);
                  expected +=
                    dnu * (q.forward[c.reaction].value() - q.backward[c.reaction].value()) * g;
                  scale += std::abs(dnu) *
                           (q.forward[c.reaction].value() + q.backward[c.reaction].value()) * g;
              }
              expect(std::abs(cv - expected) <= 1.0e-14 * scale) << s.temp << cv << expected;
          }
          const auto mix = gaslib::mixture(gas, fractions(cold_air));
          for (const auto& [temp, t_ve] : {
                 std::array {5000.0, 1000.0},
                 std::array {3000.0, 3000.0},
                 std::array {2000.0, 9000.0},
               }) {
              const auto q =
                gaslib::rates_of_progress(mech, gas, mix, Rho {0.1}, Temp {temp}, Temp {t_ve});
              expect(q.backward[0].value() == 0.0 && q.backward[1].value() == 0.0);
              const double cv =
                gaslib::chemistry_vibration_source(exchange, mech, q, Temp {temp}, Temp {t_ve})
                  .value();
              expect(cv < 0.0) << temp << t_ve << cv;
          }
      };

    "energy exchange: scaled units give the SI results after conversion"_test = [&] {
        const auto built_s = gaslib::make_air5_two_temperature<double>(reference_scales);
        const auto mech_s = gaslib::make_air5_park(
          built_s.species,
          reference_scales,
          gaslib::convert<scaled>(gaslib::Pressure<double, si> {atmosphere}, reference_scales));
        const auto exchange_s =
          gaslib::make_air5_park_exchange(built_s.species, mech_s, reference_scales);
        for (const ReferenceState& s : reference_states) {
            const auto mix = gaslib::mixture(gas, fractions(s.y));
            const auto mix_s = gaslib::mixture(built_s.model, fractions(s.y));
            const Rho rho {s.rho};
            const Temp temp {s.temp};
            const Temp t_ve {s.t_ve};
            const auto rho_s = gaslib::convert<scaled>(rho, reference_scales);
            const auto temp_s = gaslib::convert<scaled>(temp, reference_scales);
            const auto t_ve_s = gaslib::convert<scaled>(t_ve, reference_scales);
            const auto q = gaslib::rates_of_progress(mech, gas, mix, rho, temp, t_ve);
            const auto q_s =
              gaslib::rates_of_progress(mech_s, built_s.model, mix_s, rho_s, temp_s, t_ve_s);
            // As for one temperature: the Gibbs energies of the scaled thermo
            // model differ from the SI ones by up to 6e-13 relative, times
            // m / (k_B T) in K_c.
            for (std::size_t k = 0; k < reaction_count; ++k) {
                const double f = gaslib::convert<si>(q_s.forward[k], reference_scales).value();
                const double b = gaslib::convert<si>(q_s.backward[k], reference_scales).value();
                expect(relative(f, q.forward[k].value(), 1.0e-12)) << s.temp << k;
                expect(relative(b, q.backward[k].value(), 1.0e-10)) << s.temp << k;
            }
            // The SI rates converted, so that only the exchange is compared.
            gaslib::RatesOfProgress<double, scaled, reaction_count> q_converted {};
            for (std::size_t k = 0; k < reaction_count; ++k) {
                q_converted.forward[k] = gaslib::convert<scaled>(q.forward[k], reference_scales);
                q_converted.backward[k] = gaslib::convert<scaled>(q.backward[k], reference_scales);
            }
            const double cv =
              gaslib::chemistry_vibration_source(exchange, mech, q, temp, t_ve).value();
            const double cv_s =
              gaslib::convert<si>(
                gaslib::chemistry_vibration_source(exchange_s, mech_s, q_converted, temp_s, t_ve_s),
                reference_scales)
                .value();
            expect(relative(cv_s, cv, 1.0e-12)) << s.temp << cv_s << cv;
            // Q_VT is a difference of energies that cancel when T_ve is
            // close to T; states with |T - T_ve| >= 1000 K.
            if (std::abs(s.temp - s.t_ve) < 1000.0) { continue; }
            const double vt =
              gaslib::vibration_translation_source(exchange, gas, mix, rho, temp, t_ve).value();
            const double vt_s =
              gaslib::convert<si>(gaslib::vibration_translation_source(exchange_s,
                                                                       built_s.model,
                                                                       mix_s,
                                                                       rho_s,
                                                                       temp_s,
                                                                       t_ve_s),
                                  reference_scales)
                .value();
            expect(relative(vt_s, vt, 1.0e-12)) << s.temp << vt_s << vt;
        }
    };

    "energy exchange: dual-number derivatives equal central differences"_test = [&] {
        for (const std::size_t state : {std::size_t {0}, std::size_t {2}}) {
            const ReferenceState& s = reference_states[state];
            const auto sources = [&]<typename T>(T temp, T t_ve) {
                std::array<T, species_count> y {};
                for (std::size_t i = 0; i < species_count; ++i) { y[i] = T(s.y[i]); }
                const auto mix = gaslib::mixture(gas, fractions(y));
                const gaslib::Density<T, si> rho {T(s.rho)};
                const gaslib::Temperature<T, si> t {temp};
                const gaslib::Temperature<T, si> tv {t_ve};
                const auto q = gaslib::rates_of_progress(mech, gas, mix, rho, t, tv);
                // The sources, then the size of the terms each is summed
                // from.
                std::array<T, 2 * (species_count + 2)> out {};
                const auto omega = gaslib::production_rates(mech, gas, mix, rho, t, tv);
                for (std::size_t i = 0; i < species_count; ++i) { out[i] = omega[i].value(); }
                out[species_count] =
                  gaslib::vibration_translation_source(exchange, gas, mix, rho, t, tv).value();
                out[species_count + 1] =
                  gaslib::chemistry_vibration_source(exchange, mech, q, t, tv).value();
                const std::size_t scales = species_count + 2;
                for (std::size_t k = 0; k < reaction_count; ++k) {
                    const T both = q.forward[k].value() + q.backward[k].value();
                    for (std::size_t i = 0; i < species_count; ++i) {
                        const double m = built.species.mass[i].value();
                        out[scales + i] =
                          out[scales + i] + T(m * std::abs(coefficient_change(mech, k, i)) * both);
                    }
                }
                out[scales + species_count] =
                  out[species_count] < T(0.0) ? T(-out[species_count]) : out[species_count];
                for (const auto& coupling : exchange.chemistry_vibration) {
                    out[scales + species_count + 1] =
                      out[scales + species_count + 1] +
                      T(boltzmann * coupling.dissociation_temperature.value() *
                        (q.forward[coupling.reaction].value() +
                         q.backward[coupling.reaction].value()));
                }
                return out;
            };
            // Column 0: d/dT; column 1: d/dT_ve.
            for (std::size_t c = 0; c < 2; ++c) {
                const auto ad =
                  sources(D(s.temp, c == 0 ? 1.0 : 0.0), D(s.t_ve, c == 1 ? 1.0 : 0.0));
                const double h = 1.0e-5 * (c == 0 ? s.temp : s.t_ve);
                const auto plus = sources(s.temp + (c == 0 ? h : 0.0), s.t_ve + (c == 1 ? h : 0.0));
                const auto minus =
                  sources(s.temp - (c == 0 ? h : 0.0), s.t_ve - (c == 1 ? h : 0.0));
                for (std::size_t i = 0; i < species_count + 2; ++i) {
                    const double fd = (plus[i] - minus[i]) / (2.0 * h);
                    // Central differences with h = 1e-5 T have truncation
                    // error of order 1e-10 relative. Their rounding error is
                    // about eps times the size of the terms summed in the
                    // value (scale) over 2 h; 1e-13 allows for 400 ulp of
                    // the terms.
                    const double rounding = 1.0e-13 * plus[species_count + 2 + i] / h;
                    expect(std::abs(fd - ad[i].dual()) <=
                           (1.0e-6 * std::abs(ad[i].dual())) + rounding)
                      << state << c << i << ad[i].dual() << fd;
                }
            }
        }
    };

    "energy exchange: the builder rejects invalid input with the entry in the message"_test = [&] {
        using Input = gaslib::EnergyExchangeInput<15, 3>;
        const auto build = [&](const Input& input) {
            (void)gaslib::make_energy_exchange(input,
                                               gaslib::detail::air5_park_species_names,
                                               built.species,
                                               mech,
                                               si_scales);
        };
        const auto starts = [&](const Input& input, std::string_view prefix) {
            const std::string message = error_message([&] { build(input); });
            const bool ok = message.starts_with(prefix);
            expect(ok) << message;
        };
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double inf = std::numeric_limits<double>::infinity();

        Input input = gaslib::air5_park_exchange;
        input.vibration_translation[3].species = "Ar";
        starts(input, "vibration-translation pair 3 (Ar, O2): unknown species \"Ar\"");
        input = gaslib::air5_park_exchange;
        input.vibration_translation[3].partner = "Ar";
        starts(input, "vibration-translation pair 3 (N2, Ar): unknown partner \"Ar\"");
        input = gaslib::air5_park_exchange;
        input.vibration_translation[3].partner = nullptr;
        starts(input, "vibration-translation pair 3 (N2, ): partner has no name");

        for (const double a : {nan, inf}) {
            input = gaslib::air5_park_exchange;
            input.vibration_translation[4].millikan_white.a = a;
            starts(input, "vibration-translation pair 4 (N2, NO): Millikan-White a must be finite");
        }
        // Any finite a is accepted.
        for (const double a : {0.0, -27.32}) {
            input = gaslib::air5_park_exchange;
            input.vibration_translation[4].millikan_white.a = a;
            build(input);
        }
        input = gaslib::air5_park_exchange;
        input.vibration_translation[4].millikan_white.b = inf;
        starts(input, "vibration-translation pair 4 (N2, NO): Millikan-White b must be finite");

        using gaslib::RelaxationCorrection;
        for (const RelaxationCorrection model : {
               RelaxationCorrection::park_1985,
               RelaxationCorrection::park_1993,
               RelaxationCorrection::kim_2021,
             }) {
            for (const double sigma : {nan, inf, 0.0, -1.0e-21}) {
                input = gaslib::air5_park_exchange;
                input.vibration_translation[5].correction.model = model;
                input.vibration_translation[5].correction.cross_section = sigma;
                starts(input,
                       "vibration-translation pair 5 (O2, N): correction cross-section must be");
            }
        }
        for (const double t_ref : {nan, inf, 0.0, -5.0e4}) {
            input = gaslib::air5_park_exchange;
            input.vibration_translation[5].correction.reference_temperature = t_ref;
            starts(
              input,
              "vibration-translation pair 5 (O2, N): correction reference temperature must be");
        }
        for (const double exponent : {nan, inf}) {
            input = gaslib::air5_park_exchange;
            input.vibration_translation[5].correction.model = RelaxationCorrection::kim_2021;
            input.vibration_translation[5].correction.exponent = exponent;
            starts(input,
                   "vibration-translation pair 5 (O2, N): correction exponent must be finite");
        }
        input = gaslib::air5_park_exchange;
        // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange): an invalid model
        input.vibration_translation[5].correction.model = static_cast<RelaxationCorrection>(9);
        starts(input,
               "vibration-translation pair 5 (O2, N): unknown relaxation-time correction model");
        // A model ignores the constants it does not use.
        input = gaslib::air5_park_exchange;
        input.vibration_translation[5].correction = {
          .model = RelaxationCorrection::none,
          .cross_section = nan,
          .reference_temperature = nan,
          .exponent = nan,
        };
        build(input);
        input.vibration_translation[5].correction.model = RelaxationCorrection::park_1985;
        input.vibration_translation[5].correction.cross_section = 1.0e-21;
        build(input);
        input.vibration_translation[5].correction.model = RelaxationCorrection::kim_2021;
        input.vibration_translation[5].correction.exponent = 0.5;
        build(input);

        for (const int reaction : {-1, 5}) {
            input = gaslib::air5_park_exchange;
            input.chemistry_vibration[1].reaction = reaction;
            starts(input,
                   "chemistry-vibration coupling 1 (O2): reaction index must be between 0 and 4");
        }
        input = gaslib::air5_park_exchange;
        input.chemistry_vibration[1].species = "Ar";
        starts(input, "chemistry-vibration coupling 1 (Ar): unknown species \"Ar\"");
        // N2 takes no part in NO + O <=> O2 + N.
        input = gaslib::air5_park_exchange;
        input.chemistry_vibration[0].reaction = 3;
        starts(
          input,
          "chemistry-vibration coupling 0 (N2): species \"N2\" has no net change in reaction 3");
        for (const double d : {nan, inf, 0.0, -1.0}) {
            input = gaslib::air5_park_exchange;
            input.chemistry_vibration[2].dissociation_energy = d;
            starts(input, "chemistry-vibration coupling 2 (NO): dissociation energy must be");
        }
        for (const double theta : {nan, inf, 0.0, -1.0}) {
            input = gaslib::air5_park_exchange;
            input.chemistry_vibration[2].vibrational_temperature = theta;
            starts(input, "chemistry-vibration coupling 2 (NO): vibrational temperature must be");
        }

        auto species = built.species;
        species.n = 4;
        expect(error_message([&] {
                   (void)gaslib::make_air5_park_exchange(species, mech, si_scales);
               }).starts_with("species data"));
        expect(error_message([&] {
                   (void)gaslib::make_energy_exchange(
                     gaslib::air5_park_exchange,
                     std::array<std::string_view, 5> {"N2", "O2", "N", "N2", "NO"},
                     built.species,
                     mech,
                     si_scales);
               }).starts_with("species \"N2\" is named twice"));
    };
}
