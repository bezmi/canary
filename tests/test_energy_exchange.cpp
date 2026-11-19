#include "dual.hpp"
#include "nasa9_blend_offsets.hpp"
#include "oscillator_reference.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/data/air5.hpp>
#include <canary/data/air5_park.hpp>
#include <canary/data/air5_park_exchange.hpp>
#include <canary/energy_exchange.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/species.hpp>
#include <canary/status.hpp>
#include <canary/thermo.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <numbers>
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
using Air = canary::Air5TwoTemperature<double, si>;
using AirScaled = canary::Air5TwoTemperature<double, scaled>;
using AirSingle = canary::Air5<double, si>;
using Mechanism = canary::Air5Park<double, si>;
using Exchange = canary::Air5ParkExchange<double, si>;
using Rho = canary::Density<double, si>;
using Temp = canary::Temperature<double, si>;
using Extents = Kokkos::dextents<std::size_t, 1>;

constexpr std::size_t species_count = 5;
constexpr std::size_t reaction_count = 5;
using Fractions = std::array<double, species_count>;
using PerReaction = std::array<double, reaction_count>;

constexpr double atmosphere = 101325.0;
constexpr double avogadro = canary::constants::avogadro<double>;
constexpr double boltzmann = canary::constants::boltzmann<double>;
// Universal gas constant of the reference implementation [J/(mol K)].
constexpr double reference_gas_constant = 8.31451;

const canary::UnitScales<double, scaled> reference_scales {{
  .mass = 4.65e-26,
  .density = 1.0e-3,
  .temperature = 1000.0,
  .length = 1.0,
}};

template <typename T> auto fractions(const std::array<T, species_count>& y)
{
    return canary::as_kind<canary::kind::mass_fraction>(
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
    .forward = {1.6038551846505644, 3212.233354819928, 364.8961518304201, 1876.3961515915175, 27622.447324252233},
    .backward = {0.7156813177663885, 1.6643762994328262, 39.252078506703185, 42133.758577092616, 13348.95451299561},
    .source = {-399.8739442886589, -1390.9216434110062, -359.36218691421226, 523.6707879227602, 1626.486986691117},
    .total = 83082251486.96964, .vt = 83165692436.3816, .cv = -83440949.4119523,
  },
  {
    .y = mixed, .rho = 0.001, .temp = 20000.0, .t_ve = 5000.0,
    .forward = {70.20781492573974, 1771.7729507439935, 1860.8592867048633, 49.62218884836747, 942.0578300944791},
    .backward = {0.0001482170082127718, 0.0006167115915839476, 0.03215973122482218, 391.3994131968709, 52.965633125358245},
    .source = {-26.87325080115315, -67.63104960871893, 35.696879396965244, 77.70999492459175, -18.90257391168492},
    .total = 8275718148.439633, .vt = 8444383490.442057, .cv = -168665342.00242382,
  },
  {
    .y = dissociated, .rho = 0.01, .temp = 1500.0, .t_ve = 4000.0,
    .forward = {4.199262144860181e-11, 0.029209560525757138, 2.9636622312706282e-05, 0.051112163981894596, 5.4417893045302135e-05},
    .backward = {260.1042856462481, 75.49258341184914, 302.2637324998749, 26293.81998235093, 256709.92677607664},
    .source = {7198.604266620037, -838.9543139161444, -3975.4679854487217, 4520.638551957391, -6904.820519212562},
    .total = -433199166.9341856, .vt = -662418147.7794039, .cv = 229218980.84521827,
  },
  // T_ve = 6000 (1 + 1e-9).
  {
    .y = dissociated, .rho = 0.01, .temp = 6000.0, .t_ve = 6000.000006,
    .forward = {7.539116908104633, 13320.078136616828, 2472.321537415028, 855.141739678003, 2965.7585011250517},
    .backward = {28.198197571134163, 13.235054575841255, 367.34828960500874, 127876.93303414612, 74655.8056053904},
    .source = {2008.8606966408727, -4490.347905687039, -2754.392109199685, 3638.751506671321, 1597.127811574531},
    .total = -3820193514.342525, .vt = -41.28258794582087, .cv = -3820193473.059937,
  },
  {
    .y = mixed, .rho = 0.01, .temp = 8000.0, .t_ve = 8000.0,
    .forward = {592.0905010853626, 55944.841118592776, 28182.706612165155, 1153.8486379368017, 13220.54821709806},
    .backward = {1.0838985513898018, 2.0053804105740167, 38.081042398050045, 43494.89097772273, 17911.16679622102},
    .source = {114.84407014497799, -3144.9661578412865, -247.98893366539914, 2992.8790902658675, 285.23193109584054},
    .total = -22235037984.763607, .vt = 0.0, .cv = -22235037984.763607,
  },
  {
    .y = cold_air, .rho = 0.1, .temp = 5000.0, .t_ve = 1000.0,
    .forward = {3.0006878820365215e-11, 0.13262844379992447, 0.0, 0.0, 0.0},
    .backward = {0.0, 0.0, 0.0, 0.0, 0.0},
    .source = {-8.405946991464189e-13, -0.004243951047465023, 8.405946991464189e-13, 0.004243951047465023, 0.0},
    .total = 105825816246.38803, .vt = 105825816731.82959, .cv = -485.4415485772803,
  },
  // T above the upper limit 50000 K of the mechanism. The reference's
  // stand-alone rates of progress clamp T to 30000 K instead of the
  // mechanism's limits, so they give m_N2 (-(q_f,0 - q_b,0) - (q_f,4 -
  // q_b,4)) = -130 kg/(m^3 s) for N2, not the listed -316 kg/(m^3 s), which
  // uses 50000 K; the species and energy sources are compared, the rates are
  // not.
  {
    .y = mixed, .rho = 0.001, .temp = 60000.0, .t_ve = 10000.0,
    .forward = {3487.584831935751, 9609.83827637619, 45242.813088672, 68.62142173516973, 1191.0631079824614},
    .backward = {1.096568957807923e-05, 6.344260726703739e-05, 0.005035473573938065, 343.4663649334981, 32.29428595713462},
    .source = {-315.78638508101545, -466.8418568147578, 1986.9717420646678, 2375.7827922616684, -3580.126292430564},
    .total = 5863507828.2854805, .vt = 18571668428.810677, .cv = -12708160600.525196,
    .rates_listed = false,
  },
  // T below the lower limit 300 K of the mechanism.
  {
    .y = dissociated, .rho = 0.1, .temp = 250.0, .t_ve = 2000.0,
    .forward = {1.0604418849190528e-51, 2.5404098502789777e-22, 3.3671223960274183e-32, 1.5244638855451262e-22, 9.182412110315756e-47},
    .backward = {7888097.244335068, 1845602.6278664698, 472096.5915748417, 189.76248978040806, 29346583.696194474},
    .source = {1043070.0110594302, 59050.997196615615, -638636.3904774201, 402920.4356990567, -866405.0534776826},
    .total = 4203039636792.5767, .vt = -130923090.81970087, .cv = 4203170559883.3965,
  },
 },
};
// clang-format on

// Change of the coefficient of species i in reaction k.
int coefficient_change(const Mechanism& mech, std::size_t k, std::size_t i)
{
    return mech.reactions()[canary::ReactionIndex {k}].products()[canary::SpeciesIndex {i}] -
           mech.reactions()[canary::ReactionIndex {k}].reactants()[canary::SpeciesIndex {i}];
}

// Rates of progress, mass production rates and energy sources of canary in
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
    const auto mix = canary::mixture(gas, fractions(y));
    const canary::Density<double, U> rho_q {rho};
    const canary::Temperature<double, U> temp_q {temp};
    const canary::Temperature<double, U> t_ve_q {t_ve};
    const auto q = canary::rates_of_progress(mech, gas, mix, rho_q, temp_q, t_ve_q);
    const auto sources =
      canary::two_temperature_sources(mech, exchange, gas, mix, rho_q, temp_q, t_ve_q);
    Values out {};
    for (const auto k : canary::reaction_indices(0, reaction_count)) {
        out.forward[k.value()] = q.forward()[k].value();
        out.backward[k.value()] = q.backward()[k].value();
    }
    for (std::size_t i = 0; i < species_count; ++i) {
        out.production[i] = sources.species[i].value();
    }
    out.total = sources.vibrational_electronic.value();
    out.vt =
      canary::vibration_translation_source(exchange, gas, mix, rho_q, temp_q, t_ve_q).value();
    out.cv = canary::chemistry_vibration_source(exchange, mech, q, temp_q, t_ve_q).value();
    return out;
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
static_assert(canary::HasSpeciesProperties<Air, D, 5>);
static_assert(canary::ReactingModel<Air, double, double, 5>);
static_assert(canary::ReactingModel<AirScaled, D, double, 5>);
static_assert(canary::HasVibrationalElectronicEnergy<Air, double, 5>);
static_assert(canary::HasVibrationalElectronicEnergy<Air, float, 5>);
static_assert(canary::HasVibrationalElectronicEnergy<AirScaled, D, 5>);
static_assert(!canary::HasVibrationalElectronicEnergy<Air, double, 4>);
static_assert(!canary::HasVibrationalElectronicEnergy<AirSingle, double, 5>);
static_assert(canary::RelaxingModel<Air, double, double, 5>);
static_assert(canary::RelaxingModel<Air, float, double, 5>);
static_assert(canary::RelaxingModel<Air, D, double, 5>);
static_assert(!canary::RelaxingModel<AirSingle, double, double, 5>);
static_assert(std::is_trivially_copyable_v<Exchange>);
static_assert(std::is_trivially_copyable_v<canary::Air5ParkExchange<float, scaled>>);
static_assert(std::is_trivially_copyable_v<canary::EnergyExchange<double, si, 5, 0, 0>>);
static_assert(std::is_trivially_copyable_v<canary::TwoTemperatureSources<D, si, 5>>);
static_assert(canary::air5_park.reactions.size() == 5);
static_assert(canary::air5_park_exchange.vibration_translation.size() == 15);

// A built exchange is constructed only by its validated factory and is read-only:
// record lists are const references, species data are const-element views
// indexed only by SpeciesIndex, and borrowing from a temporary does not compile.
template <typename X>
concept BorrowsExchangeMass = requires { std::declval<X>().mass(); };
template <typename X>
concept BorrowsExchangeRecords = requires { std::declval<X>().vibration_translation(); };
template <typename X>
concept ReadsExchangeBoltzmann = requires { std::declval<X>().boltzmann(); };
template <typename X, typename Index>
concept ExchangeMassIndexedBy = requires(const X& x, Index i) { x.mass()[i]; };
static_assert(!std::is_default_constructible_v<Exchange>);
static_assert(std::is_copy_constructible_v<Exchange> && std::is_copy_assignable_v<Exchange>);
static_assert(BorrowsExchangeMass<const Exchange&> && BorrowsExchangeRecords<const Exchange&>);
static_assert(!BorrowsExchangeMass<Exchange&&> && !BorrowsExchangeMass<const Exchange&&>);
static_assert(!BorrowsExchangeRecords<Exchange&&> && !BorrowsExchangeRecords<const Exchange&&>);
static_assert(ReadsExchangeBoltzmann<Exchange&&>);
static_assert(ExchangeMassIndexedBy<Exchange, canary::SpeciesIndex>);
static_assert(!ExchangeMassIndexedBy<Exchange, std::size_t>);
static_assert(!ExchangeMassIndexedBy<Exchange, canary::ModeIndex>);
static_assert(std::is_const_v<decltype(std::declval<const Exchange&>().mass())::element_type>);
static_assert(std::is_same_v<decltype(std::declval<const Exchange&>().vibration_translation()),
                             const std::array<canary::VibrationTranslationPair<double, si>, 15>&>);
static_assert(std::is_same_v<decltype(std::declval<const Exchange&>().chemistry_vibration()),
                             const std::array<canary::ChemistryVibrationCoupling<double, si>, 3>&>);
static_assert(noexcept(std::declval<const Exchange&>().mass()) &&
              noexcept(std::declval<const Exchange&>().vibration_translation()) &&
              noexcept(std::declval<const Exchange&>().chemistry_vibration()) &&
              noexcept(std::declval<const Exchange&>().boltzmann()));

void test_energy_exchange()
{
    using namespace boost::ut;
    const canary::UnitScales<double, si> si_scales {};
    const auto built = canary::make_air5_two_temperature<double>(si_scales);
    const Air& gas = built.model;
    const Mechanism mech =
      canary::make_air5_park(si_scales, canary::Pressure<double, si> {atmosphere}, gas);
    const Exchange exchange = canary::make_air5_park_exchange(built.species, mech, si_scales);
    const double ru = canary::constants::universal_gas_constant<double>;

    "energy exchange: the generated mechanism and exchange hold their constants"_test = [&] {
        expect(mech.reactions()[canary::ReactionIndex {0}].forward.translational_exponent == 0.5_d);
        expect(mech.reactions()[canary::ReactionIndex {2}].forward.translational_exponent == 0.5_d);
        expect(mech.reactions()[canary::ReactionIndex {3}].forward.translational_exponent == 1.0_d);
        expect(mech.reactions()[canary::ReactionIndex {4}].forward.translational_exponent == 1.0_d);
        // N2 + M <=> N + N + M: order 2.
        expect(relative(mech.reactions()[canary::ReactionIndex {0}].forward.A,
                        7.0e15 / avogadro,
                        1.0e-15));
        expect(mech.temperature_upper.value() == 50000.0_d);

        // Pair 2 is N2-N2 and pair 12 NO-N.
        const auto& n2_n2 = exchange.vibration_translation()[2];
        expect(n2_n2.species == canary::SpeciesIndex {0} &&
               n2_n2.partner == canary::SpeciesIndex {0});
        expect(exchange.vibration_translation()[12].species == canary::SpeciesIndex {4});
        expect(exchange.vibration_translation()[12].partner == canary::SpeciesIndex {2});
        expect(n2_n2.millikan_white_a == 221.0_d);
        expect(relative(n2_n2.millikan_white_offset,
                        std::log(atmosphere) - (221.0 * 0.029) - 18.42,
                        1.0e-15));
        // 1 / K = 1 / (sigma_ref T_ref^2 sqrt(8 k_B / (pi mu))), mu = m_N2 / 2.
        const double mu = built.species.mass[0].value() / 2.0;
        const double k =
          3.0e-21 * 50000.0 * 50000.0 * std::sqrt(8.0 * boltzmann / (std::numbers::pi * mu));
        expect(relative(n2_n2.correction_coefficient, 1.0 / k, 1.0e-15));
        expect(n2_n2.correction_exponent == 1.5_d);

        const auto& n2 = exchange.chemistry_vibration()[0];
        expect(n2.reaction == canary::ReactionIndex {0} && n2.species == canary::SpeciesIndex {0});
        expect(relative(n2.dissociation_temperature.value(), 945360.0 / ru, 1.0e-15));
        expect(n2.vibrational_temperature.value() == 3393.44_d);
        expect(exchange.chemistry_vibration()[2].reaction == canary::ReactionIndex {2});
        expect(exchange.chemistry_vibration()[2].species == canary::SpeciesIndex {4});
        expect(exchange.boltzmann().value() == boltzmann);
    };

    "energy exchange: rates, species sources and energy sources equal the reference"_test = [&] {
        const double r = ru / reference_gas_constant;
        for (const ReferenceState& s : reference_states) {
            const Values got = evaluate(mech, exchange, gas, s.y, s.rho, s.temp, s.t_ve);
            // The reference joins NASA9 segments by blending their coefficients;
            // its equilibrium constants use T clamped to the mechanism limits.
            // The Gibbs energies of canary differ by the blend offsets, which
            // scale the backward rates (see nasa9_blend_offsets.hpp).
            const double t_c = std::clamp(s.temp, 300.0, 50000.0);
            const auto off = canary::test::blend_offsets(canary::air5_species, t_c);
            PerReaction forward {};
            PerReaction backward {};
            PerReaction factor {};
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                forward[k.value()] = avogadro * s.forward[k.value()];
                std::array<int, species_count> dnu {};
                for (std::size_t i = 0; i < species_count; ++i) {
                    dnu[i] = coefficient_change(mech, k.value(), i);
                }
                factor[k.value()] = canary::test::blend_rate_factor(off, dnu, t_c);
                backward[k.value()] = avogadro * s.backward[k.value()] *
                                      std::pow(r, mech.reactions()[k].coefficient_change) *
                                      factor[k.value()];
                if (!s.rates_listed) { continue; }
                // Forward rates contain no constant that differs between
                // the implementations; backward rates contain R_u through
                // the Gibbs energies, corrected above, and differ by
                // rounding amplified by the exponent of K_c (up to about
                // 50).
                expect(relative(got.forward[k.value()], forward[k.value()], 1.0e-12))
                  << s.temp << k.value() << got.forward[k.value()] << forward[k.value()];
                expect(relative(got.backward[k.value()], backward[k.value()], 1.0e-10))
                  << s.temp << k.value() << got.backward[k.value()] << backward[k.value()];
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
                    // Without listed rates, the reference backward rate is
                    // canary's divided by the blend factor.
                    expected += m * dnu * got.backward[k] * (1.0 / factor[k] - 1.0);
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
              const auto mix = canary::mixture(gas, fractions(mixed));
              const double rho = 0.01;
              const auto q =
                canary::rates_of_progress(mech, gas, mix, Rho {rho}, Temp {temp}, Temp {t_ve});
              // Below the lower limit T_l = 300 K; T_ve is not clamped.
              const double t_l = std::max(temp, 300.0);
              const double t_c = std::sqrt(t_l * t_ve);
              std::array<double, species_count> n {};
              for (std::size_t i = 0; i < species_count; ++i) {
                  n[i] = rho * mixed[i] / built.species.mass[i].value();
              }
              for (const auto k : canary::reaction_indices(0, 3)) {
                  const auto rate_records = mech.reactions();
                  const auto& rate = rate_records[k].forward;
                  double third_body = 0;
                  for (std::size_t i = 0; i < species_count; ++i) {
                      third_body +=
                        mech.reactions()[k].efficiencies()[canary::SpeciesIndex {i}] * n[i];
                  }
                  // Reactions 0 to 2 dissociate species 0, 1 and 4.
                  const std::size_t molecule = k.value() == 2 ? 4 : k.value();
                  const double expected = rate.A * std::pow(t_c, rate.n) *
                                          std::exp(-rate.activation_temperature.value() / t_c) *
                                          third_body * n[molecule];
                  // The exponent theta / T_c is at most about 410 here, so a
                  // few ulp in T_c give up to about 2e-13 in k_f.
                  expect(relative(q.forward()[k].value(), expected, 1.0e-12))
                    << temp << t_ve << k.value() << q.forward()[k].value() << expected;
              }
          }
      };

    "energy exchange: rates at T_ve = T equal the single-temperature rates"_test = [&] {
        for (const ReferenceState& s : reference_states) {
            if (s.temp < 300.0 || s.temp > 50000.0) { continue; }
            const auto mix = canary::mixture(gas, fractions(s.y));
            const auto one = canary::rates_of_progress(mech, gas, mix, Rho {s.rho}, Temp {s.temp});
            const auto two =
              canary::rates_of_progress(mech, gas, mix, Rho {s.rho}, Temp {s.temp}, Temp {s.temp});
            const auto omega_one =
              canary::production_rates(mech, gas, mix, Rho {s.rho}, Temp {s.temp});
            const auto omega_two =
              canary::production_rates(mech, gas, mix, Rho {s.rho}, Temp {s.temp}, Temp {s.temp});
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                expect(one.forward()[k] == two.forward()[k] &&
                       one.backward()[k] == two.backward()[k])
                  << s.temp << k.value();
            }
            for (std::size_t i = 0; i < species_count; ++i) {
                expect(omega_one[i] == omega_two[i]) << s.temp << i;
            }
        }
    };

    "energy exchange: rates with s = 1 and equilibrium backward rates do not depend on T_ve"_test =
      [&] {
          const auto mix = canary::mixture(gas, fractions(dissociated));
          const Rho rho {0.01};
          for (const double temp : {250.0, 1500.0, 8000.0, 60000.0}) {
              const auto a =
                canary::rates_of_progress(mech, gas, mix, rho, Temp {temp}, Temp {1000.0});
              const auto b =
                canary::rates_of_progress(mech, gas, mix, rho, Temp {temp}, Temp {9000.0});
              for (const auto k : canary::reaction_indices(0, reaction_count)) {
                  expect(a.backward()[k] == b.backward()[k]) << temp << k.value();
              }
              expect(
                a.forward()[canary::ReactionIndex {3}] == b.forward()[canary::ReactionIndex {3}] &&
                a.forward()[canary::ReactionIndex {4}] == b.forward()[canary::ReactionIndex {4}])
                << temp;
              expect(
                !(a.forward()[canary::ReactionIndex {0}] == b.forward()[canary::ReactionIndex {0}]))
                << temp;
              // Derivatives with respect to T_ve.
              const std::array<D, species_count> y_d {0.3, 0.1, 0.3, 0.25, 0.05};
              const auto q = canary::rates_of_progress(mech,
                                                       gas,
                                                       canary::mixture(gas, fractions(y_d)),
                                                       canary::Density<D, si> {D(0.01)},
                                                       canary::Temperature<D, si> {D(temp)},
                                                       canary::Temperature<D, si> {D(3000.0, 1.0)});
              for (const auto k : canary::reaction_indices(0, reaction_count)) {
                  expect(q.backward()[k].value().dual() == 0.0) << temp << k.value();
              }
              expect(q.forward()[canary::ReactionIndex {3}].value().dual() == 0.0 &&
                     q.forward()[canary::ReactionIndex {4}].value().dual() == 0.0);
              expect(q.forward()[canary::ReactionIndex {0}].value().dual() > 0.0);
          }
      };

    "energy exchange: below the lower limit T is the limit and T_ve is not clamped"_test = [&] {
        const auto mix = canary::mixture(gas, fractions(dissociated));
        const Rho rho {0.1};
        for (const double t_ve : {250.0, 2000.0}) {
            const auto below =
              canary::rates_of_progress(mech, gas, mix, rho, Temp {250.0}, Temp {t_ve});
            const auto limit =
              canary::rates_of_progress(mech, gas, mix, rho, Temp {300.0}, Temp {t_ve});
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                expect(below.forward()[k] == limit.forward()[k] &&
                       below.backward()[k] == limit.backward()[k])
                  << t_ve << k.value();
            }
        }
        const auto clamped =
          canary::rates_of_progress(mech, gas, mix, rho, Temp {250.0}, Temp {300.0});
        const auto unclamped =
          canary::rates_of_progress(mech, gas, mix, rho, Temp {250.0}, Temp {250.0});
        expect(unclamped.forward()[canary::ReactionIndex {0}].value() <
               clamped.forward()[canary::ReactionIndex {0}].value());
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
            const auto mix = canary::mixture(gas, fractions(y));
            const Rho rho {0.01};
            for (const double temp : {500.0, 3000.0, 9000.0, 25000.0}) {
                const auto same = canary::vibration_translation_source(exchange,
                                                                       gas,
                                                                       mix,
                                                                       rho,
                                                                       Temp {temp},
                                                                       Temp {temp});
                expect(same.value() == 0.0) << temp;
                for (const double factor : {0.5, 0.99, 1.01, 2.0}) {
                    const double t_ve = factor * temp;
                    const double q = canary::vibration_translation_source(exchange,
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
          canary::EnergyExchangeInput<1, 0> input {};
          input.vibration_translation[0] = canary::air5_park_exchange.vibration_translation[1];
          const auto absent = canary::make_energy_exchange(input,
                                                           canary::detail::air5_park_species_names,
                                                           built.species,
                                                           mech,
                                                           si_scales);
          const auto mix = canary::mixture(gas, fractions(cold_air));
          const double q = canary::vibration_translation_source(absent,
                                                                gas,
                                                                mix,
                                                                Rho {0.1},
                                                                Temp {5000.0},
                                                                Temp {1000.0})
                             .value();
          expect(q == 0.0) << q;

          const canary::EnergyExchangeInput<0, 0> nothing {};
          const auto empty = canary::make_energy_exchange(nothing,
                                                          canary::detail::air5_park_species_names,
                                                          built.species,
                                                          mech,
                                                          si_scales);
          const auto sources = canary::two_temperature_sources(mech,
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
          const auto mix = canary::mixture(gas, fractions(y));
          const double m = built.species.mass[0].value();
          for (const double rho : {1.0e-3, 0.1}) {
              for (const double temp : {2000.0, 6000.0, 15000.0}) {
                  const double t_ve = temp * (1.0 - 1.0e-4);
                  const double q = canary::vibration_translation_source(exchange,
                                                                        gas,
                                                                        mix,
                                                                        Rho {rho},
                                                                        Temp {temp},
                                                                        Temp {t_ve})
                                     .value();
                  const auto e_t = canary::species_vibrational_electronic_energy(gas, Temp {temp});
                  const auto e_ve = canary::species_vibrational_electronic_energy(gas, Temp {t_ve});
                  const double rate = q / (rho * (e_t[0].value() - e_ve[0].value()));
                  // tau from the formulas, with x_q = 1 and p = n k_B T.
                  const double n = rho / m;
                  const double p = n * boltzmann * temp;
                  const double tau_mw =
                    (atmosphere / p) * std::exp((221.0 * (std::cbrt(1.0 / temp) - 0.029)) - 18.42);
                  const double sigma = 3.0e-21 * (50000.0 / temp) * (50000.0 / temp);
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
        using canary::RelaxationCorrection;
        const Fractions y {1.0, 0.0, 0.0, 0.0, 0.0};
        const auto mix = canary::mixture(gas, fractions(y));
        const double m = built.species.mass[0].value();
        // sigma_0, T_ref and omega of each model; a negative Millikan-White a
        // with the last.
        const std::array<canary::RelaxationCorrectionInput, 4> corrections {
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
            canary::EnergyExchangeInput<1, 0> input {};
            input.vibration_translation[0] = {
              .species = "N2",
              .partner = "N2",
              .millikan_white = {.a = a, .b = b},
              .correction = corrections[c],
            };
            const auto single =
              canary::make_energy_exchange(input,
                                           canary::detail::air5_park_species_names,
                                           built.species,
                                           mech,
                                           si_scales);
            for (const double temp : {2000.0, 6000.0, 15000.0}) {
                const double rho = 0.01;
                const double t_ve = temp * (1.0 - 1.0e-4);
                const double q = canary::vibration_translation_source(single,
                                                                      gas,
                                                                      mix,
                                                                      Rho {rho},
                                                                      Temp {temp},
                                                                      Temp {t_ve})
                                   .value();
                const auto e_t = canary::species_vibrational_electronic_energy(gas, Temp {temp});
                const auto e_ve = canary::species_vibrational_electronic_energy(gas, Temp {t_ve});
                const double rate = q / (rho * (e_t[0].value() - e_ve[0].value()));
                const double n = rho / m;
                const double tau_mw = (atmosphere / (n * boltzmann * temp)) *
                                      std::exp((a * (std::cbrt(1.0 / temp) - b)) - 18.42);
                const canary::RelaxationCorrectionInput& k = corrections[c];
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

    "energy exchange: an absent partner gives zero when both relaxation-time terms underflow"_test =
      [&] {
          // A negative Millikan-White a without a correction: exp(a' T^(-1/3) + c)
          // underflows at T = 1e-6 K, and the correction term is zero.
          canary::EnergyExchangeInput<1, 0> input {};
          input.vibration_translation[0] = {
            .species = "N2",
            .partner = "O",
            .millikan_white = {.a = -27.32, .b = 0.1465},
            .correction = {.model = canary::RelaxationCorrection::none,
                           .cross_section = 0.0,
                           .reference_temperature = 0.0,
                           .exponent = 0.0},
          };
          const auto single = canary::make_energy_exchange(input,
                                                           canary::detail::air5_park_species_names,
                                                           built.species,
                                                           mech,
                                                           si_scales);
          const double t = 1.0e-6;
          const double a = single.vibration_translation()[0].millikan_white_a;
          const double c = single.vibration_translation()[0].millikan_white_offset;
          expect(std::exp((a / std::cbrt(t)) + c) == 0.0_d);
          const std::array<double, species_count> n {1.0e20, 0.0, 0.0, 0.0, 0.0};
          const auto inverse_tau =
            canary::detail::inverse_relaxation_times(single.vibration_translation(),
                                                     single.boltzmann(),
                                                     n,
                                                     t,
                                                     t);
          for (const double value : inverse_tau) { expect(value == 0.0_d) << value; }
      };

    "energy exchange: relaxation at fixed energy reaches the equilibrium temperature"_test = [&] {
        const Fractions y {1.0, 0.0, 0.0, 0.0, 0.0};
        const auto mix = canary::mixture(gas, fractions(y));
        const double rho = 0.01;
        const std::array t_ve_start {1000.0};
        const auto start =
          canary::thermo_state(gas,
                               mix,
                               Rho {rho},
                               Temp {8000.0},
                               canary::as_kind<canary::kind::mode_temperature, si>(
                                 Kokkos::mdspan<const double, Extents>(t_ve_start.data(), 1)));
        const double rho_e = rho * start.state.e.value();
        const auto equilibrium =
          canary::thermo_state(gas, mix, Rho {rho}, canary::EnergyDensity<double, si> {rho_e});
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
            const double q = canary::vibration_translation_source(exchange,
                                                                  gas,
                                                                  mix,
                                                                  Rho {rho},
                                                                  Temp {temp},
                                                                  Temp {t_ve})
                               .value();
            const auto e_t = canary::species_vibrational_electronic_energy(gas, Temp {temp});
            const auto e_ve = canary::species_vibrational_electronic_energy(gas, Temp {t_ve});
            const double tau = rho * (e_t[0].value() - e_ve[0].value()) / q;
            rho_e_ve[0] += 0.01 * tau * q;
            const auto next =
              canary::thermo_state(gas,
                                   mix,
                                   Rho {rho},
                                   canary::EnergyDensity<double, si> {rho_e},
                                   canary::as_kind<canary::kind::mode_energy_density, si>(
                                     Kokkos::mdspan<const double, Extents>(rho_e_ve.data(), 1)));
            converged = converged && next.status == canary::Status::Valid;
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
        const auto& references = canary::test::oscillator_references;
        expect(references.size() == exchange.chemistry_vibration().size());
        for (std::size_t i = 0; i < references.size(); ++i) {
            const auto& coupling = exchange.chemistry_vibration().at(i);
            const auto& reference = references[i];
            const double theta_v = coupling.vibrational_temperature.value();
            const double theta_d = coupling.dissociation_temperature.value();
            expect(theta_v == reference.vibrational_temperature);
            expect(theta_d == reference.dissociation_temperature);
            for (const auto& sample : reference.samples) {
                const double x = sample.inverse_temperature;
                const double g =
                  canary::detail::truncated_oscillator_energy(theta_v, theta_d, x).value;
                const double precise = sample.energy_over_boltzmann;
                expect(relative(g, precise, 1.0e-12)) << theta_d << x << g << precise;
            }
            const auto g0 = canary::detail::truncated_oscillator_energy(theta_v, theta_d, 0.0);
            expect(g0.value == (theta_d - theta_v) / 2.0);
        }
    };

    "energy exchange: G_vanish equals G_appear at T = T_ve and Q_CV < 0 without products"_test =
      [&] {
          for (const ReferenceState& s : reference_states) {
              const auto mix = canary::mixture(gas, fractions(s.y));
              const auto q = canary::rates_of_progress(mech,
                                                       gas,
                                                       mix,
                                                       Rho {s.rho},
                                                       Temp {s.temp},
                                                       Temp {s.temp});
              const double cv =
                canary::chemistry_vibration_source(exchange, mech, q, Temp {s.temp}, Temp {s.temp})
                  .value();
              // With G_vanish = G_appear = k_B (theta_D - theta_v) / 2,
              // Q_CV = sum_c dnu (q_f - q_b) G_appear.
              double expected = 0;
              double scale = 0;
              for (const auto& c : exchange.chemistry_vibration()) {
                  const double g =
                    boltzmann *
                    (c.dissociation_temperature.value() - c.vibrational_temperature.value()) / 2.0;
                  const double dnu =
                    coefficient_change(mech, c.reaction.value(), c.species.value());
                  expected +=
                    dnu * (q.forward()[c.reaction].value() - q.backward()[c.reaction].value()) * g;
                  scale += std::abs(dnu) *
                           (q.forward()[c.reaction].value() + q.backward()[c.reaction].value()) * g;
              }
              expect(std::abs(cv - expected) <= 1.0e-14 * scale) << s.temp << cv << expected;
          }
          const auto mix = canary::mixture(gas, fractions(cold_air));
          for (const auto& [temp, t_ve] : {
                 std::array {5000.0, 1000.0},
                 std::array {3000.0, 3000.0},
                 std::array {2000.0, 9000.0},
               }) {
              const auto q =
                canary::rates_of_progress(mech, gas, mix, Rho {0.1}, Temp {temp}, Temp {t_ve});
              expect(q.backward()[canary::ReactionIndex {0}].value() == 0.0 &&
                     q.backward()[canary::ReactionIndex {1}].value() == 0.0);
              const double cv =
                canary::chemistry_vibration_source(exchange, mech, q, Temp {temp}, Temp {t_ve})
                  .value();
              expect(cv < 0.0) << temp << t_ve << cv;
          }
      };

    "energy exchange: scaled units give the SI results after conversion"_test = [&] {
        const auto built_s = canary::make_air5_two_temperature<double>(reference_scales);
        const auto mech_s = canary::make_air5_park(
          reference_scales,
          canary::convert<scaled>(canary::Pressure<double, si> {atmosphere}, reference_scales),
          built_s.model);
        const auto exchange_s =
          canary::make_air5_park_exchange(built_s.species, mech_s, reference_scales);
        for (const ReferenceState& s : reference_states) {
            const auto mix = canary::mixture(gas, fractions(s.y));
            const auto mix_s = canary::mixture(built_s.model, fractions(s.y));
            const Rho rho {s.rho};
            const Temp temp {s.temp};
            const Temp t_ve {s.t_ve};
            const auto rho_s = canary::convert<scaled>(rho, reference_scales);
            const auto temp_s = canary::convert<scaled>(temp, reference_scales);
            const auto t_ve_s = canary::convert<scaled>(t_ve, reference_scales);
            const auto q = canary::rates_of_progress(mech, gas, mix, rho, temp, t_ve);
            const auto q_s =
              canary::rates_of_progress(mech_s, built_s.model, mix_s, rho_s, temp_s, t_ve_s);
            // As for one temperature: the Gibbs energies of the scaled thermo
            // model differ from the SI ones by up to 6e-13 relative, times
            // m / (k_B T) in K_c.
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                const double f = canary::convert<si>(q_s.forward()[k], reference_scales).value();
                const double b = canary::convert<si>(q_s.backward()[k], reference_scales).value();
                expect(relative(f, q.forward()[k].value(), 1.0e-12)) << s.temp << k.value();
                expect(relative(b, q.backward()[k].value(), 1.0e-10)) << s.temp << k.value();
            }
            // The SI rates converted, so that only the exchange is compared.
            canary::RatesOfProgress<double, scaled, reaction_count> q_converted {};
            for (const auto k : canary::reaction_indices(0, reaction_count)) {
                q_converted.forward()[k] =
                  canary::convert<scaled>(q.forward()[k], reference_scales);
                q_converted.backward()[k] =
                  canary::convert<scaled>(q.backward()[k], reference_scales);
            }
            const double cv =
              canary::chemistry_vibration_source(exchange, mech, q, temp, t_ve).value();
            const double cv_s =
              canary::convert<si>(
                canary::chemistry_vibration_source(exchange_s, mech_s, q_converted, temp_s, t_ve_s),
                reference_scales)
                .value();
            expect(relative(cv_s, cv, 1.0e-12)) << s.temp << cv_s << cv;
            // Q_VT is a difference of energies that cancel when T_ve is
            // close to T; states with |T - T_ve| >= 1000 K.
            if (std::abs(s.temp - s.t_ve) < 1000.0) { continue; }
            const double vt =
              canary::vibration_translation_source(exchange, gas, mix, rho, temp, t_ve).value();
            const double vt_s =
              canary::convert<si>(canary::vibration_translation_source(exchange_s,
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
                const auto mix = canary::mixture(gas, fractions(y));
                const canary::Density<T, si> rho {T(s.rho)};
                const canary::Temperature<T, si> t {temp};
                const canary::Temperature<T, si> tv {t_ve};
                const auto q = canary::rates_of_progress(mech, gas, mix, rho, t, tv);
                // The sources, then the size of the terms each is summed
                // from.
                std::array<T, 2 * (species_count + 2)> out {};
                const auto omega = canary::production_rates(mech, gas, mix, rho, t, tv);
                for (std::size_t i = 0; i < species_count; ++i) { out[i] = omega[i].value(); }
                out[species_count] =
                  canary::vibration_translation_source(exchange, gas, mix, rho, t, tv).value();
                out[species_count + 1] =
                  canary::chemistry_vibration_source(exchange, mech, q, t, tv).value();
                const std::size_t scales = species_count + 2;
                for (const auto k : canary::reaction_indices(0, reaction_count)) {
                    const T both = q.forward()[k].value() + q.backward()[k].value();
                    for (std::size_t i = 0; i < species_count; ++i) {
                        const double m = built.species.mass[i].value();
                        out[scales + i] =
                          out[scales + i] +
                          T(m * std::abs(coefficient_change(mech, k.value(), i)) * both);
                    }
                }
                out[scales + species_count] =
                  out[species_count] < T(0.0) ? T(-out[species_count]) : out[species_count];
                for (const auto& coupling : exchange.chemistry_vibration()) {
                    out[scales + species_count + 1] =
                      out[scales + species_count + 1] +
                      T(boltzmann * coupling.dissociation_temperature.value() *
                        (q.forward()[coupling.reaction].value() +
                         q.backward()[coupling.reaction].value()));
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

    "energy exchange: oscillator temperatures remain ordered after scaling and rounding"_test =
      [&] {
          for (const double theta_d : {1.0, 2.0}) {
              auto input = canary::air5_park_exchange;
              input.chemistry_vibration[0].vibrational_temperature = 2.0;
              input.chemistry_vibration[0].dissociation_energy = theta_d * ru;
              expect(
                error_message([&] {
                    (void)canary::make_energy_exchange(input,
                                                       canary::detail::air5_park_species_names,
                                                       built.species,
                                                       mech,
                                                       si_scales);
                })
                  .starts_with(
                    "chemistry-vibration coupling 0 (N2): dissociation temperature must exceed"));
          }
          const auto check_float = [&]<typename U>(const canary::UnitScales<float, U>& scales) {
              const auto gas_f = canary::make_air5<float>(scales);
              const auto mech_f = canary::make_air5_park(scales);
              auto input = canary::air5_park_exchange;
              input.chemistry_vibration[0].vibrational_temperature = 1000.0;
              input.chemistry_vibration[0].dissociation_energy = 1000.00001 * ru;
              expect(error_message([&] {
                         (void)canary::make_energy_exchange(input,
                                                            canary::detail::air5_park_species_names,
                                                            gas_f.species,
                                                            mech_f,
                                                            scales);
                     })
                       .starts_with("chemistry-vibration coupling 0 (N2): stored dissociation "
                                    "temperature must exceed"));
              input.chemistry_vibration[0].dissociation_energy = 1001.0 * ru;
              const auto accepted =
                canary::make_energy_exchange(input,
                                             canary::detail::air5_park_species_names,
                                             gas_f.species,
                                             mech_f,
                                             scales);
              expect(accepted.chemistry_vibration()[0].dissociation_temperature.value() >
                     accepted.chemistry_vibration()[0].vibrational_temperature.value());
          };
          check_float(canary::UnitScales<float, si> {});
          check_float(canary::UnitScales<float, scaled> {
            {.mass = 4.65e-26F, .density = 1.0e-3F, .temperature = 1000.0F, .length = 1.0F}});
      };

    "energy exchange: the builder rejects invalid input with the entry in the message"_test = [&] {
        using Input = canary::EnergyExchangeInput<15, 3>;
        const auto build = [&](const Input& input) {
            (void)canary::make_energy_exchange(input,
                                               canary::detail::air5_park_species_names,
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

        Input input = canary::air5_park_exchange;
        input.vibration_translation[3].species = "Ar";
        starts(input, "vibration-translation pair 3 (Ar, O2): unknown species \"Ar\"");
        input = canary::air5_park_exchange;
        input.vibration_translation[3].partner = "Ar";
        starts(input, "vibration-translation pair 3 (N2, Ar): unknown partner \"Ar\"");
        input = canary::air5_park_exchange;
        input.vibration_translation[3].partner = nullptr;
        starts(input, "vibration-translation pair 3 (N2, ): partner has no name");

        for (const double a : {nan, inf}) {
            input = canary::air5_park_exchange;
            input.vibration_translation[4].millikan_white.a = a;
            starts(input, "vibration-translation pair 4 (N2, NO): Millikan-White a must be finite");
        }
        // Any finite a is accepted.
        for (const double a : {0.0, -27.32}) {
            input = canary::air5_park_exchange;
            input.vibration_translation[4].millikan_white.a = a;
            build(input);
        }
        input = canary::air5_park_exchange;
        input.vibration_translation[4].millikan_white.b = inf;
        starts(input, "vibration-translation pair 4 (N2, NO): Millikan-White b must be finite");

        using canary::RelaxationCorrection;
        for (const RelaxationCorrection model : {
               RelaxationCorrection::park_1985,
               RelaxationCorrection::park_1993,
               RelaxationCorrection::kim_2021,
             }) {
            for (const double sigma : {nan, inf, 0.0, -1.0e-21}) {
                input = canary::air5_park_exchange;
                input.vibration_translation[5].correction.model = model;
                input.vibration_translation[5].correction.cross_section = sigma;
                starts(input,
                       "vibration-translation pair 5 (O2, N): correction cross-section must be");
            }
        }
        for (const double t_ref : {nan, inf, 0.0, -5.0e4}) {
            input = canary::air5_park_exchange;
            input.vibration_translation[5].correction.reference_temperature = t_ref;
            starts(
              input,
              "vibration-translation pair 5 (O2, N): correction reference temperature must be");
        }
        for (const double exponent : {nan, inf}) {
            input = canary::air5_park_exchange;
            input.vibration_translation[5].correction.model = RelaxationCorrection::kim_2021;
            input.vibration_translation[5].correction.exponent = exponent;
            starts(input,
                   "vibration-translation pair 5 (O2, N): correction exponent must be finite");
        }
        input = canary::air5_park_exchange;
        // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange): an invalid model
        input.vibration_translation[5].correction.model = static_cast<RelaxationCorrection>(9);
        starts(input,
               "vibration-translation pair 5 (O2, N): unknown relaxation-time correction model");
        // A model ignores the constants it does not use.
        input = canary::air5_park_exchange;
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
            input = canary::air5_park_exchange;
            input.chemistry_vibration[1].reaction = reaction;
            starts(input,
                   "chemistry-vibration coupling 1 (O2): reaction index must be between 0 and 4");
        }
        input = canary::air5_park_exchange;
        input.chemistry_vibration[1].species = "Ar";
        starts(input, "chemistry-vibration coupling 1 (Ar): unknown species \"Ar\"");
        // N2 takes no part in NO + O <=> O2 + N.
        input = canary::air5_park_exchange;
        input.chemistry_vibration[0].reaction = 3;
        starts(
          input,
          "chemistry-vibration coupling 0 (N2): species \"N2\" has no net change in reaction 3");
        for (const double d : {nan, inf, 0.0, -1.0}) {
            input = canary::air5_park_exchange;
            input.chemistry_vibration[2].dissociation_energy = d;
            starts(input, "chemistry-vibration coupling 2 (NO): dissociation energy must be");
        }
        for (const double theta : {nan, inf, 0.0, -1.0}) {
            input = canary::air5_park_exchange;
            input.chemistry_vibration[2].vibrational_temperature = theta;
            starts(input, "chemistry-vibration coupling 2 (NO): vibrational temperature must be");
        }

        auto species = built.species;
        species.n = 4;
        expect(error_message([&] {
                   (void)canary::make_air5_park_exchange(species, mech, si_scales);
               }).starts_with("species data"));
        expect(error_message([&] {
                   (void)canary::make_energy_exchange(
                     canary::air5_park_exchange,
                     std::array<std::string_view, 5> {"N2", "O2", "N", "N2", "NO"},
                     built.species,
                     mech,
                     si_scales);
               }).starts_with("species \"N2\" is named twice"));
    };
    "exchange construction rejects duplicate physical entries"_test = [&] {
        const auto build = [&](const auto& input) {
            return canary::make_energy_exchange(input,
                                                canary::detail::air5_park_species_names,
                                                built.species,
                                                mech,
                                                si_scales);
        };
        auto input = canary::air5_park_exchange;
        input.vibration_translation[1] = input.vibration_translation[0];
        expect(error_message([&] { (void)build(input); }).find("duplicates entry 0") !=
               std::string::npos);
        input.vibration_translation[1].correction.model = canary::RelaxationCorrection::none;
        expect(error_message([&] { (void)build(input); }).find("duplicates entry 0") !=
               std::string::npos);
        input = canary::air5_park_exchange;
        input.chemistry_vibration[1] = input.chemistry_vibration[0];
        expect(error_message([&] { (void)build(input); }).find("duplicates entry 0") !=
               std::string::npos);
        canary::EnergyExchangeInput<2, 0> directed {};
        directed.vibration_translation[0] = canary::air5_park_exchange.vibration_translation[3];
        directed.vibration_translation[1] = directed.vibration_translation[0];
        directed.vibration_translation[1].species = directed.vibration_translation[0].partner;
        directed.vibration_translation[1].partner = directed.vibration_translation[0].species;
        expect(nothrow([&] { (void)build(directed); }));
        canary::EnergyExchangeInput<0, 2> reactions {};
        reactions.chemistry_vibration[0] = canary::air5_park_exchange.chemistry_vibration[0];
        reactions.chemistry_vibration[1] = reactions.chemistry_vibration[0];
        reactions.chemistry_vibration[1].reaction = 4;
        expect(nothrow([&] { (void)build(reactions); }));
    };

    "exchange construction checks stored masses and correction exponents"_test = [] {
        const auto check = []<typename P, typename U>(const canary::UnitScales<P, U>& scales) {
            const auto gas = canary::make_air5<P>(scales);
            const auto mechanism = canary::make_air5_park(scales);
            for (P mass : {P(0),
                           P(-1),
                           std::numeric_limits<P>::infinity(),
                           std::numeric_limits<P>::quiet_NaN()}) {
                auto species = gas.species;
                species.mass[4] = canary::ParticleMass<P, U> {mass};
                expect(error_message([&] {
                           (void)canary::make_air5_park_exchange(species, mechanism, scales);
                       }).starts_with("species \"NO\": particle mass must be finite and positive"));
            }
            if constexpr (std::same_as<P, float>) {
                for (double exponent : {-2.0 * std::numeric_limits<float>::max(),
                                        2.0 * std::numeric_limits<float>::max()}) {
                    auto input = canary::air5_park_exchange;
                    input.vibration_translation[0].correction.model =
                      canary::RelaxationCorrection::kim_2021;
                    input.vibration_translation[0].correction.exponent = exponent;
                    expect(error_message([&] {
                               (void)canary::make_energy_exchange(
                                 input,
                                 canary::detail::air5_park_species_names,
                                 gas.species,
                                 mechanism,
                                 scales);
                           }).find("stored correction exponent must be finite") !=
                           std::string::npos);
                }
            }
        };
        check(canary::UnitScales<float, si> {});
        check(canary::UnitScales<double, si> {});
        check(canary::UnitScales<float, scaled> {
          {.mass = 4.65e-26F, .density = 1.0e-3F, .temperature = 1000.0F, .length = 1.0F}});
        check(canary::UnitScales<double, scaled> {
          {.mass = 4.65e-26, .density = 1.0e-3, .temperature = 1000.0, .length = 1.0}});
    };
}
