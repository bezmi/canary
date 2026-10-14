#include "dual.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/data/air11.hpp>
#include <canary/data/air11_species.hpp>
#include <canary/data/air5.hpp>
#include <canary/data/air5_electronic.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/status.hpp>
#include <canary/thermo.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <string>
#include <type_traits>
#include <utility>

namespace
{

using canary::scaled;
using canary::si;
using canary::Status;
using D = canary::test::Dual<double>;
using Rho = canary::Density<double, si>;
using E = canary::SpecificEnergy<double, si>;
using RhoE = canary::EnergyDensity<double, si>;
using P = canary::Pressure<double, si>;
using Temp = canary::Temperature<double, si>;
using Extents = Kokkos::dextents<std::size_t, 1>;
using Air4 = canary::Air5FourTemperature<double, si>;
using Air4Scaled = canary::Air5FourTemperature<double, scaled>;
using Air11 = canary::Air11ThreeTemperature<double, si>;
using Air11Scaled = canary::Air11ThreeTemperature<double, scaled>;

constexpr double eps = std::numeric_limits<double>::epsilon();
// Relative tolerance of the temperature solves (the model default).
constexpr double solve_tolerance = 1.0e-10;

const canary::UnitScales<double, scaled> reference_scales {{
  .mass = 4.65e-26,
  .density = 1.0e-3,
  .temperature = 1000.0,
  .length = 1.0,
}};

template <typename T, std::size_t N> auto fractions(const std::array<T, N>& y)
{
    return canary::as_kind<canary::kind::mass_fraction>(
      Kokkos::mdspan<const T, Extents>(y.data(), N));
}

template <canary::UnitSystem U, typename T, std::size_t M>
auto mode_temperatures(const std::array<T, M>& t_m)
{
    return canary::as_kind<canary::kind::mode_temperature, U>(
      Kokkos::mdspan<const T, Extents>(t_m.data(), M));
}

template <canary::UnitSystem U, typename T, std::size_t M>
auto mode_energy_densities(const std::array<T, M>& rho_e_m)
{
    return canary::as_kind<canary::kind::mode_energy_density, U>(
      Kokkos::mdspan<const T, Extents>(rho_e_m.data(), M));
}

bool relative(double a, double b, double bound) { return std::abs(a - b) <= bound * std::abs(b); }

// |a - b| <= bound * max(|b|, floor).
bool close_to(double a, double b, double bound, double floor)
{ return std::abs(a - b) <= bound * std::max(std::abs(b), floor); }

constexpr std::array<double, 5> air5_air {0.767, 0.233, 0.0, 0.0, 0.0};
constexpr std::array<double, 5> air5_dissociated {0.3, 0.1, 0.3, 0.25, 0.05};
constexpr std::array<double, 5> air5_molecules {0.6, 0.3, 0.0, 0.0, 0.1};

// Compositions of the reference states: before and after the charge-neutral
// helper (the ionised compositions), and neutral air with atoms.
constexpr std::array<double, 11> air11_weakly_ionised_input {
  0.5,
  0.05,
  0.2,
  0.2,
  0.02,
  0.01,
  0.01,
  0.004,
  0.002,
  0.004,
  0.0,
};
constexpr std::array<double, 11> air11_weakly_ionised {
  0.4999995398459613,
  0.049999953984596134,
  0.19999981593838453,
  0.19999981593838453,
  0.01999998159383845,
  0.009999990796919225,
  0.009999990796919225,
  0.003999996318767691,
  0.0019999981593838453,
  0.003999996318767691,
  9.203080774345904e-07,
};
constexpr std::array<double, 11> air11_ionised {
  0.10100951469793618,
  0.010100951469793618,
  0.4040380587917447,
  0.3030285440938085,
  0.010100951469793618,
  0.08080761175834894,
  0.0606057088187617,
  0.010100951469793618,
  0.005050475734896809,
  0.015151427204690425,
  5.8044904318440994e-06,
};
constexpr std::array<double, 11> air11_neutral {0.3, 0.1, 0.3, 0.25, 0.05, 0, 0, 0, 0, 0, 0};
// Molecules, molecular ions and electrons only: every species has all of
// its components in modes.
constexpr std::array<double, 11> air11_molecules {
  0.5,
  0.2,
  0.0,
  0.0,
  0.1,
  0.0,
  0.0,
  0.1,
  0.05,
  0.05 - 2.0e-6,
  2.0e-6,
};

// A state of an independent implementation of the same model, evaluated
// from (rho, T, T_m, y): pressure, electron pressure (nan without
// electrons), translational-rotational energy u, mode energies and the
// frozen sound speed. Its universal gas constant is 8.31451 J/(mol K), its
// level temperatures differ from ours by about 6e-6 relative, and its
// electron energy is (3/2) R_e (T_e - T_ref).
template <std::size_t N, std::size_t M> struct ReferenceState {
    std::array<double, N> y;
    double rho;
    double temp;
    std::array<double, M> t_m;
    double p;
    double p_e;
    double u;
    std::array<double, M> u_m;
    double a;
};

const double no_value = std::numeric_limits<double>::quiet_NaN();

const std::array<ReferenceState<5, 3>, 5> air5_reference {
  {
    {
      .y = air5_air,
      .rho = 0.01,
      .temp = 8000.0,
      .t_m = {3000.0, 4000.0, 5000.0},
      .p = 23055.32275962548,
      .p_e = no_value,
      .u = 5463095.6189463,
      .u_m = {385699.54313113756, 227726.69742106923, 0.0},
      .a = 1717.8634734273871,
    },
    {
      .y = air5_air,
      .rho = 1.0,
      .temp = 300.0,
      .t_m = {300.0, 300.0, 300.0},
      .p = 86457.46034859553,
      .p_e = no_value,
      .u = -84591.42008858034,
      .u_m = {1.2147275219753502, 3.800540813899158, 0.0},
      .a = 347.7192295347249,
    },
    {
      .y = air5_dissociated,
      .rho = 0.001,
      .temp = 12000.0,
      .t_m = {8000.0, 9000.0, 10000.0},
      .p = 5242.582530039553,
      .p_e = no_value,
      .u = 23222656.736808177,
      .u_m = {648226.8993315813, 451864.38030879194, 969652.0993560122},
      .a = 2667.007484403189,
    },
    {
      .y = air5_dissociated,
      .rho = 0.01,
      .temp = 5000.0,
      .t_m = {5000.0, 5000.0, 5000.0},
      .p = 21844.093875164806,
      .p_e = no_value,
      .u = 17733237.43797749,
      .u_m = {331308.82330622786, 156327.0742693314, 108170.45346900629},
      .a = 1768.7443930886536,
    },
    {
      .y = air5_dissociated,
      .rho = 0.001,
      .temp = 2000.0,
      .t_m = {6000.0, 4000.0, 3000.0},
      .p = 873.7637550065922,
      .p_e = no_value,
      .u = 15380629.167050052,
      .u_m = {428222.94909218256, 103202.51370182398, 28577.760635586914},
      .a = 1129.62576871778,
    },
  },
};

const std::array<ReferenceState<11, 2>, 4> air11_reference {
  {
    {
      .y = air11_weakly_ionised,
      .rho = 0.001,
      .temp = 10000.0,
      .t_m = {5000.0, 8000.0},
      .p = 4147.0041617847555,
      .p_e = 111.58864072219835,
      .u = 20074720.14843051,
      .u_m = {645065.0908679113, 557256.2425120978},
      .a = 2368.767243712144,
    },
    {
      .y = air11_weakly_ionised,
      .rho = 0.01,
      .temp = 5000.0,
      .t_m = {5000.0, 5000.0},
      .p = 20874.50660982653,
      .p_e = 697.4290045137396,
      .u = 16199401.867069706,
      .u_m = {645065.0908679113, 159683.98677926004},
      .a = 1702.844502280075,
    },
    {
      .y = air11_ionised,
      .rho = 0.0001,
      .temp = 20000.0,
      .t_m = {10000.0, 15000.0},
      .p = 1173.3560784614415,
      .p_e = 131.962983359068,
      .u = 52522140.579992406,
      .u_m = {479619.8345225857, 5176765.931842619},
      .a = 4001.2212620364826,
    },
    {
      .y = air11_ionised,
      .rho = 0.001,
      .temp = 2000.0,
      .t_m = {3000.0, 1000.0},
      .p = 1129.3684173417523,
      .p_e = 87.97532223937867,
      .u = 37672909.28114045,
      .u_m = {77339.45935493513, 94162.99523824808},
      .a = 1310.406297394422,
    },
  },
};

// Builds a model from inputs that make_multi_temperature_mixture must
// reject, and returns the message, or "no exception".
template <typename Build> std::string builder_error(Build build)
{
    try {
        build();
    } catch (const canary::ParameterError& error) {
        return error.what();
    }
    return "no exception";
}

std::string air11_error(const std::array<canary::Nasa9Species<4>, 11>& species,
                        const std::array<int, 11>& degrees,
                        const std::array<canary::ElectronicLevelsInput<8>, 11>& levels,
                        const std::array<canary::EnergyModeInput<11>, 2>& modes)
{
    return builder_error([&] {
        canary::make_multi_temperature_mixture(species,
                                               degrees,
                                               levels,
                                               modes,
                                               canary::UnitScales<double, si> {});
    });
}

std::string air11_modes_error(const std::array<canary::EnergyModeInput<11>, 2>& modes)
{
    return air11_error(canary::air11_species,
                       canary::air11_rotational_degrees,
                       canary::air11_electronic_levels,
                       modes);
}

std::string air11_levels_error(const std::array<canary::ElectronicLevelsInput<8>, 11>& levels)
{
    return air11_error(canary::air11_species,
                       canary::air11_rotational_degrees,
                       levels,
                       canary::air11_three_temperature_modes);
}

// Long-double evaluation of eps(T) = R sum(g theta exp(-theta / T)) /
// sum(g exp(-theta / T)) from the model's levels.
long double direct_electronic_energy(const canary::ElectronicLevels<double, si, 8>& levels,
                                     double r,
                                     long double t)
{
    long double q = 0;
    long double a = 0;
    const auto n = static_cast<std::size_t>(levels.levels);
    for (std::size_t k = 0; k < n; ++k) {
        const long double theta = levels.temperatures[k].value();
        const long double w = levels.degeneracies[k] * std::exp(-theta / t);
        q += w;
        a += w * theta;
    }
    return levels.levels == 0 ? 0.0L : static_cast<long double>(r) * a / q;
}

// Checks every member of closure c against the expected state.
template <typename C, typename S, std::size_t M>
void check_state(const C& c,
                 const S& expected,
                 const std::array<double, M>& expected_t_m,
                 double bound,
                 const char* label)
{
    using namespace boost::ut;
    expect(c.status == Status::Valid) << label;
    expect(relative(c.state.rho.value(), expected.rho.value(), bound)) << label;
    expect(relative(c.state.temp.value(), expected.temp.value(), bound))
      << label << c.state.temp.value() << expected.temp.value();
    // Energies relative to the total energy scale: e and the mode
    // energies are sums of terms of size |e| or more.
    const double e_scale = std::max(std::abs(expected.e.value()), 1.0e5);
    expect(close_to(c.state.e.value(), expected.e.value(), bound, e_scale)) << label;
    expect(relative(c.state.p.value(), expected.p.value(), bound)) << label;
    for (std::size_t m = 0; m < M; ++m) {
        expect(relative(c.state.mode_temp[m].value(), expected_t_m[m], bound))
          << label << m << c.state.mode_temp[m].value() << expected_t_m[m];
        expect(close_to(c.state.mode_e[m].value(), expected.mode_e[m].value(), bound, e_scale))
          << label << m << c.state.mode_e[m].value() << expected.mode_e[m].value();
    }
}

// Every closure family from the state of (rho, T, T_m) recovers that state,
// for each density, T and set of mode temperatures. A mode without species
// present has temperature T in the closures from mode energies.
template <typename G, std::size_t N, std::size_t M>
void check_inversion(const G& gas,
                     const std::array<double, N>& y,
                     const std::array<double, 2>& densities,
                     const std::array<double, 3>& temperatures,
                     const std::array<std::array<double, M>, 3>& mode_sets)
{
    using namespace boost::ut;
    const auto mix = canary::mixture(gas, fractions(y));
    const auto values = canary::detail::fraction_values(mix);
    for (const double rho_value : densities) {
        for (const double t : temperatures) {
            for (const auto& t_m : mode_sets) {
                const Rho rho {rho_value};
                const auto base =
                  canary::thermo_state(gas, mix, rho, Temp {t}, mode_temperatures<si>(t_m));
                expect(base.status == Status::Valid);
                std::array<double, M> rho_e_m {};
                std::array<double, M> solved_t_m {};
                for (std::size_t m = 0; m < M; ++m) {
                    rho_e_m[m] = rho_value * base.state.mode_e[m].value();
                    solved_t_m[m] = canary::detail::mode_present(gas, values, m) ? t_m[m] : t;
                }
                const auto temps = mode_temperatures<si>(t_m);
                const auto energies = mode_energy_densities<si>(rho_e_m);
                const E e = base.state.e;
                const P p = base.state.p;
                const auto& s = base.state;
                check_state(canary::thermo_state(gas, mix, rho, e, temps),
                            s,
                            t_m,
                            solve_tolerance,
                            "rho e T_m");
                check_state(canary::thermo_state(gas, mix, rho, p, temps),
                            s,
                            t_m,
                            solve_tolerance,
                            "rho p T_m");
                check_state(canary::thermo_state(gas, mix, p, Temp {t}, temps),
                            s,
                            t_m,
                            solve_tolerance,
                            "p T T_m");
                check_state(canary::thermo_state(gas, mix, rho, e, energies),
                            s,
                            solved_t_m,
                            solve_tolerance,
                            "rho e e_m");
                check_state(
                  canary::thermo_state(gas, mix, rho, RhoE {rho_value * e.value()}, energies),
                  s,
                  solved_t_m,
                  solve_tolerance,
                  "rho rho_e e_m");
                check_state(canary::thermo_state(gas, mix, rho, p, energies),
                            s,
                            solved_t_m,
                            solve_tolerance,
                            "rho p e_m");
                check_state(canary::thermo_state(gas, mix, rho, Temp {t}, energies),
                            s,
                            solved_t_m,
                            solve_tolerance,
                            "rho T e_m");
                check_state(canary::thermo_state(gas, mix, p, Temp {t}, energies),
                            s,
                            solved_t_m,
                            solve_tolerance,
                            "p T e_m");
                if (std::all_of(t_m.begin(), t_m.end(), [t](double x) { return x == t; })) {
                    check_state(canary::thermo_state(gas, mix, rho, e),
                                s,
                                t_m,
                                solve_tolerance,
                                "rho e");
                    check_state(canary::thermo_state(gas, mix, rho, p),
                                s,
                                t_m,
                                solve_tolerance,
                                "rho p");
                    check_state(canary::thermo_state(gas, mix, rho, Temp {t}),
                                s,
                                t_m,
                                solve_tolerance,
                                "rho T");
                    check_state(canary::thermo_state(gas, mix, p, Temp {t}),
                                s,
                                t_m,
                                solve_tolerance,
                                "p T");
                }
            }
        }
    }
}

// Dual-number derivatives of e, e_m, p and a from (rho, T, T_m) with
// respect to T, each T_m and the mass fraction y[k], against central
// differences of the double closure.
template <typename G, std::size_t N, std::size_t M>
void check_derivatives(const G& gas,
                       const std::array<double, N>& y,
                       double rho_value,
                       double t,
                       const std::array<double, M>& t_m,
                       std::size_t k)
{
    using namespace boost::ut;
    // Values of e, e_m, p and a at (T, T_m, y).
    const auto evaluate = [&](double temp,
                              const std::array<double, M>& modes,
                              const std::array<double, N>& fractions_value) {
        const auto mix = canary::mixture(gas, fractions(fractions_value));
        const auto c = canary::thermo_state(gas,
                                            mix,
                                            Rho {rho_value},
                                            Temp {temp},
                                            mode_temperatures<si>(modes));
        std::array<double, M + 3> out {};
        out[0] = c.state.e.value();
        for (std::size_t m = 0; m < M; ++m) { out[1 + m] = c.state.mode_e[m].value(); }
        out[M + 1] = c.state.p.value();
        out[M + 2] = canary::properties(gas, mix, c.state).sound_speed.value();
        return out;
    };
    // The same with dual numbers, seeded in T (variable 0), T_m (1 + m) or
    // y[k] (M + 1).
    const auto evaluate_dual = [&](std::size_t variable) {
        std::array<D, N> y_d {};
        for (std::size_t i = 0; i < N; ++i) {
            y_d[i] = D(y[i], variable == M + 1 && i == k ? 1.0 : 0.0);
        }
        std::array<D, M> t_m_d {};
        for (std::size_t m = 0; m < M; ++m) { t_m_d[m] = D(t_m[m], variable == 1 + m ? 1.0 : 0.0); }
        const auto mix = canary::mixture(gas, fractions(y_d));
        const auto c =
          canary::thermo_state(gas,
                               mix,
                               canary::Density<D, si> {D(rho_value)},
                               canary::Temperature<D, si> {D(t, variable == 0 ? 1.0 : 0.0)},
                               mode_temperatures<si>(t_m_d));
        std::array<double, M + 3> out {};
        out[0] = c.state.e.value().dual();
        for (std::size_t m = 0; m < M; ++m) { out[1 + m] = c.state.mode_e[m].value().dual(); }
        out[M + 1] = c.state.p.value().dual();
        out[M + 2] = canary::properties(gas, mix, c.state).sound_speed.value().dual();
        return out;
    };
    const auto base = evaluate(t, t_m, y);
    for (std::size_t variable = 0; variable < M + 2; ++variable) {
        // Steps: 1e-5 relative for the temperatures, 1e-8 for the mass
        // fraction. The truncation error is about h^2 times the third
        // derivative relative to the first: below 1e-8 for the
        // temperatures, and below 1e-7 for the mass fraction of e-, whose
        // gas constant is 4e4 times that of the mixture. The rounding error
        // is about eps |f| / h, below 1e-7 of the derivative.
        double h = 0;
        // The variable's value, 1 for the mass fraction.
        double x = 1.0;
        std::array<double, M + 3> plus {};
        std::array<double, M + 3> minus {};
        if (variable == 0) {
            x = t;
            h = 1.0e-5 * t;
            plus = evaluate(t + h, t_m, y);
            minus = evaluate(t - h, t_m, y);
        } else if (variable <= M) {
            auto up = t_m;
            auto down = t_m;
            x = t_m[variable - 1];
            h = 1.0e-5 * x;
            up[variable - 1] += h;
            down[variable - 1] -= h;
            plus = evaluate(t, up, y);
            minus = evaluate(t, down, y);
        } else {
            auto up = y;
            auto down = y;
            h = 1.0e-8;
            up[k] += h;
            down[k] -= h;
            plus = evaluate(t, t_m, up);
            minus = evaluate(t, t_m, down);
        }
        const auto dual = evaluate_dual(variable);
        for (std::size_t q = 0; q < M + 3; ++q) {
            const double fd = (plus[q] - minus[q]) / (2.0 * h);
            // 1e-6 relative to the derivative, or to the scale |f| / x of a
            // derivative that is zero or nearly so.
            expect(close_to(dual[q], fd, 1.0e-6, std::abs(base[q]) / x))
              << "variable" << variable << "quantity" << q << dual[q] << fd;
        }
    }
}

}  // namespace

static_assert(noexcept(canary::species_vibration_mode(std::declval<const Air11&>(),
                                                      canary::SpeciesIndex {0})));
static_assert(noexcept(canary::electron_species_index(std::declval<const Air11&>())));
static_assert(noexcept(canary::electron_energy_mode(std::declval<const Air11&>())));

static_assert(canary::ThermoStateOf<Air4::State<double>, double, si>);
static_assert(canary::GasModel<Air4, double>);
static_assert(canary::MixtureGasModel<Air4, double>);
static_assert(canary::MultiTemperatureGasModel<Air4, double>);
static_assert(canary::MultiTemperatureGasModel<Air4, float>);
static_assert(canary::MultiTemperatureGasModel<Air4, D>);
static_assert(canary::MultiTemperatureGasModel<Air4Scaled, double>);
static_assert(canary::MultiTemperatureGasModel<Air11, double>);
static_assert(canary::MultiTemperatureGasModel<Air11, float>);
static_assert(canary::MultiTemperatureGasModel<Air11, D>);
static_assert(canary::MultiTemperatureGasModel<Air11Scaled, D>);
static_assert(canary::AcceptsInput<Air4, canary::ModeTemperatureView<double, si>>);
static_assert(canary::AcceptsInput<Air4, canary::ModeEnergyDensityView<double, si>>);
static_assert(!canary::AcceptsInput<Air4Scaled, canary::ModeEnergyDensityView<double, si>>);
static_assert(std::is_same_v<Air4::State<double>, canary::MultiTemperatureState<double, si, 3>>);
static_assert(std::is_same_v<Air11::State<D>, canary::MultiTemperatureState<D, si, 2>>);
static_assert(std::is_trivially_copyable_v<Air4>);
static_assert(std::is_trivially_copyable_v<Air11Scaled>);
static_assert(std::is_trivially_copyable_v<Air11::State<D>>);
static_assert(std::is_trivially_copyable_v<Air11::Mixture<D>>);

void test_multi_temperature()
{
    using namespace boost::ut;
    using canary::mixture;
    using canary::thermo_state;
    const canary::UnitScales<double, si> si_scales {};
    const auto built4 = canary::make_air5_four_temperature<double>(si_scales);
    const Air4& air4 = built4.model;
    const auto built11 = canary::make_air11_three_temperature<double>(si_scales);
    const Air11& air11 = built11.model;
    const auto single5 = canary::make_air5<double>(si_scales).model;
    const auto single11 = canary::make_air11<double>(si_scales).model;

    "multi-temperature construction retains long double scale precision"_test = [] {
        const long double tau = std::nextafter(1.0L, 2.0L);
        const canary::UnitScales<long double, scaled> reference {
          {.mass = 4.65e-26L, .density = 1.0e-3L, .temperature = tau, .length = 1.0L}};
        const auto built = canary::make_air11_three_temperature<long double>(reference);
        for (std::size_t i = 0; i < 11; ++i) {
            const auto& input = canary::air11_electronic_levels[i];
            for (int k = 0; k < input.levels; ++k) {
                const auto index = static_cast<std::size_t>(k);
                expect(built.model.electronic_levels[i].temperatures[index].value() ==
                       static_cast<long double>(input.temperatures[index]) / tau);
            }
        }
    };

    "multi-temperature closures reject electron pressure overflow with finite heavy pressure"_test =
      [&] {
          const auto check = [&]<typename V> {
              std::array<V, 11> y {};
              y[0] = V(0.9);
              y[10] = V(0.1);
              const auto mix = mixture(air11, fractions(y));
              const canary::Density<V, si> rho {std::numeric_limits<V>::max() / V(1.0e9)};
              const canary::Temperature<V, si> temp {V(300)};
              const std::array<V, 2> modes {V(3000), V(10000)};
              const V heavy_pressure =
                rho.value() * V(V(0.9) * air11.thermally_perfect.gas_constant[0].value() * V(300));
              expect(std::isfinite(heavy_pressure) && heavy_pressure > V(0));
              const auto result = thermo_state(air11, mix, rho, temp, mode_temperatures<si>(modes));
              expect(result.status == Status::Invalid);
              expect(std::isinf(result.state.p.value()));
              expect(result.state.rho == rho && result.state.temp == temp);
              expect(std::isfinite(result.state.e.value()));
              for (const auto& energy : result.state.mode_e) {
                  expect(std::isfinite(energy.value()));
              }
              const auto inverse =
                thermo_state(air11, mix, rho, result.state.e, mode_temperatures<si>(modes));
              expect(inverse.status == Status::Invalid);
              expect(std::isinf(inverse.state.p.value()));
              y.fill(V(0));
              y[0] = V(1);
              const auto neutral = mixture(air11, fractions(y));
              const canary::Pressure<V, si> supplied {std::numeric_limits<V>::max() / V(4)};
              const canary::Temperature<V, si> cold {V(1) / neutral.gas_constant.value()};
              const auto retained =
                thermo_state(air11, neutral, supplied, cold, mode_temperatures<si>(modes));
              expect(retained.status == Status::Valid);
              expect(retained.state.p == supplied);
              expect(thermo_state(air11, neutral, supplied, cold).status == Status::Valid);
          };
          check.operator()<float>();
          check.operator()<double>();
      };

    "multi-temperature: the builder holds the single-temperature model, levels and modes"_test =
      [&] {
          expect(built11.species.n == 11_i);
          expect(!air4.electron);
          expect(air4.electron_mode == -1_i);
          expect(air11.electron == canary::OptionalSpeciesIndex {canary::SpeciesIndex {10}});
          expect(air11.electron_mode == 1_i);
          for (std::size_t i = 0; i < 11; ++i) {
              expect(air11.thermally_perfect.gas_constant[i] == single11.gas_constant[i]);
              const double r_i = single11.gas_constant[i].value();
              const double d = canary::air11_rotational_degrees[i];
              // One multiplication and one division of exact small numbers.
              expect(relative(air11.translational_rotational_cp[i].value(),
                              (2.5 + (0.5 * d)) * r_i,
                              2.0 * eps));
              const auto& input = canary::air11_electronic_levels[i];
              expect(air11.electronic_levels[i].levels == input.levels);
              const auto n = static_cast<std::size_t>(input.levels);
              for (std::size_t k = 0; k < n; ++k) {
                  expect(air11.electronic_levels[i].temperatures[k].value() ==
                         input.temperatures[k]);
                  expect(air11.electronic_levels[i].degeneracies[k] ==
                         static_cast<double>(input.degeneracies[k]));
              }
              const auto& vib = air11.mode_components[0][i];
              const auto& el = air11.mode_components[1][i];
              expect(vib.vibration == (canary::air11_rotational_degrees[i] > 0));
              expect(!vib.electronic && !vib.electron);
              expect(el.electronic == (i != 10));
              expect(el.electron == (i == 10));
              expect(!el.vibration);
          }
          // Mode O2 of the four-temperature model: O2 vibration and
          // electronic, O electronic.
          expect(air4.mode_components[1][1].vibration && air4.mode_components[1][1].electronic);
          expect(air4.mode_components[1][3].electronic && !air4.mode_components[1][3].vibration);
          expect(!canary::detail::holds_any(air4.mode_components[1][0]));
      };

    "multi-temperature: species roles are independent of names and ordering"_test = [] {
        auto species = canary::air11_species;
        auto degrees = canary::air11_rotational_degrees;
        auto levels = canary::air11_electronic_levels;
        auto modes = canary::air11_three_temperature_modes;
        species[10].name = "free-electron";
        species[2].name = "e-";
        const auto rename = [](const char* name) {
            if (std::string_view(name) == "e-") { return "free-electron"; }
            if (std::string_view(name) == "N") { return "e-"; }
            return name;
        };
        for (auto& entry : levels) { entry.name = rename(entry.name); }
        for (auto& mode : modes) {
            for (int i = 0; i < mode.count; ++i) {
                auto& entry = mode.components[static_cast<std::size_t>(i)];
                entry.species = rename(entry.species);
            }
        }
        for (const std::size_t electron : {std::size_t(10), std::size_t(0)}) {
            if (electron == 0) {
                std::swap(species[0], species[10]);
                std::swap(degrees[0], degrees[10]);
            }
            const auto result =
              canary::make_multi_temperature_mixture(species,
                                                     degrees,
                                                     levels,
                                                     modes,
                                                     canary::UnitScales<double, si> {});
            expect(result.model.electron ==
                   canary::OptionalSpeciesIndex {canary::SpeciesIndex {electron}});
            expect(result.model.electron_mode == 1_i);
            expect(result.model.mode_components[1][electron].electron);
            expect(!result.model.mode_components[1][2].electron);
        }
        species = canary::air11_species;
        degrees = canary::air11_rotational_degrees;
        species[2].metadata.role = canary::SpeciesRole::electron;
        species[2].charge = -1;
        const auto error = air11_error(species,
                                       degrees,
                                       canary::air11_electronic_levels,
                                       canary::air11_three_temperature_modes);
        expect(error.find("more than one electron-role species") != std::string::npos) << error;
    };

    "multi-temperature: the builder rejects null or empty species names"_test = [] {
        for (const char* name : std::array<const char*, 2> {nullptr, ""}) {
            auto species = canary::air11_species;
            species[1].name = name;
            const std::string error = air11_error(species,
                                                  canary::air11_rotational_degrees,
                                                  canary::air11_electronic_levels,
                                                  canary::air11_three_temperature_modes);
            expect(error.starts_with("species 1: name must not be null or empty")) << error;
        }
    };

    "multi-temperature: the builder rejects an unknown species in a mode"_test = [] {
        auto modes = canary::air11_three_temperature_modes;
        modes[0].components[2].species = "Ar";
        const std::string error = air11_modes_error(modes);
        expect(error.starts_with("mode 0, component 2 (Ar vibration): unknown species")) << error;
    };

    "multi-temperature: the builder rejects an unknown species in the levels"_test = [] {
        auto levels = canary::air11_electronic_levels;
        levels[3].name = "Ar";
        const std::string error = air11_levels_error(levels);
        expect(error.starts_with("electronic levels of Ar: unknown species")) << error;
    };

    "multi-temperature: the builder rejects a component given twice"_test = [] {
        auto modes = canary::air11_three_temperature_modes;
        modes[1].components[2] = {.species = "N2", .component = canary::ModeComponent::vibration};
        std::string error = air11_modes_error(modes);
        expect(error.starts_with("mode 1, component 2 (N2 vibration): already in mode 0")) << error;
        modes = canary::air11_three_temperature_modes;
        modes[0].components[1] = {.species = "N2", .component = canary::ModeComponent::vibration};
        error = air11_modes_error(modes);
        expect(error.starts_with("mode 0, component 1 (N2 vibration): already in mode 0")) << error;
    };

    "multi-temperature: the builder rejects vibration of an atom"_test = [] {
        auto modes = canary::air11_three_temperature_modes;
        modes[0].components[1] = {.species = "O+", .component = canary::ModeComponent::vibration};
        const std::string error = air11_modes_error(modes);
        expect(error.starts_with("mode 0, component 1 (O+ vibration): an atom has no vibration"))
          << error;
    };

    "multi-temperature: the builder rejects vibration of e-"_test = [] {
        auto modes = canary::air11_three_temperature_modes;
        modes[1].components[10].component = canary::ModeComponent::vibration;
        const std::string error = air11_modes_error(modes);
        expect(error.starts_with(
          "mode 1, component 10 (e- vibration): free electrons have only the electron component"))
          << error;
    };

    "multi-temperature: the builder rejects electronic energy of e-"_test = [] {
        auto modes = canary::air11_three_temperature_modes;
        modes[1].components[10].component = canary::ModeComponent::electronic;
        const std::string error = air11_modes_error(modes);
        expect(error.starts_with(
          "mode 1, component 10 (e- electronic): free electrons have only the electron component"))
          << error;
    };

    "multi-temperature: the builder rejects the electron component of a heavy species"_test = [] {
        auto modes = canary::air11_three_temperature_modes;
        modes[1].components[3].component = canary::ModeComponent::electron;
        const std::string error = air11_modes_error(modes);
        expect(error.starts_with("mode 1, component 3 (O electron): only an electron-role species "
                                 "has the electron component"))
          << error;
    };

    "multi-temperature: the builder rejects e- in no mode"_test = [] {
        auto modes = canary::air11_three_temperature_modes;
        modes[1].count = 10;
        const std::string error = air11_modes_error(modes);
        expect(error.starts_with("species e-: free electrons must be in a mode")) << error;
    };

    "multi-temperature: the builder rejects a mode without components"_test = [] {
        auto modes = canary::air11_three_temperature_modes;
        modes[0].count = 0;
        std::string error = air11_modes_error(modes);
        expect(error.starts_with("mode 0: the component count must be between 1 and 11, got 0"))
          << error;
        modes[0].count = 12;
        error = air11_modes_error(modes);
        expect(error.starts_with("mode 0: the component count must be between 1 and 11, got 12"))
          << error;
    };

    "multi-temperature: the builder rejects rotational degrees other than 0, 2 and 3"_test = [] {
        auto degrees = canary::air11_rotational_degrees;
        degrees[4] = 1;
        const std::string error = air11_error(canary::air11_species,
                                              degrees,
                                              canary::air11_electronic_levels,
                                              canary::air11_three_temperature_modes);
        expect(
          error.starts_with("species NO: rotational degrees of freedom must be 0, 2 or 3, got 1"))
          << error;
    };

    "multi-temperature: the builder rejects rotational degrees of e-"_test = [] {
        auto degrees = canary::air11_rotational_degrees;
        degrees[10] = 2;
        const std::string error = air11_error(canary::air11_species,
                                              degrees,
                                              canary::air11_electronic_levels,
                                              canary::air11_three_temperature_modes);
        expect(error.starts_with("species e-: rotational degrees of freedom must be 0, got 2"))
          << error;
    };

    "multi-temperature: the builder rejects level temperatures that are not finite"_test = [] {
        auto levels = canary::air11_electronic_levels;
        levels[1].temperatures[2] = std::numeric_limits<double>::infinity();
        const std::string error = air11_levels_error(levels);
        expect(error.starts_with(
          "electronic levels of O2: level 2: temperature must be finite and non-negative, got inf"))
          << error;
    };

    "multi-temperature: the builder rejects negative level temperatures"_test = [] {
        auto levels = canary::air11_electronic_levels;
        levels[1].temperatures[2] = -1.0;
        const std::string error = air11_levels_error(levels);
        expect(error.starts_with(
          "electronic levels of O2: level 2: temperature must be finite and non-negative, got -1"))
          << error;
    };

    "multi-temperature: the builder rejects degeneracies that are not positive"_test = [] {
        auto levels = canary::air11_electronic_levels;
        levels[2].degeneracies[1] = 0;
        const std::string error = air11_levels_error(levels);
        expect(
          error.starts_with("electronic levels of N: level 1: degeneracy must be positive, got 0"))
          << error;
    };

    "multi-temperature: the builder rejects more levels than the capacity"_test = [] {
        auto levels = canary::air11_electronic_levels;
        levels[0].levels = 9;
        const std::string error = air11_levels_error(levels);
        expect(
          error.starts_with("electronic levels of N2: level count must be between 0 and 8, got 9"))
          << error;
    };

    "multi-temperature: the builder rejects a first level above zero"_test = [] {
        auto levels = canary::air11_electronic_levels;
        levels[0].temperatures[0] = 10.0;
        const std::string error = air11_levels_error(levels);
        expect(error.starts_with(
          "electronic levels of N2: level 0: the ground level temperature must be 0, got 10"))
          << error;
    };

    "multi-temperature: the builder rejects levels given twice"_test = [] {
        auto levels = canary::air11_electronic_levels;
        levels[10].name = "N2";
        const std::string error = air11_levels_error(levels);
        expect(error.starts_with("electronic levels of N2: given twice")) << error;
    };

    "multi-temperature: the builder rejects e- without charge -1"_test = [] {
        auto species = canary::air11_species;
        species[10].charge = 1.0;
        const std::string error = air11_error(species,
                                              canary::air11_rotational_degrees,
                                              canary::air11_electronic_levels,
                                              canary::air11_three_temperature_modes);
        expect(error.starts_with("species e-: charge must be -1, got 1")) << error;
    };

    "multi-temperature: the builder rejects species names given twice"_test = [] {
        auto species = canary::air11_species;
        species[6].name = "N+";
        const std::string error = air11_error(species,
                                              canary::air11_rotational_degrees,
                                              canary::air11_electronic_levels,
                                              canary::air11_three_temperature_modes);
        expect(error.starts_with("species N+: name given twice")) << error;
    };

    "multi-temperature: the builder rejects what the single-temperature builder rejects"_test = [] {
        auto species = canary::air11_species;
        species[1].molar_mass = -1.0;
        const std::string error = air11_error(species,
                                              canary::air11_rotational_degrees,
                                              canary::air11_electronic_levels,
                                              canary::air11_three_temperature_modes);
        expect(error.starts_with("species O2: molar mass")) << error;
    };

    // Agreement bounds of the reference states: energies and pressures 2e-5
    // relative (the universal gas constants differ by 5.7e-6 relative, the
    // level temperatures by about 6e-6), relative to the larger of the
    // value and 1e-3 of the state's total energy; the sound speed 1e-5.
    "multi-temperature: (rho, T, T_m) gives the reference states of four-temperature air"_test =
      [&] {
          for (const auto& ref : air5_reference) {
              const auto mix = mixture(air4, fractions(ref.y));
              const auto c = thermo_state(air4,
                                          mix,
                                          Rho {ref.rho},
                                          Temp {ref.temp},
                                          mode_temperatures<si>(ref.t_m));
              expect(c.status == Status::Valid);
              double e_m = 0;
              double ref_e = ref.u;
              for (std::size_t m = 0; m < 3; ++m) {
                  e_m += c.state.mode_e[m].value();
                  ref_e += ref.u_m[m];
              }
              const double floor = 1.0e-3 * std::abs(ref_e);
              expect(close_to(c.state.p.value(), ref.p, 2.0e-5, 0.0)) << ref.temp;
              expect(canary::electron_pressure(air4, mix, c.state).value() == 0.0);
              expect(close_to(c.state.e.value() - e_m, ref.u, 2.0e-5, floor))
                << ref.temp << c.state.e.value() - e_m << ref.u;
              for (std::size_t m = 0; m < 3; ++m) {
                  expect(close_to(c.state.mode_e[m].value(), ref.u_m[m], 2.0e-5, floor))
                    << ref.temp << m << c.state.mode_e[m].value() << ref.u_m[m];
              }
              const auto props = canary::properties(air4, mix, c.state);
              expect(relative(props.sound_speed.value(), ref.a, 1.0e-5))
                << ref.temp << props.sound_speed.value() << ref.a;
          }
      };

    "multi-temperature: (rho, T, T_m) gives the reference states of three-temperature air"_test =
      [&] {
          const auto species = canary::species_properties(single11, Temp {298.15});
          const double r_e = single11.gas_constant[10].value();
          for (const auto& ref : air11_reference) {
              const auto mix = mixture(air11, fractions(ref.y));
              const auto c = thermo_state(air11,
                                          mix,
                                          Rho {ref.rho},
                                          Temp {ref.temp},
                                          mode_temperatures<si>(ref.t_m));
              expect(c.status == Status::Valid);
              // The reference electron energy is (3/2) R_e (T_e - T_ref);
              // ours is h_e(T_e) - R_e T_e, larger by h_e(T_ref) - R_e T_ref.
              std::array<double, 2> u_m = ref.u_m;
              u_m[1] += ref.y[10] * (species.enthalpy[10].value() - (r_e * 298.15));
              const double ref_e = ref.u + u_m[0] + u_m[1];
              const double e_m = c.state.mode_e[0].value() + c.state.mode_e[1].value();
              const double floor = 1.0e-3 * std::abs(ref_e);
              expect(close_to(c.state.p.value(), ref.p, 2.0e-5, 0.0)) << ref.temp;
              expect(close_to(canary::electron_pressure(air11, mix, c.state).value(),
                              ref.p_e,
                              2.0e-5,
                              0.0))
                << ref.temp;
              expect(close_to(c.state.e.value() - e_m, ref.u, 2.0e-5, floor))
                << ref.temp << c.state.e.value() - e_m << ref.u;
              for (std::size_t m = 0; m < 2; ++m) {
                  expect(close_to(c.state.mode_e[m].value(), u_m[m], 2.0e-5, floor))
                    << ref.temp << m << c.state.mode_e[m].value() << u_m[m];
              }
              const auto props = canary::properties(air11, mix, c.state);
              expect(relative(props.sound_speed.value(), ref.a, 1.0e-5))
                << ref.temp << props.sound_speed.value() << ref.a;
          }
      };

    "multi-temperature: equal temperatures give the single-temperature energy"_test = [&] {
        const auto compare = [](const auto& gas, const auto& single, const auto& y, auto t_m) {
            const auto mix = mixture(gas, fractions(y));
            for (const double t : {300.0, 2000.0, 8000.0, 15000.0, 30000.0}) {
                const Rho rho {0.5};
                const auto ref = thermo_state(single, mix, rho, Temp {t});
                t_m.fill(t);
                const auto c = thermo_state(gas, mix, rho, Temp {t}, mode_temperatures<si>(t_m));
                // 1e-12 relative to the size of the terms that are summed:
                // y_i (|dh_f,i| + cp_i T) with cp_i < 10 R_i.
                double scale = 0;
                for (std::size_t i = 0; i < y.size(); ++i) {
                    scale += y[i] * (std::abs(gas.formation_enthalpy[i].value()) +
                                     (10.0 * gas.thermally_perfect.gas_constant[i].value() * t));
                }
                expect(std::abs(c.state.e.value() - ref.state.e.value()) <= 1.0e-12 * scale)
                  << t << c.state.e.value() << ref.state.e.value();
                // With e- at T_e = T, p = rho R T; R is summed in another
                // order, which changes the last bits.
                expect(relative(c.state.p.value(), ref.state.p.value(), 4.0 * eps));
            }
        };
        compare(air4, single5, air5_molecules, std::array<double, 3> {});
        compare(air11, single11, air11_molecules, std::array<double, 2> {});
    };

    "multi-temperature: four-temperature air gives the two-temperature e_ve per molecule"_test =
      [&] {
          const auto two = canary::make_air5_two_temperature<double>(si_scales).model;
          const auto mix = mixture(air4, fractions(air5_molecules));
          // Mode m holds the vibration and electronic energy of molecule
          // index[m] only.
          constexpr std::array<std::size_t, 3> molecule {0, 1, 4};
          for (const std::array<double, 3> t_m : {
                 std::array {300.0, 300.0, 300.0},
                 std::array {2000.0, 9000.0, 4000.0},
                 std::array {15000.0, 700.0, 25000.0},
               }) {
              const auto c =
                thermo_state(air4, mix, Rho {0.1}, Temp {5000.0}, mode_temperatures<si>(t_m));
              for (std::size_t m = 0; m < 3; ++m) {
                  const auto e_ve =
                    canary::species_vibrational_electronic_energy(two, Temp {t_m[m]});
                  const std::size_t i = molecule[m];
                  // Both are y_i (h_i - cp_tr,i (T - T_ref) - dh_f,i) with
                  // the same terms; 4 eps of the largest term.
                  const double scale =
                    air5_molecules[i] * (std::abs(two.formation_enthalpy[i].value()) +
                                         (two.translational_rotational_cp[i].value() * t_m[m]) +
                                         std::abs(e_ve[i].value()));
                  expect(std::abs(c.state.mode_e[m].value() -
                                  (air5_molecules[i] * e_ve[i].value())) <= 4.0 * eps * scale)
                    << m << c.state.mode_e[m].value() << air5_molecules[i] * e_ve[i].value();
              }
          }
      };

    "multi-temperature: every closure family inverts (rho, T, T_m)"_test = [&] {
        const std::array<double, 2> densities {1.0e-4, 1.0};
        const std::array<double, 3> temperatures {300.0, 3000.0, 15000.0};
        // Air without NO, N and O: the NO mode has no species present.
        check_inversion(
          air4,
          air5_air,
          densities,
          temperatures,
          std::array<std::array<double, 3>, 3> {
            {{300.0, 300.0, 300.0}, {3000.0, 3000.0, 3000.0}, {12000.0, 600.0, 4000.0}},
          });
        check_inversion(
          air4,
          air5_dissociated,
          densities,
          temperatures,
          std::array<std::array<double, 3>, 3> {
            {{15000.0, 15000.0, 15000.0}, {300.0, 2000.0, 8000.0}, {20000.0, 5000.0, 700.0}},
          });
        for (const auto& y : {air11_weakly_ionised, air11_ionised, air11_neutral}) {
            check_inversion(air11,
                            y,
                            densities,
                            temperatures,
                            std::array<std::array<double, 2>, 3> {
                              {{3000.0, 3000.0}, {600.0, 20000.0}, {12000.0, 1500.0}},
                            });
        }
    };

    "multi-temperature: a mode without species present has the temperature T and is Valid"_test =
      [&] {
          const auto mix = mixture(air4, fractions(air5_air));
          const std::array<double, 3> rho_e_m {100.0, 50.0, 0.0};
          const auto c =
            thermo_state(air4, mix, Rho {1.0}, Temp {4000.0}, mode_energy_densities<si>(rho_e_m));
          expect(c.status == Status::Valid);
          expect(c.state.mode_temp[2].value() == 4000.0);
          expect(c.state.mode_e[2].value() == 0.0);
          const auto d =
            thermo_state(air4, mix, Rho {1.0}, E {1.0e6}, mode_energy_densities<si>(rho_e_m));
          expect(d.status == Status::Valid);
          expect(d.state.mode_temp[2].value() == d.state.temp.value());
      };

    "multi-temperature: the electronic energy matches a direct evaluation from 10 K to 1e5 K"_test =
      [&] {
          for (std::size_t i = 0; i < 11; ++i) {
              const auto& levels = air11.electronic_levels[i];
              const double r = air11.thermally_perfect.gas_constant[i].value();
              for (int j = 0; j <= 40; ++j) {
                  const double t = 10.0 * std::pow(10.0, 0.1 * j);
                  const auto ours =
                    canary::detail::electronic_energy(levels,
                                                      air11.thermally_perfect.gas_constant[i],
                                                      t);
                  const long double direct = direct_electronic_energy(levels, r, t);
                  // Sums of at most 8 positive terms, each with a few
                  // roundings: 1e-13 relative. Terms below the smallest
                  // double, which the long double evaluation keeps, are
                  // covered by the absolute floor of 1e-280 J/kg.
                  expect(std::abs(ours.e - static_cast<double>(direct)) <=
                         (1.0e-13 * std::abs(static_cast<double>(direct))) + 1.0e-280)
                    << i << t << ours.e << static_cast<double>(direct);
                  // Central difference of the long-double evaluation with
                  // h = 1e-7 T: truncation (h / T)^2 (theta / T)^2 / 6 below
                  // 1e-7 relative where the weights are above the smallest
                  // double, rounding about 1e-19 / 1e-7.
                  const long double h = 1.0e-7L * t;
                  const long double fd = (direct_electronic_energy(levels, r, t + h) -
                                          direct_electronic_energy(levels, r, t - h)) /
                                         (2.0L * h);
                  expect(std::abs(ours.cv - static_cast<double>(fd)) <=
                         (1.0e-6 * std::abs(static_cast<double>(fd))) + 1.0e-280)
                    << i << t << ours.cv << static_cast<double>(fd);
              }
          }
      };

    "multi-temperature: dual numbers give the derivatives of e, e_m, p and a"_test = [&] {
        // Temperatures outside the blending regions of the curves (800 to
        // 1200 K, 5500 to 6500 K, 19500 to 20500 K), where the dual-number
        // derivative treats the blending weight as constant.
        check_derivatives(air4,
                          air5_dissociated,
                          0.01,
                          9000.0,
                          std::array {3000.0, 4500.0, 12000.0},
                          3);
        check_derivatives(air4, air5_air, 1.0, 3000.0, std::array {2500.0, 700.0, 4000.0}, 1);
        check_derivatives(air11, air11_ionised, 0.001, 14000.0, std::array {9000.0, 12000.0}, 10);
        check_derivatives(air11,
                          air11_weakly_ionised,
                          0.01,
                          4000.0,
                          std::array {3000.0, 9000.0},
                          6);
    };

    "multi-temperature: dual numbers give the mode temperature derivative 1 / cv_m"_test = [&] {
        const auto mix = mixture(air11, fractions(air11_ionised));
        const std::array<D, 11> y_d = [] {
            std::array<D, 11> out {};
            for (std::size_t i = 0; i < 11; ++i) { out[i] = D(air11_ionised[i]); }
            return out;
        }();
        const auto mix_d = mixture(air11, fractions(y_d));
        const std::array t_m {9000.0, 12000.0};
        const auto base =
          thermo_state(air11, mix, Rho {1.0}, Temp {14000.0}, mode_temperatures<si>(t_m));
        const auto props = canary::detail::mode_energy(air11, mix.fractions, 1, 12000.0);
        // d T_e / d e_m at fixed rho = 1, e and the other mode energy.
        const std::array rho_e_m {
          D(base.state.mode_e[0].value()),
          D(base.state.mode_e[1].value(), 1.0),
        };
        const auto c = thermo_state(air11,
                                    mix_d,
                                    canary::Density<D, si> {D(1.0)},
                                    canary::SpecificEnergy<D, si> {D(base.state.e.value())},
                                    mode_energy_densities<si>(rho_e_m));
        expect(c.status == Status::Valid);
        // The solve converges to 1e-10 relative, and cv_m at the root
        // differs from cv_m at 12000 K by that much times d ln cv / d ln T.
        expect(relative(c.state.mode_temp[1].value().dual(), 1.0 / props.cv, 1.0e-9))
          << c.state.mode_temp[1].value().dual() << 1.0 / props.cv;
    };

    "multi-temperature: a scaled model gives the SI state after conversion"_test = [&] {
        const auto convert = [](auto q) {
            return canary::convert<si>(q, reference_scales).value();
        };
        const auto compare_models = [&](const auto& gas,
                                        const auto& gas_s,
                                        const auto& y,
                                        auto t_m_values) {
            constexpr std::size_t M = std::tuple_size_v<decltype(t_m_values)>;
            const auto mix = mixture(gas, fractions(y));
            const auto mix_s = mixture(gas_s, fractions(y));
            for (const double rho_value : {1.0e-4, 1.0}) {
                for (const double t : {1000.0, 4000.0, 15000.0}) {
                    const Rho rho {rho_value};
                    const auto rho_s = canary::convert<scaled>(rho, reference_scales);
                    const auto base =
                      thermo_state(gas, mix, rho, Temp {t}, mode_temperatures<si>(t_m_values));
                    std::array<double, M> t_m_s {};
                    std::array<double, M> rho_e_m {};
                    std::array<double, M> rho_e_m_s {};
                    for (std::size_t m = 0; m < M; ++m) {
                        t_m_s[m] = canary::convert<scaled>(
                                     canary::ModeTemperature<double, si> {t_m_values[m]},
                                     reference_scales)
                                     .value();
                        rho_e_m[m] = rho_value * base.state.mode_e[m].value();
                        rho_e_m_s[m] = canary::convert<scaled>(
                                         canary::ModeEnergyDensity<double, si> {rho_e_m[m]},
                                         reference_scales)
                                         .value();
                    }
                    const RhoE rho_e {rho_value * base.state.e.value()};
                    const auto compare = [&](const auto& a, const auto& b) {
                        expect(b.status == Status::Valid);
                        // The folded curves and level temperatures round
                        // differently: 1e-12 relative, energies relative to
                        // the total energy scale.
                        const double e_scale = std::max(std::abs(a.state.e.value()), 1.0e5);
                        expect(relative(convert(b.state.temp), a.state.temp.value(), 1.0e-12))
                          << rho_value << t;
                        expect(relative(convert(b.state.p), a.state.p.value(), 1.0e-12))
                          << rho_value << t;
                        expect(close_to(convert(b.state.e), a.state.e.value(), 1.0e-12, e_scale));
                        for (std::size_t m = 0; m < M; ++m) {
                            expect(relative(convert(b.state.mode_temp[m]),
                                            a.state.mode_temp[m].value(),
                                            1.0e-12))
                              << rho_value << t << m;
                            expect(close_to(convert(b.state.mode_e[m]),
                                            a.state.mode_e[m].value(),
                                            1.0e-12,
                                            e_scale))
                              << rho_value << t << m;
                        }
                    };
                    compare(base,
                            thermo_state(gas_s,
                                         mix_s,
                                         rho_s,
                                         canary::convert<scaled>(Temp {t}, reference_scales),
                                         mode_temperatures<scaled>(t_m_s)));
                    compare(thermo_state(gas, mix, rho, rho_e, mode_energy_densities<si>(rho_e_m)),
                            thermo_state(gas_s,
                                         mix_s,
                                         rho_s,
                                         canary::convert<scaled>(rho_e, reference_scales),
                                         mode_energy_densities<scaled>(rho_e_m_s)));
                }
            }
        };
        const auto air4_s = canary::make_air5_four_temperature<double>(reference_scales).model;
        const auto air11_s = canary::make_air11_three_temperature<double>(reference_scales).model;
        compare_models(air4, air4_s, air5_dissociated, std::array {1500.0, 7000.0, 12000.0});
        compare_models(air11, air11_s, air11_ionised, std::array {1500.0, 12000.0});
    };

    "multi-temperature: float states close with double parameters"_test = [&] {
        std::array<float, 11> y {};
        for (std::size_t i = 0; i < 11; ++i) { y[i] = static_cast<float>(air11_weakly_ionised[i]); }
        const auto mix = mixture(air11, fractions(y));
        const std::array t_m {3000.0F, 8000.0F};
        const auto c = thermo_state(air11,
                                    mix,
                                    canary::Density<float, si> {1.0e-3F},
                                    canary::Temperature<float, si> {10000.0F},
                                    mode_temperatures<si>(t_m));
        const std::array rho_e_m {
          1.0e-3F * c.state.mode_e[0].value(),
          1.0e-3F * c.state.mode_e[1].value(),
        };
        const auto back = thermo_state(air11,
                                       mix,
                                       canary::Density<float, si> {1.0e-3F},
                                       c.state.e,
                                       mode_energy_densities<si>(rho_e_m));
        expect(back.status == Status::Valid);
        // Float rounding of e (2e7 J/kg) is about 1 J/kg, which moves T by
        // about 1e-3 K; 1e-5 relative leaves a margin of 100.
        expect(std::abs(back.state.temp.value() - 10000.0F) <= 1.0e-5F * 10000.0F)
          << back.state.temp.value();
        expect(std::abs(back.state.mode_temp[0].value() - 3000.0F) <= 1.0e-5F * 3000.0F);
        expect(std::abs(back.state.mode_temp[1].value() - 8000.0F) <= 1.0e-5F * 8000.0F);
    };

    "multi-temperature: the electron pressure is rho y_e R_e T_e"_test = [&] {
        const auto mix = mixture(air11, fractions(air11_ionised));
        const std::array t_m {9000.0, 12000.0};
        const auto c =
          thermo_state(air11, mix, Rho {0.01}, Temp {14000.0}, mode_temperatures<si>(t_m));
        const double r_e = air11.thermally_perfect.gas_constant[10].value();
        const double p_e = canary::electron_pressure(air11, mix, c.state).value();
        // A product of four values: 4 eps.
        expect(relative(p_e, 0.01 * mix.fractions[10].value() * r_e * 12000.0, 4.0 * eps));
        double r_h = 0;
        for (std::size_t i = 0; i < 10; ++i) {
            r_h += mix.fractions[i].value() * air11.thermally_perfect.gas_constant[i].value();
        }
        // p - p_e = rho R_h T; p_e is 1e-4 of p, so the difference keeps
        // the relative rounding of the sums, 16 eps.
        expect(relative(c.state.p.value() - p_e, 0.01 * r_h * 14000.0, 16.0 * eps));
        const auto neutral = mixture(air11, fractions(air11_neutral));
        const auto d =
          thermo_state(air11, neutral, Rho {0.01}, Temp {14000.0}, mode_temperatures<si>(t_m));
        expect(canary::electron_pressure(air11, neutral, d.state).value() == 0.0);
    };

    "multi-temperature: the charge-neutral fractions balance the ion charge with electrons"_test =
      [&] {
          const auto& species = built11.species;
          const auto y =
            canary::charge_neutral_fractions(air11, species, fractions(air11_weakly_ionised_input));
          static_assert(std::same_as<std::remove_cv_t<decltype(y)>,
                                     canary::ChargeNeutralFractions<double, 11>>);
          expect(y.status == Status::Valid);
          double sum = 0;
          double charge = 0;
          double charge_scale = 0;
          for (std::size_t i = 0; i < 11; ++i) {
              // Reference compositions were made by the same formula: a few
              // roundings, 8 eps.
              expect(relative(y.fractions[i].value(), air11_weakly_ionised[i], 8.0 * eps))
                << i << y.fractions[i].value() << air11_weakly_ionised[i];
              sum += y.fractions[i].value();
              const double q =
                y.fractions[i].value() * species.charge[i].value() / species.mass[i].value();
              charge += q;
              charge_scale += std::abs(q);
          }
          // Sums of 11 terms: 16 eps.
          expect(std::abs(sum - 1.0) <= 16.0 * eps) << sum;
          expect(std::abs(charge) <= 16.0 * eps * charge_scale) << charge;
          const auto unchanged =
            canary::charge_neutral_fractions(air4, built4.species, fractions(air5_dissociated));
          expect(unchanged.status == Status::Valid);
          for (std::size_t i = 0; i < 5; ++i) {
              expect(unchanged.fractions[i].value() == air5_dissociated[i]);
          }
      };

    "multi-temperature: the eleven-species thermally perfect model has the five-species data"_test =
      [&] {
          for (std::size_t i = 0; i < 5; ++i) {
              expect(single11.gas_constant[i] == single5.gas_constant[i]) << i;
              expect(single11.reference_energy[i] == single5.reference_energy[i]) << i;
              for (const double t : {100.0, 300.0, 1000.0, 5000.0, 12000.0, 25000.0, 60000.0}) {
                  const auto a = canary::species_properties(single11, Temp {t});
                  const auto b = canary::species_properties(single5, Temp {t});
                  // The same data and the same arithmetic: equal.
                  expect(a.energy[i] == b.energy[i]) << i << t;
                  expect(a.cv[i] == b.cv[i]) << i << t;
                  expect(a.entropy[i] == b.entropy[i]) << i << t;
              }
          }
          for (std::size_t i = 0; i < 5; ++i) {
              const auto& a = canary::air11_electronic_levels[i];
              const auto& b = canary::air5_electronic_levels[i];
              expect(std::string(a.name) == b.name && a.levels == b.levels);
              expect(a.degeneracies == b.degeneracies && a.temperatures == b.temperatures) << i;
          }
      };

    "multi-temperature: status reports non-convergence, invalid states and composition violations"_test =
      [&] {
          const auto mix = mixture(air11, fractions(air11_ionised));
          const std::array t_m {4000.0, 9000.0};
          const auto base =
            thermo_state(air11, mix, Rho {1.0}, Temp {9000.0}, mode_temperatures<si>(t_m));
          const std::array rho_e_m {base.state.mode_e[0].value(), base.state.mode_e[1].value()};
          const auto energies = mode_energy_densities<si>(rho_e_m);
          Air11 one_iteration = air11;
          one_iteration.thermally_perfect.max_iterations = 1;
          expect(thermo_state(one_iteration, mix, Rho {1.0}, base.state.e, energies).status ==
                 Status::NotConverged);
          expect(thermo_state(one_iteration, mix, Rho {1.0}, base.state.e).status ==
                 Status::NotConverged);
          expect(thermo_state(air11, mix, Rho {0.0}, base.state.e, energies).status ==
                 Status::Invalid);
          expect(
            thermo_state(air11, mix, Rho {-1.0}, base.state.e, mode_temperatures<si>(t_m)).status ==
            Status::Invalid);
          expect(thermo_state(air11, mix, Rho {1.0}, E {no_value}, energies).status ==
                 Status::Invalid);
          const std::array negative {-5.0, 9000.0};
          expect(thermo_state(air11, mix, Rho {1.0}, Temp {9000.0}, mode_temperatures<si>(negative))
                   .status == Status::Invalid);
          // Mode energies below e_m(temperature_min) and above
          // e_m(temperature_max).
          for (const double out_of_range : {-1.0e8, 1.0e12}) {
              const std::array bad {rho_e_m[0], out_of_range};
              expect(
                thermo_state(air11, mix, Rho {1.0}, Temp {9000.0}, mode_energy_densities<si>(bad))
                  .status == Status::Invalid)
                << out_of_range;
              expect(thermo_state(air11,
                                  mix,
                                  P {base.state.p},
                                  Temp {9000.0},
                                  mode_energy_densities<si>(bad))
                       .status == Status::Invalid)
                << out_of_range;
          }
          auto violated = air11_ionised;
          violated[3] += 1.0e-10;
          expect(
            thermo_state(air11, fractions(violated), Rho {1.0}, base.state.e, energies).status ==
            Status::CompositionViolation);
          expect(thermo_state(air11, fractions(violated), Rho {1.0}, base.state.e).status ==
                 Status::CompositionViolation);
      };

    "multi-temperature: partial densities give the same closure as the equivalent fractions"_test =
      [&] {
          const Rho rho {2.0};
          std::array<double, 11> rho_i {};
          for (std::size_t i = 0; i < 11; ++i) { rho_i[i] = rho.value() * air11_ionised[i]; }
          const auto view = canary::as_kind<canary::kind::density, si>(
            Kokkos::mdspan<const double, Extents>(rho_i.data(), 11));
          const std::array t_m {2500.0, 7000.0};
          const auto a = thermo_state(air11, view, rho, Temp {3000.0}, mode_temperatures<si>(t_m));
          const auto b = thermo_state(air11,
                                      fractions(air11_ionised),
                                      rho,
                                      Temp {3000.0},
                                      mode_temperatures<si>(t_m));
          expect(a.status == Status::Valid);
          // The fractions differ by the rounding of rho_i / rho.
          expect(relative(a.state.e.value(), b.state.e.value(), 1.0e-12));
          expect(relative(a.state.p.value(), b.state.p.value(), 1.0e-12));
          expect(relative(a.state.mode_e[1].value(), b.state.mode_e[1].value(), 1.0e-12));
      };
}
