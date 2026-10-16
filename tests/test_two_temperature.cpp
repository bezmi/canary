#include "dual.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/data/air5.hpp>
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
using Air = canary::Air5TwoTemperature<double, si>;
using AirScaled = canary::Air5TwoTemperature<double, scaled>;
using AirSingle = canary::Air5<double, si>;

constexpr double eps = std::numeric_limits<double>::epsilon();
// Relative tolerance of the temperature solves (the model default).
constexpr double solve_tolerance = 1.0e-10;

const canary::UnitScales<double, scaled> reference_scales {{
  .mass = 4.65e-26,
  .density = 1.0e-3,
  .temperature = 1000.0,
  .length = 1.0,
}};

template <typename T> auto fractions(const std::array<T, 5>& y)
{
    return canary::as_kind<canary::kind::mass_fraction>(
      Kokkos::mdspan<const T, Extents>(y.data(), 5));
}

// Views of one mode over the storage of one value.
template <canary::UnitSystem U, typename T> auto mode_temperatures(const std::array<T, 1>& t_ve)
{
    return canary::as_kind<canary::kind::mode_temperature, U>(
      Kokkos::mdspan<const T, Extents>(t_ve.data(), 1));
}

template <canary::UnitSystem U, typename T>
auto mode_energy_densities(const std::array<T, 1>& rho_e_ve)
{
    return canary::as_kind<canary::kind::mode_energy_density, U>(
      Kokkos::mdspan<const T, Extents>(rho_e_ve.data(), 1));
}

bool relative(double a, double b, double bound) { return std::abs(a - b) <= bound * std::abs(b); }

// Energies are compared relative to the larger of the value and 1 J/kg.
bool energy_relative(double a, double b, double bound)
{ return std::abs(a - b) <= bound * std::max(std::abs(b), 1.0); }

constexpr std::array<double, 5> air {0.767, 0.233, 0.0, 0.0, 0.0};
constexpr std::array<double, 5> dissociated {0.3, 0.1, 0.3, 0.25, 0.05};

// A state of an independent implementation of the model with the same
// species data and the universal gas constant 8.31451 J/(mol K): GDTk's
// TwoTemperatureGas (gdtk-uq/gdtk commit bf559c3f, libgas.so through its C
// interface, five-species air). Values: mass fractions, pressure, both
// temperatures, density, internal energy per unit mass without (u) and of
// (u_ve) the vibrational-electronic mode. Sound speeds are calculated separately
// for fixed specific mode energy: a^2 = (1 + R/cv_tr) R T, with R = sum(y R_u/M)
// and cv_tr = sum(y (3 + rotational_degrees) R_u/(2 M)), using R_u = 8.31451.
struct ReferenceState {
    std::array<double, 5> y;
    double p;
    double temp;
    double t_ve;
    double rho;
    double u;
    double u_ve;
    double a;
};

constexpr std::array<ReferenceState, 5> reference_states {
  {
    {
      .y = air,
      .p = 101325.0,
      .temp = 300.0,
      .t_ve = 300.0,
      .rho = 1.1719636407483949,
      .u = -84591.42008858029,
      .u_ve = 5.015268335709959,
      .a = 347.90867262549483,
    },
    {
      .y = air,
      .p = 101325.0,
      .temp = 5000.0,
      .t_ve = 2000.0,
      .rho = 0.07031781844490369,
      .u = 3301659.110231412,
      .u_ve = 251464.28084760078,
      .a = 1420.3312083691005,
    },
    {
      .y = air,
      .p = 1000.0,
      .temp = 9000.0,
      .t_ve = 12000.0,
      .rho = 0.0003855460615999325,
      .u = 6183574.455184597,
      .u_ve = 4650385.702086919,
      .a = 1905.5742794866362,
    },
    {
      .y = dissociated,
      .p = 5000.0,
      .temp = 8000.0,
      .t_ve = 6000.0,
      .rho = 0.0014305926434206111,
      .u = 20085845.70890493,
      .u_ve = 853034.3232539449,
      .a = 2332.8440603648296,
    },
    {
      .y = dissociated,
      .p = 5000.0,
      .temp = 3000.0,
      .t_ve = 3000.0,
      .rho = 0.0038149137157882955,
      .u = 16164831.924025863,
      .u_ve = 241567.21495804077,
      .a = 1428.5693993440777,
    },
  },
};

// Every energy and pressure is proportional to the universal gas constant.
const double gas_constant_ratio = canary::constants::universal_gas_constant<double> / 8.31451;

}  // namespace

static_assert(noexcept(canary::species_vibration_mode(std::declval<const Air&>(),
                                                      canary::SpeciesIndex {0})));
static_assert(noexcept(canary::species_electronic_mode(std::declval<const Air&>(),
                                                       canary::SpeciesIndex {0})));
static_assert(noexcept(canary::electron_species_index(std::declval<const Air&>())));
static_assert(noexcept(canary::electron_energy_mode(std::declval<const Air&>())));

static_assert(canary::ThermoStateOf<Air::State<double>, double, si>);
static_assert(canary::GasModel<Air, double>);
static_assert(canary::MixtureGasModel<Air, double>);
static_assert(canary::MixtureGasModel<Air, float>);
static_assert(canary::MixtureGasModel<Air, D>);
static_assert(canary::MultiTemperatureGasModel<Air, double>);
static_assert(canary::MultiTemperatureGasModel<Air, float>);
static_assert(canary::MultiTemperatureGasModel<Air, D>);
static_assert(canary::MultiTemperatureGasModel<AirScaled, double>);
static_assert(canary::MultiTemperatureGasModel<AirScaled, float>);
static_assert(canary::MultiTemperatureGasModel<AirScaled, D>);
static_assert(!canary::MultiTemperatureGasModel<AirSingle, double>);
static_assert(canary::AcceptsInput<Air, canary::ModeTemperatureView<double, si>>);
static_assert(canary::AcceptsInput<Air, canary::ModeEnergyDensityView<double, si>>);
static_assert(!canary::AcceptsInput<AirScaled, canary::ModeEnergyDensityView<double, si>>);
static_assert(std::is_trivially_copyable_v<Air>);
static_assert(std::is_trivially_copyable_v<AirScaled>);
static_assert(std::is_trivially_copyable_v<Air::State<double>>);
static_assert(std::is_trivially_copyable_v<Air::State<D>>);
static_assert(std::is_trivially_copyable_v<Air::Mixture<double>>);
static_assert(std::is_trivially_copyable_v<Air::Mixture<D>>);

void test_two_temperature()
{
    using namespace boost::ut;
    using canary::mixture;
    using canary::thermo_state;
    const auto built = canary::make_air5_two_temperature<double>(canary::UnitScales<double, si> {});
    const Air& gas = built.model;
    const AirSingle single = canary::make_air5<double>(canary::UnitScales<double, si> {}).model;
    const double r = gas_constant_ratio;

    "bound-electronic membership in two-temperature models"_test = [] {
        const auto check =
          []<std::floating_point P, canary::UnitSystem U>(const canary::UnitScales<P, U>& scales) {
              const auto model = canary::make_air5_two_temperature<P>(scales).model;
              for (const auto i : canary::species_indices(0, 5)) {
                  expect(canary::species_electronic_mode(model, i) ==
                         canary::OptionalModeIndex {canary::ModeIndex {0}});
              }
              expect(!canary::species_electronic_mode(model, canary::SpeciesIndex {5}));
              expect(!canary::species_electronic_mode(
                model,
                canary::SpeciesIndex {std::numeric_limits<std::uint32_t>::max() - 1}));
          };
        check(canary::UnitScales<float, si> {});
        check(canary::UnitScales<double, si> {});
        check(canary::UnitScales<float, scaled> {
          {.mass = 4.65e-26, .density = 1.0e-3, .temperature = 1000.0, .length = 1.0}});
        check(canary::UnitScales<double, scaled> {
          {.mass = 4.65e-26, .density = 1.0e-3, .temperature = 1000.0, .length = 1.0}});
    };

    "typed mode views preserve two-temperature closures"_test = [&] {
        const std::array<double, 3> data {4000, -1, 5000};
        using Mapping = Kokkos::layout_stride::mapping<Extents>;
        const Kokkos::mdspan<const double, Extents, Kokkos::layout_stride> storage(
          data.data() + 2,
          Mapping(Extents(1), std::array<std::size_t, 1> {2}));
        const auto values = canary::as_kind<canary::kind::mode_temperature, si>(storage);
        const auto modes = canary::as_indexed<canary::ModeIndex>(values);
        const auto mix = mixture(gas, fractions(air));
        const auto a = thermo_state(gas, mix, modes, Temp {4000}, Rho {1});
        const auto b = thermo_state(gas, mix, Rho {1}, Temp {4000}, values);
        expect(a.status == b.status);
        expect(a.state.e == b.state.e && a.state.p == b.state.p);
        expect(a.state.mode_temp[0] == b.state.mode_temp[0]);
    };

    "two-temperature diagnostics distinguish mode and equilibrium solves"_test = [&] {
        const std::array<double, 5> y {1, 0, 0, 0, 0};
        const auto mix = mixture(gas, fractions(y));
        const std::array<double, 1> energies {1e100};
        const auto mode =
          thermo_state(gas, mix, Rho {1}, Temp {4000}, mode_energy_densities<si>(energies));
        expect(mode.status == Status::Invalid);
        expect(mode.diagnostic.reason == canary::ClosureReason::NoTemperatureBracket);
        expect(mode.diagnostic.mode == canary::OptionalModeIndex {canary::ModeIndex {0}});
        const auto equilibrium = thermo_state(gas, mix, Rho {1}, E {1e100});
        expect(equilibrium.status == Status::Invalid);
        expect(equilibrium.diagnostic.reason == canary::ClosureReason::NoTemperatureBracket);
        expect(!equilibrium.diagnostic.mode.has_value());
    };

    "two-temperature properties expose separate energy derivatives"_test = [&] {
        const std::array<D, 5> y {D(0.7), D(0.2), D(0.05), D(0.03), D(0.02)};
        const auto mix = mixture(gas, fractions(y));
        const std::array<D, 1> modes {D(5000)};
        const auto heavy = thermo_state(gas,
                                        mix,
                                        canary::Density<D, si> {D(0.01)},
                                        canary::Temperature<D, si> {D(4000, 1)},
                                        mode_temperatures<si>(modes));
        const auto props = canary::properties(gas, mix, heavy.state);
        static_assert(canary::MultiTemperaturePropertiesOf<decltype(props), D, si>);
        static_assert(noexcept(canary::properties(gas, mix, heavy.state)));
        expect(relative(props.cv_heavy.value().real(), heavy.state.e.value().dual(), 1.0e-11));
        const std::array<D, 1> seeded_modes {D(5000, 1)};
        const auto mode = thermo_state(gas,
                                       mix,
                                       canary::Density<D, si> {D(0.01)},
                                       canary::Temperature<D, si> {D(4000)},
                                       mode_temperatures<si>(seeded_modes));
        expect(
          relative(props.mode_cv[0].value().real(), mode.state.mode_e[0].value().dual(), 1.0e-11));
        expect(props.gas_constant_electron.value().real() == 0.0);
        expect(props.gas_constant_heavy == mix.gas_constant);
    };

    "two-temperature acoustics match compression with fixed specific mode energy"_test = [] {
        const auto run = []<std::floating_point V, canary::UnitSystem U>(
                           const canary::UnitScales<V, U>& scales) {
            using Dual = canary::test::Dual<V>;
            const auto model = canary::make_air5_two_temperature<V>(scales).model;
            const V tolerance = V(256) * std::numeric_limits<V>::epsilon();
            for (const auto& composition : {air, dissociated}) {
                std::array<Dual, 5> y {};
                for (const canary::SpeciesIndex i : canary::species_indices(0, y.size())) {
                    y[i.value()] = Dual(V(composition[i.value()]));
                }
                const auto mix = mixture(model, fractions(y));
                const V rho = canary::convert<U>(canary::Density<V, si> {V(0.01)}, scales).value();
                for (const V temperature : {V(600), V(4000)}) {
                    const V temp =
                      canary::convert<U>(canary::Temperature<V, si> {temperature}, scales).value();
                    V first_speed = V(0);
                    for (const V mode_temperature : {V(500), V(5000)}) {
                        const V mode_temp =
                          canary::convert<U>(canary::ModeTemperature<V, si> {mode_temperature},
                                             scales)
                            .value();
                        const std::array<Dual, 1> modes {Dual(mode_temp)};
                        const auto base = thermo_state(model,
                                                       mix,
                                                       canary::Density<Dual, U> {Dual(rho)},
                                                       canary::Temperature<Dual, U> {Dual(temp)},
                                                       mode_temperatures<U>(modes));
                        expect(fatal(base.status == Status::Valid));
                        const auto props = canary::properties(model, mix, base.state);
                        const V a = props.acoustics.sound_speed.value().real();
                        if (first_speed == V(0)) { first_speed = a; }
                        expect(a == first_speed);

                        // Continuity and pressure work give de/drho = p/rho^2.
                        // Fixed composition and T_ve keep specific mode energy fixed.
                        const Dual density(rho, V(1));
                        const Dual energy(base.state.e.value().real(),
                                          base.state.p.value().real() / (rho * rho));
                        const auto compressed =
                          thermo_state(model,
                                       mix,
                                       canary::Density<Dual, U> {density},
                                       canary::EnergyDensity<Dual, U> {density * energy},
                                       mode_temperatures<U>(modes));
                        expect(fatal(compressed.status == Status::Valid));
                        expect(relative(compressed.state.p.value().dual(), a * a, tolerance));
                        expect(relative(props.acoustics.acoustic_exponent.value().real(),
                                        rho * a * a / base.state.p.value().real(),
                                        tolerance));

                        auto changed_pressure = base.state;
                        changed_pressure.p =
                          canary::Pressure<Dual, U> {base.state.p.value() * V(2)};
                        const auto changed = canary::properties(model, mix, changed_pressure);
                        expect(changed.acoustics.sound_speed.value().real() == a);
                        expect(relative(changed.acoustics.acoustic_exponent.value().real(),
                                        props.acoustics.acoustic_exponent.value().real() / V(2),
                                        tolerance));
                    }
                }
            }
        };
        run.template operator()<float>(canary::UnitScales<float, si> {});
        run.template operator()<double>(canary::UnitScales<double, si> {});
        run.template operator()<float>(canary::UnitScales<float, scaled> {
          {.mass = 4.65e-26F, .density = 1.0e-3F, .temperature = 1000.0F, .length = 1.0F}});
        run.template operator()<double>(reference_scales);
    };

    "two-temperature closures validate pressure with and without explicit modes"_test = [&] {
        const auto check = [&]<typename V> {
            const std::array<V, 5> y {V(1), V(0), V(0), V(0), V(0)};
            const auto mix = mixture(gas, fractions(y));
            const canary::Density<V, si> rho {std::numeric_limits<V>::max() / V(16)};
            const canary::Temperature<V, si> temp {V(300)};
            const std::array<V, 1> mode {V(300)};
            const auto result = thermo_state(gas, mix, rho, temp, mode_temperatures<si>(mode));
            expect(result.status == Status::Invalid);
            expect(std::isinf(result.state.p.value()));
            expect(std::isfinite(result.state.e.value()));
            expect(std::isfinite(result.state.mode_e[0].value()));
            expect(thermo_state(gas, mix, rho, temp).status == Status::Invalid);
            const auto inverse =
              thermo_state(gas, mix, rho, result.state.e, mode_temperatures<si>(mode));
            expect(inverse.status == Status::Invalid);
            expect(std::isinf(inverse.state.p.value()));
            const canary::Pressure<V, si> supplied {std::numeric_limits<V>::max() / V(4)};
            const canary::Temperature<V, si> cold {V(1) / mix.gas_constant.value()};
            const auto retained =
              thermo_state(gas, mix, supplied, cold, mode_temperatures<si>(mode));
            expect(retained.status == Status::Valid);
            expect(retained.state.p == supplied);
            expect(thermo_state(gas, mix, supplied, cold).status == Status::Valid);
        };
        check.operator()<float>();
        check.operator()<double>();
    };

    "two-temperature: the builder holds the single-temperature model and the species data"_test =
      [&] {
          expect(built.species.n == 5_i);
          for (std::size_t i = 0; i < 5; ++i) {
              const double r_i = single.gas_constant[i].value();
              expect(gas.thermally_perfect.gas_constant[i] == single.gas_constant[i]);
              const double d = canary::air5_rotational_degrees[i];
              expect(relative(gas.translational_rotational_cp[i].value(),
                              (2.5 + (0.5 * d)) * r_i,
                              2.0 * eps));
          }
          expect(gas.thermally_perfect.temperature_min.value() == 10.0);
          expect(gas.thermally_perfect.max_iterations == 30_i);
      };

    "two-temperature: the builder rejects null or empty species names"_test = [] {
        for (const char* name : std::array<const char*, 2> {nullptr, ""}) {
            auto species = canary::air5_species;
            species[1].name = name;
            try {
                canary::make_two_temperature_mixture(species,
                                                     canary::air5_rotational_degrees,
                                                     canary::UnitScales<double, si> {});
                expect(false) << "no exception";
            } catch (const canary::ParameterError& error) {
                expect(std::string(error.what())
                         .starts_with("species 1: name must not be null or empty"))
                  << error.what();
            }
        }
    };

    "two-temperature: electron rejection follows the role rather than the name"_test = [] {
        auto species = canary::air5_species;
        species[3].name = "e-";
        expect(nothrow([&] {
            canary::make_two_temperature_mixture(species,
                                                 canary::air5_rotational_degrees,
                                                 canary::UnitScales<double, si> {});
        }));
        species[3].name = "free-electron";
        species[3].metadata.role = canary::SpeciesRole::electron;
        try {
            canary::make_two_temperature_mixture(species,
                                                 canary::air5_rotational_degrees,
                                                 canary::UnitScales<double, si> {});
            expect(false) << "electron-role species accepted";
        } catch (const canary::ParameterError& error) {
            expect(
              std::string_view(error.what()).starts_with("species free-electron: free electrons"));
        }
    };

    "two-temperature: the builder rejects rotational degrees other than 0, 2 and 3 and free electrons"_test =
      [] {
          for (const int d : {1, 4, -2}) {
              auto degrees = canary::air5_rotational_degrees;
              degrees[2] = d;
              try {
                  canary::make_two_temperature_mixture(canary::air5_species,
                                                       degrees,
                                                       canary::UnitScales<double, si> {});
                  expect(false) << "no exception";
              } catch (const canary::ParameterError& error) {
                  expect(
                    std::string(error.what())
                      .starts_with("species N: rotational degrees of freedom must be 0, 2 or 3"))
                    << error.what();
              }
          }
          auto species = canary::air5_species;
          species[3].name = "e-";
          species[3].metadata.role = canary::SpeciesRole::electron;
          try {
              canary::make_two_temperature_mixture(species,
                                                   canary::air5_rotational_degrees,
                                                   canary::UnitScales<double, si> {});
              expect(false) << "no exception";
          } catch (const canary::ParameterError& error) {
              expect(std::string(error.what()).starts_with("species e-: free electrons"))
                << error.what();
          }
          species = canary::air5_species;
          species[1].molar_mass = -1.0;
          try {
              canary::make_two_temperature_mixture(species,
                                                   canary::air5_rotational_degrees,
                                                   canary::UnitScales<double, si> {});
              expect(false) << "no exception";
          } catch (const canary::ParameterError& error) {
              expect(std::string(error.what()).starts_with("species O2: molar mass"))
                << error.what();
          }
      };

    "two-temperature: (p, T, T_ve) gives the reference density, energies and sound speed"_test =
      [&] {
          for (const ReferenceState& ref : reference_states) {
              const auto mix = mixture(gas, fractions(ref.y));
              const std::array t_ve {ref.t_ve};
              const auto c =
                thermo_state(gas, mix, P {r * ref.p}, Temp {ref.temp}, mode_temperatures<si>(t_ve));
              expect(c.status == Status::Valid);
              expect(relative(c.state.rho.value(), ref.rho, 1.0e-12))
                << ref.temp << c.state.rho.value() << ref.rho;
              expect(energy_relative(c.state.e.value(), r * (ref.u + ref.u_ve), 1.0e-12))
                << ref.temp << c.state.e.value() << r * (ref.u + ref.u_ve);
              // e_ve is a small difference of terms of size sum(y (|dh_f| + cp_tr T_ve))
              // (5 J/kg from terms near 1e5 J/kg at 300 K), so its rounding
              // in either implementation is relative to that size.
              double cancelled = 0.0;
              for (std::size_t i = 0; i < 5; ++i) {
                  cancelled += ref.y[i] * (std::abs(gas.formation_enthalpy[i].value()) +
                                           (gas.translational_rotational_cp[i].value() * ref.t_ve));
              }
              expect(std::abs(c.state.mode_e[0].value() - (r * ref.u_ve)) <= 1.0e-12 * cancelled)
                << ref.temp << ref.t_ve << c.state.mode_e[0].value() << r * ref.u_ve;
              expect(c.state.mode_temp[0].value() == ref.t_ve);
              const auto props = canary::properties(gas, mix, c.state);
              expect(relative(props.acoustics.sound_speed.value(), std::sqrt(r) * ref.a, 1.0e-10))
                << ref.temp << props.acoustics.sound_speed.value() << std::sqrt(r) * ref.a;
          }
      };

    "two-temperature: (rho, rho e, rho e_ve) gives the reference temperatures and pressure"_test =
      [&] {
          for (const ReferenceState& ref : reference_states) {
              const std::array rho_e_ve {ref.rho * r * ref.u_ve};
              const auto c = thermo_state(gas,
                                          fractions(ref.y),
                                          Rho {ref.rho},
                                          RhoE {ref.rho * r * (ref.u + ref.u_ve)},
                                          mode_energy_densities<si>(rho_e_ve));
              expect(c.status == Status::Valid);
              expect(relative(c.state.temp.value(), ref.temp, 1.0e-8))
                << c.state.temp.value() << ref.temp;
              expect(relative(c.state.mode_temp[0].value(), ref.t_ve, 1.0e-8))
                << c.state.mode_temp[0].value() << ref.t_ve;
              expect(relative(c.state.p.value(), r * ref.p, 1.0e-8))
                << c.state.p.value() << r * ref.p;
          }
      };

    "two-temperature: e_ve of every species is zero at the reference temperature"_test = [&] {
        for (std::size_t i = 0; i < 5; ++i) {
            std::array<double, 5> y {};
            y[i] = 1.0;
            const std::array t_ve {298.15};
            const auto c = thermo_state(gas,
                                        fractions(y),
                                        Rho {1.0},
                                        Temp {1000.0},
                                        mode_temperatures<si>(t_ve));
            expect(std::abs(c.state.mode_e[0].value()) <=
                   4.0 * eps * std::abs(gas.formation_enthalpy[i].value()))
              << i << c.state.mode_e[0].value();
        }
    };

    "two-temperature: equal temperatures reproduce the single-temperature model"_test = [&] {
        for (const auto& y : {air, dissociated}) {
            const auto mix = mixture(gas, fractions(y));
            for (const double t : {300.0, 2000.0, 8000.0, 15000.0}) {
                const Rho rho {0.5};
                const std::array t_ve {t};
                const auto modes = mode_temperatures<si>(t_ve);
                const auto ref = thermo_state(single, mix, rho, Temp {t});
                // Energy scale of the rounding: the terms of e_i summed.
                double e_scale = 0;
                for (std::size_t i = 0; i < 5; ++i) {
                    e_scale += y[i] * (std::abs(gas.formation_enthalpy[i].value()) +
                                       (gas.translational_rotational_cp[i].value() * t));
                }
                const double e_bound = 16.0 * eps * e_scale;

                const auto a = thermo_state(gas, mix, rho, Temp {t}, modes);
                expect(std::abs(a.state.e.value() - ref.state.e.value()) <= e_bound)
                  << t << a.state.e.value() << ref.state.e.value();
                expect(a.state.p == ref.state.p);

                const auto b = thermo_state(gas, mix, rho, ref.state.p, modes);
                expect(std::abs(b.state.e.value() - ref.state.e.value()) <= e_bound);
                expect(relative(b.state.temp.value(), t, 4.0 * eps));

                const auto c = thermo_state(gas, mix, ref.state.p, Temp {t}, modes);
                expect(std::abs(c.state.e.value() - ref.state.e.value()) <= e_bound);
                expect(relative(c.state.rho.value(), rho.value(), 4.0 * eps));

                const auto d = thermo_state(gas, mix, rho, ref.state.e, modes);
                const auto d_ref = thermo_state(single, mix, rho, ref.state.e);
                expect(relative(d.state.temp.value(), d_ref.state.temp.value(), solve_tolerance))
                  << t << d.state.temp.value() << d_ref.state.temp.value();
                expect(relative(d.state.p.value(), d_ref.state.p.value(), solve_tolerance));
                expect(d.state.e == d_ref.state.e);
            }
        }
    };

    "two-temperature: closures without mode input are the single-temperature closures"_test = [&] {
        for (const auto& y : {air, dissociated}) {
            const auto mix = mixture(gas, fractions(y));
            for (const double t : {300.0, 2000.0, 8000.0, 15000.0}) {
                const Rho rho {0.5};
                const auto base = thermo_state(single, mix, rho, Temp {t});
                const auto same = [&](const auto& two, const auto& one) {
                    return two.state.rho == one.state.rho && two.state.e == one.state.e &&
                           two.state.p == one.state.p && two.state.temp == one.state.temp &&
                           two.status == one.status &&
                           two.state.mode_temp[0].value() == one.state.temp.value();
                };
                const auto a = thermo_state(gas, mix, rho, base.state.e);
                expect(same(a, thermo_state(single, mix, rho, base.state.e)));
                expect(same(thermo_state(gas, mix, rho, base.state.p),
                            thermo_state(single, mix, rho, base.state.p)));
                const auto c = thermo_state(gas, mix, rho, Temp {t});
                expect(same(c, base));
                expect(same(thermo_state(gas, mix, base.state.p, Temp {t}),
                            thermo_state(single, mix, base.state.p, Temp {t})));
                // mode_e is e_ve(T): the (rho, T, T_ve = T) closure.
                const std::array t_ve {t};
                const auto with_modes =
                  thermo_state(gas, mix, rho, Temp {t}, mode_temperatures<si>(t_ve));
                expect(c.state.mode_e[0] == with_modes.state.mode_e[0]);
            }
        }
    };

    "two-temperature: every closure family inverts (rho, T, T_ve)"_test = [&] {
        for (const auto& y : {air, dissociated}) {
            const auto mix = mixture(gas, fractions(y));
            for (const double rho_value : {1.0e-4, 1.0}) {
                for (const double t : {300.0, 2000.0, 8000.0, 15000.0}) {
                    for (const double t_ve_value : {300.0, 2000.0, 8000.0, 15000.0}) {
                        const Rho rho {rho_value};
                        const std::array t_ve {t_ve_value};
                        const auto temps = mode_temperatures<si>(t_ve);
                        const auto base = thermo_state(gas, mix, rho, Temp {t}, temps);
                        expect(base.status == Status::Valid);
                        const E e = base.state.e;
                        const P p = base.state.p;
                        const double e_ve = base.state.mode_e[0].value();
                        const std::array rho_e_ve {rho_value * e_ve};
                        const auto energies = mode_energy_densities<si>(rho_e_ve);
                        const auto check = [&](const auto& c) {
                            expect(c.status == Status::Valid) << rho_value << t << t_ve_value;
                            expect(relative(c.state.rho.value(), rho_value, solve_tolerance));
                            expect(relative(c.state.temp.value(), t, solve_tolerance))
                              << rho_value << t << t_ve_value << c.state.temp.value();
                            expect(
                              relative(c.state.mode_temp[0].value(), t_ve_value, solve_tolerance))
                              << rho_value << t << t_ve_value << c.state.mode_temp[0].value();
                            expect(energy_relative(c.state.e.value(), e.value(), solve_tolerance));
                            expect(
                              energy_relative(c.state.mode_e[0].value(), e_ve, solve_tolerance));
                            expect(relative(c.state.p.value(), p.value(), solve_tolerance));
                        };
                        check(thermo_state(gas, mix, rho, e, temps));
                        check(thermo_state(gas, mix, rho, p, temps));
                        check(thermo_state(gas, mix, p, Temp {t}, temps));
                        check(thermo_state(gas, mix, rho, e, energies));
                        check(thermo_state(gas, mix, rho, RhoE {rho_value * e.value()}, energies));
                        check(thermo_state(gas, mix, rho, p, energies));
                        check(thermo_state(gas, mix, rho, Temp {t}, energies));
                        check(thermo_state(gas, mix, p, Temp {t}, energies));
                        if (t == t_ve_value) {
                            check(thermo_state(gas, mix, rho, e));
                            check(thermo_state(gas, mix, rho, p));
                            check(thermo_state(gas, mix, rho, Temp {t}));
                            check(thermo_state(gas, mix, p, Temp {t}));
                        }
                    }
                }
            }
        }
    };

    "two-temperature: dual numbers give the mode heat capacities"_test = [&] {
        const std::array<D, 5> y_d {0.3, 0.1, 0.3, 0.25, 0.05};
        const auto mix_d = mixture(gas, fractions(y_d));
        const auto mix = mixture(gas, fractions(dissociated));
        for (const double t : {600.0, 3000.0, 9000.0}) {
            for (const double t_ve_value : {500.0, 2500.0, 12000.0}) {
                // Mixture cv_ve and cv_tr from the single-temperature species
                // cv_i = cp_i - R_i and cp_tr,i.
                const auto sp = canary::species_properties(single, Temp {t_ve_value});
                double cv_ve = 0;
                double cv_tr = 0;
                for (std::size_t i = 0; i < 5; ++i) {
                    const double r_i = single.gas_constant[i].value();
                    const double cp_tr = gas.translational_rotational_cp[i].value();
                    cv_ve += dissociated[i] * (sp.cv[i].value() + r_i - cp_tr);
                    cv_tr += dissociated[i] * (cp_tr - r_i);
                }
                const std::array t_ve {t_ve_value};
                const auto base =
                  thermo_state(gas, mix, Rho {1.0}, Temp {t}, mode_temperatures<si>(t_ve));

                // d e_ve / d T_ve.
                const std::array t_ve_d {D(t_ve_value, 1.0)};
                const auto from_t = thermo_state(gas,
                                                 mix_d,
                                                 canary::Density<D, si> {D(1.0)},
                                                 canary::Temperature<D, si> {D(t)},
                                                 mode_temperatures<si>(t_ve_d));
                expect(relative(from_t.state.mode_e[0].value().dual(), cv_ve, 1.0e-12))
                  << t_ve_value << from_t.state.mode_e[0].value().dual() << cv_ve;

                // d T_ve / d e_ve at fixed rho and e: rho = 1, so rho e_ve
                // carries the derivative of e_ve.
                const std::array rho_e_ve {D(base.state.mode_e[0].value(), 1.0)};
                const auto from_e =
                  thermo_state(gas,
                               mix_d,
                               canary::Density<D, si> {D(1.0)},
                               canary::SpecificEnergy<D, si> {D(base.state.e.value())},
                               mode_energy_densities<si>(rho_e_ve));
                expect(from_e.status == Status::Valid);
                expect(relative(from_e.state.mode_temp[0].value().dual(), 1.0 / cv_ve, 1.0e-10))
                  << t_ve_value << from_e.state.mode_temp[0].value().dual() << 1.0 / cv_ve;

                // d T / d e at fixed e_ve.
                const std::array rho_e_ve_fixed {D(base.state.mode_e[0].value())};
                const auto from_e_t =
                  thermo_state(gas,
                               mix_d,
                               canary::Density<D, si> {D(1.0)},
                               canary::SpecificEnergy<D, si> {D(base.state.e.value(), 1.0)},
                               mode_energy_densities<si>(rho_e_ve_fixed));
                expect(relative(from_e_t.state.temp.value().dual(), 1.0 / cv_tr, 1.0e-12))
                  << from_e_t.state.temp.value().dual() << 1.0 / cv_tr;
            }
        }
    };

    "two-temperature: a scaled model gives the SI state after conversion"_test = [&] {
        // T_ve starts at 1000 K: near 300 K, e_ve is a difference of terms
        // about 1e3 to 1e7 times larger, and the different rounding of the
        // folded curves alone changes it by more than 1e-12 relative.
        const auto scaled_built = canary::make_air5_two_temperature<double>(reference_scales);
        const AirScaled& gas_s = scaled_built.model;
        for (const auto& y : {air, dissociated}) {
            const auto mix = mixture(gas, fractions(y));
            const auto mix_s = mixture(gas_s, fractions(y));
            for (const double rho_value : {1.0e-4, 1.0}) {
                for (const double t : {300.0, 2000.0, 8000.0, 15000.0}) {
                    for (const double t_ve_value : {1000.0, 2000.0, 8000.0, 15000.0}) {
                        const Rho rho {rho_value};
                        const std::array t_ve {t_ve_value};
                        const auto base =
                          thermo_state(gas, mix, rho, Temp {t}, mode_temperatures<si>(t_ve));
                        const RhoE rho_e {rho_value * base.state.e.value()};
                        const canary::ModeEnergyDensity<double, si> rho_e_ve {
                          rho_value * base.state.mode_e[0].value()};
                        const std::array rho_e_ve_si {rho_e_ve.value()};
                        const double rho_e_ve_scaled =
                          canary::convert<scaled>(rho_e_ve, reference_scales).value();
                        const double t_ve_scaled =
                          canary::convert<scaled>(canary::ModeTemperature<double, si> {t_ve_value},
                                                  reference_scales)
                            .value();
                        const std::array rho_e_ve_s {rho_e_ve_scaled};
                        const std::array t_ve_s {t_ve_scaled};
                        const auto rho_s = canary::convert<scaled>(rho, reference_scales);
                        const auto compare = [&](const auto& a, const auto& b) {
                            expect(b.status == Status::Valid);
                            const auto convert = [](auto q) {
                                return canary::convert<si>(q, reference_scales).value();
                            };
                            expect(relative(convert(b.state.temp), a.state.temp.value(), 1.0e-12))
                              << rho_value << t << t_ve_value;
                            expect(relative(convert(b.state.mode_temp[0]),
                                            a.state.mode_temp[0].value(),
                                            1.0e-12))
                              << rho_value << t << t_ve_value;
                            expect(relative(convert(b.state.p), a.state.p.value(), 1.0e-12));
                            expect(energy_relative(convert(b.state.e), a.state.e.value(), 1.0e-12));
                            expect(energy_relative(convert(b.state.mode_e[0]),
                                                   a.state.mode_e[0].value(),
                                                   1.0e-12))
                              << rho_value << t << t_ve_value;
                        };
                        compare(base,
                                thermo_state(gas_s,
                                             mix_s,
                                             rho_s,
                                             canary::convert<scaled>(Temp {t}, reference_scales),
                                             mode_temperatures<scaled>(t_ve_s)));
                        compare(thermo_state(gas,
                                             mix,
                                             rho,
                                             rho_e,
                                             mode_energy_densities<si>(rho_e_ve_si)),
                                thermo_state(gas_s,
                                             mix_s,
                                             rho_s,
                                             canary::convert<scaled>(rho_e, reference_scales),
                                             mode_energy_densities<scaled>(rho_e_ve_s)));
                    }
                }
            }
        }
    };

    "two-temperature: status reports non-convergence, invalid states and composition violations"_test =
      [&] {
          const auto mix = mixture(gas, fractions(air));
          const std::array t_ve {4000.0};
          const auto base =
            thermo_state(gas, mix, Rho {1.0}, Temp {5000.0}, mode_temperatures<si>(t_ve));
          const E e = base.state.e;
          const std::array rho_e_ve {base.state.mode_e[0].value()};
          const auto energies = mode_energy_densities<si>(rho_e_ve);
          Air one_iteration = gas;
          one_iteration.thermally_perfect.max_iterations = 1;
          expect(thermo_state(one_iteration, mix, Rho {1.0}, e, energies).status ==
                 Status::NotConverged);
          expect(thermo_state(gas, mix, Rho {0.0}, e, energies).status == Status::Invalid);
          expect(thermo_state(gas, mix, Rho {-1.0}, e, energies).status == Status::Invalid);
          expect(thermo_state(gas, mix, Rho {-1.0}, e, mode_temperatures<si>(t_ve)).status ==
                 Status::Invalid);
          const double nan = std::numeric_limits<double>::quiet_NaN();
          expect(thermo_state(gas, mix, Rho {1.0}, E {nan}, energies).status == Status::Invalid);
          expect(thermo_state(gas, mix, Rho {1.0}, E {nan}, mode_temperatures<si>(t_ve)).status ==
                 Status::Invalid);
          const std::array nan_energy {nan};
          expect(
            thermo_state(gas, mix, Rho {1.0}, e, mode_energy_densities<si>(nan_energy)).status ==
            Status::Invalid);
          expect(
            thermo_state(gas, mix, Rho {1.0}, Temp {5000.0}, mode_energy_densities<si>(nan_energy))
              .status == Status::Invalid);
          expect(thermo_state(gas, mix, Rho {1.0}, Temp {5000.0}, mode_temperatures<si>(nan_energy))
                   .status == Status::Invalid);
          const std::array negative_t_ve {-5.0};
          expect(
            thermo_state(gas, mix, Rho {1.0}, Temp {5000.0}, mode_temperatures<si>(negative_t_ve))
              .status == Status::Invalid);
          // Mode energies below e_ve(temperature_min) and above
          // e_ve(temperature_max).
          for (const double out_of_range : {-1.0e6, 1.0e12}) {
              const std::array bad {out_of_range};
              expect(thermo_state(gas, mix, Rho {1.0}, e, mode_energy_densities<si>(bad)).status ==
                     Status::Invalid)
                << out_of_range;
              expect(
                thermo_state(gas, mix, Rho {1.0}, Temp {5000.0}, mode_energy_densities<si>(bad))
                  .status == Status::Invalid)
                << out_of_range;
          }
          const std::array<double, 5> violated {0.767, 0.233 + 1.0e-10, 0.0, 0.0, 0.0};
          expect(thermo_state(gas, fractions(violated), Rho {1.0}, e, energies).status ==
                 Status::CompositionViolation);
          expect(thermo_state(gas, fractions(violated), Rho {1.0}, e, mode_temperatures<si>(t_ve))
                   .status == Status::CompositionViolation);
          expect(thermo_state(gas, fractions(violated), Rho {1.0}, e).status ==
                 Status::CompositionViolation);
          expect(thermo_state(gas, fractions(violated), Rho {-1.0}, e, energies).status ==
                 Status::Invalid);
      };

    "two-temperature: partial densities give the same closure as the equivalent fractions"_test =
      [&] {
          const Rho rho {2.0};
          std::array<double, 5> rho_i {};
          for (std::size_t i = 0; i < 5; ++i) { rho_i[i] = rho.value() * dissociated[i]; }
          const auto view = canary::as_kind<canary::kind::density, si>(
            Kokkos::mdspan<const double, Extents>(rho_i.data(), 5));
          const std::array t_ve {2500.0};
          const auto a = thermo_state(gas, view, rho, Temp {3000.0}, mode_temperatures<si>(t_ve));
          const auto b = thermo_state(gas,
                                      fractions(dissociated),
                                      rho,
                                      Temp {3000.0},
                                      mode_temperatures<si>(t_ve));
          expect(a.status == Status::Valid);
          expect(relative(a.state.e.value(), b.state.e.value(), 1.0e-12));
          expect(relative(a.state.mode_e[0].value(), b.state.mode_e[0].value(), 1.0e-12));
      };

    "two-temperature: float states close with double parameters"_test = [&] {
        const std::array<float, 5> y {0.767F, 0.233F, 0.0F, 0.0F, 0.0F};
        const auto mix = mixture(gas, fractions(y));
        const std::array t_ve {1500.0F};
        const auto c = thermo_state(gas,
                                    mix,
                                    canary::Density<float, si> {1.0F},
                                    canary::Temperature<float, si> {2000.0F},
                                    mode_temperatures<si>(t_ve));
        const std::array rho_e_ve {c.state.mode_e[0].value()};
        const auto back = thermo_state(gas,
                                       mix,
                                       canary::Density<float, si> {1.0F},
                                       c.state.e,
                                       mode_energy_densities<si>(rho_e_ve));
        expect(back.status == Status::Valid);
        expect(std::abs(back.state.temp.value() - 2000.0F) <= 1.0e-5F * 2000.0F);
        expect(std::abs(back.state.mode_temp[0].value() - 1500.0F) <= 1.0e-5F * 1500.0F);
    };

    "two-temperature: a T_ve solve that converges from one side ends Valid"_test = [&] {
        // Pure N2 with e_ve = e_ve(T_ve = 2210.65 K), solved from the guess
        // 7660.64 K: every Newton iterate lies above the root, and the last
        // Newton step rounds to the iterate itself.
        const std::array<double, 5> y {1.0, 0.0, 0.0, 0.0, 0.0};
        const auto mix = mixture(gas, fractions(y));
        const auto solved = canary::detail::solve_mode_temperature(gas,
                                                                   mix,
                                                                   E {288859.87009850127},
                                                                   7660.6364266913906);
        expect(solved.status == Status::Valid);
        expect(relative(solved.root, 2210.6513935268, 1.0e-12)) << solved.root;
    };
}
