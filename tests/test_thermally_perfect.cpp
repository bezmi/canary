#include "dual.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/data/air5.hpp>
#include <canary/models/thermally_perfect.hpp>
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
using Guess = canary::TemperatureGuess<double, si>;
using Extents = Kokkos::dextents<std::size_t, 1>;
using Air = canary::Air5<double, si>;
using AirScaled = canary::Air5<double, scaled>;

constexpr double eps = std::numeric_limits<double>::epsilon();
// Relative tolerance of the temperature solve (the model default).
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

bool relative(double a, double b, double bound) { return std::abs(a - b) <= bound * std::abs(b); }

constexpr std::array<double, 5> air {0.767, 0.233, 0.0, 0.0, 0.0};
constexpr std::array<double, 5> dissociated {0.2, 0.05, 0.3, 0.4, 0.05};
constexpr std::array<double, 5> nitrogen {1.0, 0.0, 0.0, 0.0, 0.0};

}  // namespace

static_assert(canary::GasModel<Air, double>);
static_assert(canary::MixtureGasModel<Air, double>);
static_assert(canary::MixtureGasModel<Air, float>);
static_assert(canary::MixtureGasModel<Air, D>);
static_assert(canary::MixtureGasModel<AirScaled, double>);
static_assert(canary::MixtureGasModel<AirScaled, float>);
static_assert(canary::MixtureGasModel<AirScaled, D>);
static_assert(canary::AcceptsInput<Air, canary::TemperatureGuess<double, si>>);
static_assert(canary::AcceptsInput<AirScaled, canary::TemperatureGuess<double, scaled>>);
static_assert(!canary::AcceptsInput<AirScaled, canary::TemperatureGuess<double, si>>);
static_assert(std::is_trivially_copyable_v<Air>);
static_assert(std::is_trivially_copyable_v<AirScaled>);
static_assert(std::is_trivially_copyable_v<Air::Mixture<double>>);
static_assert(std::is_trivially_copyable_v<Air::Mixture<D>>);
static_assert(std::is_trivially_copyable_v<canary::SpeciesProperties<double, si, 5>>);

void test_thermally_perfect()
{
    using namespace boost::ut;
    using canary::mixture;
    using canary::thermo_state;
    const auto built = canary::make_air5<double>(canary::UnitScales<double, si> {});
    const Air& gas = built.model;

    "thermally perfect construction validates stored species and reference parameters"_test = [] {
        const auto check = []<typename V> {
            const std::array<canary::Nasa9Species<1>, 1> good {{{
              .name = "test",
              .molar_mass = 0.028,
              .charge = 0,
              .metadata = {.role = canary::SpeciesRole::heavy},
              .fit = {.segments = 1,
                      .breaks = {200, 20000},
                      .blends = {},
                      .a = {{{0, 0, 3, 0, 0, 0, 0, 0, 0}}}},
            }}};
            const auto rejects =
              [](const auto& inputs, const auto& scales, std::string_view reason) {
                  bool rejected = false;
                  try {
                      static_cast<void>(canary::make_thermally_perfect_mixture(inputs, scales));
                  } catch (const canary::ParameterError& error) {
                      rejected = true;
                      const std::string_view message {error.what()};
                      expect(message.starts_with("species test: ")) << message;
                      expect(message.find(reason) != std::string_view::npos) << message << reason;
                  }
                  expect(rejected) << reason;
              };
            const canary::UnitScales<V, si> si_scales {};
            auto inputs = good;
            inputs[0].molar_mass = std::numeric_limits<double>::denorm_min();
            rejects(inputs, si_scales, "stored gas constant");
            inputs[0].molar_mass = std::same_as<V, float> ? 1.0e-30 : 1.0e-301;
            rejects(inputs, si_scales, "stored particle mass");
            inputs = good;
            inputs[0].charge = std::numeric_limits<double>::denorm_min();
            rejects(inputs, si_scales, "nonzero charge rounds to zero");
            if constexpr (std::same_as<V, float>) {
                inputs[0].charge = 1.0e60;
                rejects(inputs, si_scales, "stored charge");
            }
            const canary::UnitScales<V, scaled> large_mass {
              {.mass = V(1), .density = V(1), .temperature = V(1), .length = V(1)}};
            inputs = good;
            if constexpr (std::same_as<V, float>) {
                inputs[0].molar_mass = 1.0e64;
                rejects(inputs, large_mass, "stored particle mass");
            }
            inputs[0].molar_mass = std::numeric_limits<double>::max();
            const canary::UnitScales<V, scaled> small_mass {{.mass = std::numeric_limits<V>::min(),
                                                             .density = V(1),
                                                             .temperature = V(1),
                                                             .length = V(1)}};
            rejects(inputs, small_mass, "stored gas constant");

            // A one-kilogram particle makes the scaled gas constant equal the mass scale.
            inputs = good;
            inputs[0].molar_mass = canary::constants::avogadro<double>;
            inputs[0].fit.a[0][2] = std::nextafter(V(1), V(2));
            const canary::UnitScales<V, scaled> small_capacity {
              {.mass = V(1.25) / std::numeric_limits<V>::max(),
               .density = V(1),
               .temperature = V(1),
               .length = V(1)}};
            rejects(inputs, small_capacity, "reference isochoric specific heat");

            inputs = good;
            const V particle_mass = V(inputs[0].molar_mass / canary::constants::avogadro<double>);
            const canary::UnitScales<V, scaled> large_capacity {
              {.mass = (std::numeric_limits<V>::max() / V(2)) * particle_mass,
               .density = V(1),
               .temperature = V(1),
               .length = V(1)}};
            rejects(inputs, large_capacity, "reference isobaric specific heat");
            const canary::UnitScales<V, scaled> large_energy {
              {.mass = (std::numeric_limits<V>::max() / V(100)) * particle_mass,
               .density = V(1),
               .temperature = V(1),
               .length = V(1)}};
            rejects(inputs, large_energy, "reference specific energy");
            const canary::UnitScales<V, scaled> small_temperature {
              {.mass = particle_mass,
               .density = V(1),
               .temperature = std::numeric_limits<V>::min(),
               .length = V(1)}};
            bool rejected = false;
            try {
                static_cast<void>(canary::make_thermally_perfect_mixture(good, small_temperature));
            } catch (const canary::ParameterError& error) {
                rejected = true;
                expect(std::string_view(error.what()).starts_with("reference temperature"))
                  << error.what();
            }
            expect(rejected);
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "thermally perfect construction scales before narrowing species parameters"_test = [] {
        const auto check = []<typename V>(double molar_mass, V mass_scale, double charge) {
            const std::array<canary::Nasa9Species<1>, 1> inputs {{{
              .name = "test",
              .molar_mass = molar_mass,
              .charge = charge,
              .metadata = {.role = canary::SpeciesRole::heavy},
              .fit = {.segments = 1,
                      .breaks = {200, 20000},
                      .blends = {},
                      .a = {{{0, 0, 3, 0, 0, 0, 0, 0, 0}}}},
            }}};
            const canary::UnitScales<V, scaled> scales {
              {.mass = mass_scale, .density = V(1), .temperature = V(1), .length = V(1)}};
            const auto built = canary::make_thermally_perfect_mixture(inputs, scales);
            const long double expected_mass =
              static_cast<long double>(molar_mass) /
              (canary::constants::avogadro<long double> * static_cast<long double>(mass_scale));
            const auto close = [](V actual, long double expected) {
                return std::abs(static_cast<long double>(actual) - expected) <=
                       4 * std::numeric_limits<V>::epsilon() * std::abs(expected);
            };
            expect(close(built.species.mass[0].value(), expected_mass));
            expect(close(built.model.gas_constant[0].value(), 1 / expected_mass));
            expect(built.species.charge[0].value() == V(charge));
            expect(std::isfinite(built.model.reference_energy[0].value()));
            expect(std::isfinite(built.model.reference_cv[0].value()));
            expect(built.model.reference_cv[0].value() > V(0));
        };
        check(1.0e-40, 1.0e-30F, 1.0e-30);
        check(1.0e-310, 1.0e-300, std::numeric_limits<double>::denorm_min());
    };

    "thermally perfect construction validates explicit species metadata"_test = [] {
        const auto check = [](const auto& inputs, std::string_view expected) {
            try {
                canary::make_thermally_perfect_mixture(inputs, canary::UnitScales<double, si> {});
                expect(false) << "invalid metadata accepted";
            } catch (const canary::ParameterError& error) {
                expect(std::string_view(error.what()).starts_with(expected)) << error.what();
            }
        };
        auto inputs = canary::air5_species;
        inputs[0].metadata.role = canary::SpeciesRole::unspecified;
        check(inputs, "species N2: species role must be specified");
        inputs[0].metadata.role = static_cast<canary::SpeciesRole>(99);
        check(inputs, "species N2: unknown species role");
        inputs[0].metadata.role = canary::SpeciesRole::electron;
        check(inputs, "species N2: charge must be -1");
        inputs[0].charge = -1;
        inputs[0].metadata.elements = std::span<const canary::ElementCount> {};
        inputs[1].metadata.role = canary::SpeciesRole::electron;
        inputs[1].charge = -1;
        check(inputs, "species O2: more than one electron-role species");
        inputs = canary::air5_species;
        const std::array invalid_elements {
          std::array<canary::ElementCount, 2> {{{"", 1}, {"N", 1}}},
          std::array<canary::ElementCount, 2> {{{"N", 0}, {"O", 1}}},
          std::array<canary::ElementCount, 2> {{{"N", 1}, {"N", 1}}},
        };
        const std::array diagnostics {
          "species N2: element identifier must not be empty",
          "species N2: element N: count must be positive",
          "species N2: element N: given twice",
        };
        for (std::size_t i = 0; i < invalid_elements.size(); ++i) {
            inputs[0].metadata.elements =
              std::span<const canary::ElementCount>(invalid_elements[i]);
            check(inputs, diagnostics[i]);
        }
        const std::array<canary::ElementCount, 1> elements {{{"N", 2}}};
        inputs[0].metadata.elements = std::span<const canary::ElementCount>(elements);
        expect(nothrow([&] {
            canary::make_thermally_perfect_mixture(inputs, canary::UnitScales<double, si> {});
        }));
        inputs[0].metadata.role = canary::SpeciesRole::electron;
        inputs[0].charge = -1;
        check(inputs, "species N2: free-electron elemental composition must be empty");
        inputs[0].metadata.elements = std::span<const canary::ElementCount> {};
        expect(nothrow([&] {
            canary::make_thermally_perfect_mixture(inputs, canary::UnitScales<double, si> {});
        }));
        inputs[0].metadata.elements.reset();
        expect(nothrow([&] {
            canary::make_thermally_perfect_mixture(inputs, canary::UnitScales<double, si> {});
        }));
        for (const char* name : std::array<const char*, 2> {nullptr, ""}) {
            inputs[0].name = name;
            check(inputs, "species 0: name must not be null or empty");
        }
    };

    "thermally perfect closures reject pressure overflow and retain diagnostic fields"_test = [&] {
        const auto check = [&]<typename V> {
            const std::array<V, 5> y {V(1), V(0), V(0), V(0), V(0)};
            const auto mix = mixture(gas, fractions(y));
            const canary::Density<V, si> rho {std::numeric_limits<V>::max() / V(16)};
            const canary::Temperature<V, si> temp {V(300)};
            const auto result = thermo_state(gas, mix, rho, temp);
            expect(result.status == Status::Invalid);
            expect(std::isinf(result.state.p.value()));
            expect(result.state.rho == rho && result.state.temp == temp);
            expect(std::isfinite(result.state.e.value()));
            const auto inverse =
              thermo_state(gas, mix, rho, result.state.e, canary::TemperatureGuess<V, si> {temp});
            expect(inverse.status == Status::Invalid);
            expect(std::isinf(inverse.state.p.value()));
            const canary::Pressure<V, si> supplied {std::numeric_limits<V>::max() / V(4)};
            const canary::Temperature<V, si> cold {V(1) / mix.gas_constant.value()};
            const auto retained = thermo_state(gas, mix, supplied, cold);
            expect(retained.status == Status::Valid);
            expect(retained.state.p == supplied);
            expect(std::isfinite(retained.state.rho.value()));
        };
        check.operator()<float>();
        check.operator()<double>();
    };

    "the builder derives gas constants, species data and settings"_test = [&] {
        const double r_u = canary::constants::universal_gas_constant<double>;
        expect(gas.n == 5_i);
        expect(built.species.n == 5_i);
        for (std::size_t i = 0; i < 5; ++i) {
            const auto& sp = canary::air5_species[i];
            expect(relative(gas.gas_constant[i].value(), r_u / sp.molar_mass, 2.0 * eps));
            expect(relative(built.species.mass[i].value(),
                            sp.molar_mass / canary::constants::avogadro<double>,
                            2.0 * eps));
            expect(built.species.charge[i].value() == 0.0);
        }
        expect(gas.temperature_min.value() == 10.0);
        expect(gas.temperature_max.value() == 100000.0);
        expect(gas.bracket_width.value() == 1000.0);
        expect(gas.relative_tolerance == solve_tolerance);
        expect(gas.absolute_tolerance.value() == 0.0);
        expect(gas.max_iterations == 30_i);
    };

    "the temperature solver uses absolute tolerance in the model units"_test = [&] {
        const auto check = []<typename Gas>(Gas model) {
            using V = typename Gas::parameter_type;
            using U = typename Gas::units;
            using Temperature = canary::Temperature<V, U>;
            static_assert(std::same_as<decltype(model.absolute_tolerance), Temperature>);
            model.temperature_min = Temperature {V(2)};
            model.temperature_max = Temperature {V(4)};
            model.bracket_width = Temperature {V(2)};
            model.relative_tolerance = V(0);
            model.max_iterations = 1;
            const auto f = [](V x) {
                return canary::detail::Residual<V> {.value = x - V(3), .derivative = V(0)};
            };
            const auto exact_only = canary::detail::solve_temperature<V>(model, f, f, V(2));
            expect(exact_only.status == Status::NotConverged);
            model.absolute_tolerance = Temperature {V(2)};
            const auto approximate = canary::detail::solve_temperature<V>(model, f, f, V(2));
            expect(approximate.status == Status::Valid);
            expect(approximate.root == V(2));
            for (const V invalid :
                 {V(-1), std::numeric_limits<V>::quiet_NaN(), std::numeric_limits<V>::infinity()}) {
                model.absolute_tolerance = Temperature {invalid};
                expect(canary::detail::solve_temperature<V>(model, f, f, V(2)).status ==
                       Status::Invalid);
            }
        };
        check(gas);
        check(canary::make_air5<double>(reference_scales).model);
    };

    "the builder names the species of an invalid parameter"_test = [] {
        auto species = canary::air5_species;
        species[3].molar_mass = -1.0;
        try {
            canary::make_thermally_perfect_mixture(species, canary::UnitScales<double, si> {});
            expect(false) << "no exception";
        } catch (const canary::ParameterError& error) {
            expect(std::string(error.what()).starts_with("species O: molar mass")) << error.what();
        }
        species = canary::air5_species;
        species[1].fit.blends[0] = -1.0;
        try {
            canary::make_thermally_perfect_mixture(species, canary::UnitScales<double, si> {});
            expect(false) << "no exception";
        } catch (const canary::ParameterError& error) {
            expect(std::string(error.what()).starts_with("species O2: blending width"))
              << error.what();
        }
    };

    "(rho, e) inverts (rho, T) and all closures agree"_test = [&] {
        for (const auto& y : {air, dissociated}) {
            const auto mix = mixture(gas, fractions(y));
            expect(mix.status == Status::Valid);
            for (const double rho_value : {1.0e-4, 1.0, 100.0}) {
                for (const double t : {250.0, 1000.0, 5000.0, 15000.0, 30000.0}) {
                    const Rho rho {rho_value};
                    const auto from_t = thermo_state(gas, mix, rho, Temp {t});
                    const auto from_e = thermo_state(gas, mix, rho, from_t.state.e);
                    expect(from_t.status == Status::Valid);
                    expect(from_e.status == Status::Valid) << rho_value << t;
                    expect(relative(from_e.state.temp.value(), t, solve_tolerance))
                      << rho_value << t << from_e.state.temp.value();
                    expect(relative(from_e.state.p.value(), from_t.state.p.value(), 1.0e-9));
                    expect(from_e.state.e == from_t.state.e);

                    const auto from_p = thermo_state(gas, mix, rho, from_t.state.p);
                    expect(relative(from_p.state.temp.value(), t, 4.0 * eps));
                    expect(relative(from_p.state.e.value(), from_t.state.e.value(), 1.0e-12));
                    expect(from_p.status == Status::Valid);

                    const auto from_pt = thermo_state(gas, mix, from_t.state.p, Temp {t});
                    expect(relative(from_pt.state.rho.value(), rho_value, 4.0 * eps));
                    expect(from_pt.state.e == from_t.state.e);
                    expect(from_pt.status == Status::Valid);

                    // (rho, rho e) is (rho, e) with e = rho_e / rho.
                    const RhoE rho_e {canary::energy_density(rho, from_t.state.e)};
                    const auto a = thermo_state(gas, mix, rho, rho_e);
                    const auto b = thermo_state(gas, mix, rho, canary::specific_energy(rho_e, rho));
                    expect(a.state.temp == b.state.temp && a.state.e == b.state.e &&
                           a.state.p == b.state.p && a.status == b.status);
                }
            }
        }
    };

    "partial densities give the same closure as the equivalent fractions"_test = [&] {
        const Rho rho {2.0};
        std::array<double, 5> rho_i {};
        for (std::size_t i = 0; i < 5; ++i) { rho_i[i] = rho.value() * dissociated[i]; }
        const auto view = canary::as_kind<canary::kind::density, si>(
          Kokkos::mdspan<const double, Extents>(rho_i.data(), 5));
        const E e = thermo_state(gas, fractions(dissociated), rho, Temp {3000.0}).state.e;
        const auto a = thermo_state(gas, view, rho, e);
        const auto b = thermo_state(gas, fractions(dissociated), rho, e);
        expect(relative(a.state.temp.value(), b.state.temp.value(), 1.0e-12));
        expect(a.status == Status::Valid);
    };

    "nitrogen at 300 K has gamma 1.4 within 0.1 percent"_test = [&] {
        const auto mix = mixture(gas, fractions(nitrogen));
        const auto c = thermo_state(gas, mix, Rho {1.0}, Temp {300.0});
        const auto props = canary::properties(gas, mix, c.state);
        expect(relative(props.capacity_ratio.value(), 1.4, 1.0e-3)) << props.capacity_ratio.value();
        expect(relative(props.cp.value() - props.cv.value(), mix.gas_constant.value(), 1.0e-12));
        expect(relative(
          props.acoustics.sound_speed.value(),
          std::sqrt(props.capacity_ratio.value() * c.state.p.value() / c.state.rho.value()),
          4.0 * eps));
    };

    "dT/de at fixed density and composition is 1 / cv"_test = [&] {
        const auto mix_d = mixture(gas, fractions(std::array<D, 5> {0.767, 0.233, 0.0, 0.0, 0.0}));
        const auto mix = mixture(gas, fractions(air));
        for (const double t : {300.0, 2500.0, 8000.0}) {
            const auto ref = thermo_state(gas, mix, Rho {1.0}, Temp {t});
            const auto c =
              thermo_state(gas,
                           mix_d,
                           canary::Density<D, si> {D(1.0)},
                           canary::SpecificEnergy<D, si> {D(ref.state.e.value(), 1.0)});
            const auto plain = thermo_state(gas, mix, Rho {1.0}, ref.state.e);
            expect(c.status == Status::Valid && plain.status == Status::Valid);
            expect(c.state.temp.value().real() == plain.state.temp.value());
            const auto props = canary::properties(gas, mix, ref.state);
            expect(relative(c.state.temp.value().dual(), 1.0 / props.cv.value(), 1.0e-9))
              << t << c.state.temp.value().dual() << 1.0 / props.cv.value();
        }
    };

    "dT/d(epsilon) at fixed e is -(e_1 - e_2) / cv"_test = [&] {
        // y_1 = a + epsilon, y_2 = b - epsilon.
        const std::array<D, 5> y {D(0.4, 1.0), D(0.1, -1.0), D(0.2), D(0.25), D(0.05)};
        const std::array<double, 5> y_value {0.4, 0.1, 0.2, 0.25, 0.05};
        const auto mix = mixture(gas, fractions(y_value));
        for (const double t : {600.0, 3000.0, 9000.0}) {
            const auto ref = thermo_state(gas, mix, Rho {1.0}, Temp {t});
            const auto c = thermo_state(gas,
                                        fractions(y),
                                        canary::Density<D, si> {D(1.0)},
                                        canary::SpecificEnergy<D, si> {D(ref.state.e.value())});
            const auto sp = canary::species_properties(gas, ref.state.temp);
            const double cv = canary::properties(gas, mix, ref.state).cv.value();
            const double expected = -(sp.energy[0].value() - sp.energy[1].value()) / cv;
            expect(relative(c.state.temp.value().dual(), expected, 1.0e-9))
              << t << c.state.temp.value().dual() << expected;
        }
    };

    "a scaled model gives the SI state after conversion"_test = [&] {
        const auto scaled_built = canary::make_air5<double>(reference_scales);
        const AirScaled& gas_s = scaled_built.model;
        for (const auto& y : {air, dissociated}) {
            const auto mix = mixture(gas, fractions(y));
            const auto mix_s = mixture(gas_s, fractions(y));
            for (const double rho_value : {1.0e-4, 1.0, 100.0}) {
                for (const double t : {1000.0, 5000.0, 15000.0}) {
                    const Rho rho {rho_value};
                    const auto e = thermo_state(gas, mix, rho, Temp {t}).state.e;
                    const auto a = thermo_state(gas, mix, rho, e);
                    const auto b = thermo_state(gas_s,
                                                mix_s,
                                                canary::convert<scaled>(rho, reference_scales),
                                                canary::convert<scaled>(e, reference_scales));
                    expect(b.status == Status::Valid);
                    const auto temp = canary::convert<si>(b.state.temp, reference_scales);
                    const auto p = canary::convert<si>(b.state.p, reference_scales);
                    const auto e_si = canary::convert<si>(b.state.e, reference_scales);
                    expect(relative(temp.value(), a.state.temp.value(), 1.0e-12))
                      << rho_value << t << temp.value() << a.state.temp.value();
                    expect(relative(p.value(), a.state.p.value(), 1.0e-12));
                    expect(relative(e_si.value(), a.state.e.value(), 1.0e-12));
                }
            }
        }
    };

    "status reports non-convergence, invalid states and composition violations"_test = [&] {
        const auto mix = mixture(gas, fractions(air));
        const E e = thermo_state(gas, mix, Rho {1.0}, Temp {5000.0}).state.e;
        Air one_iteration = gas;
        one_iteration.max_iterations = 1;
        expect(thermo_state(one_iteration, mix, Rho {1.0}, e).status == Status::NotConverged);
        expect(thermo_state(gas, mix, Rho {0.0}, e).status == Status::Invalid);
        expect(thermo_state(gas, mix, Rho {-1.0}, e).status == Status::Invalid);
        expect(
          thermo_state(gas, mix, Rho {1.0}, E {std::numeric_limits<double>::quiet_NaN()}).status ==
          Status::Invalid);
        expect(thermo_state(gas, mix, Rho {1.0}, Temp {-5.0}).status == Status::Invalid);
        const std::array<double, 5> violated {0.767, 0.233 + 1.0e-10, 0.0, 0.0, 0.0};
        expect(thermo_state(gas, fractions(violated), Rho {1.0}, e).status ==
               Status::CompositionViolation);
        expect(thermo_state(gas, fractions(violated), Rho {1.0}, Temp {500.0}).status ==
               Status::CompositionViolation);
        expect(thermo_state(gas, fractions(violated), Rho {-1.0}, e).status == Status::Invalid);
    };

    "a temperature guess near the answer gives the same temperature"_test = [&] {
        for (const auto& y : {air, dissociated}) {
            const auto mix = mixture(gas, fractions(y));
            for (const double t : {250.0, 1000.0, 5000.0, 30000.0}) {
                const Rho rho {1.0};
                const E e = thermo_state(gas, mix, rho, Temp {t}).state.e;
                const auto without = thermo_state(gas, mix, rho, e);
                const auto with = thermo_state(gas, mix, rho, e, Guess {Temp {1.01 * t}});
                expect(with.status == Status::Valid);
                expect(relative(with.state.temp.value(),
                                without.state.temp.value(),
                                2.0 * solve_tolerance));
                // Through the front end, from a composition view and rho e.
                const auto front = thermo_state(gas,
                                                fractions(y),
                                                RhoE {canary::energy_density(rho, e)},
                                                rho,
                                                Guess {Temp {0.99 * t}});
                expect(relative(front.state.temp.value(), t, 2.0 * solve_tolerance));
                // A guess outside the temperature limits is clamped to them.
                expect(
                  relative(thermo_state(gas, mix, rho, e, Guess {Temp {1.0e9}}).state.temp.value(),
                           t,
                           2.0 * solve_tolerance));
            }
        }
    };

    "energy inversion rejects temperatures above the configured maximum"_test = [&] {
        Air limited = gas;
        limited.temperature_max = Temp {1000.0};
        for (const auto& y : {air, dissociated}) {
            const auto mix = mixture(limited, fractions(y));
            for (const double t : {1001.0, 1250.0}) {
                const E e = thermo_state(gas, mix, Rho {1.0}, Temp {t}).state.e;
                for (const double guess : {1000.0, 2000.0}) {
                    const auto c = thermo_state(limited, mix, Rho {1.0}, e, Guess {Temp {guess}});
                    expect(c.status == Status::Invalid);
                    expect(c.state.temp.value() >= limited.temperature_min.value());
                    expect(c.state.temp.value() <= limited.temperature_max.value());
                }
            }
        }
    };

    "species properties are consistent with each other and the closure"_test = [&] {
        for (const double t : {250.0, 1000.0, 5000.0, 30000.0}) {
            const auto sp = canary::species_properties(gas, Temp {t});
            expect(sp.n == 5_i);
            double e_mix = 0;
            double cv_mix = 0;
            double e_scale = 0;
            for (std::size_t i = 0; i < 5; ++i) {
                const double r = gas.gas_constant[i].value();
                const double h = sp.enthalpy[i].value();
                const double s = sp.entropy[i].value();
                const double g = sp.gibbs_energy[i].value();
                expect(std::abs(g - (h - (t * s))) <= 4.0 * eps * (std::abs(h) + std::abs(t * s)));
                expect(std::abs(sp.energy[i].value() - (h - (r * t))) <=
                       4.0 * eps * (std::abs(h) + (r * t)));
                e_mix += dissociated[i] * sp.energy[i].value();
                cv_mix += dissociated[i] * sp.cv[i].value();
                e_scale += dissociated[i] * (std::abs(h) + (r * t));
            }
            const auto mix = mixture(gas, fractions(dissociated));
            const auto c = thermo_state(gas, mix, Rho {1.0}, Temp {t});
            expect(std::abs(c.state.e.value() - e_mix) <= 16.0 * eps * e_scale)
              << t << c.state.e.value() << e_mix;
            expect(relative(canary::properties(gas, mix, c.state).cv.value(), cv_mix, 1.0e-14));
        }
    };

    "float states close with double parameters"_test = [&] {
        const std::array<float, 5> y {0.767F, 0.233F, 0.0F, 0.0F, 0.0F};
        const auto mix = mixture(gas, fractions(y));
        const auto c = thermo_state(gas,
                                    mix,
                                    canary::Density<float, si> {1.0F},
                                    canary::Temperature<float, si> {2000.0F});
        const auto back = thermo_state(gas, mix, canary::Density<float, si> {1.0F}, c.state.e);
        expect(back.status == Status::Valid);
        expect(std::abs(back.state.temp.value() - 2000.0F) <= 1.0e-5F * 2000.0F);
    };
}
