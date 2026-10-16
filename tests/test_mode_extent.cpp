#include <cstdlib>
#include <string_view>
#include <type_traits>

namespace
{
bool expect_extent_failure = false;
[[noreturn]] void fail_mode_extent() noexcept { std::_Exit(expect_extent_failure ? 0 : 2); }
}  // namespace

#define CANARY_ASSERT(condition) ((condition) ? void() : fail_mode_extent())
#include <canary/data/air11.hpp>
#include <canary/data/air11_kim_jo.hpp>
#include <canary/data/air11_kim_jo_exchange.hpp>
#include <canary/data/air5_park.hpp>
#include <canary/reactions.hpp>
#include <canary/thermo.hpp>

// First evaluate a valid view, then require the same operation to assert on a
// wrong extent. Assertion interception exits before any out-of-range access.
int main(int argc, char** argv)
{
    if (argc != 3) { return 3; }
    const std::string_view operation(argv[1]);
    const bool typed = std::string_view(argv[2]) == "typed";
    const canary::UnitScales<double, canary::si> scales {};
    const auto built = canary::make_air11_three_temperature<double>(scales);
    const auto& gas = built.model;
    const auto mechanism =
      canary::make_air11_kim_jo(scales, canary::Pressure<double, canary::si> {101325}, gas);
    const std::array<double, 11> y {0.767, 0.233, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    using Extents = Kokkos::dextents<std::size_t, 1>;
    const auto mix = canary::mixture(gas,
                                     canary::as_kind<canary::kind::mass_fraction>(
                                       Kokkos::mdspan<const double, Extents>(y.data(), y.size())));
    const std::array<double, 3> modes {5000, 8000, 9000};
    if (operation == "bound_rates" || operation == "bound_production") {
        const auto unbound = canary::make_air11_kim_jo(scales);
        const auto raw = canary::as_kind<canary::kind::mode_temperature, canary::si>(
          Kokkos::mdspan<const double, Extents>(nullptr, gas.mode_count));
        const auto check = [&](const auto& view) {
            expect_extent_failure = true;
            if (operation == "bound_rates") {
                (void)canary::rates_of_progress(unbound,
                                                gas,
                                                mix,
                                                canary::Density<double, canary::si> {0.001},
                                                canary::Temperature<double, canary::si> {10000},
                                                view);
            } else {
                (void)canary::production_rates(unbound,
                                               gas,
                                               mix,
                                               canary::Density<double, canary::si> {0.001},
                                               canary::Temperature<double, canary::si> {10000},
                                               view);
            }
        };
        if (typed) {
            check(canary::as_indexed<canary::ModeIndex>(raw));
        } else {
            check(raw);
        }
        return 1;
    }
    if (operation == "bound_scalar") {
        const auto two = canary::make_air5_two_temperature<double>(scales);
        const auto unbound = canary::make_air5_park(scales);
        const std::array<double, 5> fractions {0.767, 0.233, 0, 0, 0};
        const auto composition = canary::mixture(
          two.model,
          canary::as_kind<canary::kind::mass_fraction>(
            Kokkos::mdspan<const double, Extents>(fractions.data(), fractions.size())));
        expect_extent_failure = true;
        (void)canary::rates_of_progress(unbound,
                                        two.model,
                                        composition,
                                        canary::Density<double, canary::si> {0.001},
                                        canary::Temperature<double, canary::si> {10000},
                                        canary::Temperature<double, canary::si> {5000});
        return 1;
    }
    const auto evaluate = [&](const auto& view) {
        const canary::Density<double, canary::si> rho {0.001};
        const canary::Temperature<double, canary::si> temp {10000};
        if (operation == "closure") {
            (void)canary::thermo_state(gas, mix, rho, temp, view);
        } else if (operation == "rates") {
            (void)canary::rates_of_progress(mechanism, gas, mix, rho, temp, view);
        } else if (operation == "production") {
            (void)canary::production_rates(mechanism, gas, mix, rho, temp, view);
        } else {
            std::_Exit(3);
        }
    };
    if (operation.starts_with("state_")) {
        const auto exchange =
          canary::make_air11_kim_jo_exchange(built.species, mechanism, gas, scales);
        const bool bad_energy = std::string_view(argv[2]) == "energy";
        const auto run_state = [&](std::size_t count) {
            const auto mode_temp = canary::as_indexed<canary::ModeIndex>(
              canary::as_kind<canary::kind::mode_temperature, canary::si>(
                Kokkos::mdspan<const double, Extents>(modes.data(),
                                                      bad_energy ? gas.mode_count : count)));
            const auto mode_e = canary::as_indexed<canary::ModeIndex>(
              canary::as_kind<canary::kind::specific_energy, canary::si>(
                Kokkos::mdspan<const double, Extents>(modes.data(),
                                                      bad_energy ? count : gas.mode_count)));
            const canary::
              MultiTemperatureStateView<double, canary::si, decltype(mode_e), decltype(mode_temp)>
                state {canary::Density<double, canary::si> {0.001},
                       canary::SpecificEnergy<double, canary::si> {1e6},
                       canary::Pressure<double, canary::si> {1000},
                       canary::Temperature<double, canary::si> {10000},
                       mode_e,
                       mode_temp};
            if (operation == "state_vt") {
                (void)canary::vibration_translation_source(exchange, gas, mix, state);
            } else if (operation == "state_cv") {
                const canary::RatesOfProgress<
                  double,
                  canary::si,
                  std::remove_cvref_t<decltype(mechanism)>::reaction_count>
                  rates {};
                (void)canary::chemistry_vibration_source(exchange, mechanism, rates, state);
            } else if (operation == "state_et") {
                (void)canary::electron_translation_source(exchange, mix, state);
            } else if (operation == "state_ve") {
                (void)canary::vibration_electron_source(exchange, gas, mix, state);
            } else if (operation == "state_sources") {
                (void)canary::multi_temperature_sources(mechanism, exchange, gas, mix, state);
            } else {
                std::_Exit(3);
            }
        };
        run_state(gas.mode_count);
        expect_extent_failure = true;
        run_state(gas.mode_count + 1);
        return 1;
    }
    const auto run = [&](std::size_t count) {
        const auto raw = canary::as_kind<canary::kind::mode_temperature, canary::si>(
          Kokkos::mdspan<const double, Extents>(modes.data(), count));
        if (typed) {
            evaluate(canary::as_indexed<canary::ModeIndex>(raw));
        } else {
            evaluate(raw);
        }
    };
    run(gas.mode_count);
    expect_extent_failure = true;
    run(gas.mode_count + 1);
    return 1;
}
