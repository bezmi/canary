#include <cstdlib>
#include <string_view>

namespace
{
bool expect_extent_failure = false;
[[noreturn]] void fail_mode_extent() noexcept { std::_Exit(expect_extent_failure ? 0 : 2); }
}  // namespace

#define CANARY_ASSERT(condition) ((condition) ? void() : fail_mode_extent())
#include <canary/data/air11.hpp>
#include <canary/data/air11_kim_jo.hpp>
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
    const auto gas = canary::make_air11_three_temperature<double>(scales).model;
    const auto mechanism =
      canary::make_air11_kim_jo(scales, canary::Pressure<double, canary::si> {101325}, gas);
    const std::array<double, 11> y {0.767, 0.233, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    using Extents = Kokkos::dextents<std::size_t, 1>;
    const auto mix = canary::mixture(gas,
                                     canary::as_kind<canary::kind::mass_fraction>(
                                       Kokkos::mdspan<const double, Extents>(y.data(), y.size())));
    const std::array<double, 3> modes {5000, 8000, 9000};
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
