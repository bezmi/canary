#include "dual.hpp"
#include "test_suites.hpp"

#include <boost/ut.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>
#include <canary/status.hpp>
#include <concepts>
#include <limits>
#include <ranges>
#include <span>
#include <type_traits>

namespace
{

using canary::scaled;
using canary::si;
using D = canary::test::Dual<double>;

// A state without a temperature member.
struct NoTemperature {
    using value_type = double;
    using units = si;
    canary::Density<double, si> rho;
    canary::SpecificEnergy<double, si> e;
    canary::Pressure<double, si> p;
};

// A caller-defined result can expose both sets of properties and non-owning
// mode capacities without inheriting from a library record.
template <typename Ratio = canary::HeatCapacityRatio<double>,
          typename ModeCapacity = canary::SpecificHeat<double, si>>
struct CustomProperties {
    canary::SpecificHeat<double, si> cv;
    canary::SpecificHeat<double, si> cp;
    Ratio capacity_ratio;
    canary::SpecificGasConstant<double, si> gas_constant;
    canary::SpecificHeat<double, si> cv_heavy;
    std::span<const ModeCapacity> mode_cv;
    canary::SpecificGasConstant<double, si> gas_constant_heavy;
    canary::SpecificGasConstant<double, si> gas_constant_electron;
    canary::AcousticProperties<double, si> acoustics;
};

struct IncorrectAcousticProperties {
    canary::SoundSpeed<double, si> sound_speed;
    canary::HeatCapacityRatio<double> acoustic_exponent;
};

}  // namespace

static_assert(sizeof(canary::ModeIndex) == sizeof(std::size_t));
static_assert(std::is_trivially_copyable_v<canary::ModeIndex>);
static_assert(std::is_standard_layout_v<canary::ModeIndex>);
static_assert(std::is_nothrow_default_constructible_v<canary::ModeIndex>);
static_assert(std::is_nothrow_constructible_v<canary::ModeIndex, std::size_t>);
static_assert(!std::convertible_to<std::size_t, canary::ModeIndex>);
static_assert(!std::convertible_to<canary::ModeIndex, std::size_t>);
static_assert(!std::constructible_from<canary::ModeIndex, canary::SpeciesIndex>);
static_assert(!std::constructible_from<canary::SpeciesIndex, canary::ModeIndex>);
static_assert(canary::ModeIndex {} == canary::ModeIndex {0});
static_assert(canary::ModeIndex {std::numeric_limits<std::size_t>::max()}.value() ==
              std::numeric_limits<std::size_t>::max());
static_assert(noexcept(canary::ModeIndex {1}.value()));
static_assert(noexcept(canary::ModeIndex {1} == canary::ModeIndex {2}));

static_assert(sizeof(canary::OptionalModeIndex) == sizeof(std::size_t));
static_assert(std::is_trivially_copyable_v<canary::OptionalModeIndex>);
static_assert(std::is_nothrow_constructible_v<canary::OptionalModeIndex, canary::ModeIndex>);
static_assert(!std::constructible_from<canary::OptionalModeIndex, std::size_t>);
static_assert(!std::constructible_from<canary::OptionalModeIndex, canary::SpeciesIndex>);
static_assert(!std::constructible_from<canary::OptionalModeIndex, canary::OptionalSpeciesIndex>);
static_assert(!std::convertible_to<canary::ModeIndex, canary::OptionalModeIndex>);
static_assert(!canary::OptionalModeIndex {}.has_value());
static_assert(*canary::OptionalModeIndex {canary::ModeIndex {0}} == canary::ModeIndex {0});
static_assert(noexcept(*canary::OptionalModeIndex {canary::ModeIndex {0}}));
static_assert(std::ranges::forward_range<canary::ModeIndexRange>);
static_assert(std::ranges::sized_range<canary::ModeIndexRange>);
static_assert(std::same_as<std::ranges::range_value_t<canary::ModeIndexRange>, canary::ModeIndex>);
static_assert(noexcept(canary::mode_indices(0, 1)));
static_assert(noexcept(++std::declval<canary::ModeIndexRange::iterator&>()));
static_assert(canary::mode_indices(2, 5).size() == 3);

static_assert(canary::ThermoStateOf<canary::ThermoState<double, si>, double, si>);
static_assert(canary::ThermoStateOf<canary::ThermoState<D, scaled>, D, scaled>);
static_assert(!canary::ThermoStateOf<canary::ThermoState<double, si>, double, scaled>);
static_assert(!canary::ThermoStateOf<canary::ThermoState<double, si>, float, si>);
static_assert(!canary::ThermoStateOf<NoTemperature, double, si>);
static_assert(std::is_trivially_copyable_v<canary::ThermoState<double, si>>);
static_assert(std::is_trivially_copyable_v<canary::Closure<double, si>>);
static_assert(std::is_trivially_copyable_v<canary::SingleTemperatureProperties<float, scaled>>);
static_assert(canary::ThermoStateOf<canary::MultiTemperatureState<double, si, 1>, double, si>);
static_assert(canary::ThermoStateOf<canary::MultiTemperatureState<D, scaled, 2>, D, scaled>);
static_assert(!canary::ThermoStateOf<canary::MultiTemperatureState<double, si, 1>, double, scaled>);
static_assert(std::is_trivially_copyable_v<canary::MultiTemperatureState<double, si, 1>>);
static_assert(sizeof(canary::MultiTemperatureState<double, si, 2>) == 8 * sizeof(double));
static_assert(sizeof(canary::ThermoState<double, si>) == 4 * sizeof(double));
static_assert(
  std::same_as<canary::Closure<double, si>, canary::ClosureOf<canary::ThermoState<double, si>>>);

static_assert(canary::AcousticPropertiesOf<canary::AcousticProperties<float, si>, float, si>);
static_assert(canary::AcousticPropertiesOf<canary::AcousticProperties<D, scaled>, D, scaled>);
static_assert(!canary::AcousticPropertiesOf<IncorrectAcousticProperties, double, si>);
static_assert(!canary::AcousticPropertiesOf<canary::AcousticProperties<double, si>, float, si>);
static_assert(
  !canary::AcousticPropertiesOf<canary::AcousticProperties<double, si>, double, scaled>);
static_assert(
  canary::SingleTemperaturePropertiesOf<canary::SingleTemperatureProperties<float, si>, float, si>);
static_assert(
  canary::SingleTemperaturePropertiesOf<canary::SingleTemperatureProperties<D, scaled>, D, scaled>);
static_assert(canary::MultiTemperaturePropertiesOf<canary::MultiTemperatureProperties<float, si, 1>,
                                                   float,
                                                   si>);
static_assert(canary::MultiTemperaturePropertiesOf<canary::MultiTemperatureProperties<D, scaled, 2>,
                                                   D,
                                                   scaled>);
static_assert(
  !canary::
    SingleTemperaturePropertiesOf<canary::MultiTemperatureProperties<double, si, 2>, double, si>);
static_assert(!canary::MultiTemperaturePropertiesOf<canary::SingleTemperatureProperties<double, si>,
                                                    double,
                                                    si>);
static_assert(canary::SingleTemperaturePropertiesOf<CustomProperties<>, double, si>);
static_assert(canary::MultiTemperaturePropertiesOf<CustomProperties<>, double, si>);
static_assert(
  !canary::
    SingleTemperaturePropertiesOf<CustomProperties<canary::AcousticExponent<double>>, double, si>);
static_assert(
  !canary::MultiTemperaturePropertiesOf<
    CustomProperties<canary::HeatCapacityRatio<double>, canary::SpecificHeat<double, scaled>>,
    double,
    si>);
static_assert(!canary::AcousticPropertiesOf<int, double, si>);
static_assert(!canary::SingleTemperaturePropertiesOf<int, double, si>);
static_assert(!canary::MultiTemperaturePropertiesOf<int, double, si>);
static_assert(std::is_trivially_copyable_v<canary::AcousticProperties<double, si>>);
static_assert(std::is_trivially_copyable_v<canary::SingleTemperatureProperties<float, scaled>>);
static_assert(std::is_trivially_copyable_v<canary::MultiTemperatureProperties<D, scaled, 2>>);
static_assert(sizeof(canary::AcousticProperties<double, si>) == 2 * sizeof(double));
static_assert(sizeof(canary::SingleTemperatureProperties<float, scaled>) == 6 * sizeof(float));
static_assert(sizeof(canary::MultiTemperatureProperties<D, scaled, 2>) == 7 * sizeof(D));

void test_state()
{
    using namespace boost::ut;

    "optional modes distinguish absence from mode zero"_test = [] {
        const canary::OptionalModeIndex absent {};
        const canary::OptionalModeIndex zero {canary::ModeIndex {0}};
        const canary::OptionalModeIndex last {
          canary::ModeIndex {std::numeric_limits<std::size_t>::max() - 1}};
        expect(!absent);
        expect(zero.has_value());
        expect((*zero).value() == 0u);
        expect(zero != absent);
        expect((*last).value() == std::numeric_limits<std::size_t>::max() - 1);
        expect(last == canary::OptionalModeIndex {*last});
    };

    "mode ranges preserve bounds and independent iterator positions"_test = [] {
        std::array<std::size_t, 3> visited {};
        std::size_t count = 0;
        for (const canary::ModeIndex mode : canary::mode_indices(2, 5)) {
            visited[count++] = mode.value();
        }
        expect(count == 3u);
        expect(visited == std::array<std::size_t, 3> {2, 3, 4});
        auto current = canary::mode_indices(2, 5).begin();
        const auto previous = current++;
        expect((*previous).value() == 2u && (*current).value() == 3u);
        const auto maximum = std::numeric_limits<std::size_t>::max();
        const auto empty = canary::mode_indices(maximum, maximum);
        expect(empty.size() == 0u && empty.begin() == empty.end());
        count = 0;
        for (const canary::ModeIndex mode : canary::mode_indices(maximum - 1, maximum)) {
            expect(mode.value() == maximum - 1);
            ++count;
        }
        expect(count == 1u);
    };

    "a closure holds the state members and the status"_test = [] {
        const canary::Closure<double, si> c {
          .state =
            {
              .rho = canary::Density<double, si> {1.2},
              .e = canary::SpecificEnergy<double, si> {2.0e5},
              .p = canary::Pressure<double, si> {1.0e5},
              .temp = canary::Temperature<double, si> {290.0},
            },
          .status = canary::Status::Valid,
          .diagnostic = {},
        };
        expect(c.state.rho.value() == 1.2_d);
        expect(c.state.temp.value() == 290.0_d);
        expect(c.status == canary::Status::Valid);
        expect(c.diagnostic.reason == canary::ClosureReason::None);
        expect(!c.diagnostic.mode.has_value());
    };
}
