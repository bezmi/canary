#include "test_suites.hpp"

#include <array>
#include <bit>
#include <boost/ut.hpp>
#include <canary/data/air11.hpp>
#include <canary/data/air5.hpp>
#include <canary/layout.hpp>
#include <canary/models/calorically_perfect.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/models/two_temperature.hpp>
#include <canary/thermo.hpp>
#include <cstddef>
#include <cstdint>

namespace
{
using canary::scaled;
using canary::si;

// A 4-byte member followed by an 8-byte one: 4 bytes of padding.
struct Gapped {
    std::int32_t n;
    double x;

    static consteval bool padding_free()
    { return canary::layout::members_fill<Gapped, decltype(n), decltype(x)>(); }
};

// A member missing from the list.
struct Unlisted {
    double x;
    double y;

    static consteval bool padding_free()
    { return canary::layout::members_fill<Unlisted, decltype(x)>(); }
};

struct Listed {
    std::size_t n;
    std::array<double, 3> x;

    static consteval bool padding_free()
    { return canary::layout::members_fill<Listed, decltype(n), decltype(x)>(); }
};

// A class that does not list its members.
struct Opaque {
    double x;
};

struct HoldsOpaque {
    Opaque o;

    static consteval bool padding_free()
    { return canary::layout::members_fill<HoldsOpaque, decltype(o)>(); }
};

// No padding of its own, but its elements have padding.
struct HoldsGapped {
    std::array<Gapped, 2> g;

    static consteval bool padding_free()
    { return canary::layout::members_fill<HoldsGapped, decltype(g)>(); }
};

enum class Small : unsigned char { a };

// True when the two models are equal byte for byte.
template <typename G> bool same_bytes(const G& a, const G& b)
{
    using Bytes = std::array<std::byte, sizeof(G)>;
    return std::bit_cast<Bytes>(a) == std::bit_cast<Bytes>(b);
}

}  // namespace

static_assert(canary::PaddingFreeScalar<double> && canary::PaddingFreeScalar<float> &&
              canary::PaddingFreeScalar<std::uint32_t> && canary::PaddingFreeScalar<bool> &&
              canary::PaddingFreeScalar<Small>);
static_assert(!canary::PaddingFreeScalar<long double>);
static_assert(canary::DeclaresMembers<Gapped> && !canary::DeclaresMembers<Opaque>);
static_assert(canary::PaddingFree<Listed> && canary::PaddingFree<const Listed>);
static_assert(!canary::PaddingFree<Gapped>);
static_assert(!canary::PaddingFree<Unlisted>);
static_assert(!canary::PaddingFree<Opaque> && !canary::PaddingFree<HoldsOpaque>);
static_assert(!canary::PaddingFree<HoldsGapped>);
static_assert(canary::PaddingFree<std::array<double, 4>> &&
              canary::PaddingFree<std::array<Listed, 2>>);
static_assert(!canary::PaddingFree<std::array<double, 0>> &&
              !canary::PaddingFree<std::array<Gapped, 2>>);
static_assert(canary::PaddingFree<canary::SpecificHeat<double, si>> &&
              canary::PaddingFree<canary::OptionalSpeciesIndex> &&
              canary::PaddingFree<canary::OptionalModeIndex>);

// The double models.
static_assert(canary::PaddingFree<canary::CaloricallyPerfect<double, si>>);
static_assert(canary::PaddingFree<canary::CaloricallyPerfectMixture<double, scaled, 16>>);
static_assert(canary::PaddingFree<canary::Air5<double, si>>);
static_assert(canary::PaddingFree<canary::Air11<double, si>>);
static_assert(canary::PaddingFree<canary::Air5TwoTemperature<double, si>>);
static_assert(canary::PaddingFree<canary::Air5FourTemperature<double, si>>);
// Its mode table, 3-byte entries for 11 species and 2 modes followed by a
// 4-byte index, leaves one byte unused.
static_assert(!canary::PaddingFree<canary::Air11ThreeTemperature<double, si>>);
static_assert(canary::GasModelType<canary::Air11ThreeTemperature<double, si>>);
// With float parameters the 8-byte counts leave gaps; the models remain gas models.
static_assert(canary::PaddingFree<canary::CaloricallyPerfect<float, si>>);
static_assert(!canary::PaddingFree<canary::Air5<float, scaled>>);
static_assert(!canary::PaddingFree<canary::Air5TwoTemperature<float, scaled>>);
static_assert(canary::GasModelType<canary::Air5TwoTemperature<float, scaled>>);

void test_layout()
{
    using namespace boost::ut;
    using Mass = canary::ParticleMass<double, scaled>;
    using Q = canary::Charge<double, scaled>;
    using Gamma = canary::HeatCapacityRatio<double>;
    constexpr canary::BoltzmannConstant<double, scaled> k_b {1.0};

    "models built separately from the same inputs are equal byte for byte"_test = [&] {
        const std::array masses {Mass {1.0}, Mass {16.0}};
        const std::array charges {Q {1.0}, Q {-1.0}};
        const std::array gammas {Gamma {5.0 / 3.0}, Gamma {1.4}};
        const auto a = canary::make_calorically_perfect_mixture(masses, charges, k_b, gammas).model;
        const auto b = canary::make_calorically_perfect_mixture(masses, charges, k_b, gammas).model;
        expect(same_bytes(a, b));
        const std::array other {Gamma {5.0 / 3.0}, Gamma {1.3}};
        expect(
          !same_bytes(a,
                      canary::make_calorically_perfect_mixture(masses, charges, k_b, other).model));

        const canary::UnitScales<double, si> scales {};
        expect(same_bytes(canary::make_air5<double>(scales).model,
                          canary::make_air5<double>(scales).model));
        expect(same_bytes(canary::make_air5_four_temperature<double>(scales).model,
                          canary::make_air5_four_temperature<double>(scales).model));
    };
}
