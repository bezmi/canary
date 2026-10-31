#include "dual.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <canary/quantity.hpp>
#include <canary/species.hpp>
#include <canary/status.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using canary::scaled;
using canary::Status;
namespace kind = canary::kind;

using Extents = Kokkos::dextents<std::size_t, 1>;
using View = Kokkos::mdspan<const double, Extents>;
using StridedView = Kokkos::mdspan<const double, Extents, Kokkos::layout_stride>;

// A generated range exercises the extent guard without allocating its entries.
struct GeneratedFractions {
    std::size_t count;
    std::size_t* reads;
    static constexpr std::size_t rank() { return 1; }
    std::size_t extent(std::size_t) const { return count; }
    canary::MassFraction<double> operator[](std::size_t) const
    {
        ++*reads;
        return canary::MassFraction<double> {1.0};
    }
};

constexpr double eps = std::numeric_limits<double>::epsilon();

bool within_ulp(double a, double b, double ulps)
{ return std::abs(a - b) <= ulps * eps * std::abs(b); }

template <std::size_t N> View view_of(const std::array<double, N>& values)
{ return View(values.data(), N); }

auto fractions(const View& v) { return canary::as_kind<kind::mass_fraction>(v); }
auto densities(const View& v) { return canary::as_kind<kind::density, scaled>(v); }

using Species2 = canary::SpeciesData<double, scaled, 2>;

// Scaled particle masses (proton mass units) and charge numbers.
constexpr Species2 two_species {
  .n = 2,
  .mass = {canary::ParticleMass<double, scaled> {1.0}, canary::ParticleMass<double, scaled> {16.0}},
  .charge = {canary::Charge<double, scaled> {1.0}, canary::Charge<double, scaled> {-1.0}},
};

}  // namespace

static_assert(std::is_trivially_copyable_v<Species2>);
static_assert(
  canary::SpeciesRangeOf<decltype(fractions(std::declval<View>())), kind::mass_fraction>);
static_assert(!canary::SpeciesRangeOf<decltype(fractions(std::declval<View>())), kind::density>);
static_assert(canary::SpeciesRangeOf<decltype(densities(std::declval<View>())), kind::density>);
static_assert(!canary::SpeciesRangeOf<View, kind::mass_fraction>);
static_assert(
  std::same_as<decltype(fractions(std::declval<View>()))::layout_type, Kokkos::layout_right>);
static_assert(std::same_as<decltype(canary::as_kind<kind::mass_fraction>(
                             std::declval<StridedView>()))::layout_type,
                           Kokkos::layout_stride>);

static_assert(std::is_nothrow_default_constructible_v<canary::SpeciesIndex>);
static_assert(std::is_nothrow_constructible_v<canary::SpeciesIndex, std::size_t>);
static_assert(noexcept(std::declval<canary::SpeciesIndex>().value()));
static_assert(noexcept(canary::SpeciesIndex {0} == canary::SpeciesIndex {1}));
static_assert(std::is_nothrow_default_constructible_v<canary::OptionalSpeciesIndex>);
static_assert(std::is_nothrow_constructible_v<canary::OptionalSpeciesIndex, canary::SpeciesIndex>);
static_assert(noexcept(std::declval<canary::OptionalSpeciesIndex>().has_value()));
static_assert(noexcept(bool(std::declval<canary::OptionalSpeciesIndex>())));
static_assert(noexcept(*std::declval<canary::OptionalSpeciesIndex>()));
static_assert(noexcept(canary::OptionalSpeciesIndex {} == canary::OptionalSpeciesIndex {}));
static_assert(noexcept(canary::species_indices(0, 1)));
static_assert(noexcept(std::declval<const canary::SpeciesIndexRange&>().begin()));
static_assert(noexcept(std::declval<const canary::SpeciesIndexRange&>().end()));
static_assert(noexcept(std::declval<const canary::SpeciesIndexRange&>().size()));
static_assert(noexcept(*std::declval<canary::SpeciesIndexRange::iterator>()));
static_assert(noexcept(++std::declval<canary::SpeciesIndexRange::iterator&>()));
static_assert(noexcept(std::declval<canary::SpeciesIndexRange::iterator&>()++));
static_assert(noexcept(canary::SpeciesIndexRange::iterator {} ==
                       canary::SpeciesIndexRange::iterator {}));
static_assert(noexcept(canary::SpeciesIndexRange::iterator {} ==
                       canary::SpeciesIndexRange::sentinel {}));

static_assert(sizeof(canary::SpeciesIndex) == sizeof(std::uint32_t));
static_assert(sizeof(canary::OptionalSpeciesIndex) == sizeof(std::uint32_t));
static_assert(alignof(canary::SpeciesIndex) == alignof(std::uint32_t));
static_assert(alignof(canary::OptionalSpeciesIndex) == alignof(std::uint32_t));
static_assert(std::is_trivially_copyable_v<canary::SpeciesIndex>);
static_assert(std::is_trivially_copyable_v<canary::OptionalSpeciesIndex>);
static_assert(std::is_standard_layout_v<canary::SpeciesIndex>);
static_assert(std::is_standard_layout_v<canary::OptionalSpeciesIndex>);
static_assert(!std::is_convertible_v<std::size_t, canary::SpeciesIndex>);
static_assert(!std::is_convertible_v<canary::SpeciesIndex, std::size_t>);
static_assert(!std::is_convertible_v<canary::OptionalSpeciesIndex, canary::SpeciesIndex>);
static_assert(canary::SpeciesIndex {0}.value() == 0);
static_assert(!canary::OptionalSpeciesIndex {});
static_assert(*canary::OptionalSpeciesIndex {canary::SpeciesIndex {7}} == canary::SpeciesIndex {7});

static_assert(std::ranges::forward_range<canary::SpeciesIndexRange>);
static_assert(std::ranges::sized_range<canary::SpeciesIndexRange>);
static_assert(
  std::same_as<std::ranges::range_value_t<canary::SpeciesIndexRange>, canary::SpeciesIndex>);
static_assert(std::is_trivially_copyable_v<canary::SpeciesIndexRange>);
static_assert(std::is_trivially_copyable_v<canary::SpeciesIndexRange::iterator>);
static_assert(sizeof(canary::SpeciesIndexRange) == 2 * sizeof(std::uint32_t));
static_assert(sizeof(canary::SpeciesIndexRange::iterator) == sizeof(std::size_t));
static_assert(sizeof(canary::SpeciesIndexRange::sentinel) == sizeof(std::uint32_t));
static_assert(std::is_trivially_copyable_v<canary::SpeciesIndexRange::sentinel>);
static_assert(
  std::sentinel_for<canary::SpeciesIndexRange::sentinel, canary::SpeciesIndexRange::iterator>);
static_assert(!std::ranges::common_range<canary::SpeciesIndexRange>);
static_assert(canary::SpeciesIndexRange::iterator {} == canary::SpeciesIndexRange::sentinel {});
static_assert(!std::is_constructible_v<canary::SpeciesIndexRange, std::size_t, std::size_t>);
static_assert([] {
    std::uint32_t next = 3;
    for (const canary::SpeciesIndex i : canary::species_indices(3, 7)) {
        if (i.value() != next++) { return false; }
    }
    return next == 7 && canary::species_indices(3, 7).size() == 4;
}());
static_assert([] {
    constexpr auto end = std::numeric_limits<std::uint32_t>::max();
    const auto last = canary::species_indices(end - 1, end);
    auto it = last.begin();
    const auto previous = it++;
    return (*previous).value() == end - 1 && it == last.end() &&
           canary::species_indices(end, end).size() == 0 &&
           canary::species_indices(end, end).begin() == canary::species_indices(end, end).end();
}());

void test_species()
{
    using namespace boost::ut;
    using canary::mean_specific_charge;
    using canary::check_composition;
    using canary::mean_charge;
    using canary::mean_particle_mass;

    "composition rejects an unrepresentable species count before reading weights"_test = [] {
        if constexpr (std::numeric_limits<std::size_t>::max() >
                      std::numeric_limits<std::uint32_t>::max()) {
            std::size_t reads = 0;
            const GeneratedFractions weights {
              std::size_t(std::numeric_limits<std::uint32_t>::max()) + 1,
              &reads};
            const auto check = canary::check_composition(weights);
            expect(check.status == Status::Invalid);
            expect(check.reason == canary::CompositionReason::SpeciesCountOverflow);
            expect(reads == 0_ul);
        }
    };

    "species ranges preserve subranges and independent iterator positions"_test = [] {
        for (const std::size_t count : {std::size_t(0), std::size_t(1), std::size_t(11)}) {
            for (std::size_t first = 0; first <= count; ++first) {
                const auto indices = canary::species_indices(first, count);
                std::vector<std::uint32_t> positions;
                for (const canary::SpeciesIndex i : indices) { positions.push_back(i.value()); }
                expect(positions.size() == count - first);
                expect(indices.size() == count - first);
                expect(std::ranges::distance(indices.begin(), indices.end()) ==
                       static_cast<std::ptrdiff_t>(count - first));
                for (std::size_t k = 0; k < positions.size(); ++k) {
                    expect(positions[k] == first + k);
                }
                if (first < count) {
                    auto cursor = indices.begin();
                    const auto saved = cursor;
                    ++cursor;
                    expect((*saved).value() == first);
                    expect(cursor != saved);
                }
            }
        }
    };

    "species identities preserve zero and the largest present index"_test = [] {
        for (const std::size_t value :
             {std::size_t(0),
              std::size_t(11),
              std::size_t(std::numeric_limits<std::uint32_t>::max()) - 1}) {
            const canary::SpeciesIndex required {value};
            const canary::OptionalSpeciesIndex present {required};
            expect(required.value() == value);
            expect(present.has_value());
            expect(*present == required);
            expect(present != canary::OptionalSpeciesIndex {});
            const auto copied = present;
            expect(copied == present);
        }
        const std::array<std::string_view, 3> names {"N2", "carrier", "e-"};
        expect(canary::detail::species_index(names, "N2") ==
               canary::OptionalSpeciesIndex {canary::SpeciesIndex {0}});
        expect(canary::detail::species_index(names, "e-") ==
               canary::OptionalSpeciesIndex {canary::SpeciesIndex {2}});
        expect(!canary::detail::species_index(names, "missing"));
        expect(!canary::detail::species_index(std::array<std::string_view, 0> {}, "missing"));
    };

    "as_kind reads each element as the requested kind with the same extents"_test = [] {
        const std::array values {0.25, 0.75};
        const auto y = fractions(view_of(values));
        static_assert(std::same_as<decltype(y[0]), canary::MassFraction<double>>);
        expect(y.extent(0) == 2_u);
        expect(y[1].value() == 0.75_d);
        expect(y.data_handle() == values.data());
    };

    "composition checks retain causes for strided scalar and Dual inputs"_test = [] {
        const auto run = []<canary::Scalar T, canary::UnitSystem U>() {
            using V = canary::ScalarTraits<T>::value_type;
            using Reason = canary::CompositionReason;
            using Layout = Kokkos::layout_stride;
            const Layout::mapping<Extents> mapping(Extents(2), std::array<std::size_t, 1> {2});
            const auto check = [&](V first, V second, V total, Status status, Reason reason) {
                const std::array<T, 5> storage {T(V(99)), T(first), T(V(99)), T(second), T(V(99))};
                const auto view = canary::as_kind<kind::density, U>(
                  Kokkos::mdspan<const T, Extents, Layout>(storage.data() + 1, mapping));
                const auto result = canary::check_composition(view, canary::Density<T, U> {T(total)});
                expect(result.status == status && result.reason == reason);
            };
            const V nan = std::numeric_limits<V>::quiet_NaN();
            const V inf = std::numeric_limits<V>::infinity();
            check(V(0.25), V(0.75), V(1), Status::Valid, Reason::None);
            check(V(0.25), V(0.5), V(1), Status::CompositionViolation, Reason::SumMismatch);
            check(V(-0.25), V(1.25), V(1), Status::CompositionViolation, Reason::WeightOutOfRange);
            check(V(0), V(0), V(1), Status::Invalid, Reason::NoPositiveWeight);
            check(V(-1), V(-1), V(1), Status::Invalid, Reason::NoPositiveWeight);
            check(V(0.25), V(0.75), V(0), Status::Invalid, Reason::NonpositiveTotal);
            check(V(0.25), V(0.75), V(-1), Status::Invalid, Reason::NonpositiveTotal);
            for (const V value : {nan, inf, -inf}) {
                check(value, V(0), V(1), Status::Invalid, Reason::NonfiniteWeight);
                check(V(-1), value, V(1), Status::Invalid, Reason::NonfiniteWeight);
                check(value, V(-1), value, Status::Invalid, Reason::NonfiniteTotal);
            }
        };
        run.template operator()<float, canary::si>();
        run.template operator()<double, canary::si>();
        run.template operator()<float, scaled>();
        run.template operator()<double, scaled>();
        run.template operator()<canary::test::Dual<float>, scaled>();
        run.template operator()<canary::test::Dual<double>, canary::si>();
        static_assert(!std::is_convertible_v<canary::CompositionCheck, Status>);
        static_assert(std::is_trivially_copyable_v<canary::CompositionCheck>);
    };

    "composition clamps retain interior derivatives at exact bounds"_test = [] {
        const auto run = []<std::floating_point P, canary::UnitSystem U>() {
            using D = canary::test::Dual<P>;
            using Layout = Kokkos::layout_stride;
            const Layout::mapping<Extents> mapping(Extents(2), std::array<std::size_t, 1> {2});
            const auto check = [&](D weight, D expected) {
                const std::array<D, 5> storage {D(99), weight, D(99), D(0), D(99)};
                const auto view = canary::as_kind<kind::density, U>(
                  Kokkos::mdspan<const D, Extents, Layout>(storage.data() + 1, mapping));
                const auto result =
                  canary::detail::clamped_weight(view, canary::SpeciesIndex {0}, D(2, 3));
                expect(result.real() == expected.real());
                expect(result.dual() == expected.dual());
            };
            check(D(-1, 7), D(0));
            check(D(0, 7), D(0, 7));
            check(D(-P(0), 7), D(-P(0), 7));
            check(D(1, 7), D(1, 7));
            check(D(2, 7), D(2, 7));
            check(D(3, 7), D(2, 3));

            const canary::SpeciesData<P, U, 2> species {
              .n = 2,
              .mass = {canary::ParticleMass<P, U> {P(1)}, canary::ParticleMass<P, U> {P(16)}},
              .charge = {canary::Charge<P, U> {P(1)}, canary::Charge<P, U> {P(-1)}},
            };
            // A transfer into the absent second species keeps sum(Y) = 1.
            const std::array<D, 2> y {D(1, -1), D(0, 1)};
            const auto view = canary::as_kind<kind::mass_fraction>(
              Kokkos::mdspan<const D, Extents>(y.data(), y.size()));
            const auto mass = mean_particle_mass(species, view).value();
            expect(mass.real() == P(1));
            expect(mass.dual() == P(15) / P(16));
        };
        run.template operator()<float, canary::si>();
        run.template operator()<double, canary::si>();
        run.template operator()<float, scaled>();
        run.template operator()<double, scaled>();
    };

    "composition status of mass fractions"_test = [] {
        const std::array valid {0.25, 0.75};
        const std::array over {0.25, 0.75 + 1.0e-10};
        const std::array negative {-1.0e-10, 1.0 + 1.0e-10};
        const std::array not_finite {std::numeric_limits<double>::quiet_NaN(), 1.0};
        const std::array all_zero {0.0, 0.0};
        expect(check_composition(fractions(view_of(valid))).status == Status::Valid);
        expect(check_composition(fractions(view_of(over))).status == Status::CompositionViolation);
        expect(check_composition(fractions(view_of(negative))).status == Status::CompositionViolation);
        expect(check_composition(fractions(view_of(not_finite))).status == Status::Invalid);
        expect(check_composition(fractions(view_of(all_zero))).status == Status::Invalid);
    };

    "composition status of partial densities scales the tolerance with the density"_test = [] {
        const std::array valid {0.5, 1.5};
        const std::array over {0.5, 1.5 + 1.0e-9};
        expect(check_composition(densities(view_of(valid)),
                                  canary::Density<double, scaled> {2.0}).status == Status::Valid);
        expect(
          check_composition(densities(view_of(over)), canary::Density<double, scaled> {2.0}).status ==
          Status::CompositionViolation);
        expect(check_composition(densities(view_of(valid)),
                                  canary::Density<double, scaled> {0.0}).status == Status::Invalid);
    };

    "mean mass and charge follow the number-weighted definitions"_test = [] {
        const std::array y {0.25, 0.75};
        const double number = (0.25 / 1.0) + (0.75 / 16.0);
        const auto m = mean_particle_mass(two_species, fractions(view_of(y)));
        const auto q = mean_charge(two_species, fractions(view_of(y)));
        static_assert(std::same_as<decltype(m), const canary::ParticleMass<double, scaled>>);
        expect(m.value() == (0.25 + 0.75) / number);
        expect(q.value() == (((0.25 / 1.0) * 1.0) + ((0.75 / 16.0) * -1.0)) / number);
    };

    "charge-to-mass ratio is the number-weighted charge over the clamped weight"_test = [] {
        const std::array y {0.25, 0.75};
        const auto r = mean_specific_charge(two_species, fractions(view_of(y)));
        static_assert(std::same_as<decltype(r), const canary::SpecificCharge<double, scaled>>);
        expect(r.value() == (((0.25 / 1.0) * 1.0) + ((0.75 / 16.0) * -1.0)) / (0.25 + 0.75));
    };

    "charge-to-mass ratio equals the mean charge over the mean particle mass"_test = [] {
        using Species3 = canary::SpeciesData<double, scaled, 3>;
        Species3 s {};
        s.n = 3;
        const std::array<double, 3> masses {1.0, 3.0, 7.0};
        const std::array<double, 3> charges {1.0, 2.0, -1.0};
        for (std::size_t i = 0; i < 3; ++i) {
            s.mass[i] = canary::ParticleMass<double, scaled> {masses[i]};
            s.charge[i] = canary::Charge<double, scaled> {charges[i]};
        }
        const std::array<std::array<double, 3>, 3> compositions {
          {{0.2, 0.3, 0.5}, {0.9, 0.05, 0.05}, {1.0e-6, 0.5, 0.5 - 1.0e-6}},
        };
        for (const auto& y : compositions) {
            const double r = mean_specific_charge(s, fractions(view_of(y))).value();
            const double expected = mean_charge(s, fractions(view_of(y))).value() /
                                    mean_particle_mass(s, fractions(view_of(y))).value();
            expect(within_ulp(r, expected, 8.0)) << r << expected;
        }
        for (std::size_t i = 0; i < 3; ++i) { s.charge[i] = canary::Charge<double, scaled> {0.0}; }
        expect(mean_specific_charge(s, fractions(view_of(compositions[0]))).value() == 0.0_d);
    };

    "partial densities give the same results as the equivalent fractions"_test = [] {
        const std::array y {0.25, 0.75};
        const std::array rho_i {0.5, 1.5};
        const canary::Density<double, scaled> rho {2.0};
        expect(mean_particle_mass(two_species, densities(view_of(rho_i)), rho).value() ==
               mean_particle_mass(two_species, fractions(view_of(y))).value());
        expect(mean_charge(two_species, densities(view_of(rho_i)), rho).value() ==
               mean_charge(two_species, fractions(view_of(y))).value());
        expect(mean_specific_charge(two_species, densities(view_of(rho_i)), rho).value() ==
               mean_specific_charge(two_species, fractions(view_of(y))).value());
    };

    "weights outside [0, total] are clamped and normalised by their sum"_test = [] {
        const std::array y {1.2, -0.2};
        expect(mean_particle_mass(two_species, fractions(view_of(y))).value() == 1.0_d);
        expect(mean_charge(two_species, fractions(view_of(y))).value() == 1.0_d);
        expect(mean_specific_charge(two_species, fractions(view_of(y))).value() == 1.0_d);
        const std::array short_sum {0.2, 0.3};  // sums to 0.5
        const std::array full_sum {0.4, 0.6};
        expect(within_ulp(mean_particle_mass(two_species, fractions(view_of(short_sum))).value(),
                          mean_particle_mass(two_species, fractions(view_of(full_sum))).value(),
                          2.0));
        expect(within_ulp(mean_specific_charge(two_species, fractions(view_of(short_sum))).value(),
                          mean_specific_charge(two_species, fractions(view_of(full_sum))).value(),
                          2.0));
    };

    "results do not change when the species are permuted"_test = [] {
        using Species3 = canary::SpeciesData<double, scaled, 3>;
        const std::array<double, 3> masses {1.0, 4.0, 16.0};
        const std::array<double, 3> charges {1.0, 2.0, -1.0};
        const std::array<double, 3> y {0.2, 0.3, 0.5};
        const std::array<std::array<std::size_t, 3>, 6> orders {
          {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}},
        };
        double m0 = 0;
        double q0 = 0;
        double r0 = 0;
        for (const auto& order : orders) {
            Species3 s {};
            s.n = 3;
            std::array<double, 3> yp {};
            for (std::size_t i = 0; i < 3; ++i) {
                s.mass[i] = canary::ParticleMass<double, scaled> {masses[order[i]]};
                s.charge[i] = canary::Charge<double, scaled> {charges[order[i]]};
                yp[i] = y[order[i]];
            }
            const double m = mean_particle_mass(s, fractions(view_of(yp))).value();
            const double q = mean_charge(s, fractions(view_of(yp))).value();
            const double r = mean_specific_charge(s, fractions(view_of(yp))).value();
            if (order[0] == 0 && order[1] == 1) {
                m0 = m;
                q0 = q;
                r0 = r;
            }
            expect(within_ulp(m, m0, 4.0) && within_ulp(q, q0, 4.0)) << m << m0 << q << q0;
            expect(within_ulp(r, r0, 4.0)) << r << r0;
        }
    };

    "identical species give that species' mass and charge"_test = [] {
        using Species5 = canary::SpeciesData<double, scaled, 5>;
        Species5 s {};
        s.n = 5;
        for (std::size_t i = 0; i < 5; ++i) {
            s.mass[i] = canary::ParticleMass<double, scaled> {2.5};
            s.charge[i] = canary::Charge<double, scaled> {-1.0};
        }
        const std::array y {0.1, 0.3, 0.05, 0.25, 0.3};
        expect(within_ulp(mean_particle_mass(s, fractions(view_of(y))).value(), 2.5, 4.0));
        expect(within_ulp(mean_charge(s, fractions(view_of(y))).value(), -1.0, 4.0));
        expect(within_ulp(mean_specific_charge(s, fractions(view_of(y))).value(), -1.0 / 2.5, 4.0));
    };

    "a strided view gives results identical to a contiguous view"_test = [] {
        const std::array y {0.25, 0.75};
        std::vector<double> storage(14, -99.0);
        storage[0] = y[0];
        storage[7] = y[1];
        const StridedView strided(
          storage.data(),
          Kokkos::layout_stride::mapping<Extents>(Extents(2), std::array<std::size_t, 1> {7}));
        const auto ys = canary::as_kind<kind::mass_fraction>(strided);
        expect(ys[1].value() == 0.75_d);
        expect(mean_particle_mass(two_species, ys).value() ==
               mean_particle_mass(two_species, fractions(view_of(y))).value());
        expect(mean_charge(two_species, ys).value() ==
               mean_charge(two_species, fractions(view_of(y))).value());
        expect(mean_specific_charge(two_species, ys).value() ==
               mean_specific_charge(two_species, fractions(view_of(y))).value());
        expect(check_composition(ys).status == Status::Valid);
    };

    "dual-number fractions carry the derivative of the mean mass"_test = [] {
        // Y1 = 0.25 + t, Y2 = 0.75 - t: d m / d t = -m^2 (1/m1 - 1/m2) since sum(Y) = 1.
        using D = canary::test::Dual<double>;
        const std::array<D, 2> y {D(0.25, 1.0), D(0.75, -1.0)};
        const Kokkos::mdspan<const D, Extents> v(y.data(), 2);
        const auto m =
          mean_particle_mass(two_species, canary::as_kind<kind::mass_fraction>(v)).value();
        const double expected = -m.real() * m.real() * (1.0 - (1.0 / 16.0));
        expect(within_ulp(m.dual(), expected, 8.0)) << m.dual() << expected;
    };
}
