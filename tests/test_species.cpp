#include "dual.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <gaslib/quantity.hpp>
#include <gaslib/species.hpp>
#include <gaslib/status.hpp>
#include <limits>
#include <mdspan/mdspan.hpp>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using gaslib::scaled;
using gaslib::Status;
namespace kind = gaslib::kind;

using Extents = Kokkos::dextents<std::size_t, 1>;
using View = Kokkos::mdspan<const double, Extents>;
using StridedView = Kokkos::mdspan<const double, Extents, Kokkos::layout_stride>;

constexpr double eps = std::numeric_limits<double>::epsilon();

bool within_ulp(double a, double b, double ulps)
{ return std::abs(a - b) <= ulps * eps * std::abs(b); }

template <std::size_t N> View view_of(const std::array<double, N>& values)
{ return View(values.data(), N); }

auto fractions(const View& v) { return gaslib::as_kind<kind::mass_fraction>(v); }
auto densities(const View& v) { return gaslib::as_kind<kind::density, scaled>(v); }

using Species2 = gaslib::SpeciesData<double, scaled, 2>;

// Scaled particle masses (proton mass units) and charge numbers.
constexpr Species2 two_species {
  .n = 2,
  .mass = {gaslib::ParticleMass<double, scaled> {1.0}, gaslib::ParticleMass<double, scaled> {16.0}},
  .charge = {gaslib::Charge<double, scaled> {1.0}, gaslib::Charge<double, scaled> {-1.0}},
};

}  // namespace

static_assert(std::is_trivially_copyable_v<Species2>);
static_assert(
  gaslib::SpeciesRangeOf<decltype(fractions(std::declval<View>())), kind::mass_fraction>);
static_assert(!gaslib::SpeciesRangeOf<decltype(fractions(std::declval<View>())), kind::density>);
static_assert(gaslib::SpeciesRangeOf<decltype(densities(std::declval<View>())), kind::density>);
static_assert(!gaslib::SpeciesRangeOf<View, kind::mass_fraction>);
static_assert(
  std::same_as<decltype(fractions(std::declval<View>()))::layout_type, Kokkos::layout_right>);
static_assert(std::same_as<decltype(gaslib::as_kind<kind::mass_fraction>(
                             std::declval<StridedView>()))::layout_type,
                           Kokkos::layout_stride>);

void test_species()
{
    using namespace boost::ut;
    using gaslib::composition_status;
    using gaslib::mean_charge;
    using gaslib::mean_particle_mass;

    "as_kind reads each element as the requested kind with the same extents"_test = [] {
        const std::array values {0.25, 0.75};
        const auto y = fractions(view_of(values));
        static_assert(std::same_as<decltype(y[0]), gaslib::MassFraction<double>>);
        expect(y.extent(0) == 2_u);
        expect(y[1].value() == 0.75_d);
        expect(y.data_handle() == values.data());
    };

    "composition status of mass fractions"_test = [] {
        const std::array valid {0.25, 0.75};
        const std::array over {0.25, 0.75 + 1.0e-10};
        const std::array negative {-1.0e-10, 1.0 + 1.0e-10};
        const std::array not_finite {std::numeric_limits<double>::quiet_NaN(), 1.0};
        const std::array all_zero {0.0, 0.0};
        expect(composition_status(fractions(view_of(valid))) == Status::Valid);
        expect(composition_status(fractions(view_of(over))) == Status::CompositionViolation);
        expect(composition_status(fractions(view_of(negative))) == Status::CompositionViolation);
        expect(composition_status(fractions(view_of(not_finite))) == Status::Invalid);
        expect(composition_status(fractions(view_of(all_zero))) == Status::Invalid);
    };

    "composition status of partial densities scales the tolerance with the density"_test = [] {
        const std::array valid {0.5, 1.5};
        const std::array over {0.5, 1.5 + 1.0e-9};
        expect(composition_status(densities(view_of(valid)),
                                  gaslib::Density<double, scaled> {2.0}) == Status::Valid);
        expect(
          composition_status(densities(view_of(over)), gaslib::Density<double, scaled> {2.0}) ==
          Status::CompositionViolation);
        expect(composition_status(densities(view_of(valid)),
                                  gaslib::Density<double, scaled> {0.0}) == Status::Invalid);
    };

    "mean mass and charge follow the number-weighted definitions"_test = [] {
        const std::array y {0.25, 0.75};
        const double number = (0.25 / 1.0) + (0.75 / 16.0);
        const auto m = mean_particle_mass(two_species, fractions(view_of(y)));
        const auto q = mean_charge(two_species, fractions(view_of(y)));
        static_assert(std::same_as<decltype(m), const gaslib::ParticleMass<double, scaled>>);
        expect(m.value() == (0.25 + 0.75) / number);
        expect(q.value() == (((0.25 / 1.0) * 1.0) + ((0.75 / 16.0) * -1.0)) / number);
    };

    "partial densities give the same results as the equivalent fractions"_test = [] {
        const std::array y {0.25, 0.75};
        const std::array rho_i {0.5, 1.5};
        const gaslib::Density<double, scaled> rho {2.0};
        expect(mean_particle_mass(two_species, densities(view_of(rho_i)), rho).value() ==
               mean_particle_mass(two_species, fractions(view_of(y))).value());
        expect(mean_charge(two_species, densities(view_of(rho_i)), rho).value() ==
               mean_charge(two_species, fractions(view_of(y))).value());
    };

    "weights outside [0, total] are clamped and normalised by their sum"_test = [] {
        const std::array y {1.2, -0.2};
        expect(mean_particle_mass(two_species, fractions(view_of(y))).value() == 1.0_d);
        expect(mean_charge(two_species, fractions(view_of(y))).value() == 1.0_d);
        const std::array short_sum {0.2, 0.3};  // sums to 0.5
        const std::array full_sum {0.4, 0.6};
        expect(within_ulp(mean_particle_mass(two_species, fractions(view_of(short_sum))).value(),
                          mean_particle_mass(two_species, fractions(view_of(full_sum))).value(),
                          2.0));
    };

    "results do not change when the species are permuted"_test = [] {
        using Species3 = gaslib::SpeciesData<double, scaled, 3>;
        const std::array<double, 3> masses {1.0, 4.0, 16.0};
        const std::array<double, 3> charges {1.0, 2.0, -1.0};
        const std::array<double, 3> y {0.2, 0.3, 0.5};
        const std::array<std::array<std::size_t, 3>, 6> orders {
          {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}},
        };
        double m0 = 0;
        double q0 = 0;
        for (const auto& order : orders) {
            Species3 s {};
            s.n = 3;
            std::array<double, 3> yp {};
            for (std::size_t i = 0; i < 3; ++i) {
                s.mass[i] = gaslib::ParticleMass<double, scaled> {masses[order[i]]};
                s.charge[i] = gaslib::Charge<double, scaled> {charges[order[i]]};
                yp[i] = y[order[i]];
            }
            const double m = mean_particle_mass(s, fractions(view_of(yp))).value();
            const double q = mean_charge(s, fractions(view_of(yp))).value();
            if (order[0] == 0 && order[1] == 1) {
                m0 = m;
                q0 = q;
            }
            expect(within_ulp(m, m0, 4.0) && within_ulp(q, q0, 4.0)) << m << m0 << q << q0;
        }
    };

    "identical species give that species' mass and charge"_test = [] {
        using Species5 = gaslib::SpeciesData<double, scaled, 5>;
        Species5 s {};
        s.n = 5;
        for (std::size_t i = 0; i < 5; ++i) {
            s.mass[i] = gaslib::ParticleMass<double, scaled> {2.5};
            s.charge[i] = gaslib::Charge<double, scaled> {-1.0};
        }
        const std::array y {0.1, 0.3, 0.05, 0.25, 0.3};
        expect(within_ulp(mean_particle_mass(s, fractions(view_of(y))).value(), 2.5, 4.0));
        expect(within_ulp(mean_charge(s, fractions(view_of(y))).value(), -1.0, 4.0));
    };

    "a strided view gives results identical to a contiguous view"_test = [] {
        const std::array y {0.25, 0.75};
        std::vector<double> storage(14, -99.0);
        storage[0] = y[0];
        storage[7] = y[1];
        const StridedView strided(
          storage.data(),
          Kokkos::layout_stride::mapping<Extents>(Extents(2), std::array<std::size_t, 1> {7}));
        const auto ys = gaslib::as_kind<kind::mass_fraction>(strided);
        expect(ys[1].value() == 0.75_d);
        expect(mean_particle_mass(two_species, ys).value() ==
               mean_particle_mass(two_species, fractions(view_of(y))).value());
        expect(mean_charge(two_species, ys).value() ==
               mean_charge(two_species, fractions(view_of(y))).value());
        expect(composition_status(ys) == Status::Valid);
    };

    "dual-number fractions carry the derivative of the mean mass"_test = [] {
        // Y1 = 0.25 + t, Y2 = 0.75 - t: d m / d t = -m^2 (1/m1 - 1/m2) since sum(Y) = 1.
        using D = gaslib::test::Dual<double>;
        const std::array<D, 2> y {D(0.25, 1.0), D(0.75, -1.0)};
        const Kokkos::mdspan<const D, Extents> v(y.data(), 2);
        const auto m =
          mean_particle_mass(two_species, gaslib::as_kind<kind::mass_fraction>(v)).value();
        const double expected = -m.real() * m.real() * (1.0 - (1.0 / 16.0));
        expect(within_ulp(m.dual(), expected, 8.0)) << m.dual() << expected;
    };
}
