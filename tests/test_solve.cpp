#include "dual.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <boost/ut.hpp>
#include <cmath>
#include <canary/detail/solve.hpp>
#include <canary/status.hpp>
#include <limits>

namespace
{

using canary::Status;
using canary::detail::Bracket;
using canary::detail::Residual;
using canary::detail::SolveSettings;
using D = canary::test::Dual<double>;

constexpr double eps = std::numeric_limits<double>::epsilon();
constexpr SolveSettings<double> settings {.relative_tolerance = 1.0e-14, .max_iterations = 50};

// f(x) = x^3 + x - c, with root 2 for c = 10.
template <typename X> Residual<X> cubic(X x, X c)
{ return {.value = (x * x * x) + x - c, .derivative = (X(3.0) * x * x) + X(1.0)}; }

// f(x) = (x - 1) / (1 + (x - 1)^2): root 1; the derivative is negative for
// |x - 1| > 1, so Newton steps from there move away from the root.
Residual<double> rational(double x)
{
    const double u = x - 1.0;
    const double d = 1.0 + (u * u);
    return {.value = u / d, .derivative = (1.0 - (u * u)) / (d * d)};
}

}  // namespace

void test_solve()
{
    using namespace boost::ut;
    using canary::detail::solve;

    "the root of a cubic is found within the tolerance"_test = [] {
        const auto f = [](double x) { return cubic(x, 10.0); };
        const auto r =
          solve<double>(f, f, Bracket<double> {.lower = 0.0, .upper = 5.0}, 1.0, settings);
        expect(r.status == Status::Valid);
        expect(std::abs(r.root - 2.0) <= 4.0 * eps * 2.0) << r.root;
    };

    "the dual derivative of the root is -f_c / f_x"_test = [] {
        // dx/dc = 1 / (3 x^2 + 1) = 1/13 at x = 2.
        const D c(10.0, 1.0);
        const auto f_value = [&](double x) { return cubic(x, c.real()); };
        const auto f_full = [&](D x) { return cubic(x, c); };
        const auto r =
          solve<D>(f_value, f_full, Bracket<double> {.lower = 0.0, .upper = 5.0}, 1.0, settings);
        expect(r.status == Status::Valid);
        expect(std::abs(r.root.real() - 2.0) <= 4.0 * eps * 2.0);
        expect(std::abs(r.root.dual() - (1.0 / 13.0)) <= 8.0 * eps / 13.0) << r.root.dual();
    };

    "a start at the bracket edge converges"_test = [] {
        const auto f = [](double x) { return cubic(x, 10.0); };
        for (const double guess : {0.0, 5.0, -1.0, 9.0}) {
            const auto r =
              solve<double>(f, f, Bracket<double> {.lower = 0.0, .upper = 5.0}, guess, settings);
            expect(r.status == Status::Valid) << guess;
            expect(std::abs(r.root - 2.0) <= 4.0 * eps * 2.0) << guess;
        }
    };

    "bisection takes over where Newton steps leave the bracket"_test = [] {
        const auto r = solve<double>(rational,
                                     rational,
                                     Bracket<double> {.lower = -3.0, .upper = 4.0},
                                     3.5,
                                     settings);
        expect(r.status == Status::Valid);
        expect(std::abs(r.root - 1.0) <= 4.0 * eps) << r.root;
    };

    "a solve that converges from one side ends Valid"_test = [] {
        // exp(x) - c is convex and increasing: from the upper bracket end
        // every Newton iterate lies above the root ln(c), so the bracket's
        // lower end stays at -10, and the last Newton step often rounds to
        // the iterate itself.
        int not_valid = 0;
        double worst = 0.0;
        for (int k = 0; k < 10000; ++k) {
            const double c = 1.0 + (0.001 * k);
            const auto f = [c](double x) {
                return Residual<double> {.value = std::exp(x) - c, .derivative = std::exp(x)};
            };
            const auto r =
              solve<double>(f, f, Bracket<double> {.lower = -10.0, .upper = 10.0}, 10.0, settings);
            not_valid += r.status == Status::Valid ? 0 : 1;
            worst = std::max(worst, std::abs(r.root - std::log(c)));
        }
        expect(not_valid == 0_i);
        // The tolerance 1e-14 |x| with |x| <= ln(11), plus rounding of log.
        expect(worst <= 1.0e-13) << worst;
    };

    "one iteration on a hard case gives NotConverged"_test = [] {
        const auto r =
          solve<double>(rational,
                        rational,
                        Bracket<double> {.lower = -3.0, .upper = 4.0},
                        3.5,
                        SolveSettings<double> {.relative_tolerance = 1.0e-14, .max_iterations = 1});
        expect(r.status == Status::NotConverged);
    };

    "a bracket without a sign change gives Invalid"_test = [] {
        const auto f = [](double x) { return cubic(x, 10.0); };
        expect(
          solve<double>(f, f, Bracket<double> {.lower = 3.0, .upper = 5.0}, 4.0, settings).status ==
          Status::Invalid);
        const auto nan = [](double) {
            return Residual<double> {
              .value = std::numeric_limits<double>::quiet_NaN(),
              .derivative = 1.0,
            };
        };
        expect(solve<double>(nan, nan, Bracket<double> {.lower = 0.0, .upper = 1.0}, 0.5, settings)
                 .status == Status::Invalid);
    };
}
