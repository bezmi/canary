#include "dual.hpp"
#include "scalar_types.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <boost/ut.hpp>
#include <cmath>
#include <concepts>
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

constexpr auto accepts_solve = []<canary::Scalar T>() {
    using V = canary::ScalarTraits<T>::value_type;
    const auto f_value = [](V x) { return Residual<V> {.value = x, .derivative = V(1)}; };
    const auto f_full = [](T x) { return Residual<T> {.value = x, .derivative = T(V(1))}; };
    return requires {
        canary::detail::solve<T>(
          f_value,
          f_full,
          Bracket<V> {.lower = V(-1), .upper = V(1)},
          V(0),
          SolveSettings<V> {.relative_tolerance = V(0.1), .max_iterations = 1});
    };
};

static_assert(accepts_solve.template operator()<float>());
static_assert(accepts_solve.template operator()<double>());
static_assert(accepts_solve.template operator()<canary::test::Dual<float>>());
static_assert(accepts_solve.template operator()<D>());
static_assert(canary::Scalar<canary::test::OnlySqrt>);
static_assert(!accepts_solve.template operator()<canary::test::OnlySqrt>());

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

    "dual solves preserve the converged value and use its implicit derivative"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            const Dual coefficient(V(1), V(1));
            const auto f_value = [&](V x) {
                return Residual<V> {
                  .value = (coefficient.real() * x * x * x) + x - V(10),
                  .derivative = (V(3) * coefficient.real() * x * x) + V(1),
                };
            };
            const auto f_full = [&](Dual x) {
                return Residual<Dual> {
                  .value = (coefficient * x * x * x) + x - V(10),
                  .derivative = (V(3) * coefficient * x * x) + V(1),
                };
            };
            // A loose step tolerance leaves a nonzero residual so that a
            // separate final Newton step would change the returned value.
            const SolveSettings<V> limits {.relative_tolerance = V(0.1), .max_iterations = 50};
            const Bracket<V> bracket {.lower = V(0), .upper = V(5)};
            const auto plain = solve<V>(f_value, f_value, bracket, V(2.125), limits);
            const auto dual = solve<Dual>(f_value, f_full, bracket, V(2.125), limits);
            const V x = plain.root;
            const V expected = -(x * x * x) / ((V(3) * x * x) + V(1));
            expect(plain.status == Status::Valid && dual.status == Status::Valid);
            expect(f_value(x).value != V(0));
            expect(dual.root.real() == x);
            expect(std::abs(dual.root.dual() - expected) <=
                   V(8) * std::numeric_limits<V>::epsilon() * std::abs(expected));
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "implicit derivatives remain finite for extreme finite residual slopes"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            for (const V slope : {std::numeric_limits<V>::min(),
                                  -std::numeric_limits<V>::min(),
                                  std::numeric_limits<V>::max() / V(2),
                                  -std::numeric_limits<V>::max() / V(2)}) {
                const auto f_value = [slope](V x) {
                    return Residual<V> {.value = slope * x, .derivative = slope};
                };
                const auto f_full = [slope](Dual x) {
                    return Residual<Dual> {.value = slope * x - Dual(V(0), V(1)),
                                           .derivative = Dual(slope)};
                };
                const auto r =
                  solve<Dual>(f_value,
                              f_full,
                              Bracket<V> {.lower = V(-1), .upper = V(1)},
                              V(0),
                              SolveSettings<V> {.relative_tolerance = V(0.1), .max_iterations = 1});
                expect(r.status == Status::Valid);
                expect(r.root.real() == V(0));
                expect(std::isfinite(r.root.dual()));
                expect(r.root.dual() == V(1) / slope);
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "unusable full residual value parts return an Invalid constant root"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            const auto f_value = [](V x) {
                return Residual<V> {.value = x - V(1), .derivative = V(1)};
            };
            const Bracket<V> bracket {.lower = V(0), .upper = V(2)};
            const SolveSettings<V> limits {.relative_tolerance = V(0.1), .max_iterations = 50};
            const auto plain = solve<V>(f_value, f_value, bracket, V(0.5), limits);
            expect(plain.status == Status::Valid);
            const auto check_invalid = [&](V value, V derivative) {
                const auto f_full = [=](Dual) {
                    return Residual<Dual> {.value = Dual(value, V(1)),
                                           .derivative = Dual(derivative)};
                };
                const auto r = solve<Dual>(f_value, f_full, bracket, V(0.5), limits);
                expect(r.status == Status::Invalid);
                expect(r.root.real() == plain.root);
                expect(r.root.real() >= bracket.lower && r.root.real() <= bracket.upper);
                expect(r.root.dual() == V(0));
            };
            check_invalid(V(0), V(0));
            for (const V value : {std::numeric_limits<V>::quiet_NaN(),
                                  std::numeric_limits<V>::infinity(),
                                  -std::numeric_limits<V>::infinity()}) {
                check_invalid(value, V(1));
                check_invalid(V(0), value);
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "a singular dual root returns an Invalid constant"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            const auto f_value = [](V x) {
                return Residual<V> {.value = x * x * x, .derivative = V(3) * x * x};
            };
            const auto f_full = [](Dual x) {
                return Residual<Dual> {.value = x * x * x - Dual(V(0), V(1)),
                                       .derivative = V(3) * x * x};
            };
            const auto r =
              solve<Dual>(f_value,
                          f_full,
                          Bracket<V> {.lower = V(-1), .upper = V(1)},
                          V(0),
                          SolveSettings<V> {.relative_tolerance = V(0.1), .max_iterations = 1});
            expect(r.status == Status::Invalid);
            expect(r.root.real() == V(0));
            expect(r.root.dual() == V(0));
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "an exact endpoint root is accepted before iteration in either orientation"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            const Bracket<V> bracket {.lower = V(0), .upper = V(2)};
            const SolveSettings<V> limits {.relative_tolerance = V(0.1), .max_iterations = 0};
            for (const V endpoint : {bracket.lower, bracket.upper}) {
                for (const V slope : {V(-1), V(1)}) {
                    for (const V guess : {V(-1), V(0), V(0.5), V(1), V(2), V(3)}) {
                        int value_calls = 0;
                        int full_calls = 0;
                        const auto f_value = [&](V x) {
                            ++value_calls;
                            return Residual<V> {.value = slope * (x - endpoint),
                                                .derivative = slope};
                        };
                        const auto plain_full = [&](V x) {
                            ++full_calls;
                            return Residual<V> {.value = slope * (x - endpoint),
                                                .derivative = slope};
                        };
                        const auto plain = solve<V>(f_value, plain_full, bracket, guess, limits);
                        expect(plain.status == Status::Valid);
                        expect(plain.root == endpoint);
                        expect(value_calls == 2_i);
                        expect(full_calls == 0_i);
                        value_calls = 0;
                        const auto f_full = [&](Dual x) {
                            ++full_calls;
                            return Residual<Dual> {.value = slope * (x - Dual(endpoint, V(1))),
                                                   .derivative = Dual(slope)};
                        };
                        const auto dual = solve<Dual>(f_value, f_full, bracket, guess, limits);
                        expect(dual.status == Status::Valid);
                        expect(dual.root.real() == endpoint);
                        expect(dual.root.dual() == V(1));
                        expect(value_calls == 2_i);
                        expect(full_calls == 1_i);
                    }
                }
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "two exact endpoint roots select the nearest endpoint with a lower tie"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            const Bracket<V> bracket {.lower = V(0), .upper = V(1)};
            const SolveSettings<V> limits {.relative_tolerance = V(0.1), .max_iterations = 0};
            for (const V guess : {V(-1), V(0), V(0.25), V(0.5), V(0.75), V(1), V(2)}) {
                int value_calls = 0;
                int full_calls = 0;
                const auto f_value = [&](V x) {
                    ++value_calls;
                    return Residual<V> {.value = x * (x - V(1)), .derivative = V(2) * x - V(1)};
                };
                const auto plain = solve<V>(f_value, f_value, bracket, guess, limits);
                const V expected = guess <= V(0.5) ? bracket.lower : bracket.upper;
                expect(plain.status == Status::Valid);
                expect(plain.root == expected);
                expect(value_calls == 2_i);
                value_calls = 0;
                const auto f_full = [&](Dual x) {
                    ++full_calls;
                    return Residual<Dual> {.value = x * (x - V(1)) - Dual(V(0), V(1)),
                                           .derivative = V(2) * x - V(1)};
                };
                const auto dual = solve<Dual>(f_value, f_full, bracket, guess, limits);
                expect(dual.status == Status::Valid);
                expect(dual.root.real() == expected);
                expect(dual.root.dual() == V(1) / (V(2) * expected - V(1)));
                expect(value_calls == 2_i);
                expect(full_calls == 1_i);
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
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

    "bisection locates a root with a zero or nonfinite derivative"_test = [] {
        const auto check = []<std::floating_point V>() {
            const SolveSettings<V> limits {
              .relative_tolerance = static_cast<V>(settings.relative_tolerance),
              .max_iterations = settings.max_iterations,
            };
            for (const V derivative : {V(0),
                                       std::numeric_limits<V>::quiet_NaN(),
                                       std::numeric_limits<V>::infinity(),
                                       -std::numeric_limits<V>::infinity()}) {
                const auto f = [derivative](V x) {
                    return Residual<V> {.value = x - V(1), .derivative = derivative};
                };
                const auto r =
                  solve<V>(f, f, Bracket<V> {.lower = V(0), .upper = V(2)}, V(0.5), limits);
                expect(r.status == Status::Valid);
                expect(r.root >= V(0) && r.root <= V(2));
                expect(std::abs(f(r.root).value) <= limits.relative_tolerance * std::abs(r.root));

                const auto first =
                  solve<V>(f,
                           f,
                           Bracket<V> {.lower = V(0), .upper = V(2)},
                           V(0.5),
                           SolveSettings<V> {.relative_tolerance = limits.relative_tolerance,
                                             .max_iterations = 1});
                expect(first.status == Status::NotConverged);
                expect(first.root == V(1.25));
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
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

    "an unconverged dual solve returns the constant last bounded iterate"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            const auto f_value = [](V x) { return cubic(x, V(10)); };
            for (const int iterations : {0, 1}) {
                const SolveSettings<V> limits {
                  .relative_tolerance = static_cast<V>(settings.relative_tolerance),
                  .max_iterations = iterations,
                };
                for (const V guess : {V(-1), V(1), V(9)}) {
                    int full_calls = 0;
                    const auto f_full = [&](Dual x) {
                        ++full_calls;
                        return cubic(x, Dual(V(10), V(1)));
                    };
                    const Bracket<V> bracket {.lower = V(0), .upper = V(5)};
                    const auto plain = solve<V>(f_value, f_value, bracket, guess, limits);
                    const auto dual = solve<Dual>(f_value, f_full, bracket, guess, limits);
                    expect(plain.status == Status::NotConverged);
                    expect(dual.status == Status::NotConverged);
                    expect(dual.root.real() == plain.root);
                    expect(dual.root.real() >= bracket.lower && dual.root.real() <= bracket.upper);
                    expect(dual.root.dual() == V(0));
                    expect(full_calls == 0_i);
                }
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
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
