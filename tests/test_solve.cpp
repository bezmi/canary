// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#include "dual.hpp"
#include "scalar_types.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <canary/detail/solve.hpp>
#include <canary/status.hpp>
#include <cmath>
#include <concepts>
#include <limits>

namespace
{

using canary::Status;
using canary::detail::Bracket;
using canary::detail::Residual;
using canary::detail::SolveReason;
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
        expect(r.reason == SolveReason::None);
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
        expect(r.reason == SolveReason::None);
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
            // A loose bracket tolerance leaves a nonzero residual so that a
            // separate final Newton step would change the returned value.
            const SolveSettings<V> limits {.relative_tolerance = V(0.1), .max_iterations = 50};
            const Bracket<V> bracket {.lower = V(1.875), .upper = V(2.25)};
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
                expect(r.reason == SolveReason::None);
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
            const auto check_invalid = [&](V value, V derivative, SolveReason reason) {
                const auto f_full = [=](Dual) {
                    return Residual<Dual> {.value = Dual(value, V(1)),
                                           .derivative = Dual(derivative)};
                };
                const auto r = solve<Dual>(f_value, f_full, bracket, V(0.5), limits);
                expect(r.status == Status::Invalid);
                expect(r.reason == reason);
                expect(r.root.real() == plain.root);
                expect(r.root.real() >= bracket.lower && r.root.real() <= bracket.upper);
                expect(r.root.dual() == V(0));
            };
            check_invalid(V(0), V(0), SolveReason::ZeroFullDerivative);
            for (const V value : {std::numeric_limits<V>::quiet_NaN(),
                                  std::numeric_limits<V>::infinity(),
                                  -std::numeric_limits<V>::infinity()}) {
                check_invalid(value, V(1), SolveReason::NonfiniteFullResidual);
                check_invalid(V(0), value, SolveReason::NonfiniteFullDerivative);
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
            expect(r.reason == SolveReason::ZeroFullDerivative);
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

    "equal finite bounds accept only an exact root"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            const SolveSettings<V> limits {.relative_tolerance = V(0.1), .max_iterations = 0};
            for (const V endpoint : {V(0),
                                     V(-1),
                                     V(1),
                                     -std::numeric_limits<V>::max(),
                                     std::numeric_limits<V>::max()}) {
                const Bracket<V> bracket {.lower = endpoint, .upper = endpoint};
                for (const V slope : {V(-1), V(0), V(1)}) {
                    const auto f_value = [=](V x) {
                        return Residual<V> {.value = slope * (x - endpoint), .derivative = slope};
                    };
                    const auto f_full = [=](Dual x) {
                        return Residual<Dual> {.value = slope * (x - Dual(endpoint, V(1))),
                                               .derivative = Dual(slope)};
                    };
                    for (const V guess : {V(-2), V(0), V(2)}) {
                        const auto plain = solve<V>(f_value, f_value, bracket, guess, limits);
                        expect(plain.status == Status::Valid);
                        expect(plain.root == endpoint);
                        const auto dual = solve<Dual>(f_value, f_full, bracket, guess, limits);
                        expect(dual.status == (slope == V(0) ? Status::Invalid : Status::Valid));
                        expect(dual.root.real() == endpoint);
                        expect(dual.root.dual() == (slope == V(0) ? V(0) : V(1)));
                    }
                }
                for (const V residual : {V(-1),
                                         V(1),
                                         std::numeric_limits<V>::quiet_NaN(),
                                         -std::numeric_limits<V>::infinity(),
                                         std::numeric_limits<V>::infinity()}) {
                    int full_calls = 0;
                    const auto f_value = [=](V) {
                        return Residual<V> {.value = residual, .derivative = V(1)};
                    };
                    const auto f_full = [&](Dual x) {
                        ++full_calls;
                        return Residual<Dual> {.value = x, .derivative = Dual(V(1))};
                    };
                    const auto plain = solve<V>(f_value, f_value, bracket, endpoint, limits);
                    const auto dual = solve<Dual>(f_value, f_full, bracket, endpoint, limits);
                    expect(plain.status == Status::Invalid);
                    expect(dual.status == Status::Invalid);
                    expect(plain.root == endpoint);
                    expect(dual.root.real() == endpoint);
                    expect(dual.root.dual() == V(0));
                    expect(full_calls == 0_i);
                }
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
        expect(r.reason == SolveReason::None);
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
                expect(r.reason == SolveReason::None);
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
            for (const std::size_t iterations : {std::size_t {0}, std::size_t {1}}) {
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

    "nonfinite interior residuals stop at a constant bounded iterate"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            const Bracket<V> bracket {.lower = V(0), .upper = V(2)};
            const SolveSettings<V> limits {.relative_tolerance = V(0.01), .max_iterations = 50};
            // The invalid point is the initial guess, a Newton candidate, or a
            // bisection candidate when the derivative is zero.
            for (const V invalid_point : {V(0.5), V(1), V(1.25)}) {
                for (const V slope : {V(-1), V(1)}) {
                    for (const V invalid_value : {std::numeric_limits<V>::quiet_NaN(),
                                                  std::numeric_limits<V>::infinity(),
                                                  -std::numeric_limits<V>::infinity()}) {
                        int value_calls = 0;
                        int full_calls = 0;
                        const auto f_value = [&](V x) {
                            ++value_calls;
                            return Residual<V> {
                              .value = x == invalid_point ? invalid_value : slope * (x - V(1)),
                              .derivative = invalid_point == V(1.25) ? V(0) : slope,
                            };
                        };
                        const auto f_full = [&](Dual x) {
                            ++full_calls;
                            return Residual<Dual> {.value = x - Dual(V(1), V(1)),
                                                   .derivative = Dual(V(1))};
                        };
                        const int expected_calls = invalid_point == V(0.5) ? 3 : 4;
                        const auto plain = solve<V>(f_value, f_value, bracket, V(0.5), limits);
                        expect(plain.status == Status::Invalid);
                        expect(plain.reason == SolveReason::NonfiniteInteriorResidual);
                        expect(plain.root == invalid_point);
                        expect(plain.root >= bracket.lower && plain.root <= bracket.upper);
                        expect(value_calls == expected_calls);
                        value_calls = 0;
                        const auto dual = solve<Dual>(f_value, f_full, bracket, V(0.5), limits);
                        expect(dual.status == Status::Invalid);
                        expect(dual.reason == SolveReason::NonfiniteInteriorResidual);
                        expect(dual.root.real() == invalid_point);
                        expect(dual.root.real() >= bracket.lower &&
                               dual.root.real() <= bracket.upper);
                        expect(dual.root.dual() == V(0));
                        expect(value_calls == expected_calls);
                        expect(full_calls == 0_i);
                    }
                }
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "invalid solve inputs are rejected before residual evaluation"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            const Bracket<V> bracket {.lower = V(0), .upper = V(2)};
            const SolveSettings<V> limits {.relative_tolerance = V(0.01), .max_iterations = 50};
            const auto check_invalid =
              [&](Bracket<V> bounds, V guess, SolveSettings<V> options, SolveReason reason) {
                  int value_calls = 0;
                  int full_calls = 0;
                  const auto f_value = [&](V x) {
                      ++value_calls;
                      return Residual<V> {.value = x, .derivative = V(1)};
                  };
                  const auto f_full = [&](Dual x) {
                      ++full_calls;
                      return Residual<Dual> {.value = x, .derivative = Dual(V(1))};
                  };
                  const auto plain = solve<V>(f_value, f_value, bounds, guess, options);
                  const auto dual = solve<Dual>(f_value, f_full, bounds, guess, options);
                  expect(plain.status == Status::Invalid);
                  expect(dual.status == Status::Invalid);
                  expect(plain.reason == reason);
                  expect(dual.reason == reason);
                  expect(std::isnan(guess) ? std::isnan(plain.root) : plain.root == guess);
                  expect(std::isnan(guess) ? std::isnan(dual.root.real())
                                           : dual.root.real() == guess);
                  expect(dual.root.dual() == V(0));
                  expect(value_calls == 0_i);
                  expect(full_calls == 0_i);
              };
            check_invalid({.lower = V(2), .upper = V(0)},
                          V(1),
                          limits,
                          SolveReason::ReversedBounds);
            check_invalid(bracket,
                          V(1),
                          {.relative_tolerance = V(-1), .max_iterations = 50},
                          SolveReason::InvalidRelativeTolerance);
            check_invalid(
              bracket,
              V(1),
              {.relative_tolerance = V(0.01), .max_iterations = 50, .absolute_tolerance = V(-1)},
              SolveReason::InvalidAbsoluteTolerance);
            for (const V invalid : {std::numeric_limits<V>::quiet_NaN(),
                                    -std::numeric_limits<V>::infinity(),
                                    std::numeric_limits<V>::infinity()}) {
                check_invalid({.lower = invalid, .upper = V(2)},
                              V(1),
                              limits,
                              SolveReason::NonfiniteBounds);
                check_invalid({.lower = V(0), .upper = invalid},
                              V(1),
                              limits,
                              SolveReason::NonfiniteBounds);
                check_invalid(bracket, invalid, limits, SolveReason::NonfiniteGuess);
                check_invalid(bracket,
                              V(1),
                              {.relative_tolerance = invalid, .max_iterations = 50},
                              SolveReason::InvalidRelativeTolerance);
                check_invalid(bracket,
                              V(1),
                              {.relative_tolerance = V(0.01),
                               .max_iterations = 50,
                               .absolute_tolerance = invalid},
                              SolveReason::InvalidAbsoluteTolerance);
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "zero tolerance accepts exact roots and rejects exhausted precision"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            const SolveSettings<V> limits {.relative_tolerance = V(0), .max_iterations = 128};
            const auto linear_value = [](V x) {
                return Residual<V> {.value = x - V(1), .derivative = V(1)};
            };
            const auto linear_full = [](Dual x) {
                return Residual<Dual> {.value = x - Dual(V(1), V(1)), .derivative = Dual(V(1))};
            };
            for (const Bracket<V> bracket : {Bracket<V> {.lower = V(0), .upper = V(2)},
                                             Bracket<V> {.lower = V(0), .upper = V(1)},
                                             Bracket<V> {.lower = V(1), .upper = V(1)}}) {
                const auto plain = solve<V>(linear_value, linear_value, bracket, V(0.5), limits);
                const auto dual = solve<Dual>(linear_value, linear_full, bracket, V(0.5), limits);
                expect(plain.status == Status::Valid);
                expect(plain.root == V(1));
                expect(dual.status == Status::Valid);
                expect(dual.root.real() == V(1));
                expect(dual.root.dual() == V(1));
            }
            // A fused residual retains the nonzero error at every representable
            // approximation to sqrt(2).
            const auto quadratic_value = [](V x) {
                return Residual<V> {.value = std::fma(x, x, V(-2)), .derivative = V(2) * x};
            };
            int full_calls = 0;
            const auto quadratic_full = [&](Dual x) {
                ++full_calls;
                return Residual<Dual> {.value = x * x - Dual(V(2), V(1)), .derivative = V(2) * x};
            };
            const Bracket<V> bracket {.lower = V(1), .upper = V(2)};
            for (const std::size_t iterations :
                 {std::size_t {0}, std::size_t {1}, limits.max_iterations}) {
                const SolveSettings<V> options {.relative_tolerance = V(0),
                                                .max_iterations = iterations};
                const auto plain =
                  solve<V>(quadratic_value, quadratic_value, bracket, V(1.5), options);
                const auto dual =
                  solve<Dual>(quadratic_value, quadratic_full, bracket, V(1.5), options);
                expect(plain.status == Status::NotConverged);
                expect(plain.root >= bracket.lower && plain.root <= bracket.upper);
                expect(quadratic_value(plain.root).value != V(0));
                expect(dual.status == Status::NotConverged);
                expect(dual.root.real() == plain.root);
                expect(dual.root.dual() == V(0));
                expect(full_calls == 0_i);
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "bracket convergence uses additive tolerances and the minimum absolute endpoint"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            for (const V sign : {V(-1), V(1)}) {
                const Bracket<V> bracket = sign > V(0)
                                             ? Bracket<V> {.lower = V(2), .upper = V(4)}
                                             : Bracket<V> {.lower = V(-4), .upper = V(-2)};
                const auto f_value = [=](V x) {
                    return Residual<V> {.value = x - sign * V(3.25), .derivative = V(0)};
                };
                const auto f_full = [=](Dual x) {
                    return Residual<Dual> {.value = x - Dual(sign * V(3.25), V(1)),
                                           .derivative = Dual(V(1))};
                };
                const V guess = sign * V(4);
                // Neither tolerance alone covers the width of two; their sum does.
                const SolveSettings<V> limits {.relative_tolerance = V(0.5),
                                               .max_iterations = 1,
                                               .absolute_tolerance = V(1)};
                const auto plain = solve<V>(f_value, f_value, bracket, guess, limits);
                const auto dual = solve<Dual>(f_value, f_full, bracket, guess, limits);
                expect(plain.status == Status::Valid);
                expect(plain.root == guess);
                expect(dual.status == Status::Valid);
                expect(dual.root.real() == guess);
                expect(dual.root.dual() == V(1));
                // Using abs(guess) for the relative scale would accept this interval.
                const auto relative_only =
                  solve<V>(f_value,
                           f_value,
                           bracket,
                           guess,
                           SolveSettings<V> {.relative_tolerance = V(0.5), .max_iterations = 1});
                expect(relative_only.status == Status::NotConverged);
            }
            const auto f = [](V x) { return Residual<V> {.value = x, .derivative = V(0)}; };
            for (const V absolute : {V(0), V(1), V(2)}) {
                const auto r = solve<V>(f,
                                        f,
                                        Bracket<V> {.lower = V(-1), .upper = V(1)},
                                        V(-1),
                                        SolveSettings<V> {.relative_tolerance = V(100),
                                                          .max_iterations = 1,
                                                          .absolute_tolerance = absolute});
                expect(r.status == (absolute == V(2) ? Status::Valid : Status::NotConverged));
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "bracket tolerances handle subnormal and overflowing widths"_test = [] {
        const auto check = []<std::floating_point V>() {
            const V tiny = std::numeric_limits<V>::denorm_min();
            const auto near_zero = [=](V x) {
                return Residual<V> {.value = x / tiny - V(3), .derivative = V(0)};
            };
            const auto exact =
              solve<V>(near_zero,
                       near_zero,
                       Bracket<V> {.lower = V(0), .upper = V(8) * tiny},
                       V(0),
                       SolveSettings<V> {.relative_tolerance = V(0), .max_iterations = 8});
            expect(exact.status == Status::Valid);
            expect(exact.root == V(3) * tiny);
            const auto narrow = solve<V>(near_zero,
                                         near_zero,
                                         Bracket<V> {.lower = V(2) * tiny, .upper = V(4) * tiny},
                                         V(2) * tiny,
                                         SolveSettings<V> {.relative_tolerance = V(0),
                                                           .max_iterations = 1,
                                                           .absolute_tolerance = V(2) * tiny});
            expect(narrow.status == Status::Valid);
            expect(narrow.root == V(2) * tiny);

            const V largest = std::numeric_limits<V>::max();
            const auto wide_residual = [=](V x) {
                return Residual<V> {.value = x / largest - V(0.75), .derivative = V(0)};
            };
            const SolveSettings<V> limits {.relative_tolerance = largest,
                                           .max_iterations = 1,
                                           .absolute_tolerance = largest};
            const auto wide = solve<V>(wide_residual,
                                       wide_residual,
                                       Bracket<V> {.lower = -largest, .upper = largest},
                                       -largest,
                                       limits);
            expect(wide.status == Status::NotConverged);
            expect(std::isfinite(wide.root));
            const auto positive = solve<V>(wide_residual,
                                           wide_residual,
                                           Bracket<V> {.lower = largest / V(2), .upper = largest},
                                           largest / V(2),
                                           limits);
            expect(positive.status == Status::Valid);
            expect(positive.root == largest / V(2));
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "Newton stagnation forces bisection even when the bracket width overflows"_test = [] {
        const auto check = []<std::floating_point V>() {
            const V largest = std::numeric_limits<V>::max();
            for (const V extent : {V(6), largest}) {
                std::array<V, 5> points {};
                std::size_t calls = 0;
                const auto f = [&](V x) {
                    points.at(calls++) = x;
                    return Residual<V> {.value = x / extent - V(0.75), .derivative = largest};
                };
                const auto r =
                  solve<V>(f,
                           f,
                           Bracket<V> {.lower = -extent, .upper = extent},
                           extent,
                           SolveSettings<V> {.relative_tolerance = V(0), .max_iterations = 3});
                expect(r.status == Status::NotConverged);
                expect(calls == points.size());
                expect(points[3] == std::nextafter(extent, V(0)));
                expect(points[4] == std::midpoint(-extent, points[3]));
                expect(std::isfinite(r.root));
                expect(r.root >= -extent && r.root <= extent);
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "nonzero tolerances do not accept unresolved steps or exhausted precision"_test = [] {
        const auto check = []<std::floating_point V>() {
            const V eps = std::numeric_limits<V>::epsilon();
            const auto unresolved = [](V x) {
                return Residual<V> {.value = x - V(1), .derivative = std::numeric_limits<V>::max()};
            };
            const auto r =
              solve<V>(unresolved,
                       unresolved,
                       Bracket<V> {.lower = V(0), .upper = V(2)},
                       V(0.5),
                       SolveSettings<V> {.relative_tolerance = eps, .max_iterations = 1});
            expect(r.status == Status::NotConverged);
            expect(r.root >= V(0) && r.root <= V(2));
            const V low = std::sqrt(V(2));
            const V high = std::nextafter(low, V(2));
            const auto between = [=](V x) {
                return Residual<V> {.value = (x - low) / (high - low) - V(0.5), .derivative = V(0)};
            };
            const auto exhausted =
              solve<V>(between,
                       between,
                       Bracket<V> {.lower = low, .upper = high},
                       low,
                       SolveSettings<V> {.relative_tolerance = eps / V(8), .max_iterations = 50});
            expect(exhausted.status == Status::NotConverged);
            expect(exhausted.reason == SolveReason::PrecisionLimit);
            expect(exhausted.root >= low && exhausted.root <= high);
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "solve reasons distinguish failure causes and retain input check order"_test = [] {
        const auto check = []<std::floating_point V>() {
            using Dual = canary::test::Dual<V>;
            const V nan = std::numeric_limits<V>::quiet_NaN();
            const auto f = [](auto x) {
                using X = decltype(x);
                return Residual<X> {.value = x * x - X(V(2)), .derivative = X(V(2)) * x};
            };
            const auto check_reason = [&](auto value,
                                          Bracket<V> bounds,
                                          V guess,
                                          SolveSettings<V> limits,
                                          Status status,
                                          SolveReason reason) {
                const auto plain = solve<V>(value, value, bounds, guess, limits);
                const auto dual = solve<Dual>(value, f, bounds, guess, limits);
                expect(plain.status == status && dual.status == status);
                expect(plain.reason == reason && dual.reason == reason);
                expect(std::isnan(plain.root) ? std::isnan(dual.root.real())
                                              : plain.root == dual.root.real());
                expect(dual.root.dual() == V(0));
            };
            const Bracket<V> bounds {.lower = V(0), .upper = V(2)};
            const SolveSettings<V> limits {.relative_tolerance = V(0), .max_iterations = 0};
            check_reason(f,
                         bounds,
                         V(1),
                         limits,
                         Status::NotConverged,
                         SolveReason::IterationLimit);
            check_reason(f,
                         {.lower = V(2), .upper = V(3)},
                         V(2),
                         limits,
                         Status::Invalid,
                         SolveReason::NoBracket);
            const auto nonfinite = [=](V) {
                return Residual<V> {.value = nan, .derivative = V(1)};
            };
            check_reason(nonfinite,
                         bounds,
                         V(1),
                         limits,
                         Status::Invalid,
                         SolveReason::NonfiniteEndpointResidual);
            const V low = std::sqrt(V(2));
            const V high = std::nextafter(low, V(2));
            const auto between = [=](V x) {
                return Residual<V> {.value = (x - low) / (high - low) - V(0.5), .derivative = V(0)};
            };
            check_reason(between,
                         {.lower = low, .upper = high},
                         low,
                         {.relative_tolerance = V(0), .max_iterations = 50},
                         Status::NotConverged,
                         SolveReason::PrecisionLimit);
            const SolveSettings<V> invalid {.relative_tolerance = V(-1),
                                            .max_iterations = 0,
                                            .absolute_tolerance = V(-1)};
            check_reason(f,
                         {.lower = nan, .upper = V(0)},
                         nan,
                         invalid,
                         Status::Invalid,
                         SolveReason::NonfiniteBounds);
            check_reason(f,
                         {.lower = V(2), .upper = V(0)},
                         nan,
                         invalid,
                         Status::Invalid,
                         SolveReason::ReversedBounds);
            check_reason(f, bounds, nan, invalid, Status::Invalid, SolveReason::NonfiniteGuess);
            check_reason(f,
                         bounds,
                         V(1),
                         invalid,
                         Status::Invalid,
                         SolveReason::InvalidRelativeTolerance);
            check_reason(
              f,
              bounds,
              V(1),
              {.relative_tolerance = V(0), .max_iterations = 0, .absolute_tolerance = V(-1)},
              Status::Invalid,
              SolveReason::InvalidAbsoluteTolerance);
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
