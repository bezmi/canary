#ifndef CANARY_TESTS_SCALAR_TYPES_HPP
#define CANARY_TESTS_SCALAR_TYPES_HPP

// Minimal scalar types for concept tests.

#include <cmath>
#include <canary/scalar.hpp>

namespace canary::test
{

// Arithmetic, ScalarTraits and sqrt, but no pow, exp, log, abs, min or max.
struct OnlySqrt {
    double v;
    friend OnlySqrt operator-(OnlySqrt a) { return {-a.v}; }
    friend OnlySqrt operator+(OnlySqrt a, OnlySqrt b) { return {a.v + b.v}; }
    friend OnlySqrt operator-(OnlySqrt a, OnlySqrt b) { return {a.v - b.v}; }
    friend OnlySqrt operator*(OnlySqrt a, OnlySqrt b) { return {a.v * b.v}; }
    friend OnlySqrt operator/(OnlySqrt a, OnlySqrt b) { return {a.v / b.v}; }
    friend bool operator<(OnlySqrt a, OnlySqrt b) { return a.v < b.v; }
    friend OnlySqrt sqrt(OnlySqrt a) { return {std::sqrt(a.v)}; }
};

// Arithmetic but no ScalarTraits specialisation.
struct NoTraits {
    double v;
    friend NoTraits operator-(NoTraits a) { return {-a.v}; }
    friend NoTraits operator+(NoTraits a, NoTraits b) { return {a.v + b.v}; }
    friend NoTraits operator-(NoTraits a, NoTraits b) { return {a.v - b.v}; }
    friend NoTraits operator*(NoTraits a, NoTraits b) { return {a.v * b.v}; }
    friend NoTraits operator/(NoTraits a, NoTraits b) { return {a.v / b.v}; }
    friend bool operator<(NoTraits a, NoTraits b) { return a.v < b.v; }
};

}  // namespace canary::test

template <> struct canary::ScalarTraits<canary::test::OnlySqrt> {
    using value_type = double;

    static constexpr double value(canary::test::OnlySqrt x) { return x.v; }
};

#endif  // CANARY_TESTS_SCALAR_TYPES_HPP
