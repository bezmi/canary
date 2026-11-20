// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_DETAIL_NASA9_VALIDATION_HPP
#define CANARY_DETAIL_NASA9_VALIDATION_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace canary::detail
{

// Outward-rounded host bounds used only during curve construction.
struct CapacityInterval {
    long double lower {};
    long double upper {};
    CapacityInterval() = default;
    explicit CapacityInterval(long double value) : lower(value), upper(value) {}
    CapacityInterval(long double low, long double high) : lower(low), upper(high) {}
    bool zero() const { return lower == 0 && upper == 0; }
};

inline CapacityInterval capacity_bounds(long double low, long double high)
{
    constexpr auto infinity = std::numeric_limits<long double>::infinity();
    return {std::nextafter(low, -infinity), std::nextafter(high, infinity)};
}

inline CapacityInterval operator+(CapacityInterval a, CapacityInterval b)
{
    if (a.zero()) { return b; }
    if (b.zero()) { return a; }
    return capacity_bounds(a.lower + b.lower, a.upper + b.upper);
}

inline CapacityInterval operator-(CapacityInterval a, CapacityInterval b)
{
    if (a.lower == a.upper && b.lower == b.upper && a.lower == b.lower) {
        return CapacityInterval {};
    }
    return a + CapacityInterval {-b.upper, -b.lower};
}

inline CapacityInterval operator*(CapacityInterval a, CapacityInterval b)
{
    if (a.zero() || b.zero()) { return CapacityInterval {}; }
    const std::array products {a.lower * b.lower,
                               a.lower * b.upper,
                               a.upper * b.lower,
                               a.upper * b.upper};
    for (const auto x : products) {
        if (!std::isfinite(x)) {
            const auto infinity = std::numeric_limits<long double>::infinity();
            return {-infinity, infinity};
        }
    }
    return capacity_bounds(*std::min_element(products.begin(), products.end()),
                           *std::max_element(products.begin(), products.end()));
}

// Degree nine covers T²(cp/R-1) for a NASA9 segment or cubic smoothstep blend.
using CapacityPolynomial = std::array<CapacityInterval, 10>;

// Power coefficients in u, with T = low + (high-low)u, for T²(cp/R-1).
// Multiplication by T² preserves the sign because all admitted T are positive.
template <typename P>
CapacityPolynomial capacity_polynomial(const std::array<P, 9>& coefficients, P low, P high)
{
    CapacityPolynomial out {};
    CapacityPolynomial power {};
    power[0] = CapacityInterval {1};
    const CapacityInterval left {static_cast<long double>(low)};
    const auto width = CapacityInterval {static_cast<long double>(high)} - left;
    for (std::size_t j = 0; j < 7; ++j) {
        auto coefficient = CapacityInterval {static_cast<long double>(coefficients[j])};
        if (j == 2) { coefficient = coefficient - CapacityInterval {1}; }
        for (std::size_t k = 0; k <= j; ++k) { out[k] = out[k] + coefficient * power[k]; }
        if (j == 6) { break; }
        CapacityPolynomial next {};
        for (std::size_t k = 0; k <= j; ++k) {
            next[k] = next[k] + power[k] * left;
            next[k + 1] = next[k + 1] + power[k] * width;
        }
        power = next;
    }
    return out;
}

inline unsigned capacity_binomial(unsigned n, unsigned k)
{
    unsigned result = 1;
    for (unsigned j = 1; j <= k; ++j) { result = result * (n - j + 1) / j; }
    return result;
}

inline CapacityPolynomial capacity_bernstein(const CapacityPolynomial& power)
{
    CapacityPolynomial out {};
    for (unsigned i = 0; i < out.size(); ++i) {
        for (unsigned k = 0; k <= i; ++k) {
            const long double ratio =
              static_cast<long double>(capacity_binomial(i, k)) / capacity_binomial(9, k);
            out[i] = out[i] + power[k] * capacity_bounds(ratio, ratio);
        }
    }
    return out;
}

// Bernstein coefficients bound a polynomial on [0,1]. Subdivision tightens
// these bounds without missing an interior minimum. Unresolved positivity is
// rejected after a mantissa-sized subdivision budget; it is never inferred
// from a temperature sampling grid.
inline bool positive_capacity(const CapacityPolynomial& coefficients, unsigned remaining)
{
    bool positive = true;
    for (const auto c : coefficients) {
        if (!std::isfinite(c.lower) || !std::isfinite(c.upper)) { return false; }
        positive = positive && c.lower > 0;
    }
    if (positive) { return true; }
    if (coefficients.front().upper <= 0 || coefficients.back().upper <= 0 || remaining == 0) {
        return false;
    }
    auto work = coefficients;
    CapacityPolynomial left {}, right {};
    left.front() = work.front();
    right.back() = work.back();
    for (std::size_t level = 1; level < work.size(); ++level) {
        for (std::size_t j = 0; j < work.size() - level; ++j) {
            work[j] = (work[j] + work[j + 1]) * CapacityInterval {0.5L};
        }
        left[level] = work.front();
        right[work.size() - level - 1] = work[work.size() - level - 1];
    }
    return positive_capacity(left, remaining - 1) && positive_capacity(right, remaining - 1);
}

inline bool positive_capacity(const CapacityPolynomial& power)
{ return positive_capacity(capacity_bernstein(power), std::numeric_limits<long double>::digits); }

}  // namespace canary::detail

#endif  // CANARY_DETAIL_NASA9_VALIDATION_HPP
