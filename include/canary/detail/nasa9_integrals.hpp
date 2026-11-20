// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_DETAIL_NASA9_INTEGRALS_HPP
#define CANARY_DETAIL_NASA9_INTEGRALS_HPP

#include <array>
#include <canary/config.hpp>
#include <canary/scalar.hpp>
#include <cmath>
#include <limits>
#include <utility>

namespace canary::detail
{

// Integral from 0 to 1 of v^N (s+y*v)^Power, where s=1-y. Only the
// powers needed by the NASA9 smoothstep primitives are instantiated.
template <int N, int Power, Scalar T>
CANARY_HOST_DEVICE T nasa9_moment(T y, T s, T a3, T a4) noexcept
{
    using V = ScalarTraits<T>::value_type;
    if constexpr (Power >= 0) {
        constexpr auto coefficients = [] {
            std::array<V, Power + 1> out {};
            int numerator = 1;
            int denominator = N + 1;
            for (int k = 0; k <= Power; ++k) {
                out[k] = V(numerator) / V(denominator);
                numerator *= -(Power - k);
                denominator *= N + k + 2;
            }
            return out;
        }();
        T value(coefficients[Power]);
        for (int k = Power - 1; k >= 0; --k) { value = T(value * y + coefficients[k]); }
        return value;
    } else if constexpr (N == 2 && Power == -1) {
        return T(y / V(2) + s * s * a3);
    } else if constexpr (N == 3 && Power == -1) {
        return T(s / V(2) + y * y / V(3) - s * s * s * a4);
    } else if constexpr (N == 2 && Power == -2) {
        return T(T(V(1)) - V(2) * s * a3);
    } else if constexpr (N == 3 && Power == -2) {
        return T(y - V(0.5) + V(3) * s * s * a4);
    } else if constexpr (N == 2 && Power == -3) {
        return a3;
    } else {
        static_assert(N == 3 && Power == -3);
        return T(T(V(1)) - V(3) * s * a4);
    }
}

template <bool Entropy, Scalar T, int... J>
CANARY_HOST_DEVICE std::array<T, 7>
  nasa9_moment_weights(T u, T y, T s, T a3, T a4, std::integer_sequence<int, J...>) noexcept
{
    using V = ScalarTraits<T>::value_type;
    return {T(V(3) * nasa9_moment<2, J - (Entropy ? 3 : 2)>(y, s, a3, a4) -
              V(2) * u * nasa9_moment<3, J - (Entropy ? 3 : 2)>(y, s, a3, a4))...};
}

// Integral from start to t of w*(cp_upper-cp_lower), divided by temperature
// for entropy. With d=t-start, u=d/(end-start), y=d/t and s=start/t,
// each power t^p contributes d*u²*t^p*(3*M(2,p)-2*u*M(3,p)) to h.
// Entropy replaces d by y and p by p-1 in the moments only. This local
// representation avoids powers of the reciprocal transition width.
// inv must equal 1 / t.
template <bool Entropy, HasLog T, std::floating_point P>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE T nasa9_integral(
  const std::array<P, 9>& lower, const std::array<P, 9>& upper, T t, T inv, P start, P end) noexcept
{
    using V = ScalarTraits<T>::value_type;
    const T d = T(t - start);
    const T u = T(d / (end - start));
    const T y = T(d * inv);
    const T s = T(start * inv);
    T a4(V(0));
    const V y_value = ScalarTraits<T>::value(y);
    if (y_value <= V(0.5)) {
        // A4=sum(y^k/(k+4)) to k=N. The omitted tail is bounded by
        // y^(N+1)/((N+5)*(1-y)). With y <= 2^-m, m*N >= digits+8 bounds it
        // and its derivative below working precision, so N follows the
        // binary exponent of y: N=digits+8 for y in [1/4, 1/2], fewer terms
        // for smaller y, and the first terms only for y = 0. N+1 is rounded
        // up to a multiple of four. The terms are positive: those from k=4
        // are summed in four independent chains in y^4 (k mod 4 = 0..3), and
        // the first four by Horner's rule, which keeps the rounding of y^2
        // and y^4 out of the leading terms.
        constexpr int bits = std::numeric_limits<V>::digits + 8;
        constexpr std::size_t groups_max = (bits + 4) / 4;
        constexpr auto coefficients = [] {
            std::array<V, 4 * groups_max> out {};
            for (std::size_t k = 0; k < out.size(); ++k) { out[k] = V(1) / V(k + 4); }
            return out;
        }();
        std::size_t groups = groups_max;
        if (y_value > V(0)) {
            const int m = -std::ilogb(y_value) - 1;
            if (m > 1) { groups = (static_cast<std::size_t>((bits + m - 1) / m) + 4) / 4; }
        } else {
            groups = 1;
        }
        const T y2 = y * y;
        const T y4 = y2 * y2;
        std::array<T, 4> chain {T(V(0)), T(V(0)), T(V(0)), T(V(0))};
        for (std::size_t j = groups; j > 1; --j) {
            for (std::size_t r = 0; r < 4; ++r) {
                chain[r] = T(chain[r] * y4 + coefficients[(4 * (j - 1)) + r]);
            }
        }
        a4 = T((chain[0] + y * chain[1]) + y2 * (chain[2] + y * chain[3]));
        for (std::size_t k = 4; k > 0; --k) { a4 = T(a4 * y + coefficients[k - 1]); }
    } else {
        const T y2 = y * y;
        a4 = T((math::log(T(t / start)) - y - y2 / V(2) - y2 * y / V(3)) / (y2 * y2));
    }
    const T a3 = T(T(V(1) / V(3)) + y * a4);
    auto weighted =
      nasa9_moment_weights<Entropy>(u, y, s, a3, a4, std::make_integer_sequence<int, 7> {});
    for (std::size_t j = 0; j < 7; ++j) {
        weighted[j] = T((T(static_cast<V>(upper[j])) - T(static_cast<V>(lower[j]))) * weighted[j]);
    }
    const T sum = (weighted[0] * inv + weighted[1]) * inv + weighted[2] +
                  t * (weighted[3] + t * (weighted[4] + t * (weighted[5] + t * weighted[6])));
    return (Entropy ? y : d) * u * u * sum;
}

template <bool Entropy, HasLog T, std::floating_point P>
    requires ScalarWith<T, P>
CANARY_HOST_DEVICE T nasa9_integral(
  const std::array<P, 9>& lower, const std::array<P, 9>& upper, T t, P start, P end) noexcept
{
    using V = ScalarTraits<T>::value_type;
    return nasa9_integral<Entropy>(lower, upper, t, T(T(V(1)) / t), start, end);
}

}  // namespace canary::detail

#endif  // CANARY_DETAIL_NASA9_INTEGRALS_HPP
