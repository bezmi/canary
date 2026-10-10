#ifndef GASLIB_NASA9_HPP
#define GASLIB_NASA9_HPP

// NASA 9-coefficient polynomial fits of the ideal-gas heat capacity,
// enthalpy and entropy of one species (McBride, Zehe and Gordon,
// NASA/TP-2002-211556). These are empirical fits: they compute on plain
// scalars, in the temperature unit of the model that holds the curve.

#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>
#include <gaslib/build.hpp>
#include <gaslib/config.hpp>
#include <gaslib/scalar.hpp>

namespace gaslib
{

// Fit of one species in SI units, temperatures in kelvin. Segment j
// (0 <= j < segments) covers [breaks[j], breaks[j + 1]] with coefficients
// a[j]; with t in K,
//   cp/R = a0 t^-2 + a1 t^-1 + a2 + a3 t + a4 t^2 + a5 t^3 + a6 t^4
//   h/R  = -a0 t^-1 + a1 ln t + a2 t + a3 t^2/2 + a4 t^3/3 + a5 t^4/4
//          + a6 t^5/5 + a7
//   s/R  = -a0 t^-2/2 - a1 t^-1 + a2 ln t + a3 t + a4 t^2/2 + a5 t^3/3
//          + a6 t^4/4 + a8
// h includes the enthalpy of formation, so h/R is in kelvin; s is at the
// standard pressure of the fit. The interior boundary breaks[i]
// (1 <= i < segments) has a blending region of width blends[i - 1] centred
// on it. Entries past the used segments are ignored.
template <std::size_t S>
    requires(S > 0)
struct Nasa9Data {
    int segments;
    std::array<double, S + 1> breaks;
    std::array<double, S> blends;
    std::array<std::array<double, 9>, S> a;
};

// cp/R, h/R and s/R of a curve at one temperature.
template <Scalar T> struct Nasa9Values {
    T cp_over_r;
    T h_over_r;
    T s_over_r;
};

// cp/R and h/R of a curve at one temperature, from one coefficient
// selection.
template <Scalar T> struct Nasa9CpH {
    T cp_over_r;
    T h_over_r;
};

// A fit in the temperature unit of a model, with temperature scale tau
// (1 in SI): with t' = t / tau, the polynomial forms of Nasa9Data with
// coefficients b give cp/R, h/R in model temperature units, and s/R. The
// breaks and blending widths are in model temperature units. low and high
// hold the values at breaks[0] and breaks[segments], which continue the
// curve beyond its ends (see nasa9_cp_h).
template <std::floating_point P, std::size_t S>
    requires(S > 0)
struct Nasa9Curve {
    int segments;
    std::array<P, S + 1> breaks;
    std::array<P, S> blends;
    std::array<std::array<P, 9>, S> b;
    Nasa9Values<P> low;
    Nasa9Values<P> high;
};

namespace detail
{

// Coefficients of the curve at the value part t:
// - inside the closed blending region [b_i - w_i / 2, b_i + w_i / 2] of an
//   interior boundary b_i: lambda * (segment above) + (1 - lambda) *
//   (segment below), with lambda = (t - (b_i - w_i / 2)) / w_i;
// - otherwise those of the segment that contains t.
// lambda is computed in P from the value part only, so a dual-number
// temperature treats it as a constant. Outside blending regions lambda is 0 and the
// same formula returns the segment's coefficients exactly.
template <std::floating_point P, std::size_t S, std::floating_point V>
GASLIB_HOST_DEVICE constexpr std::array<P, 9> nasa9_coefficients(const Nasa9Curve<P, S>& curve, V t)
{
    const auto tc = static_cast<P>(t);
    int below = 0;
    int above = 0;
    P lambda = P(0);
    for (int i = 1; i < curve.segments; ++i) {
        const auto boundary = static_cast<std::size_t>(i);
        const P width = curve.blends[boundary - 1];
        const P start = curve.breaks[boundary] - (width / P(2));
        const P end = curve.breaks[boundary] + (width / P(2));
        if (tc >= start && tc <= end) {
            below = i - 1;
            above = i;
            lambda = (tc - start) / width;
        } else if (tc > end) {
            below = i;
            above = i;
        }
    }
    const auto& lower = curve.b[static_cast<std::size_t>(below)];
    const auto& upper = curve.b[static_cast<std::size_t>(above)];
    std::array<P, 9> k {};
    for (std::size_t j = 0; j < 9; ++j) {
        k[j] = (lambda * upper[j]) + ((P(1) - lambda) * lower[j]);
    }
    return k;
}

}  // namespace detail

// cp/R and h/R of the curve at temperature t (model units), with ln_t the
// natural logarithm of t, which the caller computes once for all species.
// Inside [breaks[0], breaks[segments]] the polynomials are evaluated in
// Horner form with the coefficients of detail::nasa9_coefficients. Beyond
// an end t_end, cp keeps its value at t_end and h continues as
// h(t_end) + cp(t_end) (t - t_end).
template <Scalar T, std::floating_point P, std::size_t S>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE Nasa9CpH<T> nasa9_cp_h(const Nasa9Curve<P, S>& curve, T t, T ln_t)
{
    using V = ScalarTraits<T>::value_type;
    const V tv = ScalarTraits<T>::value(t);
    const auto last = static_cast<std::size_t>(curve.segments);
    if (tv < curve.breaks[0] || tv > curve.breaks[last]) {
        const bool below = tv < curve.breaks[0];
        const Nasa9Values<P>& end = below ? curve.low : curve.high;
        const P t_end = below ? curve.breaks[0] : curve.breaks[last];
        return {
          .cp_over_r = T(static_cast<V>(end.cp_over_r)),
          .h_over_r = T(end.h_over_r + (end.cp_over_r * (t - t_end))),
        };
    }
    const std::array<P, 9> k = detail::nasa9_coefficients(curve, tv);
    const T inv_t = T(V(1)) / t;
    const auto cp = (((k[0] * inv_t) + k[1]) * inv_t) +
                    (k[2] + (t * (k[3] + (t * (k[4] + (t * (k[5] + (t * k[6]))))))));
    const auto h =
      (-(k[0] * inv_t)) + (k[1] * ln_t) + k[7] +
      (t * (k[2] + (t * ((k[3] / P(2)) +
                         (t * ((k[4] / P(3)) + (t * ((k[5] / P(4)) + (t * (k[6] / P(5)))))))))));
    return {.cp_over_r = T(cp), .h_over_r = T(h)};
}

// s/R of the curve at temperature t (model units) and the standard pressure
// of the fit, with ln_t = ln t. Beyond an end t_end, s continues as
// s(t_end) + cp(t_end) ln(t / t_end).
template <Scalar T, std::floating_point P, std::size_t S>
    requires ScalarWith<T, P>
GASLIB_HOST_DEVICE T nasa9_s(const Nasa9Curve<P, S>& curve, T t, T ln_t)
{
    using V = ScalarTraits<T>::value_type;
    const V tv = ScalarTraits<T>::value(t);
    const auto last = static_cast<std::size_t>(curve.segments);
    if (tv < curve.breaks[0] || tv > curve.breaks[last]) {
        const bool below = tv < curve.breaks[0];
        const Nasa9Values<P>& end = below ? curve.low : curve.high;
        const P ln_t_end = std::log(below ? curve.breaks[0] : curve.breaks[last]);
        return T(end.s_over_r + (end.cp_over_r * (ln_t - ln_t_end)));
    }
    const std::array<P, 9> k = detail::nasa9_coefficients(curve, tv);
    const T inv_t = T(V(1)) / t;
    const auto s =
      (((-(k[0] / P(2)) * inv_t) - k[1]) * inv_t) + (k[2] * ln_t) + k[8] +
      (t * (k[3] + (t * ((k[4] / P(2)) + (t * ((k[5] / P(3)) + (t * (k[6] / P(4)))))))));
    return T(s);
}

namespace detail
{

// Fit checks of make_nasa9_curve; throws ParameterError.
template <std::size_t S> void check_nasa9_data(const Nasa9Data<S>& data)
{
    if (data.segments < 1 || static_cast<std::size_t>(data.segments) > S) {
        throw ParameterError(
          std::format("segment count must be between 1 and {}, got {}", S, data.segments));
    }
    const auto n = static_cast<std::size_t>(data.segments);
    for (std::size_t i = 0; i <= n; ++i) {
        if (!std::isfinite(data.breaks[i]) || !(data.breaks[i] > 0.0) ||
            (i > 0 && !(data.breaks[i] > data.breaks[i - 1]))) {
            throw ParameterError(std::format(
              "segment boundaries must be finite, positive and increasing, got {} at index {}",
              data.breaks[i],
              i));
        }
    }
    // Blending regions lie inside the fit and do not overlap; adjacent
    // regions may touch.
    double previous_end = data.breaks[0];
    for (std::size_t i = 1; i < n; ++i) {
        const double width = data.blends[i - 1];
        if (!std::isfinite(width) || !(width > 0.0)) {
            throw ParameterError(
              std::format("blending width must be finite and positive, got {} at boundary {}",
                          width,
                          i));
        }
        const double start = data.breaks[i] - (0.5 * width);
        if (start < previous_end) {
            throw ParameterError(
              std::format("blending region of boundary {} starts at {}, before {}",
                          i,
                          start,
                          previous_end));
        }
        previous_end = data.breaks[i] + (0.5 * width);
    }
    if (previous_end > data.breaks[n]) {
        throw ParameterError(std::format("blending region ends at {}, after the fit end {}",
                                         previous_end,
                                         data.breaks[n]));
    }
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t k = 0; k < 9; ++k) {
            if (!std::isfinite(data.a[j][k])) {
                throw ParameterError(
                  std::format("coefficient {} of segment {} is not finite", k, j));
            }
        }
    }
}

}  // namespace detail

// Curve of a fit in a model's temperature unit, with temperature_scale the
// SI value of that unit in kelvin (1 in SI). The coefficients are folded in
// double precision and then rounded to P:
//   b_k = a_k tau^(k - 2) for k = 0..6, b_7 = (a_7 + a_1 ln tau) / tau,
//   b_8 = a_8 + a_2 ln tau,
// and the breaks and blending widths are divided by tau. With tau = 1 the
// coefficients equal a. Throws ParameterError when the scale is not finite
// and positive, the segment count is not in [1, S], the used breaks are not
// finite, positive and increasing, a blending width is not finite and
// positive, the blending regions overlap or leave the fit, or a used
// coefficient is not finite.
template <std::size_t S, std::floating_point P>
Nasa9Curve<P, S> make_nasa9_curve(const Nasa9Data<S>& data, P temperature_scale)
{
    if (!std::isfinite(temperature_scale) || !(temperature_scale > P(0))) {
        throw ParameterError(
          std::format("temperature scale must be finite and positive, got {}", temperature_scale));
    }
    detail::check_nasa9_data(data);
    const auto tau = static_cast<double>(temperature_scale);
    const double ln_tau = std::log(tau);
    const auto n = static_cast<std::size_t>(data.segments);
    Nasa9Curve<P, S> curve {};
    curve.segments = data.segments;
    for (std::size_t i = 0; i <= n; ++i) { curve.breaks[i] = static_cast<P>(data.breaks[i] / tau); }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        curve.blends[i] = static_cast<P>(data.blends[i] / tau);
    }
    for (std::size_t j = 0; j < n; ++j) {
        const auto& a = data.a[j];
        for (std::size_t k = 0; k < 7; ++k) {
            curve.b[j][k] = static_cast<P>(a[k] * std::pow(tau, static_cast<double>(k) - 2.0));
        }
        curve.b[j][7] = static_cast<P>((a[7] + (a[1] * ln_tau)) / tau);
        curve.b[j][8] = static_cast<P>(a[8] + (a[2] * ln_tau));
    }
    // The ends lie inside [breaks[0], breaks[segments]], so these calls
    // evaluate the polynomials and do not read low or high.
    const auto values_at = [&curve](P t) {
        const P ln_t = std::log(t);
        const Nasa9CpH<P> cp_h = nasa9_cp_h(curve, t, ln_t);
        return Nasa9Values<P> {
          .cp_over_r = cp_h.cp_over_r,
          .h_over_r = cp_h.h_over_r,
          .s_over_r = nasa9_s(curve, t, ln_t),
        };
    };
    curve.low = values_at(curve.breaks[0]);
    curve.high = values_at(curve.breaks[n]);
    return curve;
}

}  // namespace gaslib

#endif  // GASLIB_NASA9_HPP
