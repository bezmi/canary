#ifndef CANARY_NASA9_HPP
#define CANARY_NASA9_HPP

// NASA 9-coefficient polynomial fits of the ideal-gas heat capacity,
// enthalpy and entropy of one species (McBride, Zehe and Gordon,
// NASA/TP-2002-211556). These are empirical fits: they compute on plain
// scalars, in the temperature unit of the model that holds the curve.

#include <algorithm>
#include <array>
#include <canary/build.hpp>
#include <canary/config.hpp>
#include <canary/detail/nasa9_integrals.hpp>
#include <canary/detail/nasa9_validation.hpp>
#include <canary/layout.hpp>
#include <canary/scalar.hpp>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>
#include <limits>

namespace canary
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

    static consteval bool padding_free()
    {
        return layout::
          members_fill<Nasa9Values, decltype(cp_over_r), decltype(h_over_r), decltype(s_over_r)>();
    }
};

// cp/R and h/R of a curve at one temperature, from one coefficient
// selection.
template <Scalar T> struct Nasa9CpH {
    T cp_over_r;
    T h_over_r;
};

namespace detail
{

// Ordinary NASA9 polynomials; integration constants are stored in k[7:9].
// inv must equal 1 / t, so that callers evaluating several polynomials at
// one temperature divide once. The fixed divisors of the integrated
// coefficients are applied as multiplications by their reciprocals.
template <Scalar T, typename P>
CANARY_HOST_DEVICE T nasa9_cp(const std::array<P, 9>& k, T t, T inv) noexcept
{ return T(((k[0] * inv) + k[1]) * inv + k[2] + t * (k[3] + t * (k[4] + t * (k[5] + t * k[6])))); }

template <Scalar T, typename P>
CANARY_HOST_DEVICE T nasa9_h(const std::array<P, 9>& k, T t, T ln_t, T inv) noexcept
{
    using V = ScalarTraits<T>::value_type;
    constexpr V third = V(1) / V(3);
    constexpr V fifth = V(1) / V(5);
    return T(-k[0] * inv + k[1] * ln_t + k[7] +
             t * (k[2] + t * (k[3] * V(0.5) +
                              t * (k[4] * third + t * (k[5] * V(0.25) + t * (k[6] * fifth))))));
}

template <Scalar T, typename P>
CANARY_HOST_DEVICE T nasa9_entropy(const std::array<P, 9>& k, T t, T ln_t, T inv) noexcept
{
    using V = ScalarTraits<T>::value_type;
    constexpr V third = V(1) / V(3);
    return T(((-k[0] * V(0.5)) * inv - k[1]) * inv + k[2] * ln_t + k[8] +
             t * (k[3] + t * (k[4] * V(0.5) + t * (k[5] * third + t * (k[6] * V(0.25))))));
}

template <Scalar T, typename P>
CANARY_HOST_DEVICE T nasa9_cp(const std::array<P, 9>& k, T t) noexcept
{
    using V = ScalarTraits<T>::value_type;
    return nasa9_cp(k, t, T(T(V(1)) / t));
}

template <Scalar T, typename P>
CANARY_HOST_DEVICE T nasa9_h(const std::array<P, 9>& k, T t, T ln_t) noexcept
{
    using V = ScalarTraits<T>::value_type;
    return nasa9_h(k, t, ln_t, T(T(V(1)) / t));
}

template <Scalar T, typename P>
CANARY_HOST_DEVICE T nasa9_entropy(const std::array<P, 9>& k, T t, T ln_t) noexcept
{
    using V = ScalarTraits<T>::value_type;
    return nasa9_entropy(k, t, ln_t, T(T(V(1)) / t));
}

}  // namespace detail

// Validated fixed-size fit in model temperature units. Construction is through
// make_nasa9_curve. Copies retain the validation; evaluation is const and allocates
// no storage dynamically. The heat capacity uses cubic smoothstep transitions.
// Enthalpy and entropy integrate that capacity and retain their 298.15 K anchors.
template <std::floating_point P, std::size_t S>
    requires(S > 0)
class Nasa9Curve
{
    std::size_t segments {};
    std::array<P, S + 1> breaks {};
    std::array<P, S> blends {};
    std::array<std::array<P, 9>, S> b {};
    Nasa9Values<P> low {};
    Nasa9Values<P> high {};

  public:
    static consteval bool padding_free()
    {
        return layout::members_fill<Nasa9Curve,
                                    decltype(segments),
                                    decltype(breaks),
                                    decltype(blends),
                                    decltype(b),
                                    decltype(low),
                                    decltype(high)>();
    }

  private:
    Nasa9Curve() = default;
    void check() const;

    // Number of fitted segments. make_nasa9_curve guarantees 1 <= segments <= S;
    // the bound by S also lets compilers see that indices stay inside the arrays.
    [[nodiscard]] CANARY_HOST_DEVICE constexpr std::size_t segment_count() const noexcept
    { return std::min(segments, S); }

    // Selection uses the scalar value; all transition arithmetic retains T.
    struct Region {
        std::size_t lower;
        P start;
        P end;
        bool transition;
    };

    template <std::floating_point V> CANARY_HOST_DEVICE Region region(V t) const noexcept
    {
        const std::size_t n = segment_count();
        for (std::size_t i = 1; i < n; ++i) {
            const P start = breaks[i] - blends[i - 1] / P(2);
            const P end = breaks[i] + blends[i - 1] / P(2);
            if (t < start) { return {i - 1, start, end, false}; }
            if (t <= end) { return {i - 1, start, end, true}; }
        }
        return {n - 1, P(0), P(0), false};
    }

    template <HasLog T>
        requires ScalarWith<T, P>
    CANARY_HOST_DEVICE Nasa9CpH<T> inside_cp_h(T t, T ln_t) const noexcept
    {
        using V = ScalarTraits<T>::value_type;
        const auto r = region(ScalarTraits<T>::value(t));
        const auto& lower = b[r.lower];
        const T inv = T(T(V(1)) / t);
        const T cp = detail::nasa9_cp(lower, t, inv);
        const T h = detail::nasa9_h(lower, t, ln_t, inv);
        if (!r.transition) { return {cp, h}; }
        const auto& upper = b[r.lower + 1];
        const T u = T((t - r.start) / (r.end - r.start));
        const T w = u * u * (T(V(3)) - T(V(2)) * u);
        return {T((T(V(1)) - w) * cp + w * detail::nasa9_cp(upper, t, inv)),
                T(h + detail::nasa9_integral<false>(lower, upper, t, inv, r.start, r.end))};
    }

    template <HasLog T>
        requires ScalarWith<T, P>
    CANARY_HOST_DEVICE T inside_s(T t, T ln_t) const noexcept
    {
        using V = ScalarTraits<T>::value_type;
        const auto r = region(ScalarTraits<T>::value(t));
        const T inv = T(T(V(1)) / t);
        const T s = detail::nasa9_entropy(b[r.lower], t, ln_t, inv);
        if (!r.transition) { return s; }
        return T(s + detail::nasa9_integral<true>(b[r.lower], b[r.lower + 1], t, inv, r.start, r.end));
    }

    template <std::size_t M, std::floating_point Q>
    friend Nasa9Curve<Q, M> make_nasa9_curve(const Nasa9Data<M>&, Q);
    template <Scalar T, std::floating_point Q, std::size_t M>
        requires(ScalarWith<T, Q> && HasLog<T>)
    friend CANARY_HOST_DEVICE Nasa9CpH<T> nasa9_cp_h(const Nasa9Curve<Q, M>&, T, T) noexcept;
    template <Scalar T, std::floating_point Q, std::size_t M>
        requires(ScalarWith<T, Q> && HasLog<T>)
    friend CANARY_HOST_DEVICE T nasa9_s(const Nasa9Curve<Q, M>&, T, T) noexcept;
};

// cp/R and h/R at positive temperature t in model units; ln_t must equal ln(t).
// The caller may share ln_t across species. Wide transitions additionally use
// log(t/start), so custom scalars must provide HasLog. Scalar operations,
// ScalarTraits access and log must not throw. Outside the fit cp is
// constant and h continues as h(end) + cp(end)*(t-end).
template <Scalar T, std::floating_point P, std::size_t S>
    requires(ScalarWith<T, P> && HasLog<T>)
CANARY_HOST_DEVICE Nasa9CpH<T> nasa9_cp_h(const Nasa9Curve<P, S>& curve, T t, T ln_t) noexcept
{
    using V = ScalarTraits<T>::value_type;
    const V tv = ScalarTraits<T>::value(t);
    const auto last = curve.segment_count();
    if (tv < curve.breaks[0] || tv > curve.breaks[last]) {
        const bool below = tv < curve.breaks[0];
        const Nasa9Values<P>& end = below ? curve.low : curve.high;
        const P t_end = below ? curve.breaks[0] : curve.breaks[last];
        return {T(static_cast<V>(end.cp_over_r)), T(end.h_over_r + end.cp_over_r * (t - t_end))};
    }
    return curve.inside_cp_h(t, ln_t);
}

// s/R at positive temperature t and the fit's standard pressure, with ln_t=ln(t).
// Requires HasLog for wide transitions; scalar operations, ScalarTraits access
// and log must not throw. Outside the fit s continues as
// s(end) + cp(end)*ln(t/end). No exceptions or allocation occur during evaluation.
template <Scalar T, std::floating_point P, std::size_t S>
    requires(ScalarWith<T, P> && HasLog<T>)
CANARY_HOST_DEVICE T nasa9_s(const Nasa9Curve<P, S>& curve, T t, T ln_t) noexcept
{
    const auto tv = ScalarTraits<T>::value(t);
    const auto last = curve.segment_count();
    if (tv < curve.breaks[0] || tv > curve.breaks[last]) {
        const bool below = tv < curve.breaks[0];
        const Nasa9Values<P>& end = below ? curve.low : curve.high;
        const P ln_t_end = std::log(below ? curve.breaks[0] : curve.breaks[last]);
        return T(end.s_over_r + end.cp_over_r * (ln_t - ln_t_end));
    }
    return curve.inside_s(t, ln_t);
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

// Checked narrowing of host-folded parameters. Coefficients may round to zero;
// positive bounds and widths are checked after storing them in P.
template <std::floating_point P> P nasa9_parameter(long double value, const char* parameter)
{
    if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<P>::max()) {
        throw ParameterError(
          std::format("converted {} is outside the finite parameter range", parameter));
    }
    return static_cast<P>(value);
}

}  // namespace detail

template <std::floating_point P, std::size_t S>
    requires(S > 0)
void Nasa9Curve<P, S>::check() const
{
    using namespace detail;
    const auto n = segment_count();
    for (std::size_t i = 0; i <= n; ++i) {
        if (!(breaks[i] > P(0)) || (i > 0 && !(breaks[i] > breaks[i - 1]))) {
            throw ParameterError(
              std::format("converted boundary {} must be positive and strictly increasing", i));
        }
    }
    P previous_end = breaks[0];
    const auto require_capacity =
      [](const CapacityPolynomial& polynomial, std::size_t region, const char* kind) {
          if (!positive_capacity(polynomial)) {
              throw ParameterError(
                std::format("cp/R > 1 could not be certified throughout {} {}", kind, region));
          }
      };
    for (std::size_t i = 1; i < n; ++i) {
        const P width = blends[i - 1];
        const P start = breaks[i] - width / P(2);
        const P end = breaks[i] + width / P(2);
        if (!(width > P(0)) || !std::isfinite(start) || !std::isfinite(end) ||
            !(start < breaks[i]) || !(end > breaks[i]) || start < previous_end || end > breaks[n]) {
            throw ParameterError(
              std::format("converted transition {} collapses, overlaps or leaves the fit", i));
        }
        require_capacity(capacity_polynomial(b[i - 1], previous_end, start), i - 1, "segment");
        auto lower = capacity_polynomial(b[i - 1], start, end);
        const auto upper = capacity_polynomial(b[i], start, end);
        const auto extent = CapacityInterval {static_cast<long double>(end)} -
                            CapacityInterval {static_cast<long double>(start)};
        const P span = end - start;
        const long double reciprocal = 1.0L / static_cast<long double>(span);
        const auto ratio = extent * capacity_bounds(reciprocal, reciprocal);
        auto blend = lower;
        const auto quadratic = CapacityInterval {3} * ratio * ratio;
        const auto cubic = CapacityInterval {-2} * ratio * ratio * ratio;
        for (std::size_t k = 0; k < 7; ++k) {
            const auto delta = upper[k] - lower[k];
            blend[k + 2] = blend[k + 2] + delta * quadratic;
            blend[k + 3] = blend[k + 3] + delta * cubic;
        }
        require_capacity(blend, i, "transition");
        previous_end = end;
    }
    require_capacity(capacity_polynomial(b[n - 1], previous_end, breaks[n]), n - 1, "segment");
}

// Curve of a fit in a model's temperature unit, with temperature_scale the
// SI value of that unit in kelvin (1 in SI). The coefficients are folded in
// long double precision and then rounded to P:
//   b_k = a_k tau^(k - 2) for k = 0..6, b_7 = (a_7 + a_1 ln tau) / tau,
//   b_8 = a_8 + a_2 ln tau,
// and the breaks and blending widths are divided by tau. Integration constants
// are then shifted to integrate the blended cp continuously while preserving
// the original evaluator at 298.15 K. Throws ParameterError when the scale is not finite
// and positive, the segment count is not in [1, S], the used breaks are not
// finite, positive and increasing, a blending width is not finite and
// positive, the blending regions overlap or leave the fit, or a used
// coefficient is not finite. Converted parameters must remain finite; bounds
// and transitions must remain positive, ordered and noncollapsed. Construction
// rejects a curve unless interval bounds establish cp/R > 1 throughout the fit.
// The endpoint continuation then retains positive ideal-gas cv outside the fit.
// Stored endpoint cp/R, h/R and s/R must all be finite.
template <std::size_t S, std::floating_point P>
Nasa9Curve<P, S> make_nasa9_curve(const Nasa9Data<S>& data, P temperature_scale)
{
    if (!std::isfinite(temperature_scale) || !(temperature_scale > P(0))) {
        throw ParameterError(
          std::format("temperature scale must be finite and positive, got {}", temperature_scale));
    }
    detail::check_nasa9_data(data);
    const auto tau = static_cast<long double>(temperature_scale);
    const long double ln_tau = std::log(tau);
    const auto n = static_cast<std::size_t>(data.segments);
    Nasa9Curve<P, S> curve {};
    curve.segments = n;
    for (std::size_t i = 0; i <= n; ++i) {
        curve.breaks[i] = detail::nasa9_parameter<P>(data.breaks[i] / tau, "boundary");
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        curve.blends[i] = detail::nasa9_parameter<P>(data.blends[i] / tau, "transition width");
    }
    for (std::size_t j = 0; j < n; ++j) {
        const auto& a = data.a[j];
        for (std::size_t k = 0; k < 7; ++k) {
            curve.b[j][k] =
              detail::nasa9_parameter<P>(a[k] * std::pow(tau, static_cast<int>(k) - 2),
                                         "coefficient");
        }
        curve.b[j][7] =
          detail::nasa9_parameter<P>((a[7] + (a[1] * ln_tau)) / tau, "enthalpy constant");
        curve.b[j][8] = detail::nasa9_parameter<P>(a[8] + (a[2] * ln_tau), "entropy constant");
    }
    curve.check();
    // Host evaluation in widened arithmetic avoids cancellation when deriving
    // regional integration constants. The reference retains the original linear
    // coefficient blend and constant-cp continuation before constants are shifted.
    const long double anchor_t = 298.15L / tau;
    const auto anchor = [&](bool original) {
        const long double t = std::clamp(anchor_t,
                                         static_cast<long double>(curve.breaks[0]),
                                         static_cast<long double>(curve.breaks[n]));
        Nasa9Values<long double> out {};
        if (original) {
            const auto r = curve.region(t);
            std::array<long double, 9> k {};
            const long double u = r.transition ? (t - r.start) / P(r.end - r.start) : 0;
            for (std::size_t j = 0; j < 9; ++j) {
                k[j] = (1 - u) * curve.b[r.lower][j];
                if (r.transition) { k[j] += u * curve.b[r.lower + 1][j]; }
            }
            out = {detail::nasa9_cp(k, t),
                   detail::nasa9_h(k, t, std::log(t)),
                   detail::nasa9_entropy(k, t, std::log(t))};
        } else {
            const auto cp_h = curve.inside_cp_h(t, std::log(t));
            out = {cp_h.cp_over_r, cp_h.h_over_r, curve.inside_s(t, std::log(t))};
        }
        out.h_over_r += out.cp_over_r * (anchor_t - t);
        out.s_over_r += out.cp_over_r * std::log(anchor_t / t);
        return out;
    };
    const auto original_anchor = anchor(true);
    for (std::size_t i = 1; i < n; ++i) {
        const P start = curve.breaks[i] - curve.blends[i - 1] / P(2);
        const P end = curve.breaks[i] + curve.blends[i - 1] / P(2);
        const long double t = end;
        const auto& lower = curve.b[i - 1];
        auto& upper = curve.b[i];
        const long double h = detail::nasa9_h(lower, t, std::log(t)) +
                              detail::nasa9_integral<false>(lower, upper, t, start, end);
        const long double s = detail::nasa9_entropy(lower, t, std::log(t)) +
                              detail::nasa9_integral<true>(lower, upper, t, start, end);
        upper[7] = detail::nasa9_parameter<P>(upper[7] + h - detail::nasa9_h(upper, t, std::log(t)),
                                              "enthalpy constant");
        upper[8] =
          detail::nasa9_parameter<P>(upper[8] + s - detail::nasa9_entropy(upper, t, std::log(t)),
                                     "entropy constant");
    }
    const auto blended_anchor = anchor(false);
    for (std::size_t i = 0; i < n; ++i) {
        curve.b[i][7] = detail::nasa9_parameter<P>(
          curve.b[i][7] + (original_anchor.h_over_r - blended_anchor.h_over_r),
          "enthalpy constant");
        curve.b[i][8] = detail::nasa9_parameter<P>(
          curve.b[i][8] + (original_anchor.s_over_r - blended_anchor.s_over_r),
          "entropy constant");
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
    for (const auto& end : {curve.low, curve.high}) {
        if (!std::isfinite(end.cp_over_r) || !(end.cp_over_r > P(1)) ||
            !std::isfinite(end.h_over_r) || !std::isfinite(end.s_over_r)) {
            throw ParameterError("converted endpoint values must be finite with cp/R > 1");
        }
    }
    return curve;
}

}  // namespace canary

#endif  // CANARY_NASA9_HPP
