#ifndef CANARY_TEST_NASA9_BLEND_OFFSETS_HPP
#define CANARY_TEST_NASA9_BLEND_OFFSETS_HPP

// Difference between Canary's NASA9 thermodynamics and the reference
// implementation's at one temperature, for tests that compare against
// values captured from that implementation. Both codes use the same NASA9
// coefficients and differ only in how they join two segments at a segment
// boundary:
// - the reference implementation takes, inside the blending region of width
//   w centred on the boundary, the linear combination (1 - u) a_lower +
//   u a_upper of the nine coefficients, u = (t - start) / w, and the plain
//   segment coefficients elsewhere;
// - Canary blends cp with the cubic smoothstep u^2 (3 - 2u), integrates that
//   cp for h and s, shifts the integration constants of each segment above a
//   blend so that h and s are continuous, and shifts all constants so that
//   the 298.15 K values equal the reference rule's.
// Beyond the fit both continue with constant cp. Above a blend the two
// rules therefore differ by two constants per species in h/R and s/R.
//
// Both rules are evaluated here from the raw coefficients in long double,
// the smoothstep integrals by Gauss-Legendre quadrature; nothing is shared
// with canary/nasa9.hpp or canary/detail/nasa9_integrals.hpp, so a test
// using these offsets still checks the library's h and s. Temperatures are
// in kelvin and fits are evaluated in SI (temperature scale 1).

#include <algorithm>
#include <array>
#include <canary/nasa9.hpp>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <utility>

namespace canary::test
{

// Canary's h/R (K) and s/R minus the reference implementation's.
struct BlendOffset {
    long double h_over_r;
    long double s_over_r;
};

namespace detail
{

using Coefficients = std::array<long double, 9>;

struct HS {
    long double h_over_r;
    long double s_over_r;
};

inline long double cp_of(const Coefficients& k, long double t)
{
    return k[0] / (t * t) + k[1] / t + k[2] + k[3] * t + k[4] * t * t + k[5] * t * t * t +
           k[6] * t * t * t * t;
}

inline long double h_of(const Coefficients& k, long double t)
{
    const long double lt = std::log(t);
    return -k[0] / t + k[1] * lt + k[2] * t + k[3] * t * t / 2.0L + k[4] * t * t * t / 3.0L +
           k[5] * t * t * t * t / 4.0L + k[6] * t * t * t * t * t / 5.0L + k[7];
}

inline long double s_of(const Coefficients& k, long double t)
{
    const long double lt = std::log(t);
    return -k[0] / (2.0L * t * t) - k[1] / t + k[2] * lt + k[3] * t + k[4] * t * t / 2.0L +
           k[5] * t * t * t / 3.0L + k[6] * t * t * t * t / 4.0L + k[8];
}

template <std::size_t S> Coefficients segment(const Nasa9Data<S>& fit, std::size_t j)
{
    Coefficients k {};
    for (std::size_t m = 0; m < 9; ++m) { k[m] = fit.a[j][m]; }
    return k;
}

// Start and end of the blending region of interior boundary i.
template <std::size_t S>
std::pair<long double, long double> blend_bounds(const Nasa9Data<S>& fit, std::size_t i)
{ return {fit.breaks[i] - fit.blends[i - 1] / 2.0L, fit.breaks[i] + fit.blends[i - 1] / 2.0L}; }

// Reference rule inside the fit.
template <std::size_t S> HS reference_inside(const Nasa9Data<S>& fit, long double t)
{
    const auto n = static_cast<std::size_t>(fit.segments);
    for (std::size_t i = 1; i < n; ++i) {
        const auto [start, end] = blend_bounds(fit, i);
        if (t < start) {
            const Coefficients k = segment(fit, i - 1);
            return {h_of(k, t), s_of(k, t)};
        }
        if (t <= end) {
            const long double u = (t - start) / (end - start);
            const Coefficients lo = segment(fit, i - 1);
            const Coefficients hi = segment(fit, i);
            Coefficients k {};
            for (std::size_t m = 0; m < 9; ++m) { k[m] = (1.0L - u) * lo[m] + u * hi[m]; }
            return {h_of(k, t), s_of(k, t)};
        }
    }
    const Coefficients k = segment(fit, n - 1);
    return {h_of(k, t), s_of(k, t)};
}

template <std::size_t S> long double reference_cp(const Nasa9Data<S>& fit, long double t)
{
    const auto n = static_cast<std::size_t>(fit.segments);
    for (std::size_t i = 1; i < n; ++i) {
        const auto [start, end] = blend_bounds(fit, i);
        if (t < start) { return cp_of(segment(fit, i - 1), t); }
        if (t <= end) {
            const long double u = (t - start) / (end - start);
            return (1.0L - u) * cp_of(segment(fit, i - 1), t) + u * cp_of(segment(fit, i), t);
        }
    }
    return cp_of(segment(fit, n - 1), t);
}

// Reference rule at any t, constant cp beyond the fit.
template <std::size_t S> HS reference(const Nasa9Data<S>& fit, long double t)
{
    const auto n = static_cast<std::size_t>(fit.segments);
    const long double lo = fit.breaks[0];
    const long double hi = fit.breaks[n];
    if (t < lo || t > hi) {
        const long double te = t < lo ? lo : hi;
        const HS e = reference_inside(fit, te);
        const long double cp = reference_cp(fit, te);
        return {e.h_over_r + cp * (t - te), e.s_over_r + cp * std::log(t / te)};
    }
    return reference_inside(fit, t);
}

// Gauss-Legendre nodes and weights on [-1, 1], by Newton iteration on the
// Legendre polynomial of degree Q.
template <std::size_t Q> struct GaussLegendre {
    std::array<long double, Q> node {};
    std::array<long double, Q> weight {};

    GaussLegendre()
    {
        for (std::size_t i = 0; i < Q; ++i) {
            long double x = std::cos(std::numbers::pi_v<long double> * (i + 0.75L) / (Q + 0.5L));
            long double dp = 0.0L;
            for (int it = 0; it < 100; ++it) {
                long double p0 = 1.0L;
                long double p1 = x;
                for (std::size_t m = 2; m <= Q; ++m) {
                    const long double p2 = ((2.0L * m - 1.0L) * x * p1 - (m - 1.0L) * p0) / m;
                    p0 = p1;
                    p1 = p2;
                }
                dp = Q * (x * p1 - p0) / (x * x - 1.0L);
                const long double dx = p1 / dp;
                x -= dx;
                if (std::abs(dx) < 1.0e-19L) { break; }
            }
            node[i] = x;
            weight[i] = 2.0L / ((1.0L - x * x) * dp * dp);
        }
    }
};

inline const GaussLegendre<64>& quadrature()
{
    static const GaussLegendre<64> q;
    return q;
}

// Integral of f over [a, b]. The integrands here are polynomials in t and
// 1/t times the smoothstep weight, smooth on [a, b] with a > 0; 64 points
// integrate them to long double precision.
template <typename F> long double integrate(F&& f, long double a, long double b)
{
    const auto& q = quadrature();
    const long double mid = (a + b) / 2.0L;
    const long double half = (b - a) / 2.0L;
    long double sum = 0.0L;
    for (std::size_t i = 0; i < q.node.size(); ++i) {
        sum += q.weight[i] * f(mid + half * q.node[i]);
    }
    return half * sum;
}

// Canary's rule: the segments with their integration constants shifted for
// continuity across every blend, in ascending boundary order, then anchored.
template <std::size_t S> struct CanaryRule {
    std::size_t segments;
    std::array<Coefficients, S> k;
    const Nasa9Data<S>* fit;

    long double blend_cp(std::size_t i, long double t) const
    {
        const auto [start, end] = blend_bounds(*fit, i);
        const long double u = (t - start) / (end - start);
        const long double w = u * u * (3.0L - 2.0L * u);
        return (1.0L - w) * cp_of(k[i - 1], t) + w * cp_of(k[i], t);
    }

    HS inside(long double t) const
    {
        for (std::size_t i = 1; i < segments; ++i) {
            const auto [start, end] = blend_bounds(*fit, i);
            if (t < start) { return {h_of(k[i - 1], t), s_of(k[i - 1], t)}; }
            if (t <= end) {
                const auto dcp = [&](long double x) { return blend_cp(i, x) - cp_of(k[i - 1], x); };
                const long double h = h_of(k[i - 1], t) + integrate(dcp, start, t);
                const long double s =
                  s_of(k[i - 1], t) +
                  integrate([&](long double x) { return dcp(x) / x; }, start, t);
                return {h, s};
            }
        }
        return {h_of(k[segments - 1], t), s_of(k[segments - 1], t)};
    }

    long double cp(long double t) const
    {
        for (std::size_t i = 1; i < segments; ++i) {
            const auto [start, end] = blend_bounds(*fit, i);
            if (t < start) { return cp_of(k[i - 1], t); }
            if (t <= end) { return blend_cp(i, t); }
        }
        return cp_of(k[segments - 1], t);
    }

    HS at(long double t) const
    {
        const long double lo = fit->breaks[0];
        const long double hi = fit->breaks[segments];
        if (t < lo || t > hi) {
            const long double te = t < lo ? lo : hi;
            const HS e = inside(te);
            const long double c = cp(te);
            return {e.h_over_r + c * (t - te), e.s_over_r + c * std::log(t / te)};
        }
        return inside(t);
    }
};

template <std::size_t S> CanaryRule<S> canary_rule(const Nasa9Data<S>& fit)
{
    CanaryRule<S> c {};
    c.segments = static_cast<std::size_t>(fit.segments);
    c.fit = &fit;
    for (std::size_t j = 0; j < c.segments; ++j) { c.k[j] = segment(fit, j); }
    // Continuity: the upper segment's constants take the values the blended
    // integrals reach at the end of the blend.
    for (std::size_t i = 1; i < c.segments; ++i) {
        const long double end = blend_bounds(fit, i).second;
        const HS at_end = c.inside(end);
        c.k[i][7] += at_end.h_over_r - h_of(c.k[i], end);
        c.k[i][8] += at_end.s_over_r - s_of(c.k[i], end);
    }
    // Anchor: the reference rule's values at 298.15 K, clamped into the fit
    // and continued with constant cp outside it.
    const long double anchor_t = 298.15L;
    const long double ta = std::clamp(anchor_t,
                                      static_cast<long double>(fit.breaks[0]),
                                      static_cast<long double>(fit.breaks[c.segments]));
    HS original = reference_inside(fit, ta);
    const long double cp_original = reference_cp(fit, ta);
    original.h_over_r += cp_original * (anchor_t - ta);
    original.s_over_r += cp_original * std::log(anchor_t / ta);
    HS blended = c.inside(ta);
    const long double cp_blended = c.cp(ta);
    blended.h_over_r += cp_blended * (anchor_t - ta);
    blended.s_over_r += cp_blended * std::log(anchor_t / ta);
    for (std::size_t j = 0; j < c.segments; ++j) {
        c.k[j][7] += original.h_over_r - blended.h_over_r;
        c.k[j][8] += original.s_over_r - blended.s_over_r;
    }
    return c;
}

}  // namespace detail

// Offset of one fit at temperature t (K).
template <std::size_t S> BlendOffset blend_offset(const Nasa9Data<S>& fit, double t)
{
    const detail::HS c = detail::canary_rule(fit).at(t);
    const detail::HS r = detail::reference(fit, t);
    return {c.h_over_r - r.h_over_r, c.s_over_r - r.s_over_r};
}

// Offsets of every species, in species order; Species is any record with a
// Nasa9Data member named fit.
template <typename Species, std::size_t N>
std::array<BlendOffset, N> blend_offsets(const std::array<Species, N>& species, double t)
{
    std::array<BlendOffset, N> out {};
    for (std::size_t i = 0; i < N; ++i) { out[i] = blend_offset(species[i].fit, t); }
    return out;
}

// Factor by which an equilibrium backward rate under Canary's rule exceeds
// the reference implementation's, for a reaction with coefficient changes
// dnu_i and Gibbs energies at t (K). With g_i = h_i - t s_i, ln K_c
// contains -sum_i dnu_i g_i / (R t) and the backward rate is k_f / K_c, so
// the ratio is exp(sum_i dnu_i (dh_i/R / t - ds_i/R)).
template <std::size_t N>
double
  blend_rate_factor(const std::array<BlendOffset, N>& off, const std::array<int, N>& dnu, double t)
{
    long double sum = 0.0L;
    for (std::size_t i = 0; i < N; ++i) { sum += dnu[i] * (off[i].h_over_r / t - off[i].s_over_r); }
    return static_cast<double>(std::exp(sum));
}

}  // namespace canary::test

#endif  // CANARY_TEST_NASA9_BLEND_OFFSETS_HPP
