#include "dual.hpp"
#include "test_suites.hpp"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/data/air11_species.hpp>
#include <canary/data/air5.hpp>
#include <canary/nasa9.hpp>
#include <canary/scalar.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

using D = canary::test::Dual<double>;
using Curve = canary::Nasa9Curve<double, 4>;
using Data = canary::Nasa9Data<4>;

constexpr double eps = std::numeric_limits<double>::epsilon();
constexpr double kilo = 1000.0;

// Values of the polynomial forms in long double, evaluated term by term
// with std::pow, and the sums of the absolute values of the terms. A
// rounding error of a few eps in each term of a double evaluation gives an
// error of a few eps times that sum, which is the scale the tolerances use.
struct Reference {
    long double cp;
    long double h;
    long double s;
    long double cp_scale;
    long double h_scale;
    long double s_scale;
};

Reference reference(const std::array<long double, 9>& a, long double t)
{
    const long double ln_t = std::log(t);
    const std::array<long double, 7> cp_terms {
      a[0] / (t * t),
      a[1] / t,
      a[2],
      a[3] * t,
      a[4] * t * t,
      a[5] * std::pow(t, 3.0L),
      a[6] * std::pow(t, 4.0L),
    };
    const std::array<long double, 8> h_terms {
      -a[0] / t,
      a[1] * ln_t,
      a[2] * t,
      a[3] * t * t / 2,
      a[4] * std::pow(t, 3.0L) / 3,
      a[5] * std::pow(t, 4.0L) / 4,
      a[6] * std::pow(t, 5.0L) / 5,
      a[7],
    };
    const std::array<long double, 8> s_terms {
      -a[0] / (2 * t * t),
      -a[1] / t,
      a[2] * ln_t,
      a[3] * t,
      a[4] * t * t / 2,
      a[5] * std::pow(t, 3.0L) / 3,
      a[6] * std::pow(t, 4.0L) / 4,
      a[8],
    };
    Reference r {};
    for (const long double x : cp_terms) {
        r.cp += x;
        r.cp_scale += std::abs(x);
    }
    for (const long double x : h_terms) {
        r.h += x;
        r.h_scale += std::abs(x);
    }
    for (const long double x : s_terms) {
        r.s += x;
        r.s_scale += std::abs(x);
    }
    return r;
}

std::array<long double, 9> segment(const Data& data, std::size_t j)
{
    std::array<long double, 9> a {};
    for (std::size_t k = 0; k < 9; ++k) { a[k] = data.a[j][k]; }
    return a;
}

// Original evaluator retained only for the reference anchor.
Reference original_reference_at(const Data& data, double t)
{
    const auto n = static_cast<std::size_t>(data.segments);
    if (t < data.breaks[0] || t > data.breaks[n]) {
        const auto i = t < data.breaks[0] ? std::size_t(0) : n - 1;
        const double end = t < data.breaks[0] ? data.breaks[0] : data.breaks[n];
        auto out = reference(segment(data, i), end);
        out.h += out.cp * (t - end);
        out.s += out.cp * std::log(static_cast<long double>(t) / end);
        return out;
    }
    std::size_t seg = 0;
    for (std::size_t i = 1; i < n; ++i) {
        const double start = data.breaks[i] - (data.blends[i - 1] / 2);
        const double end = data.breaks[i] + (data.blends[i - 1] / 2);
        if (t >= start && t <= end) {
            const long double lambda = (t - start) / (end - start);
            const auto lower = segment(data, i - 1);
            const auto upper = segment(data, i);
            std::array<long double, 9> a {};
            for (std::size_t k = 0; k < 9; ++k) {
                a[k] = (lambda * upper[k]) + ((1 - lambda) * lower[k]);
            }
            return reference(a, t);
        }
        if (t > end) { seg = i; }
    }
    return reference(segment(data, seg), t);
}

// Independent quadrature integrates the prescribed cp, without using the
// production analytic moments. Sixteen-point Gauss-Legendre quadrature resolves
// the smooth rational NASA9 terms over the shipped transition intervals.
template <typename F> long double integrate(F f, long double a, long double b)
{
    static const auto rule = [] {
        std::array<std::pair<long double, long double>, 16> out {};
        for (int i = 0; i < 16; ++i) {
            long double x = std::cos(std::acos(-1.0L) * (i + 0.75L) / 16.5L);
            long double derivative = 0;
            for (int iteration = 0; iteration < 20; ++iteration) {
                long double p0 = 1, p1 = x;
                for (int n = 2; n <= 16; ++n) {
                    const long double p = ((2 * n - 1) * x * p1 - (n - 1) * p0) / n;
                    p0 = p1;
                    p1 = p;
                }
                derivative = 16 * (x * p1 - p0) / (x * x - 1);
                const long double next = x - p1 / derivative;
                if (next == x) { break; }
                x = next;
            }
            out[i] = {x, 2 / ((1 - x * x) * derivative * derivative)};
        }
        return out;
    }();
    long double result = 0;
    for (const auto [x, w] : rule) { result += w * f((a + b) / 2 + (b - a) * x / 2); }
    return result * (b - a) / 2;
}

Reference smooth_capacity(const Data& data, long double t)
{
    t = std::clamp(t,
                   static_cast<long double>(data.breaks[0]),
                   static_cast<long double>(data.breaks[data.segments]));
    for (std::size_t i = 1; i < static_cast<std::size_t>(data.segments); ++i) {
        const long double start = data.breaks[i] - data.blends[i - 1] / 2;
        const long double end = data.breaks[i] + data.blends[i - 1] / 2;
        if (t < start) { return reference(segment(data, i - 1), t); }
        if (t <= end) {
            const long double u = (t - start) / (end - start);
            const long double w = u * u * (3 - 2 * u);
            auto lower = reference(segment(data, i - 1), t);
            const auto upper = reference(segment(data, i), t);
            lower.cp = (1 - w) * lower.cp + w * upper.cp;
            lower.cp_scale = (1 - w) * lower.cp_scale + w * upper.cp_scale;
            return lower;
        }
    }
    return reference(segment(data, static_cast<std::size_t>(data.segments - 1)), t);
}

Reference reference_at(const Data& data, double t)
{
    auto out = smooth_capacity(data, t);
    const auto anchor = original_reference_at(data, 298.15);
    out.h = anchor.h;
    out.s = anchor.s;
    out.h_scale += anchor.h_scale;
    out.s_scale += anchor.s_scale;
    std::vector<long double> points {298.15L, t};
    const long double lo = std::min(298.15L, static_cast<long double>(t));
    const long double hi = std::max(298.15L, static_cast<long double>(t));
    for (std::size_t i = 1; i < static_cast<std::size_t>(data.segments); ++i) {
        for (const double edge :
             {data.breaks[i] - data.blends[i - 1] / 2, data.breaks[i] + data.blends[i - 1] / 2}) {
            if (edge > lo && edge < hi) { points.push_back(edge); }
        }
    }
    for (const double edge : {data.breaks[0], data.breaks[data.segments]}) {
        if (edge > lo && edge < hi) { points.push_back(edge); }
    }
    std::sort(points.begin(), points.end());
    const long double sign = t < 298.15 ? -1 : 1;
    for (std::size_t i = 1; i < points.size(); ++i) {
        const long double a = points[i - 1], b = points[i];
        // Split long ordinary intervals as well to resolve inverse powers.
        const int panels = static_cast<int>(std::ceil(b / a));
        for (int panel = 0; panel < panels; ++panel) {
            const long double left = a + (b - a) * panel / panels;
            const long double right = a + (b - a) * (panel + 1) / panels;
            const long double h =
              integrate([&](long double x) { return smooth_capacity(data, x).cp; }, left, right);
            const long double entropy =
              integrate([&](long double x) { return smooth_capacity(data, x).cp / x; },
                        left,
                        right);
            out.h += sign * h;
            out.s += sign * entropy;
            out.h_scale += std::abs(h);
            out.s_scale += std::abs(entropy);
        }
    }
    return out;
}

struct Values {
    double cp;
    double h;
    double s;
};

Values evaluate(const Curve& curve, double t)
{
    const double ln_t = std::log(t);
    const auto cp_h = canary::nasa9_cp_h(curve, t, ln_t);
    return {.cp = cp_h.cp_over_r, .h = cp_h.h_over_r, .s = canary::nasa9_s(curve, t, ln_t)};
}

bool close_to(double value, long double expected, long double scale, double ulps)
{ return std::abs(value - expected) <= ulps * eps * scale; }

bool relative_jump(double a, double b, double bound)
{ return std::abs(a - b) <= bound * std::abs(b); }

// JANAF reference row: Cp and S in J/(mol K), H - H(298.15) in kJ/mol,
// at 1 bar.
struct Janaf {
    canary::SpeciesIndex species;
    double temp;
    double cp;
    double s;
    double dh;
};

constexpr std::array<Janaf, 10> janaf {
  {
    {.species = canary::SpeciesIndex {0}, .temp = 1000, .cp = 32.697, .s = 228.170, .dh = 21.463},
    {.species = canary::SpeciesIndex {0}, .temp = 5000, .cp = 37.912, .s = 286.041, .dh = 167.763},
    {.species = canary::SpeciesIndex {1}, .temp = 1000, .cp = 34.870, .s = 243.578, .dh = 22.703},
    {.species = canary::SpeciesIndex {1}, .temp = 5000, .cp = 42.675, .s = 305.538, .dh = 180.749},
    {.species = canary::SpeciesIndex {4}, .temp = 1000, .cp = 33.987, .s = 248.536, .dh = 22.229},
    {.species = canary::SpeciesIndex {4}, .temp = 5000, .cp = 38.208, .s = 307.499, .dh = 170.730},
    {.species = canary::SpeciesIndex {2}, .temp = 1000, .cp = 20.786, .s = 178.454, .dh = 14.589},
    {.species = canary::SpeciesIndex {2}, .temp = 5000, .cp = 23.461, .s = 212.469, .dh = 100.111},
    {.species = canary::SpeciesIndex {3}, .temp = 1000, .cp = 20.915, .s = 186.790, .dh = 14.860},
    {.species = canary::SpeciesIndex {3}, .temp = 5000, .cp = 21.799, .s = 220.578, .dh = 99.222},
  },
};

constexpr const auto& species = canary::air5_species;

}  // namespace

void test_nasa9()
{
    using namespace boost::ut;
    // The fits themselves differ from JANAF by up to 1.04 % in Cp (NO at
    // 5000 K), 0.06 % in S (O2 at 5000 K) and 1.83 % in H - H(298.15) (O at
    // 1000 K, whose fit has cp/R = 2.5 below 1000 K, without the electronic
    // contribution that JANAF includes). The tolerances bound these
    // differences, so the test checks the data and their evaluation against
    // an independent table; exact evaluation is checked separately.
    "curves reproduce the JANAF tables within the accuracy of the fits"_test = [] {
        constexpr double cp_tolerance = 1.5e-2;
        constexpr double s_tolerance = 1.0e-3;
        constexpr double dh_tolerance = 2.0e-2;
        const double r_u = canary::constants::universal_gas_constant<double>;
        for (const Janaf& row : janaf) {
            const auto& sp = species[row.species.value()];
            const Curve curve = canary::make_nasa9_curve(sp.fit, 1.0);
            const Values v = evaluate(curve, row.temp);
            const Values v298 = evaluate(curve, 298.15);
            const double cp = r_u * v.cp;
            const double s = r_u * v.s;
            const double dh = r_u * (v.h - v298.h) / kilo;
            const std::string where = std::string(sp.name) + " at " + std::to_string(row.temp);
            expect(relative_jump(cp, row.cp, cp_tolerance)) << "Cp of" << where << cp << row.cp;
            expect(relative_jump(s, row.s, s_tolerance)) << "S of" << where << s << row.s;
            expect(relative_jump(dh, row.dh, dh_tolerance))
              << "H - H298 of" << where << dh << row.dh;
        }
    };

    "validated curves are copyable values with preserved reference anchors"_test = [] {
        static_assert(!std::is_default_constructible_v<Curve>);
        static_assert(!std::is_aggregate_v<Curve>);
        static_assert(std::is_trivially_copyable_v<Curve>);
        static_assert(std::is_nothrow_copy_constructible_v<Curve>);
        static_assert(noexcept(canary::nasa9_cp_h(std::declval<const Curve&>(), 1.0, 0.0)));
        static_assert(noexcept(canary::nasa9_s(std::declval<const Curve&>(), 1.0, 0.0)));
        for (const auto& sp : species) {
            const Curve curve = canary::make_nasa9_curve(sp.fit, 1.0);
            const Curve copy = curve;
            const auto v = evaluate(copy, 298.15);
            const auto r = original_reference_at(sp.fit, 298.15);
            expect(close_to(v.cp, r.cp, r.cp_scale, 4.0)) << sp.name;
            expect(close_to(v.h, r.h, r.h_scale, 4.0)) << sp.name;
            expect(close_to(v.s, r.s, r.s_scale, 4.0)) << sp.name;
        }
    };

    "a scaled curve gives the SI values at the same temperature"_test = [] {
        // tau = 1000 K. Folding rounds each coefficient once or twice
        // (tau^(k-2) and the product; ln tau, the sum and the division for
        // b7 and b8), and both evaluations round each term a few times, so
        // the results agree to a few eps times the sum of the absolute
        // values of the terms. Temperatures are chosen so that t = tau t'
        // is exact.
        constexpr double tau = 1000.0;
        for (const auto& sp : species) {
            const Curve si = canary::make_nasa9_curve(sp.fit, 1.0);
            const Curve scaled = canary::make_nasa9_curve(sp.fit, tau);
            for (const double t : {0.25, 0.75, 1.0, 1.5, 3.0, 6.25, 12.5, 19.5, 30.0, 45.0}) {
                if (t * tau > sp.fit.breaks[static_cast<std::size_t>(sp.fit.segments)]) {
                    continue;
                }
                const Values a = evaluate(scaled, t);
                const Values b = evaluate(si, t * tau);
                const Reference r = reference_at(sp.fit, t * tau);
                expect(close_to(a.cp, b.cp, r.cp_scale, 4.0)) << sp.name << t << a.cp << b.cp;
                expect(close_to(a.s, b.s, r.s_scale, 4.0)) << sp.name << t << a.s << b.s;
                expect(close_to(a.h * tau, b.h, r.h_scale, 4.0)) << sp.name << t << a.h << b.h;
            }
        }
    };

    "values are continuous at the edges of every blending region"_test = [] {
        for (const auto& sp : species) {
            const Curve curve = canary::make_nasa9_curve(sp.fit, 1.0);
            const auto n = static_cast<std::size_t>(sp.fit.segments);
            for (std::size_t i = 1; i < n; ++i) {
                const double half = sp.fit.blends[i - 1] / 2;
                const double start = sp.fit.breaks[i] - half;
                const double end = sp.fit.breaks[i] + half;
                const Values in_start = evaluate(curve, start);
                const Values out_start = evaluate(curve, std::nextafter(start, 0.0));
                const Values in_end = evaluate(curve, end);
                const Values out_end = evaluate(curve, std::nextafter(end, 1.0e300));
                expect(relative_jump(in_start.cp, out_start.cp, 1.0e-12)) << sp.name << i;
                expect(relative_jump(in_start.h, out_start.h, 1.0e-12)) << sp.name << i;
                expect(relative_jump(in_start.s, out_start.s, 1.0e-12)) << sp.name << i;
                expect(relative_jump(in_end.cp, out_end.cp, 1.0e-12)) << sp.name << i;
                expect(relative_jump(in_end.h, out_end.h, 1.0e-12)) << sp.name << i;
                expect(relative_jump(in_end.s, out_end.s, 1.0e-12)) << sp.name << i;
            }
        }
    };

    "transition values integrate the smoothstep heat capacity"_test = [] {
        for (const auto& sp : species) {
            const Curve curve = canary::make_nasa9_curve(sp.fit, 1.0);
            const auto n = static_cast<std::size_t>(sp.fit.segments);
            for (std::size_t i = 1; i < n; ++i) {
                const double w = sp.fit.blends[i - 1];
                for (const double lambda : {0.0, 0.125, 0.5, 0.75, 1.0}) {
                    const double t = sp.fit.breaks[i] - (w / 2) + (lambda * w);
                    const Values v = evaluate(curve, t);
                    const Reference r = reference_at(sp.fit, t);
                    expect(close_to(v.cp, r.cp, r.cp_scale, 8.0)) << sp.name << t << v.cp;
                    expect(close_to(v.h, r.h, r.h_scale, 8.0)) << sp.name << t << v.h;
                    expect(close_to(v.s, r.s, r.s_scale, 8.0)) << sp.name << t << v.s;
                }
            }
        }
    };

    "segment capacities retain their polynomials and primitives remain continuous"_test = [] {
        for (const auto& sp : species) {
            const Curve curve = canary::make_nasa9_curve(sp.fit, 1.0);
            for (const double t : {200.0, 300.0, 700.0, 2000.0, 5000.0, 10000.0, 20000.0, 4.0e4}) {
                if (t > sp.fit.breaks[static_cast<std::size_t>(sp.fit.segments)]) { continue; }
                const Values v = evaluate(curve, t);
                const Reference r = reference_at(sp.fit, t);
                expect(close_to(v.cp, r.cp, r.cp_scale, 8.0)) << sp.name << t << v.cp;
                expect(close_to(v.h, r.h, r.h_scale, 8.0)) << sp.name << t << v.h;
                expect(close_to(v.s, r.s, r.s_scale, 8.0)) << sp.name << t << v.s;
            }
        }
    };

    "beyond the fit cp is constant and h and s continue from the ends"_test = [] {
        for (const auto& sp : species) {
            const Curve curve = canary::make_nasa9_curve(sp.fit, 1.0);
            const auto n = static_cast<std::size_t>(sp.fit.segments);
            for (const double t_end : {sp.fit.breaks[0], sp.fit.breaks[n]}) {
                const double out = t_end == sp.fit.breaks[0] ? std::nextafter(t_end, 0.0)
                                                             : std::nextafter(t_end, 1.0e300);
                const Values at = evaluate(curve, t_end);
                const Values beyond = evaluate(curve, out);
                expect(relative_jump(beyond.cp, at.cp, 1.0e-12)) << sp.name << t_end;
                expect(relative_jump(beyond.h, at.h, 1.0e-12)) << sp.name << t_end;
                expect(relative_jump(beyond.s, at.s, 1.0e-12)) << sp.name << t_end;
                for (const double factor : {0.05, 0.5, 2.0, 3.0}) {
                    if ((factor < 1.0) != (t_end == sp.fit.breaks[0])) { continue; }
                    const double t = factor * t_end;
                    const Values v = evaluate(curve, t);
                    expect(v.cp == at.cp) << sp.name << t;
                    const double h = at.h + (at.cp * (t - t_end));
                    const double s = at.s + (at.cp * std::log(t / t_end));
                    expect(relative_jump(v.h, h, 4.0 * eps)) << sp.name << t << v.h << h;
                    expect(std::abs(v.s - s) <= 8.0 * eps * (std::abs(at.s) + std::abs(s)))
                      << sp.name << t << v.s << s;
                }
            }
        }
    };

    "dual-number derivatives of h and s give cp inside segments"_test = [] {
        // Polynomial derivatives are summed in a different order from cp,
        // so they agree to a few eps times the sum of the absolute values
        // of the cp terms.
        for (const auto& sp : species) {
            const Curve curve = canary::make_nasa9_curve(sp.fit, 1.0);
            for (const double t : {300.0, 700.0, 2000.0, 5000.0, 10000.0, 15000.0, 4.0e4}) {
                if (t > sp.fit.breaks[static_cast<std::size_t>(sp.fit.segments)]) { continue; }
                const D td(t, 1.0);
                const D ln_t = canary::math::log(td);
                const auto cp_h = canary::nasa9_cp_h(curve, td, ln_t);
                const D s = canary::nasa9_s(curve, td, ln_t);
                const Reference r = reference_at(sp.fit, t);
                expect(close_to(cp_h.h_over_r.dual(), r.cp, r.cp_scale, 16.0))
                  << sp.name << t << cp_h.h_over_r.dual() << cp_h.cp_over_r.real();
                expect(close_to(t * s.dual(), r.cp, r.cp_scale, 16.0))
                  << sp.name << t << t * s.dual() << cp_h.cp_over_r.real();
                expect(cp_h.cp_over_r.real() == evaluate(curve, t).cp);
            }
        }
    };

    "smoothstep derivatives agree with capacity in rounded SI and scaled curves"_test = [] {
        const auto check = []<typename P> {
            using Dual = canary::test::Dual<P>;
            const P epsilon = std::numeric_limits<P>::epsilon();
            for (const auto& sp : canary::air11_species) {
                for (const P scale : {P(1), P(3000)}) {
                    const auto curve = canary::make_nasa9_curve(sp.fit, scale);
                    for (int i = 1; i < sp.fit.segments; ++i) {
                        const P boundary = P(static_cast<long double>(sp.fit.breaks[i]) / scale);
                        const P width = P(static_cast<long double>(sp.fit.blends[i - 1]) / scale);
                        const P start = boundary - width / P(2), end = boundary + width / P(2);
                        for (const P u : {P(0), P(0.125), P(0.5), P(0.875), P(1)}) {
                            const P t = start + u * (end - start);
                            const Dual td(t, P(1));
                            const auto values =
                              canary::nasa9_cp_h(curve, td, canary::math::log(td));
                            const auto entropy = canary::nasa9_s(curve, td, canary::math::log(td));
                            const auto ref =
                              smooth_capacity(sp.fit, static_cast<long double>(t) * scale);
                            if (u == P(0) || u == P(1)) {
                                const P outside =
                                  std::nextafter(t,
                                                 u == P(0) ? P(0) : std::numeric_limits<P>::max());
                                const Dual adjacent(outside, P(1));
                                const auto other =
                                  canary::nasa9_cp_h(curve, adjacent, canary::math::log(adjacent));
                                expect(std::abs(values.cp_over_r.dual() - other.cp_over_r.dual()) <=
                                       64 * epsilon * ref.cp_scale / (end - start))
                                  << sp.name << scale << t;
                            }
                            // Cancellation in the fitted polynomial is measured by
                            // the sum of absolute cp terms, not by the final cp alone.
                            const long double bound = 64 * epsilon * ref.cp_scale;
                            expect(std::abs(values.h_over_r.dual() - values.cp_over_r.real()) <=
                                   bound)
                              << sp.name << scale << t;
                            expect(std::abs(t * entropy.dual() - values.cp_over_r.real()) <= bound)
                              << sp.name << scale << t;
                            if constexpr (std::same_as<P, double>) {
                                const P delta = std::cbrt(epsilon) * std::min(t, end - start);
                                const auto plus =
                                  canary::nasa9_cp_h(curve, t + delta, std::log(t + delta));
                                const auto minus =
                                  canary::nasa9_cp_h(curve, t - delta, std::log(t - delta));
                                const P sp = canary::nasa9_s(curve, t + delta, std::log(t + delta));
                                const P sm = canary::nasa9_s(curve, t - delta, std::log(t - delta));
                                const long double truncation =
                                  64 * std::pow(epsilon, 2.0 / 3.0) * ref.cp_scale;
                                const long double h_bound =
                                  truncation +
                                  4 * epsilon *
                                    (std::abs(plus.h_over_r) + std::abs(minus.h_over_r)) / delta;
                                const long double s_bound =
                                  truncation +
                                  4 * epsilon * t * (std::abs(sp) + std::abs(sm)) / delta;
                                expect(std::abs((plus.h_over_r - minus.h_over_r) / (2 * delta) -
                                                values.cp_over_r.real()) <= h_bound);
                                expect(std::abs(t * (sp - sm) / (2 * delta) -
                                                values.cp_over_r.real()) <= s_bound);
                            }
                        }
                    }
                }
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "wide transitions integrate each NASA9 capacity power with a reference inside the blend"_test =
      [] {
          const auto check = []<typename P> {
              using Dual = canary::test::Dual<P>;
              const P epsilon = std::numeric_limits<P>::epsilon();
              for (std::size_t power = 0; power < 7; ++power) {
                  Data data {};
                  data.segments = 2;
                  data.breaks = {50, 1000, 3000};
                  data.blends[0] = 1800;
                  data.a[0][2] = data.a[1][2] = 3;
                  data.a[1][power] += std::pow(1000.0, 2 - static_cast<int>(power));
                  data.a[0][7] = 100;
                  data.a[1][7] = 200;
                  data.a[0][8] = 2;
                  data.a[1][8] = 5;
                  for (const P scale : {P(1), P(1000)}) {
                      const auto curve = canary::make_nasa9_curve(data, scale);
                      for (const double physical_t :
                           {100.0, 199.0, 201.0, 298.15, 800.0, 1000.0, 1800.0, 1900.0, 2000.0}) {
                          const P t = P(physical_t) / scale;
                          const auto ref = reference_at(data, static_cast<double>(t) * scale);
                          const Dual td(t, P(1));
                          const auto values = canary::nasa9_cp_h(curve, td, canary::math::log(td));
                          const auto entropy = canary::nasa9_s(curve, td, canary::math::log(td));
                          expect(std::abs(values.cp_over_r.real() - ref.cp) <=
                                 64 * epsilon * ref.cp_scale)
                            << power << t;
                          expect(std::abs(values.h_over_r.real() * scale - ref.h) <=
                                 64 * epsilon * ref.h_scale)
                            << power << t;
                          expect(std::abs(entropy.real() - ref.s) <= 64 * epsilon * ref.s_scale)
                            << power << t;
                          expect(std::abs(values.h_over_r.dual() - values.cp_over_r.real()) <=
                                 64 * epsilon * ref.cp_scale)
                            << power << t;
                          expect(std::abs(t * entropy.dual() - values.cp_over_r.real()) <=
                                 64 * epsilon * ref.cp_scale)
                            << power << t;
                      }
                  }
              }
          };
          check.template operator()<float>();
          check.template operator()<double>();
      };

    "reference anchors outside the fit retain constant capacity continuation"_test = [] {
        const auto check = []<typename P> {
            using Dual = canary::test::Dual<P>;
            const P epsilon = std::numeric_limits<P>::epsilon();
            for (const double offset : {0.0, 400.0}) {
                Data data {};
                data.segments = 2;
                data.breaks = {50 + offset, 150 + offset, 250 + offset};
                data.blends[0] = 100;
                data.a[0][2] = 3;
                data.a[1][2] = 5;
                data.a[0][7] = 100;
                data.a[1][7] = 200;
                data.a[0][8] = 2;
                data.a[1][8] = 5;
                for (const P scale : {P(1), P(1000)}) {
                    const auto curve = canary::make_nasa9_curve(data, scale);
                    for (const double physical_t : {298.15,
                                                    25 + offset,
                                                    50 + offset,
                                                    125 + offset,
                                                    175 + offset,
                                                    250 + offset,
                                                    275 + offset}) {
                        const P t = P(physical_t) / scale;
                        const auto ref = reference_at(data, static_cast<double>(t) * scale);
                        const Dual td(t, P(1));
                        const auto values = canary::nasa9_cp_h(curve, td, canary::math::log(td));
                        const auto entropy = canary::nasa9_s(curve, td, canary::math::log(td));
                        expect(std::abs(values.cp_over_r.real() - ref.cp) <=
                               64 * epsilon * ref.cp_scale)
                          << offset << scale << t;
                        expect(std::abs(values.h_over_r.real() * scale - ref.h) <=
                               64 * epsilon * ref.h_scale)
                          << offset << scale << t;
                        expect(std::abs(entropy.real() - ref.s) <= 64 * epsilon * ref.s_scale)
                          << offset << scale << t;
                        expect(std::abs(values.h_over_r.dual() - values.cp_over_r.real()) <=
                               64 * epsilon * ref.cp_scale)
                          << offset << scale << t;
                        expect(std::abs(t * entropy.dual() - values.cp_over_r.real()) <=
                               64 * epsilon * ref.cp_scale)
                          << offset << scale << t;
                    }
                }
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "rounded transition endpoints retain adjacent heat capacities"_test = [] {
        const auto check = []<typename P> {
            Data data {};
            data.segments = 2;
            data.breaks = {0.75, 1.5, 3.0};
            data.blends[0] = 1.5 * (std::nextafter(P(1.5), P(2)) - P(1.5));
            // The rounded endpoint separation is larger than the stored width.
            // Heat capacity remains between the adjacent values throughout.
            data.a[0][2] = 5;
            data.a[1][2] = 2;
            const auto curve = canary::make_nasa9_curve(data, P(1));
            const P start = P(data.breaks[1]) - P(data.blends[0]) / P(2);
            const P end = P(data.breaks[1]) + P(data.blends[0]) / P(2);
            expect(end - start != P(data.blends[0]));
            for (const auto [t, cp] : std::array {std::pair {start, P(5)}, std::pair {end, P(2)}}) {
                const auto values = canary::nasa9_cp_h(curve, t, std::log(t));
                expect(values.cp_over_r == cp);
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "NASA9 construction checks capacity throughout segments and transitions"_test = [] {
        const auto check = []<typename P> {
            Data data {};
            data.segments = 1;
            data.breaks = {1.0, 3.0};
            for (const double cp : {-1.0, 0.5, 1.0}) {
                data.a[0][2] = cp;
                expect(
                  throws<canary::ParameterError>([&] { canary::make_nasa9_curve(data, P(1)); }));
            }
            data.a[0][2] = std::nextafter(P(1), P(2));
            expect(nothrow([&] { canary::make_nasa9_curve(data, P(1)); }));

            // cp/R = minimum + (T-2)^2 has positive endpoint cv in all cases.
            data.a[0][3] = -4.0;
            data.a[0][4] = 1.0;
            for (const double minimum : {0.5, 1.0}) {
                data.a[0][2] = 4.0 + minimum;
                expect(
                  throws<canary::ParameterError>([&] { canary::make_nasa9_curve(data, P(1)); }));
            }
            data.a[0][2] = 5.25;
            expect(nothrow([&] { canary::make_nasa9_curve(data, P(1)); }));

            // Negative cv occupies only |T-1.5| < 1/1024 inside the fit.
            data.a[0][2] = 3.25 - 0x1p-20;
            data.a[0][3] = -3.0;
            expect(throws<canary::ParameterError>([&] { canary::make_nasa9_curve(data, P(1)); }));

            data = {};
            data.segments = 2;
            data.breaks = {1.0, 2.0, 3.0};
            data.blends[0] = 1.0;
            data.a[0][2] = 14.0;
            data.a[0][3] = -8.0;
            data.a[1][2] = -18.0;
            data.a[1][3] = 8.0;
            // Both unblended regions and both transition endpoints have cp/R >= 2.
            // The blended capacity has cp/R = -2 at T=2.
            expect(throws<canary::ParameterError>([&] { canary::make_nasa9_curve(data, P(1)); }));
            data.a[0][2] += 3.5;
            data.a[1][2] += 3.5;
            expect(nothrow([&] { canary::make_nasa9_curve(data, P(1)); }));

            // T²(cp/R-1) = (T-2)^6 + offset in the lower segment;
            // the upper segment doubles it, giving a degree-nine transition.
            for (const double offset : {-0.015625, 0.0, 0.015625}) {
                const std::array<double, 7>
                  coefficients {64.0 + offset, -192.0, 240.0, -160.0, 60.0, -12.0, 1.0};
                for (std::size_t k = 0; k < coefficients.size(); ++k) {
                    data.a[0][k] = coefficients[k];
                    data.a[1][k] = 2.0 * coefficients[k];
                }
                data.a[0][2] += 1.0;
                data.a[1][2] += 1.0;
                if (offset > 0) {
                    expect(nothrow([&] { canary::make_nasa9_curve(data, P(1)); }));
                } else {
                    expect(throws<canary::ParameterError>(
                      [&] { canary::make_nasa9_curve(data, P(1)); }));
                }
            }
        };
        check.template operator()<float>();
        check.template operator()<double>();
    };

    "NASA9 construction rejects unrepresentable converted parameters"_test = [] {
        const auto check = []<typename P> {
            Data good {};
            good.segments = 2;
            good.breaks = {1.0, 2.0, 3.0};
            good.blends[0] = 1.0;
            good.a[0][2] = 3.0;
            good.a[1][2] = 3.0;
            const auto rejects = [](const Data& data, P scale, const std::string& reason) {
                bool rejected = false;
                try {
                    static_cast<void>(canary::make_nasa9_curve(data, scale));
                } catch (const canary::ParameterError& error) {
                    rejected = true;
                    expect(std::string(error.what()).find(reason) != std::string::npos)
                      << error.what();
                }
                expect(rejected) << reason;
            };
            Data data = good;
            data.breaks = {1.0e-100, 2.0e-100, 3.0e-100};
            data.blends[0] = 1.0e-100;
            rejects(data, std::numeric_limits<P>::max(), "converted boundary");

            data = good;
            data.blends[0] = std::numeric_limits<double>::denorm_min();
            rejects(data, P(1), "converted transition");
            data.blends[0] = std::numeric_limits<P>::epsilon() / 2;
            rejects(data, P(1), "converted transition");

            data = good;
            data.a[0][0] = static_cast<double>(std::numeric_limits<P>::max());
            rejects(data, P(0.5), "converted coefficient");
            data = good;
            data.a[0][7] = static_cast<double>(std::numeric_limits<P>::max());
            rejects(data, P(0.5), "converted enthalpy constant");
            data = good;
            data.a[0][2] = static_cast<double>(std::numeric_limits<P>::max());
            data.a[1][2] = data.a[0][2];
            rejects(data, P(1), "converted endpoint");
            // The logarithmic unit shift overflows the stored entropy constant.
            data.a[0][8] = data.a[0][2];
            rejects(data, P(2), "converted entropy constant");

            data = good;
            data.breaks = {1.0, 2.0, 10.0};
            rejects(data, std::numeric_limits<P>::min(), "converted boundary");
        };
        check.template operator()<float>();
        check.template operator()<double>();

        Data data {};
        data.segments = 1;
        data.breaks = {1.0, std::nextafter(1.0, 2.0)};
        data.a[0][2] = 3.0;
        expect(nothrow([&] { canary::make_nasa9_curve(data, 1.0); }));
        expect(throws<canary::ParameterError>([&] { canary::make_nasa9_curve(data, 1.0F); }));
    };

    "rounded SI and scaled air curves retain positive capacity"_test = [] {
        const auto check = []<typename P>(const auto& inputs) {
            for (const auto& sp : inputs) {
                for (const P scale : {P(1), P(3000)}) {
                    const auto curve = canary::make_nasa9_curve(sp.fit, scale);
                    for (int region = 0; region < sp.fit.segments; ++region) {
                        const auto j = static_cast<std::size_t>(region);
                        for (const P fraction : {P(0), P(0.25), P(0.5), P(0.75), P(1)}) {
                            const P t = P(sp.fit.breaks[j] / scale) +
                                        fraction * (P(sp.fit.breaks[j + 1] / scale) -
                                                    P(sp.fit.breaks[j] / scale));
                            const auto values = canary::nasa9_cp_h(curve, t, std::log(t));
                            expect(std::isfinite(values.cp_over_r) && values.cp_over_r > P(1))
                              << sp.name;
                            expect(std::isfinite(values.h_over_r)) << sp.name;
                            expect(std::isfinite(canary::nasa9_s(curve, t, std::log(t))))
                              << sp.name;
                        }
                    }
                }
            }
        };
        check.template operator()<float>(species);
        check.template operator()<double>(species);
        check.template operator()<float>(canary::air11_species);
        check.template operator()<double>(canary::air11_species);
    };

    "make_nasa9_curve rejects invalid fits"_test = [] {
        const Data good = species[0].fit;
        expect(nothrow([&] { canary::make_nasa9_curve(good, 1.0); }));
        std::vector<Data> bad(8, good);
        bad[0].segments = 0;
        bad[1].segments = 5;
        bad[2].breaks[2] = bad[2].breaks[1];
        bad[3].blends[1] = 0.0;
        bad[4].blends[1] = 10000.0;  // overlaps the region of the first boundary
        bad[5].a[1][3] = std::numeric_limits<double>::infinity();
        bad[6].breaks[0] = 0.0;
        bad[7].blends[0] = 2000.0;  // starts before the fit
        for (const Data& data : bad) {
            expect(throws<canary::ParameterError>([&] { canary::make_nasa9_curve(data, 1.0); }));
        }
        expect(throws<canary::ParameterError>([&] { canary::make_nasa9_curve(good, -1.0); }));
    };
}
