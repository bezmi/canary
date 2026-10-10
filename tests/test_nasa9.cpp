#include "dual.hpp"
#include "test_suites.hpp"

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <cstddef>
#include <gaslib/build.hpp>
#include <gaslib/data/air5.hpp>
#include <gaslib/nasa9.hpp>
#include <gaslib/scalar.hpp>
#include <gaslib/units.hpp>
#include <limits>
#include <string>
#include <vector>

namespace
{

using D = gaslib::test::Dual<double>;
using Curve = gaslib::Nasa9Curve<double, 4>;
using Data = gaslib::Nasa9Data<4>;

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

// Reference at t inside the fit, with lambda-weighted coefficients inside a
// blending region.
Reference reference_at(const Data& data, double t)
{
    const auto n = static_cast<std::size_t>(data.segments);
    std::size_t seg = 0;
    for (std::size_t i = 1; i < n; ++i) {
        const double start = data.breaks[i] - (data.blends[i - 1] / 2);
        const double end = data.breaks[i] + (data.blends[i - 1] / 2);
        if (t >= start && t <= end) {
            const long double lambda = (t - start) / data.blends[i - 1];
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

struct Values {
    double cp;
    double h;
    double s;
};

Values evaluate(const Curve& curve, double t)
{
    const double ln_t = std::log(t);
    const auto cp_h = gaslib::nasa9_cp_h(curve, t, ln_t);
    return {.cp = cp_h.cp_over_r, .h = cp_h.h_over_r, .s = gaslib::nasa9_s(curve, t, ln_t)};
}

bool close_to(double value, long double expected, long double scale, double ulps)
{ return std::abs(value - expected) <= ulps * eps * scale; }

bool relative_jump(double a, double b, double bound)
{ return std::abs(a - b) <= bound * std::abs(b); }

// JANAF reference row: Cp and S in J/(mol K), H - H(298.15) in kJ/mol,
// at 1 bar.
struct Janaf {
    std::size_t species;
    double temp;
    double cp;
    double s;
    double dh;
};

constexpr std::array<Janaf, 10> janaf {
  {
    {.species = 0, .temp = 1000, .cp = 32.697, .s = 228.170, .dh = 21.463},
    {.species = 0, .temp = 5000, .cp = 37.912, .s = 286.041, .dh = 167.763},
    {.species = 1, .temp = 1000, .cp = 34.870, .s = 243.578, .dh = 22.703},
    {.species = 1, .temp = 5000, .cp = 42.675, .s = 305.538, .dh = 180.749},
    {.species = 4, .temp = 1000, .cp = 33.987, .s = 248.536, .dh = 22.229},
    {.species = 4, .temp = 5000, .cp = 38.208, .s = 307.499, .dh = 170.730},
    {.species = 2, .temp = 1000, .cp = 20.786, .s = 178.454, .dh = 14.589},
    {.species = 2, .temp = 5000, .cp = 23.461, .s = 212.469, .dh = 100.111},
    {.species = 3, .temp = 1000, .cp = 20.915, .s = 186.790, .dh = 14.860},
    {.species = 3, .temp = 5000, .cp = 21.799, .s = 220.578, .dh = 99.222},
  },
};

constexpr const auto& species = gaslib::air5_species;

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
        const double r_u = gaslib::constants::universal_gas_constant<double>;
        for (const Janaf& row : janaf) {
            const auto& sp = species[row.species];
            const Curve curve = gaslib::make_nasa9_curve(sp.fit, 1.0);
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

    "a curve with temperature scale 1 has the fit's coefficients"_test = [] {
        for (const auto& sp : species) {
            const Curve curve = gaslib::make_nasa9_curve(sp.fit, 1.0);
            const auto n = static_cast<std::size_t>(sp.fit.segments);
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t k = 0; k < 9; ++k) {
                    expect(curve.b[j][k] == sp.fit.a[j][k]) << sp.name << j << k;
                }
            }
            for (std::size_t i = 0; i <= n; ++i) { expect(curve.breaks[i] == sp.fit.breaks[i]); }
            for (std::size_t i = 0; i + 1 < n; ++i) { expect(curve.blends[i] == sp.fit.blends[i]); }
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
            const Curve si = gaslib::make_nasa9_curve(sp.fit, 1.0);
            const Curve scaled = gaslib::make_nasa9_curve(sp.fit, tau);
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
            const Curve curve = gaslib::make_nasa9_curve(sp.fit, 1.0);
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

    "values inside a blending region use the lambda-weighted coefficients"_test = [] {
        for (const auto& sp : species) {
            const Curve curve = gaslib::make_nasa9_curve(sp.fit, 1.0);
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

    "values inside segments equal the polynomial forms"_test = [] {
        for (const auto& sp : species) {
            const Curve curve = gaslib::make_nasa9_curve(sp.fit, 1.0);
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
            const Curve curve = gaslib::make_nasa9_curve(sp.fit, 1.0);
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
            const Curve curve = gaslib::make_nasa9_curve(sp.fit, 1.0);
            for (const double t : {300.0, 700.0, 2000.0, 5000.0, 10000.0, 15000.0, 4.0e4}) {
                if (t > sp.fit.breaks[static_cast<std::size_t>(sp.fit.segments)]) { continue; }
                const D td(t, 1.0);
                const D ln_t = gaslib::math::log(td);
                const auto cp_h = gaslib::nasa9_cp_h(curve, td, ln_t);
                const D s = gaslib::nasa9_s(curve, td, ln_t);
                const Reference r = reference_at(sp.fit, t);
                expect(close_to(cp_h.h_over_r.dual(), r.cp, r.cp_scale, 16.0))
                  << sp.name << t << cp_h.h_over_r.dual() << cp_h.cp_over_r.real();
                expect(close_to(t * s.dual(), r.cp, r.cp_scale, 16.0))
                  << sp.name << t << t * s.dual() << cp_h.cp_over_r.real();
                expect(cp_h.cp_over_r.real() == evaluate(curve, t).cp);
            }
        }
    };

    "make_nasa9_curve rejects invalid fits"_test = [] {
        const Data good = species[0].fit;
        expect(nothrow([&] { gaslib::make_nasa9_curve(good, 1.0); }));
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
            expect(throws<gaslib::ParameterError>([&] { gaslib::make_nasa9_curve(data, 1.0); }));
        }
        expect(throws<gaslib::ParameterError>([&] { gaslib::make_nasa9_curve(good, -1.0); }));
    };
}
