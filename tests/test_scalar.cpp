#include "dual.hpp"
#include "scalar_types.hpp"
#include "test_suites.hpp"

#include <boost/ut.hpp>
#include <cmath>
#include <concepts>
#include <gaslib/scalar.hpp>
#include <limits>
#include <numbers>

namespace
{

using gaslib::test::Dual;
using D = Dual<double>;

using gaslib::test::NoTraits;
using gaslib::test::OnlySqrt;

}  // namespace

static_assert(gaslib::Scalar<double>);
static_assert(gaslib::Scalar<float>);
static_assert(gaslib::Scalar<D>);
static_assert(gaslib::Scalar<OnlySqrt>);
static_assert(!gaslib::Scalar<int>);
static_assert(!gaslib::Scalar<int*>);
static_assert(!gaslib::Scalar<NoTraits>);

static_assert(gaslib::ScalarWith<double, float>);
static_assert(gaslib::ScalarWith<D, double>);
static_assert(gaslib::ScalarWith<D, float>);
static_assert(gaslib::ScalarWith<float, float>);
static_assert(gaslib::ScalarWith<float, double>);
static_assert(!gaslib::ScalarWith<double, int>);
static_assert(!gaslib::ScalarWith<OnlySqrt, double>);

static_assert(gaslib::HasSqrt<D> && gaslib::HasPow<D> && gaslib::HasExp<D> && gaslib::HasLog<D>);
static_assert(gaslib::HasAbs<D> && gaslib::HasMin<D> && gaslib::HasMax<D>);
static_assert(gaslib::HasSqrt<OnlySqrt>);
static_assert(!gaslib::HasPow<OnlySqrt> && !gaslib::HasExp<OnlySqrt> && !gaslib::HasLog<OnlySqrt>);
static_assert(!gaslib::HasAbs<OnlySqrt> && !gaslib::HasMin<OnlySqrt> && !gaslib::HasMax<OnlySqrt>);

static_assert(std::same_as<decltype(gaslib::math::sqrt(4.0F)), float>);
static_assert(std::same_as<decltype(gaslib::math::pow(2.0F, 3.0F)), float>);
static_assert(std::same_as<decltype(gaslib::math::sqrt(D(4.0))), D>);

void test_scalar()
{
    using namespace boost::ut;
    namespace math = gaslib::math;

    // Each case seeds the derivative of x with 1, so dual() is d f / d x. The
    // dual-number formulas evaluate the analytic derivative, so the results
    // are compared exactly.
    "sqrt of a dual number gives the analytic derivative"_test = [] {
        const D r = math::sqrt(D(4.0, 1.0));
        expect(r.real() == 2.0_d);
        expect(r.dual() == 0.25_d);
    };

    "exp of a dual number gives the analytic derivative"_test = [] {
        const D r = math::exp(D(1.0, 1.0));
        expect(r.real() == std::numbers::e);
        expect(r.dual() == std::numbers::e);
    };

    "log of a dual number gives the analytic derivative"_test = [] {
        const D r = math::log(D(2.0, 1.0));
        expect(r.real() == std::numbers::ln2);
        expect(r.dual() == 0.5_d);
    };

    "pow of a dual number gives the analytic derivative in the base and the exponent"_test = [] {
        const D in_base = math::pow(D(2.0, 1.0), D(3.0));
        expect(in_base.real() == 8.0_d);
        expect(in_base.dual() == 12.0_d);  // 3 x^2 at x = 2
        const D in_exponent = math::pow(D(2.0), D(3.0, 1.0));
        expect(in_exponent.real() == 8.0_d);
        expect(in_exponent.dual() == 8.0 * std::numbers::ln2);  // 2^y ln 2 at y = 3
    };

    "abs of a dual number flips the derivative for a negative value"_test = [] {
        const D r = math::abs(D(-3.0, 1.0));
        expect(r.real() == 3.0_d);
        expect(r.dual() == -1.0_d);
        const D s = math::abs(D(3.0, 1.0));
        expect(s.dual() == 1.0_d);
    };

    "min and max of dual numbers return the selected operand with its derivative"_test = [] {
        const D a(1.0, 1.0);
        const D b(2.0, 0.0);
        expect(math::min(a, b).real() == 1.0_d && math::min(a, b).dual() == 1.0_d);
        expect(math::min(b, a).real() == 1.0_d && math::min(b, a).dual() == 1.0_d);
        expect(math::max(a, b).real() == 2.0_d && math::max(a, b).dual() == 0.0_d);
        expect(math::max(b, a).real() == 2.0_d && math::max(b, a).dual() == 0.0_d);
    };

    "isfinite checks only the value part"_test = [] {
        constexpr double inf = std::numeric_limits<double>::infinity();
        constexpr double nan = std::numeric_limits<double>::quiet_NaN();
        expect(math::isfinite(D(1.0, inf)));
        expect(!math::isfinite(D(inf, 0.0)));
        expect(!math::isfinite(D(nan, 0.0)));
        expect(math::isfinite(1.0));
        expect(!math::isfinite(nan));
        expect(!math::isfinite(-inf));
    };

    "floating-point arguments use the standard functions"_test = [] {
        expect(math::sqrt(5.0) == std::sqrt(5.0));
        expect(math::pow(2.0, 0.3) == std::pow(2.0, 0.3));
        expect(math::exp(0.7) == std::exp(0.7));
        expect(math::log(0.7) == std::log(0.7));
        expect(math::abs(-0.7) == 0.7_d);
        expect(math::min(1.0, 2.0) == 1.0_d && math::max(1.0, 2.0) == 2.0_d);
        expect(math::sqrt(6.0F) == std::sqrt(6.0F));
    };

    "a type with only sqrt can use math::sqrt"_test = [] {
        expect(math::sqrt(OnlySqrt {9.0}).v == 3.0_d);
    };

    "mixed operations with a parameter carry the derivative"_test = [] {
        const D x(2.0, 1.0);
        const D r = 3.0F * x + 1.0 - x / 4.0;  // d/dx = 3 - 1/4
        expect(r.real() == 6.5_d);
        expect(r.dual() == 2.75_d);
    };
}
