#ifndef CANARY_TESTS_DUAL_HPP
#define CANARY_TESTS_DUAL_HPP

// Forward-mode dual number shared by host and device tests. It has the operators and the free
// functions sqrt, pow, abs, min, max, isnan and isinf of Ibis::Dual (ibis
// src/util/util/dual.h), with the same overloads and formulas, plus exp and
// log. Mixed operations with the value type are non-template friends taking
// const V&, as in Ibis::Dual, so a parameter of another floating-point type
// converts to V on the call.

#include <canary/config.hpp>
#include <canary/scalar.hpp>
#include <cmath>
#include <concepts>

namespace canary::test
{

template <std::floating_point V> class Dual
{
  public:
    constexpr Dual() = default;
    CANARY_HOST_DEVICE constexpr Dual(V real, V dual) : real_(real), dual_(dual) {}
    // Implicit, as in Ibis::Dual.
    CANARY_HOST_DEVICE constexpr Dual(V real) : real_(real) {}  // NOLINT(*-explicit-constructor)

    CANARY_HOST_DEVICE constexpr Dual& operator=(const V& real)
    {
        real_ = real;
        dual_ = V(0);
        return *this;
    }

    friend CANARY_HOST_DEVICE constexpr bool operator<(const Dual& lhs, const Dual& rhs)
    { return lhs.real_ < rhs.real_; }
    friend CANARY_HOST_DEVICE constexpr bool operator>(const Dual& lhs, const Dual& rhs)
    { return lhs.real_ > rhs.real_; }
    friend CANARY_HOST_DEVICE constexpr bool operator<=(const Dual& lhs, const Dual& rhs)
    { return lhs.real_ <= rhs.real_; }
    friend CANARY_HOST_DEVICE constexpr bool operator>=(const Dual& lhs, const Dual& rhs)
    { return lhs.real_ >= rhs.real_; }
    friend CANARY_HOST_DEVICE constexpr bool operator==(const Dual& lhs, const Dual& rhs)
    { return lhs.real_ == rhs.real_; }
    friend CANARY_HOST_DEVICE constexpr bool operator!=(const Dual& lhs, const Dual& rhs)
    { return lhs.real_ != rhs.real_; }

    friend CANARY_HOST_DEVICE constexpr Dual operator+(const Dual& lhs, const Dual& rhs)
    { return {lhs.real_ + rhs.real_, lhs.dual_ + rhs.dual_}; }
    friend CANARY_HOST_DEVICE constexpr Dual operator+(const Dual& lhs, const V& rhs)
    { return {lhs.real_ + rhs, lhs.dual_}; }
    friend CANARY_HOST_DEVICE constexpr Dual operator+(const V& lhs, const Dual& rhs)
    { return {lhs + rhs.real_, rhs.dual_}; }

    friend CANARY_HOST_DEVICE constexpr Dual operator-(const Dual& lhs, const Dual& rhs)
    { return {lhs.real_ - rhs.real_, lhs.dual_ - rhs.dual_}; }
    friend CANARY_HOST_DEVICE constexpr Dual operator-(const Dual& lhs, const V& rhs)
    { return {lhs.real_ - rhs, lhs.dual_}; }
    friend CANARY_HOST_DEVICE constexpr Dual operator-(const V& lhs, const Dual& rhs)
    { return {lhs - rhs.real_, -rhs.dual_}; }

    CANARY_HOST_DEVICE constexpr Dual operator-() const { return {-real_, -dual_}; }

    friend CANARY_HOST_DEVICE constexpr Dual operator*(const Dual& lhs, const Dual& rhs)
    { return {lhs.real_ * rhs.real_, (lhs.real_ * rhs.dual_) + (lhs.dual_ * rhs.real_)}; }
    friend CANARY_HOST_DEVICE constexpr Dual operator*(const Dual& lhs, const V& rhs)
    { return {lhs.real_ * rhs, lhs.dual_ * rhs}; }
    friend CANARY_HOST_DEVICE constexpr Dual operator*(const V& lhs, const Dual& rhs)
    { return {rhs.real_ * lhs, rhs.dual_ * lhs}; }

    friend CANARY_HOST_DEVICE constexpr Dual operator/(const Dual& num, const Dual& den)
    {
        return {
          num.real_ / den.real_,
          ((num.dual_ * den.real_) - (num.real_ * den.dual_)) / (den.real_ * den.real_),
        };
    }
    friend CANARY_HOST_DEVICE constexpr Dual operator/(const Dual& num, const V& den)
    { return {num.real_ / den, num.dual_ / den}; }
    friend CANARY_HOST_DEVICE constexpr Dual operator/(const V& num, const Dual& den)
    { return {num / den.real_, -den.dual_ * num / (den.real_ * den.real_)}; }

    CANARY_HOST_DEVICE constexpr Dual& operator+=(const Dual& other)
    { return *this = *this + other; }
    CANARY_HOST_DEVICE constexpr Dual& operator+=(const V& other) { return *this = *this + other; }
    CANARY_HOST_DEVICE constexpr Dual& operator-=(const Dual& other)
    { return *this = *this - other; }
    CANARY_HOST_DEVICE constexpr Dual& operator-=(const V& other) { return *this = *this - other; }
    CANARY_HOST_DEVICE constexpr Dual& operator*=(const Dual& other)
    { return *this = *this * other; }
    CANARY_HOST_DEVICE constexpr Dual& operator*=(const V& other) { return *this = *this * other; }
    CANARY_HOST_DEVICE constexpr Dual& operator/=(const Dual& other)
    { return *this = *this / other; }
    CANARY_HOST_DEVICE constexpr Dual& operator/=(const V& other) { return *this = *this / other; }

    [[nodiscard]] CANARY_HOST_DEVICE constexpr V real() const { return real_; }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr V dual() const { return dual_; }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr Dual abs() const
    { return real_ < V(0) ? -*this : *this; }
    [[nodiscard]] CANARY_HOST_DEVICE constexpr Dual conjugate() const { return {real_, -dual_}; }

  private:
    V real_ {};
    V dual_ {};
};

template <std::floating_point V> CANARY_HOST_DEVICE Dual<V> sqrt(const Dual<V>& d)
{
    const V real = std::sqrt(d.real());
    return {real, d.dual() / (V(2) * real)};
}

template <std::floating_point V>
CANARY_HOST_DEVICE Dual<V> pow(const Dual<V>& base, const Dual<V>& power)
{
    const V real = std::pow(base.real(), power.real());
    const V dual = (base.dual() * power.real() * std::pow(base.real(), power.real() - V(1))) +
                   (power.dual() * real * std::log(base.real()));
    return {real, dual};
}

template <std::floating_point V> CANARY_HOST_DEVICE Dual<V> exp(const Dual<V>& d)
{
    const V real = std::exp(d.real());
    return {real, d.dual() * real};
}

template <std::floating_point V> CANARY_HOST_DEVICE Dual<V> log(const Dual<V>& d)
{ return {std::log(d.real()), d.dual() / d.real()}; }

template <std::floating_point V> CANARY_HOST_DEVICE Dual<V> abs(const Dual<V>& d)
{ return d.abs(); }

template <std::floating_point V>
CANARY_HOST_DEVICE Dual<V> max(const Dual<V>& d1, const Dual<V>& d2)
{ return d1.real() > d2.real() ? d1 : d2; }
template <std::floating_point V> CANARY_HOST_DEVICE Dual<V> max(const V& d1, const Dual<V>& d2)
{ return d1 > d2.real() ? Dual<V>(d1) : d2; }
template <std::floating_point V> CANARY_HOST_DEVICE Dual<V> max(const Dual<V>& d1, const V& d2)
{ return d1.real() > d2 ? d1 : Dual<V>(d2); }

template <std::floating_point V>
CANARY_HOST_DEVICE Dual<V> min(const Dual<V>& d1, const Dual<V>& d2)
{ return d1.real() < d2.real() ? d1 : d2; }
template <std::floating_point V> CANARY_HOST_DEVICE Dual<V> min(const V& d1, const Dual<V>& d2)
{ return d1 < d2.real() ? Dual<V>(d1) : d2; }
template <std::floating_point V> CANARY_HOST_DEVICE Dual<V> min(const Dual<V>& d1, const V& d2)
{ return d1.real() < d2 ? d1 : Dual<V>(d2); }

template <std::floating_point V> CANARY_HOST_DEVICE bool isnan(const Dual<V>& d)
{ return std::isnan(d.real()) || std::isnan(d.dual()); }

template <std::floating_point V> CANARY_HOST_DEVICE bool isinf(const Dual<V>& d)
{ return std::isinf(d.real()) || std::isinf(d.dual()); }

}  // namespace canary::test

template <std::floating_point V> struct canary::ScalarTraits<canary::test::Dual<V>> {
    using value_type = V;

    static CANARY_HOST_DEVICE constexpr V value(const canary::test::Dual<V>& x) { return x.real(); }

    template <std::floating_point W>
    static CANARY_HOST_DEVICE constexpr canary::test::Dual<W> cast(const canary::test::Dual<V>& x)
    { return {static_cast<W>(x.real()), static_cast<W>(x.dual())}; }
};

#endif  // CANARY_TESTS_DUAL_HPP
