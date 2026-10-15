#ifndef CANARY_SCALAR_HPP
#define CANARY_SCALAR_HPP

#include <bit>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <canary/config.hpp>

namespace canary
{

// Access to the value part of a state scalar: the number itself for a
// floating-point type, the real part for a dual number. A consumer
// specialises it for its own scalar type, giving value_type (a
// floating-point type) and a static function value(x) that returns it.
// Unit conversion that changes precision also requires cast<W>(x), returning
// the same scalar family with value_type W. It must convert every component,
// including derivatives, and support conversion back to the original type.
template <typename T> struct ScalarTraits;

template <std::floating_point T> struct ScalarTraits<T> {
    using value_type = T;

    CANARY_HOST_DEVICE static constexpr T value(T x) { return x; }

    template <std::floating_point W> CANARY_HOST_DEVICE static constexpr W cast(T x)
    { return static_cast<W>(x); }
};

// Type of the per-cell values that canary functions compute with: double,
// float, or a dual number that carries derivatives. Comparisons use the
// value part. A constant (a value with zero derivative) is built from a
// plain value of the value type. The ScalarTraits requirement is checked
// first, so a type without it is rejected before any operator lookup.
template <typename T>
concept HasScalarTraits = requires { typename ScalarTraits<T>::value_type; } &&
                          std::floating_point<typename ScalarTraits<T>::value_type>;

template <typename T>
concept Scalar = HasScalarTraits<T> && std::copyable<T> &&
                 std::constructible_from<T, typename ScalarTraits<T>::value_type> &&
                 requires(const T a, const T b) {
                     {
                         ScalarTraits<T>::value(a)
                     } -> std::same_as<typename ScalarTraits<T>::value_type>;
                     { -a } -> std::convertible_to<T>;
                     { a + b } -> std::convertible_to<T>;
                     { a - b } -> std::convertible_to<T>;
                     { a * b } -> std::convertible_to<T>;
                     { a / b } -> std::convertible_to<T>;
                     { a < b } -> std::convertible_to<bool>;
                 };

// A state scalar T that combines directly with a model parameter of
// floating-point type P, in either operand position, without first
// converting the parameter to T. The result of a mixed operation may have a
// wider type than T (float state with a double parameter gives double), but
// must convert to T.
template <typename T, typename P>
concept ScalarWith = Scalar<T> && std::floating_point<P> && requires(const T t, const P p) {
    { t + p } -> std::convertible_to<T>;
    { p + t } -> std::convertible_to<T>;
    { t - p } -> std::convertible_to<T>;
    { p - t } -> std::convertible_to<T>;
    { t * p } -> std::convertible_to<T>;
    { p * t } -> std::convertible_to<T>;
    { t / p } -> std::convertible_to<T>;
    { p / t } -> std::convertible_to<T>;
};

// Two scalar types that combine in either operand position, each result
// being a scalar type. Used where both operands may be state values or
// parameters, for example a gas constant that is a model parameter for one
// species and a per-cell value for a mixture.
template <typename A, typename B>
concept ScalarPair = Scalar<A> && Scalar<B> && requires(const A a, const B b) {
    { a + b } -> Scalar;
    { b + a } -> Scalar;
    { a - b } -> Scalar;
    { b - a } -> Scalar;
    { a * b } -> Scalar;
    { b * a } -> Scalar;
    { a / b } -> Scalar;
    { b / a } -> Scalar;
};

namespace detail::scalar_lookup
{

// Deleted declarations stop ordinary lookup here, so that an unqualified call
// below finds only the overloads declared next to the argument's type
// (argument-dependent lookup).
void sqrt() = delete;
void pow() = delete;
void exp() = delete;
void log() = delete;
void abs() = delete;
void min() = delete;
void max() = delete;

template <typename T>
concept HasOwnSqrt = requires(const T x) {
    { sqrt(x) } -> std::convertible_to<T>;
};
template <typename T>
concept HasOwnPow = requires(const T x, const T y) {
    { pow(x, y) } -> std::convertible_to<T>;
};
template <typename T>
concept HasOwnExp = requires(const T x) {
    { exp(x) } -> std::convertible_to<T>;
};
template <typename T>
concept HasOwnLog = requires(const T x) {
    { log(x) } -> std::convertible_to<T>;
};
template <typename T>
concept HasOwnAbs = requires(const T x) {
    { abs(x) } -> std::convertible_to<T>;
};
template <typename T>
concept HasOwnMin = requires(const T x, const T y) {
    { min(x, y) } -> std::convertible_to<T>;
};
template <typename T>
concept HasOwnMax = requires(const T x, const T y) {
    { max(x, y) } -> std::convertible_to<T>;
};

template <HasOwnSqrt T> CANARY_HOST_DEVICE T call_sqrt(const T& x) { return sqrt(x); }
template <HasOwnPow T> CANARY_HOST_DEVICE T call_pow(const T& x, const T& y) { return pow(x, y); }
template <HasOwnExp T> CANARY_HOST_DEVICE T call_exp(const T& x) { return exp(x); }
template <HasOwnLog T> CANARY_HOST_DEVICE T call_log(const T& x) { return log(x); }
template <HasOwnAbs T> CANARY_HOST_DEVICE T call_abs(const T& x) { return abs(x); }
template <HasOwnMin T> CANARY_HOST_DEVICE T call_min(const T& x, const T& y) { return min(x, y); }
template <HasOwnMax T> CANARY_HOST_DEVICE T call_max(const T& x, const T& y) { return max(x, y); }

}  // namespace detail::scalar_lookup

// Each math function below has its own concept. A floating-point type uses
// the standard library function; any other scalar type must declare the
// function next to its type. A scalar type that lacks one of them can still
// be used with every model that does not call it.
template <typename T>
concept HasSqrt = Scalar<T> && (std::floating_point<T> || detail::scalar_lookup::HasOwnSqrt<T>);
template <typename T>
concept HasPow = Scalar<T> && (std::floating_point<T> || detail::scalar_lookup::HasOwnPow<T>);
template <typename T>
concept HasExp = Scalar<T> && (std::floating_point<T> || detail::scalar_lookup::HasOwnExp<T>);
template <typename T>
concept HasLog = Scalar<T> && (std::floating_point<T> || detail::scalar_lookup::HasOwnLog<T>);
template <typename T>
concept HasAbs = Scalar<T> && (std::floating_point<T> || detail::scalar_lookup::HasOwnAbs<T>);
template <typename T>
concept HasMin = Scalar<T> && (std::floating_point<T> || detail::scalar_lookup::HasOwnMin<T>);
template <typename T>
concept HasMax = Scalar<T> && (std::floating_point<T> || detail::scalar_lookup::HasOwnMax<T>);

namespace math
{

template <HasSqrt T> CANARY_HOST_DEVICE T sqrt(const T& x)
{
    if constexpr (std::floating_point<T>) {
        return std::sqrt(x);
    } else {
        return detail::scalar_lookup::call_sqrt(x);
    }
}

// The exponent has the scalar type, so a dual-number exponent carries its
// derivative.
template <HasPow T> CANARY_HOST_DEVICE T pow(const T& base, const T& exponent)
{
    if constexpr (std::floating_point<T>) {
        return std::pow(base, exponent);
    } else {
        return detail::scalar_lookup::call_pow(base, exponent);
    }
}

template <HasExp T> CANARY_HOST_DEVICE T exp(const T& x)
{
    if constexpr (std::floating_point<T>) {
        return std::exp(x);
    } else {
        return detail::scalar_lookup::call_exp(x);
    }
}

template <HasLog T> CANARY_HOST_DEVICE T log(const T& x)
{
    if constexpr (std::floating_point<T>) {
        return std::log(x);
    } else {
        return detail::scalar_lookup::call_log(x);
    }
}

template <HasAbs T> CANARY_HOST_DEVICE T abs(const T& x)
{
    if constexpr (std::floating_point<T>) {
        return std::abs(x);
    } else {
        return detail::scalar_lookup::call_abs(x);
    }
}

template <HasMin T> CANARY_HOST_DEVICE T min(const T& x, const T& y)
{
    if constexpr (std::floating_point<T>) {
        return y < x ? y : x;
    } else {
        return detail::scalar_lookup::call_min(x, y);
    }
}

template <HasMax T> CANARY_HOST_DEVICE T max(const T& x, const T& y)
{
    if constexpr (std::floating_point<T>) {
        return x < y ? y : x;
    } else {
        return detail::scalar_lookup::call_max(x, y);
    }
}

// True when the value part is finite. The derivative part of a dual number
// is not checked.
template <Scalar T> CANARY_HOST_DEVICE bool isfinite(const T& x)
{ return std::isfinite(ScalarTraits<T>::value(x)); }

}  // namespace math

namespace detail
{

// True when v is finite and greater than zero. For IEEE binary32 and
// binary64 the test is one unsigned comparison of the bit pattern: the
// patterns of the positive finite values are exactly 1 up to the pattern of
// the largest finite value, so zero, values with the sign bit set, infinity
// and NaN all fail. Other floating-point types use isfinite and v > 0.
template <std::floating_point V> CANARY_HOST_DEVICE constexpr bool finite_positive(V v)
{
    if constexpr (std::numeric_limits<V>::is_iec559 && (sizeof(V) == 4 || sizeof(V) == 8)) {
        using Bits = std::conditional_t<sizeof(V) == 4, std::uint32_t, std::uint64_t>;
        constexpr Bits largest = std::bit_cast<Bits>(std::numeric_limits<V>::max());
        return std::bit_cast<Bits>(v) - Bits(1) < largest;
    } else {
        return std::isfinite(v) && v > V(0);
    }
}

// True when the value part of every argument is finite and greater than
// zero. Every argument is tested, without early exit, so that the tests
// need no branches.
template <Scalar... S> CANARY_HOST_DEVICE constexpr bool all_finite_positive(const S&... x)
{ return (finite_positive(ScalarTraits<S>::value(x)) & ...); }

}  // namespace detail

}  // namespace canary

#endif  // CANARY_SCALAR_HPP
