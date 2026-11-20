// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#ifndef CANARY_DETAIL_SOLVE_HPP
#define CANARY_DETAIL_SOLVE_HPP

// Root finding for closures that have no closed form, such as temperature
// from internal energy. Internal to the models.

#include <algorithm>
#include <canary/config.hpp>
#include <canary/scalar.hpp>
#include <canary/status.hpp>
#include <cmath>
#include <concepts>
#include <numeric>
#include <type_traits>

namespace canary::detail
{

// Value and derivative of a residual f(x).
template <Scalar X> struct Residual {
    X value;
    X derivative;
};

// A residual function of x with value type X.
template <typename F, typename X>
concept ResidualFunction = Scalar<X> && std::regular_invocable<const F&, X> &&
                           std::same_as<std::invoke_result_t<const F&, X>, Residual<X>>;

// Interval containing a root: an endpoint residual is zero, or the residuals differ in sign.
template <std::floating_point V> struct Bracket {
    V lower;
    V upper;
};

template <std::floating_point V> struct SolveSettings {
    V relative_tolerance;  // finite and nonnegative
    std::size_t max_iterations;
    V absolute_tolerance = V(0);  // finite and nonnegative, in the units of x
};

// Causes are retained for the enclosing closure to report without repeating the solve.
// Invalid input checks use the order shown here; residual and convergence reasons
// identify the stage at which the solve stopped.
enum class SolveReason : unsigned char {
    None,
    NonfiniteBounds,
    ReversedBounds,
    NonfiniteGuess,
    InvalidRelativeTolerance,
    InvalidAbsoluteTolerance,
    NonfiniteEndpointResidual,
    NoBracket,
    NonfiniteInteriorResidual,
    IterationLimit,
    PrecisionLimit,
    NonfiniteFullResidual,
    NonfiniteFullDerivative,
    ZeroFullDerivative,
};

template <Scalar T> struct SolveResult {
    T root;
    Status status;       // Valid, NotConverged, or Invalid for bracket or full-residual failures
    SolveReason reason;  // None exactly when status is Valid
};

// Root of f in the bracket, starting from guess (clamped into the bracket).
// Bounds and guess must be finite, lower <= upper, and both tolerances finite
// and nonnegative. Invalid inputs return the guess as a constant with Invalid
// before either residual function is called.
// Exact endpoint roots are accepted before iteration. If both endpoints are
// roots, select the one nearest the clamped guess, with the lower endpoint
// selected on a tie. Equal finite bounds are accepted only when the endpoint
// residual is finite and exactly zero; Dual derivative checks still apply.
//
// The iteration runs on value parts with f_value. Convergence requires an
// exact finite root or a finite sign-changing bracket whose width is at most
// absolute_tolerance + relative_tolerance * scale. The scale is the minimum
// absolute endpoint when the interval excludes zero, and zero otherwise.
// With both tolerances zero, only an exact finite root converges.
// Newton proposals use the endpoint with smaller absolute residual (lower on
// a tie). Candidates must lie strictly inside the bracket. A Newton step that
// fails to halve the bracket is followed by bisection. On bracket convergence,
// return the endpoint with smaller absolute residual.
// Exhausted iterations or representable bracket points give NotConverged with
// the last bounded iterate as a constant; f_full is not called. A nonfinite
// interior residual immediately returns Invalid with the current bounded
// iterate as a constant; f_full is not called.
//
// When T is not the value type (a dual number), f_full is evaluated at
// the converged value as a constant. The returned value is unchanged and
// its input derivatives are -(df/dinput) / (df/dx). Subtracting the value
// part of the full residual retains only its input derivatives; division
// uses the plain value of df/dx. T must support division by V with a result
// convertible to T. Nonfinite full-residual value parts or a zero df/dx
// give Invalid with a constant root. For floating-point T, f_full is not called.
template <Scalar T, std::floating_point V, typename FValue, typename FFull>
    requires std::same_as<typename ScalarTraits<T>::value_type, V> && ResidualFunction<FValue, V> &&
             ResidualFunction<FFull, T> && requires(const T residual, const V derivative) {
                 { residual / derivative } -> std::convertible_to<T>;
             }
CANARY_HOST_DEVICE SolveResult<T> solve(const FValue& f_value,
                                        const FFull& f_full,
                                        Bracket<V> bracket,
                                        V guess,
                                        SolveSettings<V> settings)
{
    const V lower = bracket.lower;
    const V upper = bracket.upper;
    if (!std::isfinite(lower) || !std::isfinite(upper)) {
        return {.root = T(guess),
                .status = Status::Invalid,
                .reason = SolveReason::NonfiniteBounds};
    }
    if (lower > upper) {
        return {.root = T(guess), .status = Status::Invalid, .reason = SolveReason::ReversedBounds};
    }
    if (!std::isfinite(guess)) {
        return {.root = T(guess), .status = Status::Invalid, .reason = SolveReason::NonfiniteGuess};
    }
    if (!std::isfinite(settings.relative_tolerance) || settings.relative_tolerance < V(0)) {
        return {.root = T(guess),
                .status = Status::Invalid,
                .reason = SolveReason::InvalidRelativeTolerance};
    }
    if (!std::isfinite(settings.absolute_tolerance) || settings.absolute_tolerance < V(0)) {
        return {.root = T(guess),
                .status = Status::Invalid,
                .reason = SolveReason::InvalidAbsoluteTolerance};
    }
    const Residual<V> r_lower = f_value(lower);
    const Residual<V> r_upper = f_value(upper);
    const V f_lower = r_lower.value;
    const V f_upper = r_upper.value;
    if (!std::isfinite(f_lower) || !std::isfinite(f_upper)) {
        return {.root = T(guess),
                .status = Status::Invalid,
                .reason = SolveReason::NonfiniteEndpointResidual};
    }
    V x = std::clamp(guess, lower, upper);
    Status status = Status::NotConverged;
    SolveReason reason = SolveReason::IterationLimit;
    if (f_lower == V(0) || f_upper == V(0)) {
        if (f_lower == V(0) && f_upper == V(0)) {
            x = x <= std::midpoint(lower, upper) ? lower : upper;
        } else {
            x = f_lower == V(0) ? lower : upper;
        }
        status = Status::Valid;
        reason = SolveReason::None;
    } else {
        if ((f_lower > V(0)) == (f_upper > V(0))) {
            return {.root = T(guess), .status = Status::Invalid, .reason = SolveReason::NoBracket};
        }
        V low = lower;
        V high = upper;
        Residual<V> left = r_lower;
        Residual<V> right = r_upper;
        V previous_half_width = high / V(2) - low / V(2);
        bool previous_newton = false;

        for (std::size_t i = 0; i < settings.max_iterations; ++i) {
            const Residual<V> r = f_value(x);
            if (!std::isfinite(r.value)) {
                return {.root = T(x),
                        .status = Status::Invalid,
                        .reason = SolveReason::NonfiniteInteriorResidual};
            }
            if (r.value == V(0)) {
                status = Status::Valid;
                reason = SolveReason::None;
                break;
            }
            if ((r.value > V(0)) == (left.value > V(0))) {
                low = x;
                left = r;
            } else {
                high = x;
                right = r;
            }
            const bool prefer_left = std::abs(left.value) <= std::abs(right.value);
            const V best = prefer_left ? low : high;
            const Residual<V>& best_residual = prefer_left ? left : right;
            const V width = high - low;
            const V scale =
              low <= V(0) && high >= V(0) ? V(0) : std::min(std::abs(low), std::abs(high));
            // Width overflow can occur only across zero, where scale is zero
            // and the tolerance remains finite. An overflowing tolerance with
            // a nonzero scale exceeds every possible same-sign bracket width.
            if (width <= settings.absolute_tolerance + settings.relative_tolerance * scale) {
                x = best;
                status = Status::Valid;
                reason = SolveReason::None;
                break;
            }
            V next = std::midpoint(low, high);
            bool use_newton = false;
            if ((!previous_newton || width <= previous_half_width) &&
                std::isfinite(best_residual.derivative) && best_residual.derivative != V(0)) {
                const V newton = best_residual.value / best_residual.derivative;
                V candidate = best - newton;
                // An adjacent-value probe still requires a finite residual
                // and an exact root or a sufficiently narrow bracket.
                if (candidate == best) {
                    candidate = std::nextafter(best, std::signbit(newton) ? high : low);
                }
                if (candidate > low && candidate < high) {
                    next = candidate;
                    use_newton = true;
                }
            }
            if (!(next > low && next < high)) {
                reason = SolveReason::PrecisionLimit;
                break;
            }
            // Half the width stays representable even for opposite extreme bounds.
            previous_half_width = high / V(2) - low / V(2);
            previous_newton = use_newton;
            x = next;
        }
    }

    if constexpr (std::same_as<T, V>) {
        return {.root = x, .status = status, .reason = reason};
    } else {
        if (status == Status::NotConverged) {
            return {.root = T(x), .status = status, .reason = reason};
        }
        const Residual<T> r = f_full(T(x));
        const V value = ScalarTraits<T>::value(r.value);
        const V derivative = ScalarTraits<T>::value(r.derivative);
        if (!std::isfinite(value)) {
            return {.root = T(x),
                    .status = Status::Invalid,
                    .reason = SolveReason::NonfiniteFullResidual};
        }
        if (!std::isfinite(derivative)) {
            return {.root = T(x),
                    .status = Status::Invalid,
                    .reason = SolveReason::NonfiniteFullDerivative};
        }
        if (derivative == V(0)) {
            return {.root = T(x),
                    .status = Status::Invalid,
                    .reason = SolveReason::ZeroFullDerivative};
        }
        const T input_residual = r.value - T(value);
        return {.root = T(x) - T(input_residual / derivative), .status = status, .reason = reason};
    }
}

}  // namespace canary::detail

#endif  // CANARY_DETAIL_SOLVE_HPP
