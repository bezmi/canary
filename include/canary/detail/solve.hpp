#ifndef CANARY_DETAIL_SOLVE_HPP
#define CANARY_DETAIL_SOLVE_HPP

// Root finding for closures that have no closed form, such as temperature
// from internal energy. Internal to the models.

#include <algorithm>
#include <cmath>
#include <concepts>
#include <canary/config.hpp>
#include <canary/scalar.hpp>
#include <canary/status.hpp>
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
    V relative_tolerance;  // finite and nonnegative; zero requires an exact root
    int max_iterations;
};

template <Scalar T> struct SolveResult {
    T root;
    Status status;  // Valid, NotConverged, or Invalid for bracket or full-residual failures
};

// Root of f in the bracket, starting from guess (clamped into the bracket).
// Bounds and guess must be finite, lower <= upper, relative_tolerance finite
// and nonnegative, and max_iterations nonnegative. Invalid inputs return the
// guess as a constant with Invalid before either residual function is called.
// Exact endpoint roots are accepted before iteration. If both endpoints are
// roots, select the one nearest the clamped guess, with the lower endpoint
// selected on a tie. Equal finite bounds are accepted only when the endpoint
// residual is finite and exactly zero; Dual derivative checks still apply.
//
// With zero relative_tolerance, only an exact finite root converges; exhausted
// iterations or representable bracket points give NotConverged.
// The iteration runs on value parts with f_value. When the Newton step
// f(x) / f'(x) is at most relative_tolerance * |x|, it takes that step and
// stops; otherwise it takes the Newton step, replaced by bisection when the
// step leaves the bracket or the derivative is zero or nonfinite, and also
// stops when that step is at most
// relative_tolerance * |x|. The bracket shrinks to the side where f changes sign. After
// max_iterations it stops with status NotConverged and returns the last
// bounded iterate as a constant with zero derivative; f_full is not called.
// A nonfinite interior residual immediately returns Invalid with the current
// bounded iterate as a constant; f_full is not called.
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
    if (!std::isfinite(lower) || !std::isfinite(upper) || lower > upper || !std::isfinite(guess) ||
        !std::isfinite(settings.relative_tolerance) || settings.relative_tolerance < V(0) ||
        settings.max_iterations < 0) {
        return {.root = T(guess), .status = Status::Invalid};
    }
    const V f_lower = f_value(lower).value;
    const V f_upper = f_value(upper).value;
    if (!std::isfinite(f_lower) || !std::isfinite(f_upper)) {
        return {.root = T(guess), .status = Status::Invalid};
    }
    V x = std::clamp(guess, lower, upper);
    Status status = Status::NotConverged;
    if (f_lower == V(0) || f_upper == V(0)) {
        if (f_lower == V(0) && f_upper == V(0)) {
            x = x <= std::midpoint(lower, upper) ? lower : upper;
        } else {
            x = f_lower == V(0) ? lower : upper;
        }
        status = Status::Valid;
    } else {
        if ((f_lower > V(0)) == (f_upper > V(0))) {
            return {.root = T(guess), .status = Status::Invalid};
        }
        // Orient the bracket so that f(negative_side) <= 0 <= f(positive_side).
        V negative_side = f_lower > V(0) ? upper : lower;
        V positive_side = f_lower > V(0) ? lower : upper;

        for (int i = 0; i < settings.max_iterations; ++i) {
            const Residual<V> r = f_value(x);
            if (!std::isfinite(r.value)) { return {.root = T(x), .status = Status::Invalid}; }
            if (r.value == V(0)) {
                status = Status::Valid;
                break;
            }
            if (r.value < V(0)) {
                negative_side = x;
            } else {
                positive_side = x;
            }
            const V low = negative_side < positive_side ? negative_side : positive_side;
            const V high = negative_side < positive_side ? positive_side : negative_side;
            V next = std::midpoint(low, high);
            if (std::isfinite(r.derivative) && r.derivative != V(0)) {
                // A Newton step below the tolerance is taken without the bracket
                // test and ends the iteration, also when it rounds to x itself,
                // which the bracket test would reject.
                const V newton = r.value / r.derivative;
                if (settings.relative_tolerance > V(0) &&
                    std::abs(newton) <= settings.relative_tolerance * std::abs(x)) {
                    x -= newton;
                    status = Status::Valid;
                    break;
                }
                const V candidate = x - newton;
                if (candidate > low && candidate < high) { next = candidate; }
            }
            if (settings.relative_tolerance == V(0) && !(next > low && next < high)) { break; }
            const V step = next - x;
            x = next;
            if (settings.relative_tolerance > V(0) &&
                std::abs(step) <= settings.relative_tolerance * std::abs(x)) {
                status = Status::Valid;
                break;
            }
        }
    }

    if constexpr (std::same_as<T, V>) {
        return {.root = x, .status = status};
    } else {
        if (status == Status::NotConverged) { return {.root = T(x), .status = status}; }
        const Residual<T> r = f_full(T(x));
        const V value = ScalarTraits<T>::value(r.value);
        const V derivative = ScalarTraits<T>::value(r.derivative);
        if (!std::isfinite(value) || !std::isfinite(derivative) || derivative == V(0)) {
            return {.root = T(x), .status = Status::Invalid};
        }
        const T input_residual = r.value - T(value);
        return {.root = T(x) - T(input_residual / derivative), .status = status};
    }
}

}  // namespace canary::detail

#endif  // CANARY_DETAIL_SOLVE_HPP
