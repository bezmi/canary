#ifndef GASLIB_DETAIL_SOLVE_HPP
#define GASLIB_DETAIL_SOLVE_HPP

// Root finding for closures that have no closed form, such as temperature
// from internal energy. Internal to the models.

#include <algorithm>
#include <cmath>
#include <concepts>
#include <gaslib/config.hpp>
#include <gaslib/scalar.hpp>
#include <gaslib/status.hpp>
#include <type_traits>

namespace gaslib::detail
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

// Interval known to contain the root: f(lower) and f(upper) differ in sign.
template <std::floating_point V> struct Bracket {
    V lower;
    V upper;
};

template <std::floating_point V> struct SolveSettings {
    V relative_tolerance;  // on the step, relative to |x|
    int max_iterations;
};

template <Scalar T> struct SolveResult {
    T root;
    Status status;  // Valid, NotConverged, or Invalid for a bracket without a sign change
};

// Root of f in the bracket, starting from guess (clamped into the bracket).
//
// The iteration runs on value parts with f_value. When the Newton step
// f(x) / f'(x) is at most relative_tolerance * |x|, it takes that step and
// stops; otherwise it takes the Newton step, replaced by bisection when the
// step leaves the bracket, and also stops when that step is at most
// relative_tolerance * |x|. The bracket shrinks to the side where f changes sign. After
// max_iterations it stops with status NotConverged.
//
// When T is not the value type (a dual number), one final Newton step is
// then taken in T with f_full, starting from the converged value as a
// constant: x = x0 - f(x0) / f'(x0). For a floating-point T the converged
// value is returned and f_full is not called.
// With dual-number inputs captured in f_full, the derivative of the root is
// -(df/dinput) / (df/dx), the implicit-function-theorem derivative, and the
// value changes by at most one more Newton step.
template <Scalar T, std::floating_point V, typename FValue, typename FFull>
    requires std::same_as<typename ScalarTraits<T>::value_type, V> && ResidualFunction<FValue, V> &&
             ResidualFunction<FFull, T>
GASLIB_HOST_DEVICE SolveResult<T> solve(const FValue& f_value,
                                        const FFull& f_full,
                                        Bracket<V> bracket,
                                        V guess,
                                        SolveSettings<V> settings)
{
    V lower = bracket.lower;
    V upper = bracket.upper;
    const V f_lower = f_value(lower).value;
    const V f_upper = f_value(upper).value;
    if (!std::isfinite(f_lower) || !std::isfinite(f_upper) ||
        (f_lower > V(0)) == (f_upper > V(0))) {
        return {.root = T(guess), .status = Status::Invalid};
    }
    // Orient the bracket so that f(negative_side) <= 0 <= f(positive_side).
    V negative_side = f_lower > V(0) ? upper : lower;
    V positive_side = f_lower > V(0) ? lower : upper;

    V x = std::clamp(guess, lower, upper);
    Status status = Status::NotConverged;
    for (int i = 0; i < settings.max_iterations; ++i) {
        const Residual<V> r = f_value(x);
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
        // A Newton step below the tolerance is taken without the bracket
        // test and ends the iteration, also when it rounds to x itself,
        // which the bracket test would reject.
        const V newton = r.value / r.derivative;
        if (std::abs(newton) <= settings.relative_tolerance * std::abs(x)) {
            x -= newton;
            status = Status::Valid;
            break;
        }
        V next = x - newton;
        if (!(next > low && next < high)) { next = V(0.5) * (low + high); }
        const V step = next - x;
        x = next;
        if (std::abs(step) <= settings.relative_tolerance * std::abs(x)) {
            status = Status::Valid;
            break;
        }
    }

    if constexpr (std::same_as<T, V>) {
        return {.root = x, .status = status};
    } else {
        const T x0(x);
        const Residual<T> r = f_full(x0);
        return {.root = x0 - (r.value / r.derivative), .status = status};
    }
}

}  // namespace gaslib::detail

#endif  // GASLIB_DETAIL_SOLVE_HPP
