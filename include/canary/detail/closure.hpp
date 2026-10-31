#ifndef CANARY_DETAIL_CLOSURE_HPP
#define CANARY_DETAIL_CLOSURE_HPP

#include <canary/detail/solve.hpp>
#include <canary/species.hpp>
#include <canary/state.hpp>

namespace canary::detail
{

struct ClosureCheck {
    Status status {Status::Valid};
    ClosureDiagnostic diagnostic {};
};

// Preserve the first diagnostic when severities are equal.
CANARY_HOST_DEVICE constexpr ClosureCheck combine(ClosureCheck a, ClosureCheck b) noexcept
{ return a.status < b.status ? b : a; }

CANARY_HOST_DEVICE constexpr ClosureCheck invalid_closure(ClosureReason reason) noexcept
{ return {Status::Invalid, {reason, {}}}; }

CANARY_HOST_DEVICE constexpr ClosureCheck for_mode(ClosureCheck check, ModeIndex mode) noexcept
{
    if (check.status != Status::Valid) { check.diagnostic.mode = OptionalModeIndex {mode}; }
    return check;
}

template <typename Mixture>
CANARY_HOST_DEVICE constexpr ClosureCheck composition_check(const Mixture& mix) noexcept
{
    ClosureReason reason = ClosureReason::None;
    switch (mix.reason) {
    case CompositionReason::None: break;
    case CompositionReason::SpeciesCountOverflow:
        reason = ClosureReason::CompositionSpeciesCountOverflow;
        break;
    case CompositionReason::NonfiniteTotal:
        reason = ClosureReason::CompositionNonfiniteTotal;
        break;
    case CompositionReason::NonpositiveTotal:
        reason = ClosureReason::CompositionNonpositiveTotal;
        break;
    case CompositionReason::NonfiniteWeight:
        reason = ClosureReason::CompositionNonfiniteWeight;
        break;
    case CompositionReason::NoPositiveWeight:
        reason = ClosureReason::CompositionNoPositiveWeight;
        break;
    case CompositionReason::WeightOutOfRange:
        reason = ClosureReason::CompositionWeightOutOfRange;
        break;
    case CompositionReason::SumMismatch: reason = ClosureReason::CompositionSumMismatch; break;
    }
    return {mix.status, {reason, {}}};
}

template <Scalar T>
CANARY_HOST_DEVICE constexpr ClosureCheck solve_check(const SolveResult<T>& result) noexcept
{
    ClosureReason reason = ClosureReason::None;
    switch (result.reason) {
    case SolveReason::None: break;
    case SolveReason::NonfiniteBounds: reason = ClosureReason::NonfiniteTemperatureBounds; break;
    case SolveReason::ReversedBounds: reason = ClosureReason::ReversedTemperatureBounds; break;
    case SolveReason::NonfiniteGuess: reason = ClosureReason::NonfiniteTemperatureGuess; break;
    case SolveReason::InvalidRelativeTolerance:
        reason = ClosureReason::InvalidRelativeTolerance;
        break;
    case SolveReason::InvalidAbsoluteTolerance:
        reason = ClosureReason::InvalidAbsoluteTolerance;
        break;
    case SolveReason::NonfiniteEndpointResidual:
        reason = ClosureReason::NonfiniteEndpointResidual;
        break;
    case SolveReason::NoBracket: reason = ClosureReason::NoTemperatureBracket; break;
    case SolveReason::NonfiniteInteriorResidual:
        reason = ClosureReason::NonfiniteInteriorResidual;
        break;
    case SolveReason::IterationLimit: reason = ClosureReason::IterationLimit; break;
    case SolveReason::PrecisionLimit: reason = ClosureReason::PrecisionLimit; break;
    case SolveReason::NonfiniteFullResidual: reason = ClosureReason::NonfiniteFullResidual; break;
    case SolveReason::NonfiniteFullDerivative:
        reason = ClosureReason::NonfiniteFullDerivative;
        break;
    case SolveReason::ZeroFullDerivative: reason = ClosureReason::ZeroFullDerivative; break;
    }
    return {result.status, {reason, {}}};
}

template <Scalar T>
CANARY_HOST_DEVICE constexpr ClosureCheck finite_check(T value, ClosureReason nonfinite) noexcept
{ return math::isfinite(value) ? ClosureCheck {} : invalid_closure(nonfinite); }

template <Scalar T>
CANARY_HOST_DEVICE constexpr ClosureCheck
  positive_check(T value, ClosureReason nonfinite, ClosureReason nonpositive) noexcept
{
    if (!math::isfinite(value)) { return invalid_closure(nonfinite); }
    if (!(ScalarTraits<T>::value(value) > typename ScalarTraits<T>::value_type(0))) {
        return invalid_closure(nonpositive);
    }
    return {};
}

}  // namespace canary::detail

#endif  // CANARY_DETAIL_CLOSURE_HPP
