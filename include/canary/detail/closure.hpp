// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
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

// The check with the more severe status; a when the severities are equal, so
// the first diagnostic is preserved. Every field is selected on its own: a
// choice between the two whole records compiles to a copy through local
// memory in device code, which this keeps out of the kernels.
CANARY_HOST_DEVICE constexpr ClosureCheck combine(ClosureCheck a, ClosureCheck b) noexcept
{
    const bool take_b = a.status < b.status;
    return {
      .status = take_b ? b.status : a.status,
      .diagnostic =
        {
          .reason = take_b ? b.diagnostic.reason : a.diagnostic.reason,
          .mode = take_b ? b.diagnostic.mode : a.diagnostic.mode,
        },
    };
}

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

// Valid for a finite value part, otherwise Invalid with the reason nonfinite.
// Field-wise selection, as in combine.
template <Scalar T>
CANARY_HOST_DEVICE constexpr ClosureCheck finite_check(T value, ClosureReason nonfinite) noexcept
{
    const bool finite = math::isfinite(value);
    return {
      .status = finite ? Status::Valid : Status::Invalid,
      .diagnostic = {.reason = finite ? ClosureReason::None : nonfinite, .mode = {}},
    };
}

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
