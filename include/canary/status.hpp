#ifndef CANARY_STATUS_HPP
#define CANARY_STATUS_HPP

#include <canary/config.hpp>

namespace canary
{

// Outcome of a per-cell calculation. Enumerators are in order of increasing
// severity:
// - Valid: inputs in range, result computed.
// - CompositionViolation: species weights outside [0, total] or not summing
//   to the total within tolerance; the result uses clamped, normalised
//   weights.
// - NotConverged: an iterative solve stopped at its iteration limit; the
//   result is the last iterate.
// - Invalid: an input or result is not finite or not positive where it must
//   be; the result must not be used.
enum class Status : unsigned char {
    Valid,
    CompositionViolation,
    NotConverged,
    Invalid,
};

// The more severe of two statuses. Commutative and associative, with Valid
// as identity, so statuses of several steps combine in any order.
CANARY_HOST_DEVICE constexpr Status combine(Status a, Status b) noexcept { return a < b ? b : a; }

}  // namespace canary

#endif  // CANARY_STATUS_HPP
