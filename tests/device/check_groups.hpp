#ifndef CANARY_TESTS_DEVICE_CHECK_GROUPS_HPP
#define CANARY_TESTS_DEVICE_CHECK_GROUPS_HPP

// The groups of numerical checks. Each group evaluates its checks for one
// parameter type P, state scalar T and unit system U and returns their report.
// A generated unit defines and explicitly instantiates one group (see
// check.in); the runtime executable calls the instances through these
// declarations. Define CANARY_HOST_DEVICE before including this header.

#include "../dual.hpp"
#include "runtime_checks.hpp"

#include <canary/units.hpp>
#include <concepts>

namespace canary::test
{

template <std::floating_point P, typename T, UnitSystem U> Report thermodynamic_single_checks();
template <std::floating_point P, typename T, UnitSystem U> Report thermodynamic_modal_checks();
template <std::floating_point P, typename T, UnitSystem U> Report thermodynamic_mixed_checks();
template <std::floating_point P, typename T, UnitSystem U> Report reaction_checks();
template <std::floating_point P, typename T, UnitSystem U> Report reaction_mixed_checks();
template <std::floating_point P, typename T, UnitSystem U> Report exchange_checks();
template <std::floating_point P, typename T, UnitSystem U> Report transport_checks();
template <std::floating_point P, typename T, UnitSystem U> Report transport_structural_checks();

}  // namespace canary::test

#endif  // CANARY_TESTS_DEVICE_CHECK_GROUPS_HPP
