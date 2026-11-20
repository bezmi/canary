// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
#include "test_suites.hpp"

#include <boost/ut.hpp>
#include <canary/config.hpp>
#include <string_view>

#define CANARY_TEST_STRING(x)        #x
#define CANARY_TEST_EXPAND_STRING(x) CANARY_TEST_STRING(x)

namespace
{
CANARY_HOST_DEVICE int annotated() { return 1; }
CANARY_HOST_DEVICE CANARY_FORCE_INLINE constexpr int force_inlined(int x) { return x + 1; }
}  // namespace

void test_config_default()
{
    using namespace boost::ut;

    "CANARY_HOST_DEVICE is empty by default"_test = [] {
        expect(std::string_view(CANARY_TEST_EXPAND_STRING(CANARY_HOST_DEVICE)).empty());
        expect(annotated() == 1_i);
    };

    "CANARY_FORCE_INLINE requests always-inline by default"_test = [] {
        static_assert(force_inlined(1) == 2);
        const std::string_view expansion = CANARY_TEST_EXPAND_STRING(CANARY_FORCE_INLINE);
        expect(expansion.find("inline") != std::string_view::npos);
#if defined(__GNUC__)
        expect(expansion.find("always_inline") != std::string_view::npos);
#endif
        expect(force_inlined(annotated()) == 2_i);
    };

    "CANARY_ASSERT accepts a true condition by default"_test = [] {
        CANARY_ASSERT(annotated() == 1);  // NOLINT(misc-static-assert): checks the run-time macro
        expect(true);
    };
}
