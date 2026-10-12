#include "test_suites.hpp"

#include <boost/ut.hpp>
#include <gaslib/config.hpp>
#include <string_view>

#define GASLIB_TEST_STRING(x)        #x
#define GASLIB_TEST_EXPAND_STRING(x) GASLIB_TEST_STRING(x)

namespace
{
GASLIB_HOST_DEVICE int annotated() { return 1; }
GASLIB_HOST_DEVICE GASLIB_FORCE_INLINE constexpr int force_inlined(int x) { return x + 1; }
}  // namespace

void test_config_default()
{
    using namespace boost::ut;

    "GASLIB_HOST_DEVICE is empty by default"_test = [] {
        expect(std::string_view(GASLIB_TEST_EXPAND_STRING(GASLIB_HOST_DEVICE)).empty());
        expect(annotated() == 1_i);
    };

    "GASLIB_FORCE_INLINE requests always-inline by default"_test = [] {
        static_assert(force_inlined(1) == 2);
        const std::string_view expansion = GASLIB_TEST_EXPAND_STRING(GASLIB_FORCE_INLINE);
        expect(expansion.find("inline") != std::string_view::npos);
#if defined(__GNUC__)
        expect(expansion.find("always_inline") != std::string_view::npos);
#endif
        expect(force_inlined(annotated()) == 2_i);
    };

    "GASLIB_ASSERT accepts a true condition by default"_test = [] {
        GASLIB_ASSERT(annotated() == 1);  // NOLINT(misc-static-assert): checks the run-time macro
        expect(true);
    };
}
