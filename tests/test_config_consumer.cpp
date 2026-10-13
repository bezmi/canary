// The consumer's definitions come before the first canary include.
namespace
{
int failed_assertions = 0;
}  // namespace

#define CANARY_HOST_DEVICE  inline
#define CANARY_FORCE_INLINE inline
#define CANARY_ASSERT(cond) ((cond) ? void() : void(++failed_assertions))

#include "test_suites.hpp"

#include <boost/ut.hpp>
#include <canary/config.hpp>  // IWYU pragma: keep
#include <string_view>

#define CANARY_TEST_STRING(x)        #x
#define CANARY_TEST_EXPAND_STRING(x) CANARY_TEST_STRING(x)

void test_config_consumer()
{
    using namespace boost::ut;

    "a consumer definition of CANARY_HOST_DEVICE is kept"_test = [] {
        expect(std::string_view(CANARY_TEST_EXPAND_STRING(CANARY_HOST_DEVICE)) == "inline");
    };

    "a consumer definition of CANARY_FORCE_INLINE is kept"_test = [] {
        expect(std::string_view(CANARY_TEST_EXPAND_STRING(CANARY_FORCE_INLINE)) == "inline");
    };

    "a consumer definition of CANARY_ASSERT is used"_test = [] {
        failed_assertions = 0;
        CANARY_ASSERT(1 + 1 == 2);
        CANARY_ASSERT(1 + 1 == 3);
        expect(failed_assertions == 1_i);
    };
}
