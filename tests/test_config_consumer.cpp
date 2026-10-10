// The consumer's definitions come before the first gaslib include.
namespace
{
int failed_assertions = 0;
}  // namespace

#define GASLIB_HOST_DEVICE  inline
#define GASLIB_ASSERT(cond) ((cond) ? void() : void(++failed_assertions))

#include "test_suites.hpp"

#include <boost/ut.hpp>
#include <gaslib/config.hpp>  // IWYU pragma: keep
#include <string_view>

#define GASLIB_TEST_STRING(x)        #x
#define GASLIB_TEST_EXPAND_STRING(x) GASLIB_TEST_STRING(x)

void test_config_consumer()
{
    using namespace boost::ut;

    "a consumer definition of GASLIB_HOST_DEVICE is kept"_test = [] {
        expect(std::string_view(GASLIB_TEST_EXPAND_STRING(GASLIB_HOST_DEVICE)) == "inline");
    };

    "a consumer definition of GASLIB_ASSERT is used"_test = [] {
        failed_assertions = 0;
        GASLIB_ASSERT(1 + 1 == 2);
        GASLIB_ASSERT(1 + 1 == 3);
        expect(failed_assertions == 1_i);
    };
}
