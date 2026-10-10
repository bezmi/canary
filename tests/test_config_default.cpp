#include <boost/ut.hpp>
#include <gaslib/config.hpp>
#include <string_view>

#define GASLIB_TEST_STRING(x)        #x
#define GASLIB_TEST_EXPAND_STRING(x) GASLIB_TEST_STRING(x)

namespace
{
GASLIB_HOST_DEVICE int annotated() { return 1; }
}  // namespace

void test_config_default()
{
    using namespace boost::ut;

    "GASLIB_HOST_DEVICE is empty by default"_test = [] {
        expect(std::string_view(GASLIB_TEST_EXPAND_STRING(GASLIB_HOST_DEVICE)).empty());
        expect(annotated() == 1_i);
    };

    "GASLIB_ASSERT accepts a true condition by default"_test = [] {
        GASLIB_ASSERT(annotated() == 1);
        expect(true);
    };
}
