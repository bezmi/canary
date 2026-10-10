#include "test_suites.hpp"

#include <boost/ut.hpp>

int main()
{
    using namespace boost::ut;

    "compiled as C++20"_test = [] { expect(__cplusplus == 202002L); };

    test_config_default();
    test_config_consumer();
    test_scalar();
    test_quantity();
    test_units();
    test_status();
}
