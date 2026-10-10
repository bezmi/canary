#include <boost/ut.hpp>

void test_config_default();
void test_config_consumer();

int main()
{
    using namespace boost::ut;

    "compiled as C++20"_test = [] { expect(__cplusplus == 202002L); };

    test_config_default();
    test_config_consumer();
}
