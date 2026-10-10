#include <boost/ut.hpp>

int main()
{
    using namespace boost::ut;

    "compiled as C++20"_test = [] { expect(__cplusplus == 202002L); };
}
