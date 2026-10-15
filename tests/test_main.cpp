#include "test_suites.hpp"

#include <boost/ut.hpp>

int main()
{
    using namespace boost::ut;

    "compiled as C++20 or later"_test = [] { expect(__cplusplus >= 202002L); };

    test_config_default();
    test_config_consumer();
    test_scalar();
    test_quantity();
    test_units();
    test_status();
    test_indexed_view();
    test_species();
    test_state();
    test_thermo();
    test_calorically_perfect();
    test_calorically_perfect_mixture();
    test_cerberus_parity();
    test_transport();
    test_solve();
    test_nasa9();
    test_thermally_perfect();
    test_two_temperature();
    test_multi_temperature();
    test_reactions();
    test_energy_exchange();
    test_multi_temperature_kinetics();
    test_mixture_transport();
}
