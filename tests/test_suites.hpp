#ifndef CANARY_TESTS_TEST_SUITES_HPP
#define CANARY_TESTS_TEST_SUITES_HPP

// Each test source defines one function that registers its Boost.UT tests;
// test_main.cpp calls them in order.

void test_config_default();
void test_config_consumer();
void test_scalar();
void test_quantity();
void test_units();
void test_status();
void test_species();
void test_state();
void test_thermo();
void test_calorically_perfect();
void test_calorically_perfect_mixture();
void test_cerberus_parity();
void test_transport();
void test_solve();
void test_nasa9();
void test_thermally_perfect();
void test_two_temperature();
void test_multi_temperature();
void test_reactions();
void test_energy_exchange();
void test_multi_temperature_kinetics();
void test_mixture_transport();

#endif  // CANARY_TESTS_TEST_SUITES_HPP
