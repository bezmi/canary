#ifndef GASLIB_TESTS_TEST_SUITES_HPP
#define GASLIB_TESTS_TEST_SUITES_HPP

// Each test source defines one function that registers its Boost.UT tests;
// test_main.cpp calls them in order.

void test_config_default();
void test_config_consumer();
void test_scalar();
void test_quantity();
void test_units();

#endif  // GASLIB_TESTS_TEST_SUITES_HPP
