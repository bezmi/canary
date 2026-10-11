// Gupta-Yos transport is a transport model of the two-temperature mixture,
// with one conductivity per mode, but not of the calorically perfect
// mixture.
#include <gaslib/data/air5.hpp>
#include <gaslib/models/calorically_perfect.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/transport/transport.hpp>

#ifdef GASLIB_COMPILE_FAIL
static_assert(gaslib::TransportModel<gaslib::Air5GuptaYosTransport<double, gaslib::si>,
                                     gaslib::CaloricallyPerfectMixture<double, gaslib::si, 5>,
                                     double>);
#endif
static_assert(gaslib::TransportModel<gaslib::Air5GuptaYosTransport<double, gaslib::si>,
                                     gaslib::Air5TwoTemperature<double, gaslib::si>,
                                     double>);
