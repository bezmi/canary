# Test targets are built as strict C++20 so that a header needing a later
# standard or a compiler extension fails here.
function(gaslib_test_options target)
  set_target_properties(${target} PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
endfunction()

add_executable(gaslib_tests test_main.cpp test_config_default.cpp test_config_consumer.cpp test_scalar.cpp test_quantity.cpp test_units.cpp test_status.cpp test_species.cpp test_state.cpp test_thermo.cpp test_calorically_perfect.cpp test_calorically_perfect_mixture.cpp)
target_link_libraries(gaslib_tests PRIVATE gaslib::gaslib Boost::ut)
gaslib_test_options(gaslib_tests)
add_test(NAME gaslib.unit COMMAND gaslib_tests)

# Code that must not compile. compile_fail/<name>.cpp is compiled twice: as
# written, as part of the normal build, where it must compile; and with
# GASLIB_COMPILE_FAIL defined, by the test, which passes only if the build
# output matches <regex>.
function(gaslib_compile_fail_test name regex)
  add_library(compile_ok_${name} OBJECT compile_fail/${name}.cpp)
  add_library(compile_fail_${name} OBJECT EXCLUDE_FROM_ALL compile_fail/${name}.cpp)
  foreach(target compile_ok_${name} compile_fail_${name})
    target_link_libraries(${target} PRIVATE gaslib::gaslib)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_LIST_DIR})
    gaslib_test_options(${target})
  endforeach()
  target_compile_definitions(compile_fail_${name} PRIVATE GASLIB_COMPILE_FAIL)
  add_test(NAME gaslib.compile_fail.${name}
           COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target compile_fail_${name})
  set_tests_properties(gaslib.compile_fail.${name} PROPERTIES PASS_REGULAR_EXPRESSION "${regex}")
endfunction()

gaslib_compile_fail_test(self_test "compile-fail self test")
gaslib_compile_fail_test(scalar_pointer "satisf(y|action of) [^ ]*Scalar")
gaslib_compile_fail_test(scalar_without_traits "satisf(y|action of) [^ ]*Scalar")
gaslib_compile_fail_test(math_function_missing "satisf(y|action of) [^ ]*HasExp")
gaslib_compile_fail_test(quantity_add_other_kind "(no match for [^ ]*operator\\+|invalid operands to binary expression).*Pressure.*Temperature")
gaslib_compile_fail_test(quantity_add_other_units "(no match for [^ ]*operator\\+|invalid operands to binary expression).*Density<double, gaslib::si>.*Density<double, gaslib::scaled>")
gaslib_compile_fail_test(quantity_product_without_relation "(no match for [^ ]*operator\\*|invalid operands to binary expression).*Density.*Temperature")
gaslib_compile_fail_test(convert_other_kind "conversion from .*energy_density.* to .*pressure|no viable conversion from .*energy_density.* to .*pressure")
gaslib_compile_fail_test(species_density_as_fractions "SpeciesRangeOf<.*kind::mass_fraction")
gaslib_compile_fail_test(thermo_density_density "ClosesInAnyOrder<.*IdealGas")
gaslib_compile_fail_test(thermo_gas_model_missing_closure "GasModel<.*NoPressureTemperatureGas")
gaslib_compile_fail_test(thermo_missing_extra_input "ClosesInAnyOrder<.*GuessRequiredGas")
gaslib_compile_fail_test(quantity_implicit_precision "conversion from .*double.* to .*float|no viable conversion from .*double.* to .*float")
gaslib_compile_fail_test(calorically_perfect_si_in_scaled "ClosesInAnyOrder<.*CaloricallyPerfect<double, gaslib::scaled>")
gaslib_compile_fail_test(calorically_perfect_swapped_kinds "no matching function for call to [^ ]*make_calorically_perfect")
gaslib_compile_fail_test(quantity_molar_scaled "KindUnits<(K, U|gaslib::kind::molar_mass, gaslib::scaled)>")

# Each public header is compiled alone, as strict C++20, with only the
# include directories of gaslib::gaslib (include/ and mdspan). This checks
# that every header includes what it uses and needs nothing newer than C++20.
file(GLOB_RECURSE gaslib_headers CONFIGURE_DEPENDS RELATIVE ${PROJECT_SOURCE_DIR}/include
     ${PROJECT_SOURCE_DIR}/include/gaslib/*.hpp)
set(header_check_sources)
foreach(header IN LISTS gaslib_headers)
  string(MAKE_C_IDENTIFIER ${header} source_name)
  set(source ${CMAKE_CURRENT_BINARY_DIR}/header_check/${source_name}.cpp)
  file(CONFIGURE OUTPUT ${source} CONTENT "#include <${header}>\n")
  list(APPEND header_check_sources ${source})
endforeach()
add_library(gaslib_header_check OBJECT ${header_check_sources})
target_link_libraries(gaslib_header_check PRIVATE gaslib::gaslib)
gaslib_test_options(gaslib_header_check)
