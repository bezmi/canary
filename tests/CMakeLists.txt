# Test targets are built as strict C++20 so that a header needing a later
# standard or a compiler extension fails here.
function(gaslib_test_options target)
  set_target_properties(${target} PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
endfunction()

add_executable(gaslib_tests test_main.cpp test_config_default.cpp test_config_consumer.cpp test_scalar.cpp test_quantity.cpp test_units.cpp test_status.cpp test_species.cpp test_state.cpp test_thermo.cpp test_calorically_perfect.cpp test_calorically_perfect_mixture.cpp
                            test_cerberus_parity.cpp test_transport.cpp test_solve.cpp test_nasa9.cpp test_thermally_perfect.cpp
                            test_two_temperature.cpp test_reactions.cpp test_energy_exchange.cpp test_mixture_transport.cpp)
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
gaslib_compile_fail_test(thermally_perfect_si_in_scaled "ClosesInAnyOrder<.*ThermallyPerfectMixture<double, gaslib::scaled")
gaslib_compile_fail_test(two_temperature_si_modes_in_scaled "ClosesInAnyOrder<.*TwoTemperatureMixture<double, gaslib::scaled")
gaslib_compile_fail_test(calorically_perfect_swapped_kinds "no matching function for call to [^ ]*make_calorically_perfect")
gaslib_compile_fail_test(transport_other_units "no matching function for call to [^ ]*viscosity")
gaslib_compile_fail_test(quantity_molar_scaled "KindUnits<(K, U|gaslib::kind::molar_mass, gaslib::scaled)>")
gaslib_compile_fail_test(energy_exchange_single_temperature "RelaxingModel<.*ThermallyPerfectMixture")
gaslib_compile_fail_test(energy_exchange_mode_energy_missing "HasVibrationalElectronicEnergy<.*ThermallyPerfectMixture")
gaslib_compile_fail_test(mixture_transport_composition "HasNormalisedComposition<.*CaloricallyPerfectMixture")
gaslib_compile_fail_test(mixture_transport_heat_capacity "HasSpeciesIsobaricSpecificHeat<.*CaloricallyPerfectMixture")
gaslib_compile_fail_test(mixture_transport_model "TransportModel<.*GuptaYosTransport")

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

# Species database library: C++ tests (the public header is C++20) and the
# Python converter and schema checks, when the interpreter and its yaml and
# jsonschema modules are available.
if(TARGET gaslib_database)
  add_executable(gaslib_database_tests test_database.cpp)
  target_link_libraries(gaslib_database_tests PRIVATE gaslib::database Boost::ut)
  target_compile_definitions(gaslib_database_tests PRIVATE GASLIB_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  gaslib_test_options(gaslib_database_tests)
  add_test(NAME gaslib.database COMMAND gaslib_database_tests)

  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_FOUND)
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import yaml, jsonschema" RESULT_VARIABLE python_modules_missing
                    OUTPUT_QUIET ERROR_QUIET)
    if(NOT python_modules_missing)
      set(species_tools ${PROJECT_SOURCE_DIR}/tools/species_database)
      add_test(NAME gaslib.database.converters
               COMMAND ${Python3_EXECUTABLE} ${species_tools}/test_converters.py ${PROJECT_SOURCE_DIR})
      add_test(NAME gaslib.database.schema
               COMMAND ${Python3_EXECUTABLE} ${species_tools}/check_species_files.py
                       ${PROJECT_SOURCE_DIR}/database/species.schema.json ${PROJECT_SOURCE_DIR}/database/data/air5.json --
                       ${PROJECT_SOURCE_DIR}/database/reactions.schema.json
                       ${PROJECT_SOURCE_DIR}/database/data/air5_kim_jo.json
                       ${PROJECT_SOURCE_DIR}/database/data/air5_park.json --
                       ${PROJECT_SOURCE_DIR}/database/energy_exchange.schema.json
                       ${PROJECT_SOURCE_DIR}/database/data/air5_park_exchange.json --
                       ${PROJECT_SOURCE_DIR}/database/collision_integrals.schema.json
                       ${PROJECT_SOURCE_DIR}/database/data/air5_collision_integrals.json)
    else()
      message(STATUS "gaslib: Python modules yaml or jsonschema missing; species converter tests skipped")
    endif()
  endif()
endif()
