# Numerical tests use the requested language standard, with C++20 as the minimum.
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
elseif(CMAKE_CXX_STANDARD LESS 20)
  message(FATAL_ERROR "canary tests require CMAKE_CXX_STANDARD >= 20")
endif()

function(canary_test_options target)
  set_target_properties(${target} PROPERTIES CXX_STANDARD ${CMAKE_CXX_STANDARD}
                        CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
  if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC" OR MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX /permissive-)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic -Werror)
  else()
    message(WARNING "canary: strict warning flags are not configured for ${CMAKE_CXX_COMPILER_ID}")
  endif()
endfunction()

add_executable(canary_tests test_main.cpp test_config_default.cpp test_config_consumer.cpp test_scalar.cpp test_quantity.cpp test_units.cpp test_status.cpp test_indexed_view.cpp test_species.cpp test_state.cpp test_thermo.cpp test_calorically_perfect.cpp test_calorically_perfect_mixture.cpp
                            test_cerberus_parity.cpp test_transport.cpp test_solve.cpp test_nasa9.cpp test_thermally_perfect.cpp
                            test_two_temperature.cpp test_multi_temperature.cpp test_reactions.cpp test_energy_exchange.cpp
                            test_multi_temperature_kinetics.cpp test_mixture_transport.cpp test_layout.cpp)
target_link_libraries(canary_tests PRIVATE canary::canary Boost::ut)
canary_test_options(canary_tests)
add_test(NAME canary.unit COMMAND canary_tests)

# Extent preconditions are checked in isolated processes before storage access.
add_executable(canary_mode_extent_tests test_mode_extent.cpp)
target_link_libraries(canary_mode_extent_tests PRIVATE canary::canary)
canary_test_options(canary_mode_extent_tests)
foreach(operation closure rates production bound_rates bound_production)
  foreach(view raw typed)
    add_test(NAME canary.mode_extent.${operation}.${view}
             COMMAND canary_mode_extent_tests ${operation} ${view})
  endforeach()
endforeach()

add_test(NAME canary.mode_extent.bound_scalar
         COMMAND canary_mode_extent_tests bound_scalar raw)

foreach(operation state_vt state_cv state_et state_ve state_sources)
  foreach(field energy temperature)
    add_test(NAME canary.mode_extent.${operation}.${field}
             COMMAND canary_mode_extent_tests ${operation} ${field})
  endforeach()
endforeach()

# Code that must not compile. compile_fail/<name>.cpp is compiled twice: as
# written, as part of the normal build, where it must compile; and with
# CANARY_COMPILE_FAIL defined, by the test, which passes only if the build
# fails and its diagnostic matches <regex>. Shared build operations are serialized.
function(canary_compile_fail_test name regex)
  add_library(compile_ok_${name} OBJECT compile_fail/${name}.cpp)
  add_library(compile_fail_${name} OBJECT EXCLUDE_FROM_ALL compile_fail/${name}.cpp)
  foreach(target compile_ok_${name} compile_fail_${name})
    target_link_libraries(${target} PRIVATE canary::canary)
    target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_LIST_DIR})
    canary_test_options(${target})
  endforeach()
  target_compile_definitions(compile_fail_${name} PRIVATE CANARY_COMPILE_FAIL)
  add_test(NAME canary.compile_fail.${name}
           COMMAND ${CMAKE_COMMAND}
                   "-DBUILD_DIR=${CMAKE_BINARY_DIR}"
                   "-DBUILD_TARGET=compile_fail_${name}"
                   "-DBUILD_CONFIG=$<CONFIG>"
                   "-DEXPECTED_DIAGNOSTIC=${regex}"
                   -P ${PROJECT_SOURCE_DIR}/cmake/CanaryExpectCompileFailure.cmake)
  set_tests_properties(canary.compile_fail.${name} PROPERTIES RESOURCE_LOCK canary_compile_fail_build)
endfunction()

canary_compile_fail_test(self_test "compile-fail self test")
canary_compile_fail_test(scalar_pointer "satisf(y|action of) [^ ]*Scalar")
canary_compile_fail_test(scalar_without_traits "satisf(y|action of) [^ ]*Scalar")
canary_compile_fail_test(math_function_missing "satisf(y|action of) [^ ]*HasExp")
canary_compile_fail_test(quantity_add_other_kind "(no match for [^ ]*operator\\+|invalid operands to binary expression).*Pressure.*Temperature")
canary_compile_fail_test(quantity_add_other_units "(no match for [^ ]*operator\\+|invalid operands to binary expression).*Density<double, canary::si>.*Density<double, canary::scaled>")
canary_compile_fail_test(quantity_product_without_relation "(no match for [^ ]*operator\\*|invalid operands to binary expression).*Density.*Temperature")
canary_compile_fail_test(convert_other_kind "conversion from .*energy_density.* to .*pressure|no viable conversion from .*energy_density.* to .*pressure")
canary_compile_fail_test(species_density_as_fractions "SpeciesRangeOf<.*kind::mass_fraction")
canary_compile_fail_test(thermo_density_density "ClosesInAnyOrder<.*IdealGas")
canary_compile_fail_test(thermo_gas_model_missing_closure "GasModel<.*NoPressureTemperatureGas")
canary_compile_fail_test(thermo_missing_extra_input "ClosesInAnyOrder<.*GuessRequiredGas")
canary_compile_fail_test(quantity_implicit_precision "conversion from .*double.* to .*float|no viable conversion from .*double.* to .*float")
canary_compile_fail_test(calorically_perfect_si_in_scaled "ClosesInAnyOrder<.*CaloricallyPerfect<double, canary::scaled>")
canary_compile_fail_test(thermally_perfect_si_in_scaled "ClosesInAnyOrder<.*ThermallyPerfectMixture<double, canary::scaled")
canary_compile_fail_test(two_temperature_si_modes_in_scaled "ClosesInAnyOrder<.*TwoTemperatureMixture<double, canary::scaled")
canary_compile_fail_test(multi_temperature_si_modes_in_scaled "ClosesInAnyOrder<.*MultiTemperatureMixture<double, canary::scaled")
canary_compile_fail_test(padding_free_gas_model_with_padding "satisf(y|action of) [^ ]*PaddingFree.*ThermallyPerfectMixture<float")
canary_compile_fail_test(calorically_perfect_swapped_kinds "no matching function for call to [^ ]*make_calorically_perfect")
canary_compile_fail_test(transport_other_units "no matching function for call to [^ ]*viscosity")
canary_compile_fail_test(quantity_molar_scaled "KindUnits<(K, U|canary::kind::molar_mass, canary::scaled)>")
canary_compile_fail_test(quantity_charge_density_other_kind "no matching function for call to [^ ]*charge_density")
canary_compile_fail_test(quantity_current_density_other_kind "no matching function for call to [^ ]*current_density")
canary_compile_fail_test(energy_exchange_single_temperature "RelaxingModel<.*ThermallyPerfectMixture")
canary_compile_fail_test(energy_exchange_mode_energy_missing "HasVibrationalElectronicEnergy<.*ThermallyPerfectMixture")
canary_compile_fail_test(reaction_mechanism_without_modes "HasEnergyModes<.*ThermallyPerfectMixture")
canary_compile_fail_test(multi_temperature_exchange_mode_energy_missing "HasModeEnergies<.*ThermallyPerfectMixture")
canary_compile_fail_test(multi_temperature_exchange_single_temperature "ModeRelaxingModel<.*ThermallyPerfectMixture")
canary_compile_fail_test(mixture_transport_composition "HasNormalisedComposition<.*CaloricallyPerfectMixture")
canary_compile_fail_test(mixture_transport_heat_capacity "HasSpeciesIsobaricSpecificHeat<.*CaloricallyPerfectMixture")
canary_compile_fail_test(mixture_transport_model "TransportModel<.*GuptaYosTransport")
canary_compile_fail_test(mixture_transport_vibrational_electronic_mode "HasVibrationalElectronicMode<.*MultiTemperatureMixture")

# Each public header is compiled alone, as strict C++20, with only the
# include directories of canary::canary (include/ and mdspan). This checks
# that every header includes what it uses and needs nothing newer than C++20.
file(GLOB_RECURSE canary_headers CONFIGURE_DEPENDS RELATIVE ${PROJECT_SOURCE_DIR}/include
     ${PROJECT_SOURCE_DIR}/include/canary/*.hpp)
set(header_check_sources)
foreach(header IN LISTS canary_headers)
  string(MAKE_C_IDENTIFIER ${header} source_name)
  set(source ${CMAKE_CURRENT_BINARY_DIR}/header_check/${source_name}.cpp)
  file(CONFIGURE OUTPUT ${source} CONTENT "#include <${header}>\n")
  list(APPEND header_check_sources ${source})
endforeach()
add_library(canary_header_check OBJECT ${header_check_sources})
target_link_libraries(canary_header_check PRIVATE canary::canary)
canary_test_options(canary_header_check)
set_target_properties(canary_header_check PROPERTIES CXX_STANDARD 20)

# Species database library: C++ tests (the public header is C++20) and the
# Python converter and schema checks, when the interpreter and its yaml and
# jsonschema modules are available.
if(TARGET canary_database)
  add_executable(canary_database_tests test_database.cpp)
  target_link_libraries(canary_database_tests PRIVATE canary::database Boost::ut)
  target_compile_definitions(canary_database_tests PRIVATE CANARY_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  canary_test_options(canary_database_tests)
  add_test(NAME canary.database COMMAND canary_database_tests)

  if(CANARY_REQUIRE_VALIDATION)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
  else()
    find_package(Python3 COMPONENTS Interpreter)
  endif()
  if(Python3_FOUND)
    add_test(NAME canary.database.output
             COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/species_database/test_json_output.py)
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import yaml, jsonschema" RESULT_VARIABLE python_modules_missing
                    OUTPUT_QUIET ERROR_QUIET)
    if(NOT python_modules_missing)
      set(species_tools ${PROJECT_SOURCE_DIR}/tools/species_database)
      set(generated_database_dir ${CMAKE_CURRENT_BINARY_DIR}/generated_database)
      set(generated_database_headers generated_species.hpp generated_transport.hpp
          generated_reactions.hpp generated_exchange.hpp negative_capacity.hpp small_mass.hpp)
      list(TRANSFORM generated_database_headers PREPEND "${generated_database_dir}/")
      file(GLOB database_generator_dependencies CONFIGURE_DEPENDS ${species_tools}/*.py
           ${PROJECT_SOURCE_DIR}/database/*.schema.json)
      add_custom_command(OUTPUT ${generated_database_headers}
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/generate_database_headers.py
                ${PROJECT_SOURCE_DIR} ${generated_database_dir}
        DEPENDS generate_database_headers.py ${database_generator_dependencies}
                ${PROJECT_SOURCE_DIR}/database/data/air5.json
                ${PROJECT_SOURCE_DIR}/database/data/air5_park.json
                ${PROJECT_SOURCE_DIR}/database/data/air5_park_exchange.json
                ${PROJECT_SOURCE_DIR}/database/data/air5_collision_integrals.json
        BYPRODUCTS ${generated_database_dir}/negative_capacity.json ${generated_database_dir}/small_mass.json
        VERBATIM)
      add_executable(canary_generated_database_tests test_generated_database.cpp ${generated_database_headers})
      target_link_libraries(canary_generated_database_tests PRIVATE canary::canary Boost::ut)
      target_include_directories(canary_generated_database_tests PRIVATE ${generated_database_dir})
      canary_test_options(canary_generated_database_tests)
      add_test(NAME canary.database.generated COMMAND canary_generated_database_tests)
      add_test(NAME canary.database.converters
               COMMAND ${Python3_EXECUTABLE} ${species_tools}/test_converters.py ${PROJECT_SOURCE_DIR})
      add_test(NAME canary.database.reaction_validation
               COMMAND ${Python3_EXECUTABLE} ${species_tools}/test_reaction_json.py)
      add_test(NAME canary.database.exchange_validation
               COMMAND ${Python3_EXECUTABLE} ${species_tools}/test_exchange_json.py)
      add_test(NAME canary.database.species_validation
               COMMAND ${Python3_EXECUTABLE} ${species_tools}/test_species_validation.py)
      add_test(NAME canary.database.schema_rejection
               COMMAND ${Python3_EXECUTABLE} ${species_tools}/test_schemas.py ${PROJECT_SOURCE_DIR})
      add_test(NAME canary.database.schema
               COMMAND ${Python3_EXECUTABLE} ${species_tools}/check_species_files.py
                       ${PROJECT_SOURCE_DIR}/database/species.schema.json ${PROJECT_SOURCE_DIR}/database/data/air5.json
                       ${PROJECT_SOURCE_DIR}/database/data/air11.json --
                       ${PROJECT_SOURCE_DIR}/database/reactions.schema.json
                       ${PROJECT_SOURCE_DIR}/database/data/air5_kim_jo.json
                       ${PROJECT_SOURCE_DIR}/database/data/air5_park.json
                       ${PROJECT_SOURCE_DIR}/database/data/air11_kim_jo.json --
                       ${PROJECT_SOURCE_DIR}/database/energy_exchange.schema.json
                       ${PROJECT_SOURCE_DIR}/database/data/air5_park_exchange.json
                       ${PROJECT_SOURCE_DIR}/database/data/air11_kim_jo_exchange.json --
                       ${PROJECT_SOURCE_DIR}/database/collision_integrals.schema.json
                       ${PROJECT_SOURCE_DIR}/database/data/air5_collision_integrals.json)
    elseif(CANARY_REQUIRE_VALIDATION)
      message(FATAL_ERROR "canary validation requires Python modules yaml and jsonschema")
    else()
      message(STATUS "canary: Python modules yaml or jsonschema missing; species converter tests skipped")
    endif()
  else()
    message(STATUS "canary: Python interpreter missing; species converter tests skipped")
  endif()
endif()
