cmake_minimum_required(VERSION 3.29)

project(
  canary
  DESCRIPTION "Header-only gas thermodynamics, transport and rate functions for one cell"
  LANGUAGES CXX)

option(CANARY_BUILD_TESTS "Build the canary tests" ${PROJECT_IS_TOP_LEVEL})
option(CANARY_BUILD_BENCHMARKS "Build the canary benchmarks" ${PROJECT_IS_TOP_LEVEL})
option(CANARY_BUILD_DATABASE "Build the species database library (C++23, glaze)" ${PROJECT_IS_TOP_LEVEL})
option(CANARY_REQUIRE_VALIDATION "Require all test targets and validation dependencies" OFF)
option(CANARY_FIND_MDSPAN "Use an installed mdspan found with find_package instead of fetching it" OFF)

if(CANARY_REQUIRE_VALIDATION AND (NOT CANARY_BUILD_TESTS OR NOT CANARY_BUILD_DATABASE))
  message(FATAL_ERROR "CANARY_REQUIRE_VALIDATION requires CANARY_BUILD_TESTS and CANARY_BUILD_DATABASE")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# Keep managed dependency versions in separate source-cache entries. Explicit
# CPM_SOURCE_CACHE (including OFF) and environment settings still take precedence.
if(NOT DEFINED CPM_SOURCE_CACHE AND NOT DEFINED ENV{CPM_SOURCE_CACHE})
  set(CPM_SOURCE_CACHE "${CMAKE_BINARY_DIR}/_deps/cpm" CACHE PATH "CPM dependency source cache")
endif()
include(CPM)

if(PROJECT_IS_TOP_LEVEL)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  # canary has no C++ modules.
  set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CanarySetupMdspan)

add_library(canary INTERFACE)
add_library(canary::canary ALIAS canary)
target_include_directories(canary INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(canary INTERFACE cxx_std_20)
target_link_libraries(canary INTERFACE mdspan::mdspan)

if(CANARY_BUILD_DATABASE)
  include(CanarySetupGlaze)
  add_subdirectory(database)
endif()

if(CANARY_BUILD_TESTS)
  enable_testing()
  include(CanarySetupBoostUT)
  add_subdirectory(tests)
endif()

if(CANARY_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

install(DIRECTORY include/canary DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS canary EXPORT canaryTargets)
if(CANARY_BUILD_DATABASE)
  install(TARGETS canary_database EXPORT canaryTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

set(canary_cmake_install_dir "${CMAKE_INSTALL_LIBDIR}/cmake/canary")
configure_package_config_file(cmake/canaryConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/canaryConfig.cmake"
  INSTALL_DESTINATION "${canary_cmake_install_dir}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/canaryConfig.cmake"
  DESTINATION "${canary_cmake_install_dir}")
install(EXPORT canaryTargets NAMESPACE canary::
  DESTINATION "${canary_cmake_install_dir}")
