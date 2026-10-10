cmake_minimum_required(VERSION 3.29)

project(
  gaslib
  DESCRIPTION "Header-only gas thermodynamics, transport and rate functions for one cell"
  LANGUAGES CXX)

option(GASLIB_BUILD_TESTS "Build the gaslib tests" ${PROJECT_IS_TOP_LEVEL})
option(GASLIB_BUILD_BENCHMARKS "Build the gaslib benchmarks" ${PROJECT_IS_TOP_LEVEL})
option(GASLIB_BUILD_DATABASE "Build the species database library (C++23, glaze)" ${PROJECT_IS_TOP_LEVEL})
option(GASLIB_FIND_MDSPAN "Use an installed mdspan found with find_package instead of fetching it" OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# Keep managed dependency versions in separate source-cache entries. Explicit
# CPM_SOURCE_CACHE (including OFF) and environment settings still take precedence.
if(NOT DEFINED CPM_SOURCE_CACHE AND NOT DEFINED ENV{CPM_SOURCE_CACHE})
  set(CPM_SOURCE_CACHE "${CMAKE_BINARY_DIR}/_deps/cpm" CACHE PATH "CPM dependency source cache")
endif()
include(CPM)

if(PROJECT_IS_TOP_LEVEL)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  # gaslib has no C++ modules.
  set(CMAKE_CXX_SCAN_FOR_MODULES OFF)
endif()

include(GaslibSetupMdspan)

add_library(gaslib INTERFACE)
add_library(gaslib::gaslib ALIAS gaslib)
target_include_directories(gaslib INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>)
target_compile_features(gaslib INTERFACE cxx_std_20)
target_link_libraries(gaslib INTERFACE mdspan::mdspan)

if(GASLIB_BUILD_DATABASE)
  include(GaslibSetupGlaze)
  add_subdirectory(database)
endif()

if(GASLIB_BUILD_TESTS)
  enable_testing()
  include(GaslibSetupBoostUT)
  add_subdirectory(tests)
endif()

if(GASLIB_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()
