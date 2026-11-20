# SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
# Copyright (C) 2026 Shahzeb Imran and the canary contributors
cmake_minimum_required(VERSION 3.29)

# Calendar version YYYY.MM of the release.
project(
  canary
  VERSION 2026.11
  DESCRIPTION "Header-only gas thermodynamics, transport and rate functions for one cell"
  LANGUAGES CXX)

option(CANARY_BUILD_TESTS "Build the canary tests" ${PROJECT_IS_TOP_LEVEL})
option(CANARY_BUILD_BENCHMARKS "Build the canary benchmarks" ${PROJECT_IS_TOP_LEVEL})
option(CANARY_BUILD_DATABASE "Build the species database library (C++23, glaze)" ${PROJECT_IS_TOP_LEVEL})
option(CANARY_REQUIRE_VALIDATION "Require all test targets and validation dependencies" OFF)
option(CANARY_FIND_MDSPAN "Use an installed mdspan found with find_package instead of fetching it" OFF)
option(CANARY_BUILD_CUDA_CHECKS "Compile the numerical device checks with CUDA" OFF)
option(CANARY_BUILD_HIP_CHECKS "Compile the numerical device checks with HIP" OFF)
option(CANARY_BUILD_CUDA_RUNTIME_TESTS "Run the numerical device checks on a CUDA device" OFF)
option(CANARY_BUILD_HIP_RUNTIME_TESTS "Run the numerical device checks on a HIP device" OFF)

if(CANARY_REQUIRE_VALIDATION AND (NOT CANARY_BUILD_TESTS OR NOT CANARY_BUILD_DATABASE))
  message(FATAL_ERROR "CANARY_REQUIRE_VALIDATION requires CANARY_BUILD_TESTS and CANARY_BUILD_DATABASE")
endif()
if(CANARY_BUILD_CUDA_RUNTIME_TESTS AND NOT CANARY_BUILD_CUDA_CHECKS)
  message(FATAL_ERROR "CANARY_BUILD_CUDA_RUNTIME_TESTS requires CANARY_BUILD_CUDA_CHECKS")
endif()
if(CANARY_BUILD_HIP_RUNTIME_TESTS AND NOT CANARY_BUILD_HIP_CHECKS)
  message(FATAL_ERROR "CANARY_BUILD_HIP_RUNTIME_TESTS requires CANARY_BUILD_HIP_CHECKS")
endif()

# Each requested native language is required: a missing compiler fails
# configuration rather than disabling the request.
if(CANARY_BUILD_CUDA_CHECKS)
  enable_language(CUDA)
endif()
if(CANARY_BUILD_HIP_CHECKS)
  enable_language(HIP)
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

if(CANARY_BUILD_TESTS OR CANARY_BUILD_CUDA_RUNTIME_TESTS OR CANARY_BUILD_HIP_RUNTIME_TESTS)
  enable_testing()
endif()
if(CANARY_BUILD_TESTS)
  include(CanarySetupBoostUT)
  add_subdirectory(tests)
endif()

# Numerical device checks are independent of the host tests, database and
# benchmarks.
if(CANARY_BUILD_CUDA_CHECKS OR CANARY_BUILD_HIP_CHECKS)
  add_subdirectory(tests/device)
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
# A consumer that requests a version receives exactly that release: releases
# make no compatibility promise to each other.
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/canaryConfigVersion.cmake"
  COMPATIBILITY ExactVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/canaryConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/canaryConfigVersion.cmake"
  DESTINATION "${canary_cmake_install_dir}")
install(EXPORT canaryTargets NAMESPACE canary::
  DESTINATION "${canary_cmake_install_dir}")
