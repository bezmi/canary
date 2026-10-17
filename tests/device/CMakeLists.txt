# Numerical device checks: the public numerical operations compiled as native
# kernels (CUDA, HIP) and, with the same sources, as ordinary host functions.
# Compilation only; no target here runs a device. See README.md.

# Required instances, as header|group|P|T|units. Each becomes its own
# translation unit per language, so a parallel build compiles them concurrently.
set(canary_device_same_precision "float|float" "float|canary::test::Dual<float>"
                                 "double|double" "double|canary::test::Dual<double>")
set(canary_device_mixed_precision "double|float" "double|canary::test::Dual<float>"
                                  "float|double" "float|canary::test::Dual<double>")
set(canary_device_check_instances)
foreach(units IN ITEMS si scaled)
  foreach(scalars IN LISTS canary_device_same_precision)
    foreach(group IN ITEMS
            "thermodynamics_cases.hpp|thermodynamic_single_checks"
            "thermodynamics_cases.hpp|thermodynamic_modal_checks"
            "reactions_cases.hpp|reaction_checks"
            "exchange_cases.hpp|exchange_checks"
            "transport_cases.hpp|transport_checks")
      list(APPEND canary_device_check_instances "${group}|${scalars}|${units}")
    endforeach()
  endforeach()
  foreach(scalars IN LISTS canary_device_mixed_precision)
    foreach(group IN ITEMS
            "thermodynamics_cases.hpp|thermodynamic_mixed_checks"
            "reactions_cases.hpp|reaction_mixed_checks"
            "transport_cases.hpp|transport_structural_checks")
      list(APPEND canary_device_check_instances "${group}|${scalars}|${units}")
    endforeach()
  endforeach()
endforeach()

# Generates one source per required instance with the given extension and
# returns their paths in out_var.
function(canary_device_check_sources out_var language extension annotation entry)
  set(sources)
  set(index 0)
  foreach(instance IN LISTS canary_device_check_instances)
    string(REPLACE "|" ";" parts "${instance}")
    list(GET parts 0 CANARY_CHECK_HEADER)
    list(GET parts 1 CANARY_CHECK_GROUP)
    list(GET parts 2 CANARY_CHECK_P)
    list(GET parts 3 CANARY_CHECK_T)
    list(GET parts 4 CANARY_CHECK_UNITS)
    set(CANARY_CHECK_LANGUAGE "${language}")
    set(CANARY_CHECK_ANNOTATION "${annotation}")
    set(CANARY_CHECK_ENTRY "${entry}")
    set(source "${CMAKE_CURRENT_BINARY_DIR}/${language}/${CANARY_CHECK_GROUP}_${index}.${extension}")
    configure_file(check.in "${source}" @ONLY)
    list(APPEND sources "${source}")
    math(EXPR index "${index} + 1")
  endforeach()
  set(${out_var} "${sources}" PARENT_SCOPE)
endfunction()

# Host companion: the same entries as ordinary functions, with the default empty
# annotation, built whenever a native check is selected.
canary_device_check_sources(host_sources host cpp "" "")
add_library(canary_device_host_check OBJECT ${host_sources})
target_include_directories(canary_device_host_check PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(canary_device_host_check PRIVATE canary::canary)
set_target_properties(canary_device_host_check PROPERTIES
  CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
target_compile_options(canary_device_host_check PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic -Werror>)

set(native_annotation "#define CANARY_HOST_DEVICE __host__ __device__")

if(CANARY_BUILD_CUDA_CHECKS)
  canary_device_check_sources(cuda_sources CUDA cu "${native_annotation}" __global__)
  add_library(canary_cuda_check OBJECT ${cuda_sources})
  target_include_directories(canary_cuda_check PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(canary_cuda_check PRIVATE canary::canary)
  set_target_properties(canary_cuda_check PROPERTIES
    CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON CUDA_EXTENSIONS OFF)
  # Canary's constexpr standard-library and mdspan calls need relaxed constexpr
  # under NVIDIA's compiler (see include/canary/config.hpp).
  target_compile_options(canary_cuda_check PRIVATE
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--expt-relaxed-constexpr>)
endif()

if(CANARY_BUILD_HIP_CHECKS)
  canary_device_check_sources(hip_sources HIP hip "${native_annotation}" __global__)
  set_source_files_properties(${hip_sources} PROPERTIES LANGUAGE HIP)
  add_library(canary_hip_check OBJECT ${hip_sources})
  target_include_directories(canary_hip_check PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(canary_hip_check PRIVATE canary::canary)
  set_target_properties(canary_hip_check PROPERTIES
    HIP_STANDARD 20 HIP_STANDARD_REQUIRED ON HIP_EXTENSIONS OFF)
endif()
