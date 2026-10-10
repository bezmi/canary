# Closure throughput. Run the Release build; results are printed per cell.
add_executable(gaslib_benchmark_closures closures.cpp)
target_link_libraries(gaslib_benchmark_closures PRIVATE gaslib::gaslib)
set_target_properties(gaslib_benchmark_closures PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

# Functions whose generated code is compared: the one-species closure through
# gaslib and the same formulas on plain doubles. Inspect with
#   objdump -d --no-show-raw-insn <object> | c++filt
add_library(gaslib_codegen OBJECT codegen.cpp)
target_link_libraries(gaslib_codegen PRIVATE gaslib::gaslib)
set_target_properties(gaslib_codegen PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
