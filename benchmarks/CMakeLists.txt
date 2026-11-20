# SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
# Copyright (C) 2026 Shahzeb Imran and the canary contributors
# Closure throughput. Run the Release build; results are printed per cell.
add_executable(canary_benchmark_closures closures.cpp)
target_link_libraries(canary_benchmark_closures PRIVATE canary::canary)
set_target_properties(canary_benchmark_closures PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

# Functions whose generated code is compared: the one-species closure through
# canary and the same formulas on plain doubles. Inspect with
#   objdump -d --no-show-raw-insn <object> | c++filt
add_library(canary_codegen OBJECT codegen.cpp)
target_link_libraries(canary_codegen PRIVATE canary::canary)
set_target_properties(canary_codegen PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
