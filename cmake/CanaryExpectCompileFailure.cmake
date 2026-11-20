# SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
# Copyright (C) 2026 Shahzeb Imran and the canary contributors
# Require both a failed build and the diagnostic for the rejected expression.
foreach(required BUILD_DIR BUILD_TARGET EXPECTED_DIAGNOSTIC)
  if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
    message(FATAL_ERROR "Missing ${required} for compile-failure check")
  endif()
endforeach()

set(build_command "${CMAKE_COMMAND}" --build "${BUILD_DIR}" --target "${BUILD_TARGET}" --parallel 1)
if(DEFINED BUILD_CONFIG AND NOT BUILD_CONFIG STREQUAL "")
  list(APPEND build_command --config "${BUILD_CONFIG}")
endif()
execute_process(COMMAND ${build_command}
                RESULT_VARIABLE build_result
                OUTPUT_VARIABLE build_stdout
                ERROR_VARIABLE build_stderr)
set(build_output "${build_stdout}\n${build_stderr}")
if(NOT build_result MATCHES "^[0-9]+$")
  message(FATAL_ERROR "Could not run compile-failure build: ${build_result}\n${build_output}")
elseif(build_result EQUAL 0)
  message(FATAL_ERROR "Expected compilation to fail, but the build succeeded.\n${build_output}")
elseif(NOT build_output MATCHES "${EXPECTED_DIAGNOSTIC}")
  message(FATAL_ERROR "Compilation failed without the expected diagnostic: ${EXPECTED_DIAGNOSTIC}\n${build_output}")
endif()
