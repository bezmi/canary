# Tests are built as strict C++20 so that a header needing a later standard or
# a compiler extension fails here.
add_executable(gaslib_tests test_main.cpp)
target_link_libraries(gaslib_tests PRIVATE gaslib::gaslib Boost::ut)
set_target_properties(gaslib_tests PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)
target_compile_options(gaslib_tests PRIVATE -Wall -Wextra -Wpedantic -Werror)
add_test(NAME gaslib.unit COMMAND gaslib_tests)
