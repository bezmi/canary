# Species database files (JSON). The public header is C++20; the sources are
# C++23 and use glaze, which consumers do not see.
add_library(gaslib_database STATIC species_file.cpp)
add_library(gaslib::database ALIAS gaslib_database)
target_link_libraries(gaslib_database PUBLIC gaslib::gaslib PRIVATE glaze::glaze)
set_target_properties(gaslib_database PROPERTIES CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

# Writes the JSON Schema of species files; database/species.schema.json is
# its output.
add_executable(gaslib_species_schema species_schema_main.cpp)
target_link_libraries(gaslib_species_schema PRIVATE gaslib::database)
