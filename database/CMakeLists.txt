# Species database files (JSON). The public header is C++20; the sources are
# C++23 and use glaze, which consumers do not see.
add_library(gaslib_database STATIC species_file.cpp reactions_file.cpp exchange_file.cpp)
add_library(gaslib::database ALIAS gaslib_database)
target_link_libraries(gaslib_database PUBLIC gaslib::gaslib PRIVATE glaze::glaze)
set_target_properties(gaslib_database PROPERTIES CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

# Write the JSON Schemas of species, reactions and energy-exchange files;
# database/species.schema.json, database/reactions.schema.json and
# database/energy_exchange.schema.json are their outputs.
add_executable(gaslib_species_schema species_schema_main.cpp)
target_link_libraries(gaslib_species_schema PRIVATE gaslib::database)
add_executable(gaslib_reactions_schema reactions_schema_main.cpp)
target_link_libraries(gaslib_reactions_schema PRIVATE gaslib::database)
add_executable(gaslib_energy_exchange_schema exchange_schema_main.cpp)
target_link_libraries(gaslib_energy_exchange_schema PRIVATE gaslib::database)
