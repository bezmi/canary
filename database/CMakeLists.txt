# SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
# Copyright (C) 2026 Shahzeb Imran and the canary contributors
# Species database files (JSON). The public header is C++20; the sources are
# C++23 and use glaze, which consumers do not see.
add_library(canary_database STATIC species_file.cpp reactions_file.cpp exchange_file.cpp collision_file.cpp file_io.cpp)
add_library(canary::database ALIAS canary_database)
target_link_libraries(canary_database PUBLIC canary::canary PRIVATE $<BUILD_INTERFACE:glaze::glaze>)
set_target_properties(canary_database PROPERTIES EXPORT_NAME database CXX_STANDARD 23 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF)

# Write the JSON Schemas of species, reactions, energy-exchange and
# collision-integral files; database/species.schema.json,
# database/reactions.schema.json, database/energy_exchange.schema.json and
# database/collision_integrals.schema.json are their outputs.
add_executable(canary_species_schema species_schema_main.cpp)
target_link_libraries(canary_species_schema PRIVATE canary::database)
add_executable(canary_reactions_schema reactions_schema_main.cpp)
target_link_libraries(canary_reactions_schema PRIVATE canary::database)
add_executable(canary_energy_exchange_schema exchange_schema_main.cpp)
target_link_libraries(canary_energy_exchange_schema PRIVATE canary::database)
add_executable(canary_collision_integrals_schema collision_schema_main.cpp)
target_link_libraries(canary_collision_integrals_schema PRIVATE canary::database)
