// Reading and writing species files with glaze. Built as C++23; glaze does
// not appear in the public header.

#include "json_io.hpp"

#include <canary/database.hpp>
#include <filesystem>
#include <format>
#include <string>
#include <string_view>

namespace canary::database
{

namespace
{

void check_header(const SpeciesFile& file)
{
    if (file.format != species_format) {
        throw DatabaseError(
          std::format(R"(format must be "{}", got "{}")", species_format, file.format));
    }
    if (file.units != Units {}) {
        throw DatabaseError(std::format("units must be molar_mass \"kg/mol\", energy \"J/mol\", "
                                        "temperature \"K\" and length \"m\", got \"{}\", \"{}\", "
                                        "\"{}\" and \"{}\"",
                                        file.units.molar_mass,
                                        file.units.energy,
                                        file.units.temperature,
                                        file.units.length));
    }
}

}  // namespace

SpeciesFile parse_species(std::string_view json)
{ return detail::parse_json<SpeciesFile>(json, "species", check_header); }

SpeciesFile read_species_file(const std::filesystem::path& path)
{ return detail::read_json<SpeciesFile>(path, "species", check_header); }

std::string write_species(const SpeciesFile& file) { return detail::write_json(file, "species"); }

void write_species_file(const std::filesystem::path& path, const SpeciesFile& file)
{ detail::write_text(path, write_species(file)); }

std::string species_schema() { return detail::json_schema<SpeciesFile>("species"); }

}  // namespace canary::database
