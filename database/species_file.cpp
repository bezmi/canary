// Reading and writing species files with glaze. Built as C++23; glaze does
// not appear in the public header.

#include "json_io.hpp"

#include <filesystem>
#include <format>
#include <canary/database.hpp>
#include <glaze/glaze.hpp>
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
{
    SpeciesFile file {};
    // glaze reads from a null-terminated buffer.
    const std::string buffer(json);
    if (const auto error = glz::read<detail::read_options>(file, buffer)) {
        throw DatabaseError(
          std::format("invalid species file: {}", glz::format_error(error, buffer)));
    }
    check_header(file);
    return file;
}

SpeciesFile read_species_file(const std::filesystem::path& path)
{
    const std::string text = detail::read_text(path, "species");
    try {
        return parse_species(text);
    } catch (const DatabaseError& error) {
        throw DatabaseError(std::format("{}: {}", path.string(), error.what()));
    }
}

std::string write_species(const SpeciesFile& file)
{
    std::string out;
    if (const auto error = glz::write<detail::write_options>(file, out)) {
        throw DatabaseError(std::format("cannot write species file: {}", glz::format_error(error)));
    }
    out += '\n';
    return out;
}

void write_species_file(const std::filesystem::path& path, const SpeciesFile& file)
{ detail::write_text(path, write_species(file)); }

std::string species_schema()
{
    auto schema = glz::write_json_schema<SpeciesFile, detail::schema_options>();
    if (!schema) {
        throw DatabaseError(
          std::format("cannot write the species schema: {}", glz::format_error(schema.error())));
    }
    return *schema + '\n';
}

}  // namespace canary::database
