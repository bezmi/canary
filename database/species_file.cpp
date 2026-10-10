// Reading and writing species files with glaze. Built as C++23; glaze does
// not appear in the public header.

#include <filesystem>
#include <format>
#include <fstream>
#include <gaslib/database.hpp>
#include <glaze/glaze.hpp>
#include <ios>
#include <iterator>
#include <string>
#include <string_view>

namespace gaslib::database
{

namespace
{

// Unknown keys and missing required keys are errors; optional members may be
// absent.
constexpr glz::opts read_options {.error_on_unknown_keys = true, .error_on_missing_keys = true};
constexpr glz::opts write_options {.prettify = true};
// The schema marks a key required exactly when reading requires it.
constexpr glz::opts schema_options {.prettify = true, .error_on_missing_keys = true};

void check_header(const SpeciesFile& file)
{
    if (file.format != species_format) {
        throw DatabaseError(
          std::format(R"(format must be "{}", got "{}")", species_format, file.format));
    }
    if (file.units != Units {}) {
        throw DatabaseError(std::format("units must be molar_mass \"kg/mol\", energy \"J/mol\" and "
                                        "temperature \"K\", got \"{}\", \"{}\" and \"{}\"",
                                        file.units.molar_mass,
                                        file.units.energy,
                                        file.units.temperature));
    }
}

}  // namespace

SpeciesFile parse_species(std::string_view json)
{
    SpeciesFile file {};
    // glaze reads from a null-terminated buffer.
    const std::string buffer(json);
    if (const auto error = glz::read<read_options>(file, buffer)) {
        throw DatabaseError(
          std::format("invalid species file: {}", glz::format_error(error, buffer)));
    }
    check_header(file);
    return file;
}

SpeciesFile read_species_file(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) { throw DatabaseError(std::format("cannot open species file {}", path.string())); }
    const std::string text {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
    try {
        return parse_species(text);
    } catch (const DatabaseError& error) {
        throw DatabaseError(std::format("{}: {}", path.string(), error.what()));
    }
}

std::string write_species(const SpeciesFile& file)
{
    std::string out;
    if (const auto error = glz::write<write_options>(file, out)) {
        throw DatabaseError(std::format("cannot write species file: {}", glz::format_error(error)));
    }
    out += '\n';
    return out;
}

void write_species_file(const std::filesystem::path& path, const SpeciesFile& file)
{
    const std::string text = write_species(file);
    std::ofstream out(path, std::ios::binary);
    if (!out) { throw DatabaseError(std::format("cannot open {} for writing", path.string())); }
    out << text;
    if (!out) { throw DatabaseError(std::format("cannot write {}", path.string())); }
}

std::string species_schema()
{
    auto schema = glz::write_json_schema<SpeciesFile, schema_options>();
    if (!schema) {
        throw DatabaseError(
          std::format("cannot write the species schema: {}", glz::format_error(schema.error())));
    }
    return *schema + '\n';
}

}  // namespace gaslib::database
