// Reading and writing species files with glaze. Built as C++23; glaze does
// not appear in the public header.

#include "json_io.hpp"

#include <canary/database.hpp>
#include <filesystem>
#include <format>
#include <string>
#include <string_view>

template <> struct glz::meta<canary::SpeciesRole> {
    using enum canary::SpeciesRole;
    static constexpr auto value = glz::enumerate(heavy, electron);
};

namespace canary::database::detail
{

// Read the format before required record fields so older files receive the
// migration diagnostic even when their records have no role field.
// The reflected type requires external linkage with Clang.
struct SpeciesFormat {
    std::string format;
};

}  // namespace canary::database::detail

namespace canary::database
{

namespace
{

void check_format(std::string_view format)
{
    if (format != species_format) {
        throw DatabaseError(std::format(
          R"(format must be "{}", got "{}"; regenerate species files with an explicit heavy/electron role for every species)",
          species_format,
          format));
    }
}

void check_header(const SpeciesFile& file)
{
    check_format(file.format);
    if (file.units != Units {}) {
        throw DatabaseError(std::format("units must be molar_mass \"kg/mol\", energy \"J/mol\", "
                                        "temperature \"K\" and length \"m\", got \"{}\", \"{}\", "
                                        "\"{}\" and \"{}\"",
                                        file.units.molar_mass,
                                        file.units.energy,
                                        file.units.temperature,
                                        file.units.length));
    }
    for (const auto& [name, record] : file.species) {
        if (record.role != SpeciesRole::heavy && record.role != SpeciesRole::electron) {
            throw DatabaseError(std::format("species {}: role must be heavy or electron", name));
        }
        if (record.role == SpeciesRole::electron && record.charge != -1.0) {
            throw DatabaseError(std::format("species {}: electron role requires charge -1", name));
        }
    }
}

}  // namespace

SpeciesFile parse_species(std::string_view json)
{
    detail::SpeciesFormat header {};
    const std::string buffer(json);
    constexpr glz::opts options {.error_on_unknown_keys = false, .error_on_missing_keys = true};
    if (const auto error = glz::read<options>(header, buffer)) {
        throw DatabaseError(
          std::format("invalid species file: {}", glz::format_error(error, buffer)));
    }
    check_format(header.format);
    return detail::parse_json<SpeciesFile>(json, "species", check_header);
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

std::string write_species(const SpeciesFile& file) { return detail::write_json(file, "species"); }

void write_species_file(const std::filesystem::path& path, const SpeciesFile& file)
{ detail::write_text(path, write_species(file)); }

std::string species_schema() { return detail::json_schema<SpeciesFile>("species"); }

}  // namespace canary::database
