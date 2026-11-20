// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
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

template <> struct glz::meta<canary::database::ElementName> {
    using mimic = std::string;
    static constexpr auto value = &canary::database::ElementName::value;
};

template <> struct glz::meta<canary::database::ElementComposition> {
    static constexpr auto value = &canary::database::ElementComposition::entries;
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
        if (record.elements) {
            const auto elements = detail::element_views(*record.elements);
            try {
                canary::detail::check_species_metadata(
                  {.role = record.role, .elements = std::span<const ElementCount> {elements}},
                  record.charge);
            } catch (const ParameterError& error) {
                throw DatabaseError(std::format("species {}: {}", name, error.what()));
            }
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

std::string write_species(const SpeciesFile& file)
{
    check_header(file);
    return detail::write_json(file, "species");
}

void write_species_file(const std::filesystem::path& path, const SpeciesFile& file)
{ detail::write_text(path, write_species(file)); }

std::string species_schema()
{
    auto schema = glz::read_json<glz::generic>(detail::json_schema<SpeciesFile>("species"));
    if (!schema) {
        throw DatabaseError(std::format("cannot read generated species schema: {}",
                                        glz::format_error(schema.error())));
    }
    auto& elements = (*schema)["properties"]["species"]["additionalProperties"]["properties"]
                              ["elements"]["anyOf"][0];
    elements["propertyNames"] = glz::generic {{"minLength", 1}};
    elements["additionalProperties"]["minimum"] = 1;
    return detail::write_json(*schema, "species schema");
}

}  // namespace canary::database
