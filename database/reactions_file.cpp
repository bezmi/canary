// Reading and writing reactions files with glaze. Built as C++23; glaze does
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

void check_header(const ReactionsFile& file)
{
    if (file.format != reactions_format) {
        throw DatabaseError(
          std::format(R"(format must be "{}", got "{}")", reactions_format, file.format));
    }
    if (file.units != ReactionUnits {}) {
        throw DatabaseError(
          std::format("units must be quantity \"mol\", length \"m\", time \"s\" and temperature "
                      "\"K\", got \"{}\", \"{}\", \"{}\" and \"{}\"",
                      file.units.quantity,
                      file.units.length,
                      file.units.time,
                      file.units.temperature));
    }
}

}  // namespace

ReactionsFile parse_reactions(std::string_view json)
{
    ReactionsFile file {};
    // glaze reads from a null-terminated buffer.
    const std::string buffer(json);
    if (const auto error = glz::read<detail::read_options>(file, buffer)) {
        throw DatabaseError(
          std::format("invalid reactions file: {}", glz::format_error(error, buffer)));
    }
    check_header(file);
    return file;
}

ReactionsFile read_reactions_file(const std::filesystem::path& path)
{
    const std::string text = detail::read_text(path, "reactions");
    try {
        return parse_reactions(text);
    } catch (const DatabaseError& error) {
        throw DatabaseError(std::format("{}: {}", path.string(), error.what()));
    }
}

std::string write_reactions(const ReactionsFile& file)
{
    std::string out;
    if (const auto error = glz::write<detail::write_options>(file, out)) {
        throw DatabaseError(
          std::format("cannot write reactions file: {}", glz::format_error(error)));
    }
    out += '\n';
    return out;
}

void write_reactions_file(const std::filesystem::path& path, const ReactionsFile& file)
{ detail::write_text(path, write_reactions(file)); }

std::string reactions_schema()
{
    auto schema = glz::write_json_schema<ReactionsFile, detail::schema_options>();
    if (!schema) {
        throw DatabaseError(
          std::format("cannot write the reactions schema: {}", glz::format_error(schema.error())));
    }
    return *schema + '\n';
}

}  // namespace canary::database
