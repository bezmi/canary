#ifndef CANARY_DATABASE_JSON_IO_HPP
#define CANARY_DATABASE_JSON_IO_HPP

// Typed JSON parsing, serialization, schemas and file I/O shared by the
// species, reactions, energy-exchange and collision-integral files.
// Private to the database library.

#include <canary/database.hpp>
#include <filesystem>
#include <format>
#include <fstream>
#include <glaze/glaze.hpp>
#include <ios>
#include <iterator>
#include <string>
#include <string_view>

namespace canary::database::detail
{

// Unknown keys and missing required keys are errors; optional members may be
// absent.
inline constexpr glz::opts read_options {
  .error_on_unknown_keys = true,
  .error_on_missing_keys = true,
};
inline constexpr glz::opts write_options {.prettify = true};
// The schema marks a key required exactly when reading requires it.
inline constexpr glz::opts schema_options {.prettify = true, .error_on_missing_keys = true};

// Contents of the file at path. Throws DatabaseError naming the kind of file
// when it cannot be opened.
inline std::string read_text(const std::filesystem::path& path, std::string_view kind)
{
    std::ifstream in(path, std::ios::binary);
    if (!in) { throw DatabaseError(std::format("cannot open {} file {}", kind, path.string())); }
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

// Writes text to path. Throws DatabaseError when it cannot.
inline void write_text(const std::filesystem::path& path, const std::string& text)
{
    std::ofstream out(path, std::ios::binary);
    if (!out) { throw DatabaseError(std::format("cannot open {} for writing", path.string())); }
    out << text;
    if (!out) { throw DatabaseError(std::format("cannot write {}", path.string())); }
}

template <typename File, typename Validate>
File parse_json(std::string_view json, std::string_view kind, Validate validate)
{
    File file {};
    // glaze reads from a null-terminated buffer.
    const std::string buffer(json);
    if (const auto error = glz::read<read_options>(file, buffer)) {
        throw DatabaseError(
          std::format("invalid {} file: {}", kind, glz::format_error(error, buffer)));
    }
    validate(file);
    return file;
}

template <typename File, typename Validate>
File read_json(const std::filesystem::path& path, std::string_view kind, Validate validate)
{
    const std::string text = read_text(path, kind);
    try {
        return parse_json<File>(text, kind, validate);
    } catch (const DatabaseError& error) {
        throw DatabaseError(std::format("{}: {}", path.string(), error.what()));
    }
}

template <typename File> std::string write_json(const File& file, std::string_view kind)
{
    std::string out;
    if (const auto error = glz::write<write_options>(file, out)) {
        throw DatabaseError(
          std::format("cannot write {} file: {}", kind, glz::format_error(error)));
    }
    out += '\n';
    return out;
}

template <typename File> std::string json_schema(std::string_view kind)
{
    auto schema = glz::write_json_schema<File, schema_options>();
    if (!schema) {
        throw DatabaseError(
          std::format("cannot write the {} schema: {}", kind, glz::format_error(schema.error())));
    }
    return *schema + '\n';
}

}  // namespace canary::database::detail

#endif  // CANARY_DATABASE_JSON_IO_HPP
