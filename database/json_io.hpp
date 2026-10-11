#ifndef GASLIB_DATABASE_JSON_IO_HPP
#define GASLIB_DATABASE_JSON_IO_HPP

// glaze options and file input and output shared by the species,
// reactions, energy-exchange and collision-integral files. Private to the database library.

#include <filesystem>
#include <format>
#include <fstream>
#include <gaslib/database.hpp>
#include <glaze/glaze.hpp>
#include <ios>
#include <iterator>
#include <string>
#include <string_view>

namespace gaslib::database::detail
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

}  // namespace gaslib::database::detail

#endif  // GASLIB_DATABASE_JSON_IO_HPP
