// Reading and writing collision-integral files with glaze. Built as C++23;
// glaze does not appear in the public header.

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

void check_header(const CollisionIntegralsFile& file)
{
    if (file.format != collision_integrals_format) {
        throw DatabaseError(
          std::format(R"(format must be "{}", got "{}")", collision_integrals_format, file.format));
    }
    if (file.units != CollisionUnits {}) {
        throw DatabaseError(
          std::format(R"(units must be temperature "K" and area "m^2", got "{}" and "{}")",
                      file.units.temperature,
                      file.units.area));
    }
}

}  // namespace

CollisionIntegralsFile parse_collision_integrals(std::string_view json)
{
    CollisionIntegralsFile file {};
    // glaze reads from a null-terminated buffer.
    const std::string buffer(json);
    if (const auto error = glz::read<detail::read_options>(file, buffer)) {
        throw DatabaseError(
          std::format("invalid collision-integral file: {}", glz::format_error(error, buffer)));
    }
    check_header(file);
    return file;
}

CollisionIntegralsFile read_collision_integrals_file(const std::filesystem::path& path)
{
    const std::string text = detail::read_text(path, "collision-integral");
    try {
        return parse_collision_integrals(text);
    } catch (const DatabaseError& error) {
        throw DatabaseError(std::format("{}: {}", path.string(), error.what()));
    }
}

std::string write_collision_integrals(const CollisionIntegralsFile& file)
{
    std::string out;
    if (const auto error = glz::write<detail::write_options>(file, out)) {
        throw DatabaseError(
          std::format("cannot write collision-integral file: {}", glz::format_error(error)));
    }
    out += '\n';
    return out;
}

void write_collision_integrals_file(const std::filesystem::path& path,
                                    const CollisionIntegralsFile& file)
{ detail::write_text(path, write_collision_integrals(file)); }

std::string collision_integrals_schema()
{
    auto schema = glz::write_json_schema<CollisionIntegralsFile, detail::schema_options>();
    if (!schema) {
        throw DatabaseError(std::format("cannot write the collision-integral schema: {}",
                                        glz::format_error(schema.error())));
    }
    return *schema + '\n';
}

}  // namespace canary::database
