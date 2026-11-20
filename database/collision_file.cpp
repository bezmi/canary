// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// Reading and writing collision-integral files with glaze. Built as C++23;
// glaze does not appear in the public header.

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
{ return detail::parse_json<CollisionIntegralsFile>(json, "collision-integral", check_header); }

CollisionIntegralsFile read_collision_integrals_file(const std::filesystem::path& path)
{ return detail::read_json<CollisionIntegralsFile>(path, "collision-integral", check_header); }

std::string write_collision_integrals(const CollisionIntegralsFile& file)
{ return detail::write_json(file, "collision-integral"); }

void write_collision_integrals_file(const std::filesystem::path& path,
                                    const CollisionIntegralsFile& file)
{ detail::write_text(path, write_collision_integrals(file)); }

std::string collision_integrals_schema()
{ return detail::json_schema<CollisionIntegralsFile>("collision-integral"); }

}  // namespace canary::database
