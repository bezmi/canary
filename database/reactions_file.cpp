// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// Reading and writing reactions files with glaze. Built as C++23; glaze does
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
{ return detail::parse_json<ReactionsFile>(json, "reactions", check_header); }

ReactionsFile read_reactions_file(const std::filesystem::path& path)
{ return detail::read_json<ReactionsFile>(path, "reactions", check_header); }

std::string write_reactions(const ReactionsFile& file)
{ return detail::write_json(file, "reactions"); }

void write_reactions_file(const std::filesystem::path& path, const ReactionsFile& file)
{ detail::write_text(path, write_reactions(file)); }

std::string reactions_schema() { return detail::json_schema<ReactionsFile>("reactions"); }

}  // namespace canary::database
