// SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
// Copyright (C) 2026 Shahzeb Imran and the canary contributors
// Reading and writing energy-exchange files with glaze. Built as C++23;
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

void check_header(const EnergyExchangeFile& file)
{
    if (file.format != energy_exchange_format) {
        throw DatabaseError(
          std::format(R"(format must be "{}", got "{}")", energy_exchange_format, file.format));
    }
    if (file.units != ExchangeUnits {}) {
        throw DatabaseError(std::format(
          "units must be energy \"J/mol\", temperature \"K\" and area \"m^2\", got \"{}\", \"{}\" "
          "and \"{}\"",
          file.units.energy,
          file.units.temperature,
          file.units.area));
    }
}

}  // namespace

EnergyExchangeFile parse_energy_exchange(std::string_view json)
{ return detail::parse_json<EnergyExchangeFile>(json, "energy-exchange", check_header); }

EnergyExchangeFile read_energy_exchange_file(const std::filesystem::path& path)
{ return detail::read_json<EnergyExchangeFile>(path, "energy-exchange", check_header); }

std::string write_energy_exchange(const EnergyExchangeFile& file)
{ return detail::write_json(file, "energy-exchange"); }

void write_energy_exchange_file(const std::filesystem::path& path, const EnergyExchangeFile& file)
{ detail::write_text(path, write_energy_exchange(file)); }

std::string energy_exchange_schema()
{ return detail::json_schema<EnergyExchangeFile>("energy-exchange"); }

}  // namespace canary::database
