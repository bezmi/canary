// Reading and writing energy-exchange files with glaze. Built as C++23;
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
{
    EnergyExchangeFile file {};
    // glaze reads from a null-terminated buffer.
    const std::string buffer(json);
    if (const auto error = glz::read<detail::read_options>(file, buffer)) {
        throw DatabaseError(
          std::format("invalid energy-exchange file: {}", glz::format_error(error, buffer)));
    }
    check_header(file);
    return file;
}

EnergyExchangeFile read_energy_exchange_file(const std::filesystem::path& path)
{
    const std::string text = detail::read_text(path, "energy-exchange");
    try {
        return parse_energy_exchange(text);
    } catch (const DatabaseError& error) {
        throw DatabaseError(std::format("{}: {}", path.string(), error.what()));
    }
}

std::string write_energy_exchange(const EnergyExchangeFile& file)
{
    std::string out;
    if (const auto error = glz::write<detail::write_options>(file, out)) {
        throw DatabaseError(
          std::format("cannot write energy-exchange file: {}", glz::format_error(error)));
    }
    out += '\n';
    return out;
}

void write_energy_exchange_file(const std::filesystem::path& path, const EnergyExchangeFile& file)
{ detail::write_text(path, write_energy_exchange(file)); }

std::string energy_exchange_schema()
{
    auto schema = glz::write_json_schema<EnergyExchangeFile, detail::schema_options>();
    if (!schema) {
        throw DatabaseError(std::format("cannot write the energy-exchange schema: {}",
                                        glz::format_error(schema.error())));
    }
    return *schema + '\n';
}

}  // namespace canary::database
