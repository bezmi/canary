#ifndef GASLIB_DATABASE_HPP
#define GASLIB_DATABASE_HPP

// Species database files: JSON files of species data in SI units, read and
// written by the compiled library gaslib_database (link gaslib::database).
// This header is C++20 and does not depend on the JSON library; the
// conversions to model inputs below are header-only.
//
// File layout (strict JSON):
//   {
//     "format": "gaslib-species-1",
//     "units": {"molar_mass": "kg/mol", "energy": "J/mol", "temperature": "K"},
//     "species": {
//       "N2": {
//         "molar_mass": 0.0280134,
//         "charge": 0,
//         "thermo": {
//           "model": "nasa9",
//           "temperature_ranges": [200, 1000, 6000, 20000],
//           "blend_widths": [400, 1000],
//           "coefficients": [[9 numbers], [9 numbers], [9 numbers]],
//           "reference": "..."
//         },
//         "reference": "...", "notes": "..."
//       }
//     },
//     "reference": "...", "notes": "..."
//   }
// "reference" and "notes" are optional everywhere; every other key is
// required, and unknown keys are errors. Numbers are written in the
// shortest form that reads back to the same double.

#include <array>
#include <cstddef>
#include <filesystem>
#include <format>
#include <gaslib/build.hpp>
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/nasa9.hpp>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace gaslib::database
{

// Value of the "format" key of a species file.
inline constexpr std::string_view species_format = "gaslib-species-1";

// Units of the values in a file. Only these SI units are accepted.
struct Units {
    std::string molar_mass = "kg/mol";
    std::string energy = "J/mol";
    std::string temperature = "K";

    friend bool operator==(const Units&, const Units&) = default;
};

// Thermodynamic data of one species. For model "nasa9": temperature_ranges
// are the segment boundaries in K, coefficients holds the 9 coefficients of
// each segment (see nasa9.hpp), and blend_widths the width in K of the
// blending region around each interior boundary.
struct ThermoRecord {
    std::string model;
    std::vector<double> temperature_ranges;
    std::vector<double> blend_widths;
    std::vector<std::array<double, 9>> coefficients;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const ThermoRecord&, const ThermoRecord&) = default;
};

// One species: molar mass in kg/mol, charge in elementary charges.
struct SpeciesRecord {
    double molar_mass {};
    double charge {};
    ThermoRecord thermo;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const SpeciesRecord&, const SpeciesRecord&) = default;
};

// A species file. Species are keyed by name.
struct SpeciesFile {
    std::string format {species_format};
    Units units;
    std::map<std::string, SpeciesRecord> species;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const SpeciesFile&, const SpeciesFile&) = default;
};

// Thrown when a file cannot be read, is not valid JSON, has an unknown or a
// missing key, has another format or other units, or lacks a requested
// species. The message says which.
class DatabaseError : public std::runtime_error
{
  public:
    using std::runtime_error::runtime_error;
};

// The species file in json. Throws DatabaseError.
SpeciesFile parse_species(std::string_view json);

// The species file at path. Throws DatabaseError.
SpeciesFile read_species_file(const std::filesystem::path& path);

// file as indented JSON.
std::string write_species(const SpeciesFile& file);

// Writes file as indented JSON to path. Throws DatabaseError when the file
// cannot be written.
void write_species_file(const std::filesystem::path& path, const SpeciesFile& file);

// JSON Schema of species files.
std::string species_schema();

// The NASA-9 fit of a thermo record with at most S segments. Throws
// ParameterError when the model is not "nasa9", the record has more than S
// segments, or the numbers of ranges, blending widths and coefficient sets
// do not match. The values themselves are checked by make_nasa9_curve.
template <std::size_t S> Nasa9Data<S> nasa9_data(const ThermoRecord& thermo)
{
    if (thermo.model != "nasa9") {
        throw ParameterError(
          std::format(R"(thermo model must be "nasa9", got "{}")", thermo.model));
    }
    const std::size_t segments = thermo.coefficients.size();
    if (segments == 0 || segments > S) {
        throw ParameterError(
          std::format("number of segments must be between 1 and {}, got {}", S, segments));
    }
    if (thermo.temperature_ranges.size() != segments + 1 ||
        thermo.blend_widths.size() != segments - 1) {
        throw ParameterError(std::format(
          "{} segments need {} temperature ranges and {} blending widths, got {} and {}",
          segments,
          segments + 1,
          segments - 1,
          thermo.temperature_ranges.size(),
          thermo.blend_widths.size()));
    }
    Nasa9Data<S> data {};
    data.segments = static_cast<int>(segments);
    for (std::size_t i = 0; i <= segments; ++i) { data.breaks[i] = thermo.temperature_ranges[i]; }
    for (std::size_t i = 0; i + 1 < segments; ++i) { data.blends[i] = thermo.blend_widths[i]; }
    for (std::size_t j = 0; j < segments; ++j) { data.a[j] = thermo.coefficients[j]; }
    return data;
}

// Model inputs for the named species of file, in the given order, for
// make_thermally_perfect_mixture. Each name points to the species' key in
// file, so the result is valid while file exists. Throws DatabaseError for
// a name that file lacks, and ParameterError (with the species name) for a
// thermo record that nasa9_data rejects.
template <std::size_t S, std::size_t N>
std::array<Nasa9Species<S>, N> nasa9_species(const SpeciesFile& file,
                                             const std::array<std::string_view, N>& names)
{
    std::array<Nasa9Species<S>, N> out {};
    for (std::size_t i = 0; i < N; ++i) {
        const auto entry = file.species.find(std::string(names[i]));
        if (entry == file.species.end()) {
            throw DatabaseError(std::format("species \"{}\" is not in the file", names[i]));
        }
        try {
            out[i] = Nasa9Species<S> {
              .name = entry->first.c_str(),
              .molar_mass = entry->second.molar_mass,
              .charge = entry->second.charge,
              .fit = nasa9_data<S>(entry->second.thermo),
            };
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("species {}: {}", entry->first, error.what()));
        }
    }
    return out;
}

}  // namespace gaslib::database

#endif  // GASLIB_DATABASE_HPP
