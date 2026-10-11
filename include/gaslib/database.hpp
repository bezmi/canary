#ifndef GASLIB_DATABASE_HPP
#define GASLIB_DATABASE_HPP

// Species database files: JSON files of species data, of reaction
// mechanisms and of energy exchange in SI units, read and written by the compiled library
// gaslib_database (link gaslib::database). This header is C++20 and does
// not depend on the JSON library; the conversions to model inputs below are
// header-only.
//
// Species file layout (strict JSON):
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
//
// Reactions files hold a mechanism of reversible reactions (see
// reactions.hpp), in SI units per mole:
//   {
//     "format": "gaslib-reactions-1",
//     "units": {"quantity": "mol", "length": "m", "time": "s", "temperature": "K"},
//     "temperature_limits": {"lower": 300.0, "upper": 30000.0},
//     "reactions": [
//       {
//         "equation": "NO + M <=> N + O + M",
//         "reactants": {"NO": 1.0},
//         "products": {"N": 1.0, "O": 1.0},
//         "forward": {"model": "arrhenius", "A": 1450000000.0, "n": 0.0,
//                     "activation_temperature": 75200.0,
//                     "translational_exponent": 0.5},
//         "backward": {"model": "equilibrium"},
//         "third_body": {"efficiencies": {"NO": 0.664827586}},
//         "label": "...", "reference": "...", "notes": "..."
//       }
//     ],
//     "reference": "...", "notes": "..."
//   }
// "third_body", "label", "reference" and "notes" are optional, and so are
// the keys A, n, activation_temperature and translational_exponent of a
// rate; every other key is required, and unknown keys are errors. A rate is
// either "arrhenius" with A, n and activation_temperature and optionally
// translational_exponent (default 1), or (backward only) "equilibrium" with
// none of them; mechanism_input checks this.
//
// Energy-exchange files hold the vibration-translation exchange and the
// chemistry-vibration coupling of a two-temperature model (see
// energy_exchange.hpp), in SI units per mole:
//   {
//     "format": "gaslib-energy-exchange-1",
//     "units": {"energy": "J/mol", "temperature": "K", "area": "m^2"},
//     "vibration_translation": [
//       {
//         "species": "N2", "partner": "O", "model": "landau_teller",
//         "relaxation_time": {
//           "model": "millikan_white", "a": 72.4, "b": 0.015,
//           "high_temperature_correction": {"model": "park_1993",
//             "cross_section": 1e-21, "reference_temperature": 50000.0}
//         },
//         "reference": "...", "notes": "..."
//       }
//     ],
//     "chemistry_vibration": [
//       {
//         "reaction": 0, "species": "N2", "model": "impartial_dissociation",
//         "dissociation_energy": 945360.0, "vibrational_temperature": 3393.44,
//         "reference": "...", "notes": "..."
//       }
//     ],
//     "reference": "...", "notes": "..."
//   }
// "high_temperature_correction", "reference" and "notes" are optional;
// every other key is required, and unknown keys are errors. "reaction" is
// the 0-based index of the reaction in the mechanism used with the file.

#include <array>
#include <cstddef>
#include <filesystem>
#include <format>
#include <gaslib/build.hpp>
#include <gaslib/energy_exchange.hpp>
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/nasa9.hpp>
#include <gaslib/reactions.hpp>
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

// Value of the "format" key of a reactions file.
inline constexpr std::string_view reactions_format = "gaslib-reactions-1";

// Units of the values in a reactions file: amounts of substance in mol,
// lengths in m, times in s and temperatures in K. Only these are accepted.
struct ReactionUnits {
    std::string quantity = "mol";
    std::string length = "m";
    std::string time = "s";
    std::string temperature = "K";

    friend bool operator==(const ReactionUnits&, const ReactionUnits&) = default;
};

// Temperature limits of rate evaluation in K.
struct TemperatureLimits {
    double lower {};
    double upper {};

    friend bool operator==(const TemperatureLimits&, const TemperatureLimits&) = default;
};

// A rate constant: model "arrhenius" with A, n and activation_temperature
// and optionally translational_exponent (see ArrheniusInput), or model
// "equilibrium" with none of them.
struct RateRecord {
    std::string model;
    std::optional<double> A;
    std::optional<double> n;
    std::optional<double> activation_temperature;
    std::optional<double> translational_exponent;

    friend bool operator==(const RateRecord&, const RateRecord&) = default;
};

// Third-body efficiencies by species name; species not listed have 1.
struct ThirdBodyRecord {
    std::map<std::string, double> efficiencies;

    friend bool operator==(const ThirdBodyRecord&, const ThirdBodyRecord&) = default;
};

// One reversible reaction. reactants and products map species names to
// stoichiometric coefficients; equation is informational.
struct ReactionRecord {
    std::string equation;
    std::map<std::string, double> reactants;
    std::map<std::string, double> products;
    RateRecord forward;
    RateRecord backward;
    std::optional<ThirdBodyRecord> third_body;
    std::optional<std::string> label;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const ReactionRecord&, const ReactionRecord&) = default;
};

// A reactions file. Reactions keep their order.
struct ReactionsFile {
    std::string format {reactions_format};
    ReactionUnits units;
    TemperatureLimits temperature_limits;
    std::vector<ReactionRecord> reactions;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const ReactionsFile&, const ReactionsFile&) = default;
};

// Value of the "format" key of an energy-exchange file.
inline constexpr std::string_view energy_exchange_format = "gaslib-energy-exchange-1";

// Units of the values in an energy-exchange file: energies per mole in
// J/mol, temperatures in K and cross-sections in m^2. Only these are
// accepted.
struct ExchangeUnits {
    std::string energy = "J/mol";
    std::string temperature = "K";
    std::string area = "m^2";

    friend bool operator==(const ExchangeUnits&, const ExchangeUnits&) = default;
};

// High-temperature correction of a relaxation time: model "park_1993" with
// the collision cross-section in m^2 at the reference temperature in K (see
// ParkCorrectionInput).
struct HighTemperatureCorrectionRecord {
    std::string model;
    double cross_section {};
    double reference_temperature {};

    friend bool operator==(const HighTemperatureCorrectionRecord&,
                           const HighTemperatureCorrectionRecord&) = default;
};

// Relaxation time of a vibration-translation pair: model "millikan_white"
// with a in K^(1/3) and b in K^(-1/3) (see MillikanWhiteInput), optionally
// with a high-temperature correction.
struct RelaxationTimeRecord {
    std::string model;
    double a {};
    double b {};
    std::optional<HighTemperatureCorrectionRecord> high_temperature_correction;

    friend bool operator==(const RelaxationTimeRecord&, const RelaxationTimeRecord&) = default;
};

// Vibration-translation exchange of species with collision partner: model
// "landau_teller".
struct VibrationTranslationRecord {
    std::string species;
    std::string partner;
    std::string model;
    RelaxationTimeRecord relaxation_time;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const VibrationTranslationRecord&,
                           const VibrationTranslationRecord&) = default;
};

// Chemistry-vibration coupling of species in the reaction with 0-based
// index reaction: model "impartial_dissociation" with the dissociation
// energy in J/mol and the vibrational temperature in K (see
// ChemistryVibrationInput).
struct ChemistryVibrationRecord {
    int reaction {};
    std::string species;
    std::string model;
    double dissociation_energy {};
    double vibrational_temperature {};
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const ChemistryVibrationRecord&,
                           const ChemistryVibrationRecord&) = default;
};

// An energy-exchange file. Entries keep their order.
struct EnergyExchangeFile {
    std::string format {energy_exchange_format};
    ExchangeUnits units;
    std::vector<VibrationTranslationRecord> vibration_translation;
    std::vector<ChemistryVibrationRecord> chemistry_vibration;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const EnergyExchangeFile&, const EnergyExchangeFile&) = default;
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

// The reactions file in json. Throws DatabaseError.
ReactionsFile parse_reactions(std::string_view json);

// The reactions file at path. Throws DatabaseError.
ReactionsFile read_reactions_file(const std::filesystem::path& path);

// file as indented JSON.
std::string write_reactions(const ReactionsFile& file);

// Writes file as indented JSON to path. Throws DatabaseError when the file
// cannot be written.
void write_reactions_file(const std::filesystem::path& path, const ReactionsFile& file);

// JSON Schema of reactions files.
std::string reactions_schema();

// The energy-exchange file in json. Throws DatabaseError.
EnergyExchangeFile parse_energy_exchange(std::string_view json);

// The energy-exchange file at path. Throws DatabaseError.
EnergyExchangeFile read_energy_exchange_file(const std::filesystem::path& path);

// file as indented JSON.
std::string write_energy_exchange(const EnergyExchangeFile& file);

// Writes file as indented JSON to path. Throws DatabaseError when the file
// cannot be written.
void write_energy_exchange_file(const std::filesystem::path& path, const EnergyExchangeFile& file);

// JSON Schema of energy-exchange files.
std::string energy_exchange_schema();

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

namespace detail
{

// Up to K species coefficients from a map; each name points to the map's
// key. Throws ParameterError when the map has more than K entries.
template <std::size_t K>
SpeciesCoefficients<K> species_coefficients(const std::map<std::string, double>& values,
                                            std::string_view role)
{
    if (values.size() > K) {
        throw ParameterError(
          std::format("{} list {} species, more than {}", role, values.size(), K));
    }
    SpeciesCoefficients<K> out {};
    out.count = static_cast<int>(values.size());
    std::size_t e = 0;
    for (const auto& [name, value] : values) {
        out.entries[e++] = SpeciesCoefficient {.name = name.c_str(), .value = value};
    }
    return out;
}

// The Arrhenius parameters of a rate record. Throws ParameterError when the
// model is not "arrhenius" or a parameter is missing.
inline ArrheniusInput arrhenius_input(const RateRecord& rate, std::string_view role)
{
    if (rate.model != "arrhenius") {
        throw ParameterError(
          std::format(R"({} rate model must be "arrhenius", got "{}")", role, rate.model));
    }
    if (!rate.A || !rate.n || !rate.activation_temperature) {
        throw ParameterError(std::format("{} rate needs A, n and activation_temperature", role));
    }
    return {
      .A = *rate.A,
      .n = *rate.n,
      .activation_temperature = *rate.activation_temperature,
      .translational_exponent = rate.translational_exponent.value_or(1.0),
    };
}

}  // namespace detail

// Model inputs of the R reactions of file, in file order, with at most K
// species in each list of species coefficients, for
// make_reaction_mechanism. Names and equations point into file, so the
// result is valid while file exists. Throws ParameterError when file does
// not have R reactions, and (with the reaction index and equation) when a
// list has more than K species, the forward rate is not an Arrhenius rate
// with A, n and activation_temperature, or the backward rate is neither
// such a rate nor "equilibrium" without parameters. The values themselves
// are checked by make_reaction_mechanism.
template <std::size_t R, std::size_t K>
    requires(R > 0 && K > 0)
MechanismInput<R, K> mechanism_input(const ReactionsFile& file)
{
    if (file.reactions.size() != R) {
        throw ParameterError(
          std::format("expected {} reactions, the file has {}", R, file.reactions.size()));
    }
    MechanismInput<R, K> out {};
    out.temperature_lower = file.temperature_limits.lower;
    out.temperature_upper = file.temperature_limits.upper;
    for (std::size_t k = 0; k < R; ++k) {
        const ReactionRecord& record = file.reactions[k];
        try {
            ReactionInput<K>& reaction = out.reactions[k];
            reaction.equation = record.equation.c_str();
            reaction.reactants = detail::species_coefficients<K>(record.reactants, "reactants");
            reaction.products = detail::species_coefficients<K>(record.products, "products");
            reaction.forward = detail::arrhenius_input(record.forward, "forward");
            if (record.backward.model == "equilibrium") {
                if (record.backward.A || record.backward.n ||
                    record.backward.activation_temperature ||
                    record.backward.translational_exponent) {
                    throw ParameterError("an equilibrium backward rate has no A, n, "
                                         "activation_temperature or translational_exponent");
                }
                reaction.backward_model = BackwardRateModel::equilibrium;
            } else {
                reaction.backward_model = BackwardRateModel::arrhenius;
                reaction.backward = detail::arrhenius_input(record.backward, "backward");
            }
            reaction.third_body = record.third_body.has_value();
            if (record.third_body) {
                reaction.efficiencies =
                  detail::species_coefficients<K>(record.third_body->efficiencies,
                                                  "third-body efficiencies");
            }
        } catch (const ParameterError& error) {
            throw ParameterError(
              std::format("reaction {} ({}): {}", k, record.equation, error.what()));
        }
    }
    return out;
}

// Model inputs of the V vibration-translation pairs and C chemistry-vibration
// couplings of file, in file order, for make_energy_exchange. Names point
// into file, so the result is valid while file exists. Throws ParameterError
// when file does not have V pairs and C couplings, and (with the entry's
// list, index and species) when a model is not "landau_teller",
// "millikan_white", "park_1993" or "impartial_dissociation" where each is
// expected. The values themselves are checked by make_energy_exchange.
template <std::size_t V, std::size_t C>
EnergyExchangeInput<V, C> energy_exchange_input(const EnergyExchangeFile& file)
{
    if (file.vibration_translation.size() != V || file.chemistry_vibration.size() != C) {
        throw ParameterError(std::format(
          "expected {} vibration-translation pairs and {} chemistry-vibration couplings, the "
          "file has {} and {}",
          V,
          C,
          file.vibration_translation.size(),
          file.chemistry_vibration.size()));
    }
    const auto expect_model =
      [](const std::string& got, std::string_view expected, std::string_view entry) {
          if (got != expected) {
              throw ParameterError(
                std::format(R"({}: model must be "{}", got "{}")", entry, expected, got));
          }
      };
    EnergyExchangeInput<V, C> out {};
    for (std::size_t k = 0; k < V; ++k) {
        const VibrationTranslationRecord& record = file.vibration_translation[k];
        const std::string entry =
          std::format("vibration_translation {} ({}-{})", k, record.species, record.partner);
        expect_model(record.model, "landau_teller", entry);
        expect_model(record.relaxation_time.model, "millikan_white", entry);
        VibrationTranslationInput& pair = out.vibration_translation[k];
        pair.species = record.species.c_str();
        pair.partner = record.partner.c_str();
        pair.millikan_white = {.a = record.relaxation_time.a, .b = record.relaxation_time.b};
        const auto& correction = record.relaxation_time.high_temperature_correction;
        pair.park_correction = correction.has_value();
        if (correction) {
            expect_model(correction->model, "park_1993", entry);
            pair.park = {
              .cross_section = correction->cross_section,
              .reference_temperature = correction->reference_temperature,
            };
        }
    }
    for (std::size_t k = 0; k < C; ++k) {
        const ChemistryVibrationRecord& record = file.chemistry_vibration[k];
        expect_model(record.model,
                     "impartial_dissociation",
                     std::format("chemistry_vibration {} ({})", k, record.species));
        out.chemistry_vibration[k] = ChemistryVibrationInput {
          .reaction = record.reaction,
          .species = record.species.c_str(),
          .dissociation_energy = record.dissociation_energy,
          .vibrational_temperature = record.vibrational_temperature,
        };
    }
    return out;
}

}  // namespace gaslib::database

#endif  // GASLIB_DATABASE_HPP
