#ifndef CANARY_DATABASE_HPP
#define CANARY_DATABASE_HPP

// Species database files: JSON files of species data, of reaction
// mechanisms, of energy exchange and of collision integrals in SI units,
// read and written by the compiled library canary_database (link
// canary::database). This header is C++20 and does
// not depend on the JSON library; the conversions to model inputs below are
// header-only.
//
// Species file layout (strict JSON):
//   {
//     "format": "canary-species-1",
//     "units": {"molar_mass": "kg/mol", "energy": "J/mol", "temperature": "K",
//               "length": "m"},
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
//         "transport": {
//           "viscosity": {"model": "cea", "temperature_ranges": [200, 1000, 5000, 15000],
//                         "coefficients": [[4 numbers], [4 numbers], [4 numbers]]},
//           "thermal_conductivity": {"model": "cea", ...},
//           "lennard_jones": {"sigma": 3.621e-10, "epsilon": 97.53},
//           "reference": "...", "notes": "..."
//         },
//         "electronic_levels": {"degeneracies": [1, 3, 6], "temperatures": [0, 72231.9, 85778.6],
//                               "reference": "...", "notes": "..."},
//         "reference": "...", "notes": "..."
//       }
//     },
//     "reference": "...", "notes": "..."
//   }
// "reference" and "notes" are optional everywhere, and so are "transport",
// each of its three records, and "electronic_levels"; every other key is required, and unknown
// keys are errors. Numbers are written in the shortest form that reads back
// to the same double.
//
// Collision-integral files hold curve fits of the collision integrals
// pi Omega(1,1) and pi Omega(2,2) of species pairs (Gupta et al. 1990 form)
// in SI units:
//   {
//     "format": "canary-collision-integrals-1",
//     "units": {"temperature": "K", "area": "m^2"},
//     "pairs": [
//       {"species": ["O2", "N2"],
//        "pi_omega_11": {"A": 0.0, "B": -0.0465, "C": 0.5729, "D": -44.4332},
//        "pi_omega_22": {"A": 0.0, "B": -0.0558, "C": 0.7590, "D": -45.1562},
//        "reference": "...", "notes": "..."}
//     ],
//     "reference": "...", "notes": "..."
//   }
// with ln(pi Omega / m^2) = A (ln T)^3 + B (ln T)^2 + C ln T + D, T in K. A
// pair applies in either order; "reference" and "notes" are optional, every
// other key is required, and unknown keys are errors.
//
// Reactions files hold a mechanism of reversible reactions (see
// reactions.hpp), in SI units per mole:
//   {
//     "format": "canary-reactions-1",
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
// Energy-exchange files hold the vibration-translation exchange, the
// chemistry-vibration coupling and the exchanges of free electrons of a
// multi-temperature model (see energy_exchange.hpp), in SI units per mole:
//   {
//     "format": "canary-energy-exchange-1",
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
//     "electron_translation": [
//       {"partner": "N2", "model": "elastic",
//        "cross_section": {"model": "gnoffo_neutral", "a": 7.5e-20, "b": 5.5e-24,
//                          "c": -1e-28, "maximum_temperature": 30000.0}},
//       {"partner": "N+", "model": "elastic", "cross_section": {"model": "coulomb"}}
//     ],
//     "vibration_electron": [
//       {"species": "N2", "model": "landau_teller",
//        "relaxation_time": {"model": "candler", "a_low": 7.5, "b_low": -57.0,
//          "c_low": 98.7, "a_high": 2.36, "b_high": -17.9, "c_high": 24.35,
//          "blend_lower": 6950.0, "blend_upper": 7050.0}}
//     ],
//     "electron_impact_ionisation": [{"reaction": 38, "species": "O"}],
//     "reference": "...", "notes": "..."
//   }
// "high_temperature_correction", the cross-section parameters, the
// correction's "reference_temperature" and "exponent", "reference" and
// "notes" are optional; every other key is required (the five lists may be
// empty), and unknown keys are errors. "reaction" is the 0-based index of
// the reaction in the mechanism used with the file; the ionisation energy of
// an electron_impact_ionisation species comes from the species file.

#include <array>
#include <cstddef>
#include <filesystem>
#include <format>
#include <canary/build.hpp>
#include <canary/energy_exchange.hpp>
#include <canary/models/multi_temperature.hpp>
#include <canary/models/thermally_perfect.hpp>
#include <canary/nasa9.hpp>
#include <canary/reactions.hpp>
#include <canary/transport/mixture.hpp>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace canary::database
{

// Value of the "format" key of a species file.
inline constexpr std::string_view species_format = "canary-species-1";

// Units of the values in a file. Only these SI units are accepted.
struct Units {
    std::string molar_mass = "kg/mol";
    std::string energy = "J/mol";
    std::string temperature = "K";
    std::string length = "m";

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

// A transport property of one species as a function of temperature. For
// model "cea" (the NASA CEA transport fits): temperature_ranges are the
// segment boundaries in K and coefficients holds A, B, C, D of each segment,
// with ln(property / SI unit) = A ln T + B / T + C / T^2 + D; the property is
// a viscosity in Pa s or a thermal conductivity in W/(m K).
struct TransportCurveRecord {
    std::string model;
    std::vector<double> temperature_ranges;
    std::vector<std::array<double, 4>> coefficients;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const TransportCurveRecord&, const TransportCurveRecord&) = default;
};

// Lennard-Jones potential of one species: collision diameter sigma in m and
// well depth epsilon / k_B in K.
struct LennardJonesRecord {
    double sigma {};
    double epsilon {};
    std::optional<std::string> reference;

    friend bool operator==(const LennardJonesRecord&, const LennardJonesRecord&) = default;
};

// Transport data of one species; each record is optional.
struct TransportRecord {
    std::optional<TransportCurveRecord> viscosity;
    std::optional<TransportCurveRecord> thermal_conductivity;
    std::optional<LennardJonesRecord> lennard_jones;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const TransportRecord&, const TransportRecord&) = default;
};

// Electronic energy levels of one species: the degeneracy g_k and the
// level energy as a temperature theta_k = E_k / k_B in K of every level.
struct ElectronicLevelsRecord {
    std::vector<int> degeneracies;
    std::vector<double> temperatures;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const ElectronicLevelsRecord&, const ElectronicLevelsRecord&) = default;
};

// One species: molar mass in kg/mol, charge in elementary charges, thermo
// data and optionally transport data, electronic levels and the ionisation
// energy in J/mol.
struct SpeciesRecord {
    double molar_mass {};
    double charge {};
    ThermoRecord thermo;
    std::optional<TransportRecord> transport;
    std::optional<ElectronicLevelsRecord> electronic_levels;
    std::optional<double> ionisation_energy;
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
inline constexpr std::string_view reactions_format = "canary-reactions-1";

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
    std::optional<std::string> other_temperature;

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
inline constexpr std::string_view energy_exchange_format = "canary-energy-exchange-1";

// Units of the values in an energy-exchange file: energies per mole in
// J/mol, temperatures in K and cross-sections in m^2. Only these are
// accepted.
struct ExchangeUnits {
    std::string energy = "J/mol";
    std::string temperature = "K";
    std::string area = "m^2";

    friend bool operator==(const ExchangeUnits&, const ExchangeUnits&) = default;
};

// High-temperature correction of a relaxation time, by its collision
// cross-section sigma in m^2: model "park_1985" (sigma constant),
// "park_1993" (sigma (reference_temperature / T)^2, reference_temperature in
// K) or "kim_2021" (sigma T^exponent, T in K).
struct HighTemperatureCorrectionRecord {
    std::string model;
    double cross_section {};
    std::optional<double> reference_temperature;
    std::optional<double> exponent;

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

// Momentum-transfer cross-section of free electrons with a heavy species:
// model "gnoffo_neutral" with sigma = a + b T_e + c T_e^2 in m^2, T_e in K
// limited to maximum_temperature; or model "coulomb" (charged partner) with
// no parameters.
struct ElectronCrossSectionRecord {
    std::string model;
    std::optional<double> a;
    std::optional<double> b;
    std::optional<double> c;
    std::optional<double> maximum_temperature;

    friend bool operator==(const ElectronCrossSectionRecord&,
                           const ElectronCrossSectionRecord&) = default;
};

// Elastic energy exchange between free electrons and a heavy partner: model
// "elastic".
struct ElectronTranslationRecord {
    std::string partner;
    std::string model;
    ElectronCrossSectionRecord cross_section;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const ElectronTranslationRecord&,
                           const ElectronTranslationRecord&) = default;
};

// Electron-vibration relaxation time of Candler's form: log10(p_e tau) with
// p_e in atm and tau in s is a_low x^2 + b_low x + c_low below blend_lower
// and a_high x^2 + b_high x + c_high above blend_upper, x = log10(T_e / K),
// blended in between.
struct CandlerRelaxationRecord {
    std::string model;
    double a_low {};
    double b_low {};
    double c_low {};
    double a_high {};
    double b_high {};
    double c_high {};
    double blend_lower {};
    double blend_upper {};

    friend bool operator==(const CandlerRelaxationRecord&,
                           const CandlerRelaxationRecord&) = default;
};

// Vibration-electron exchange of species: model "landau_teller" with a
// relaxation time of model "candler".
struct VibrationElectronRecord {
    std::string species;
    std::string model;
    CandlerRelaxationRecord relaxation_time;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const VibrationElectronRecord&,
                           const VibrationElectronRecord&) = default;
};

// Electron-impact ionisation of species in the reaction with 0-based index
// reaction: the electrons lose the species' ionisation energy per event.
struct ElectronImpactIonisationRecord {
    int reaction {};
    std::string species;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const ElectronImpactIonisationRecord&,
                           const ElectronImpactIonisationRecord&) = default;
};

// An energy-exchange file. Entries keep their order; any list may be empty.
struct EnergyExchangeFile {
    std::string format {energy_exchange_format};
    ExchangeUnits units;
    std::vector<VibrationTranslationRecord> vibration_translation;
    std::vector<ChemistryVibrationRecord> chemistry_vibration;
    std::vector<ElectronTranslationRecord> electron_translation;
    std::vector<VibrationElectronRecord> vibration_electron;
    std::vector<ElectronImpactIonisationRecord> electron_impact_ionisation;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const EnergyExchangeFile&, const EnergyExchangeFile&) = default;
};

// Value of the "format" key of a collision-integral file.
inline constexpr std::string_view collision_integrals_format = "canary-collision-integrals-1";

// Units of the values in a collision-integral file: temperatures in K and
// collision integrals in m^2. Only these are accepted.
struct CollisionUnits {
    std::string temperature = "K";
    std::string area = "m^2";

    friend bool operator==(const CollisionUnits&, const CollisionUnits&) = default;
};

// Curve fit ln(pi Omega / m^2) = A (ln T)^3 + B (ln T)^2 + C ln T + D, T in K.
struct CollisionFitRecord {
    double A {};
    double B {};
    double C {};
    double D {};

    friend bool operator==(const CollisionFitRecord&, const CollisionFitRecord&) = default;
};

// The collision integrals pi Omega(1,1) and pi Omega(2,2) of a pair of
// species, named in either order.
struct CollisionPairRecord {
    std::array<std::string, 2> species;
    CollisionFitRecord pi_omega_11;
    CollisionFitRecord pi_omega_22;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const CollisionPairRecord&, const CollisionPairRecord&) = default;
};

// A collision-integral file. Pairs keep their order.
struct CollisionIntegralsFile {
    std::string format {collision_integrals_format};
    CollisionUnits units;
    std::vector<CollisionPairRecord> pairs;
    std::optional<std::string> reference;
    std::optional<std::string> notes;

    friend bool operator==(const CollisionIntegralsFile&, const CollisionIntegralsFile&) = default;
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

// The collision-integral file in json. Throws DatabaseError.
CollisionIntegralsFile parse_collision_integrals(std::string_view json);

// The collision-integral file at path. Throws DatabaseError.
CollisionIntegralsFile read_collision_integrals_file(const std::filesystem::path& path);

// file as indented JSON.
std::string write_collision_integrals(const CollisionIntegralsFile& file);

// Writes file as indented JSON to path. Throws DatabaseError when the file
// cannot be written.
void write_collision_integrals_file(const std::filesystem::path& path,
                                    const CollisionIntegralsFile& file);

// JSON Schema of collision-integral files.
std::string collision_integrals_schema();

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

// Converted model inputs can contain pointers into file-owned strings.
// Keep each source file alive and its records unmodified until model
// construction from those inputs finishes.
//
// Model inputs for the named species of file, in the given order, for
// make_thermally_perfect_mixture. Each name points to the species' key in
// file, so the result is valid while file exists and its records remain unmodified. Throws
// DatabaseError for a name that file lacks, and ParameterError (with the species name) for a thermo
// record that nasa9_data rejects.
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
      .other_temperature = rate.other_temperature ? rate.other_temperature->c_str() : nullptr,
    };
}

}  // namespace detail

// Model inputs of the R reactions of file, in file order, with at most K
// species in each list of species coefficients, for
// make_reaction_mechanism. Names and equations point into file, so the
// result is valid while file exists and its records remain unmodified. Throws ParameterError when
// file does not have R reactions, and (with the reaction index and equation) when a list has more
// than K species, the forward rate is not an Arrhenius rate with A, n and activation_temperature,
// or the backward rate is neither such a rate nor "equilibrium" without parameters. The values
// themselves are checked by make_reaction_mechanism.
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
                if (record.backward.other_temperature) {
                    reaction.backward.other_temperature =
                      record.backward.other_temperature->c_str();
                }
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
// into file, so the result is valid while file exists and its records remain unmodified. Throws
// ParameterError when file does not have V pairs and C couplings, and (with the entry's list, index
// and species) when a model is not "landau_teller", "millikan_white", one of the correction models
// or "impartial_dissociation" where each is expected. The values themselves are checked by
// make_energy_exchange.
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
        pair.correction = {.model = RelaxationCorrection::none};
        if (correction) {
            pair.correction.cross_section = correction->cross_section;
            pair.correction.reference_temperature = correction->reference_temperature.value_or(0.0);
            pair.correction.exponent = correction->exponent.value_or(0.0);
            if (correction->model == "park_1985") {
                pair.correction.model = RelaxationCorrection::park_1985;
            } else if (correction->model == "park_1993") {
                pair.correction.model = RelaxationCorrection::park_1993;
            } else if (correction->model == "kim_2021") {
                pair.correction.model = RelaxationCorrection::kim_2021;
            } else {
                throw ParameterError(std::format(
                  R"({}: correction model must be "park_1985", "park_1993" or "kim_2021", got "{}")",
                  entry,
                  correction->model));
            }
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

// Model inputs of the E elastic electron-heavy pairs, VE vibration-electron
// exchanges and I electron-impact ionisations of file, in file order, for
// make_multi_temperature_exchange, with the ionisation energies from
// species. Names point into file, so the result is valid while file exists and its records remain
// unmodified. Throws ParameterError when the counts differ, a model is not "elastic" with
// "gnoffo_neutral" (a, b, c and maximum_temperature given) or "coulomb", or "landau_teller" with
// "candler", or (with the species) when an ionised species has no ionisation_energy in species;
// DatabaseError when that species is not in species.
template <std::size_t E, std::size_t VE, std::size_t I>
ElectronExchangeInput<E, VE, I> electron_exchange_input(const EnergyExchangeFile& file,
                                                        const SpeciesFile& species)
{
    if (file.electron_translation.size() != E || file.vibration_electron.size() != VE ||
        file.electron_impact_ionisation.size() != I) {
        throw ParameterError(
          std::format("expected {}, {} and {} electron exchanges, the file has {}, {} and {}",
                      E,
                      VE,
                      I,
                      file.electron_translation.size(),
                      file.vibration_electron.size(),
                      file.electron_impact_ionisation.size()));
    }
    ElectronExchangeInput<E, VE, I> out {};
    for (std::size_t k = 0; k < E; ++k) {
        const ElectronTranslationRecord& record = file.electron_translation[k];
        const ElectronCrossSectionRecord& section = record.cross_section;
        ElectronTranslationInput& entry = out.electron_translation[k];
        entry.partner = record.partner.c_str();
        if (record.model == "elastic" && section.model == "gnoffo_neutral" && section.a &&
            section.b && section.c && section.maximum_temperature) {
            entry.cross_section = ElectronCrossSection::gnoffo_neutral;
            entry.a = *section.a;
            entry.b = *section.b;
            entry.c = *section.c;
            entry.maximum_temperature = *section.maximum_temperature;
        } else if (record.model == "elastic" && section.model == "coulomb") {
            entry.cross_section = ElectronCrossSection::coulomb;
        } else {
            throw ParameterError(std::format(
              "electron_translation {} ({}): needs model \"elastic\" with a \"gnoffo_neutral\" "
              "cross-section (a, b, c, maximum_temperature) or a \"coulomb\" one",
              k,
              record.partner));
        }
    }
    for (std::size_t k = 0; k < VE; ++k) {
        const VibrationElectronRecord& record = file.vibration_electron[k];
        const CandlerRelaxationRecord& fit = record.relaxation_time;
        if (record.model != "landau_teller" || fit.model != "candler") {
            throw ParameterError(std::format(
              "vibration_electron {} ({}): needs model \"landau_teller\" with a \"candler\" "
              "relaxation time",
              k,
              record.species));
        }
        out.vibration_electron[k] = VibrationElectronInput {
          .species = record.species.c_str(),
          .relaxation_time =
            {
              .a_low = fit.a_low,
              .b_low = fit.b_low,
              .c_low = fit.c_low,
              .a_high = fit.a_high,
              .b_high = fit.b_high,
              .c_high = fit.c_high,
              .blend_lower = fit.blend_lower,
              .blend_upper = fit.blend_upper,
            },
        };
    }
    for (std::size_t k = 0; k < I; ++k) {
        const ElectronImpactIonisationRecord& record = file.electron_impact_ionisation[k];
        const auto found = species.species.find(record.species);
        if (found == species.species.end()) {
            throw DatabaseError(std::format("species \"{}\" is not in the file", record.species));
        }
        if (!found->second.ionisation_energy) {
            throw ParameterError(std::format(
              "electron_impact_ionisation {} ({}): the species has no ionisation_energy",
              k,
              record.species));
        }
        out.electron_impact_ionisation[k] = ElectronImpactIonisationInput {
          .reaction = record.reaction,
          .species = record.species.c_str(),
          .ionisation_energy = *found->second.ionisation_energy,
        };
    }
    return out;
}

namespace detail
{

// A CEA transport curve of a species record with at most S segments.
// Throws ParameterError when the model is not "cea", there are more than S
// segments or none, or the numbers of ranges and coefficient sets do not
// match. The values themselves are checked by the transport builders.
template <std::size_t S>
TransportCurveInput<S> transport_curve(const TransportCurveRecord& curve, std::string_view role)
{
    if (curve.model != "cea") {
        throw ParameterError(std::format(R"({} model must be "cea", got "{}")", role, curve.model));
    }
    const std::size_t segments = curve.coefficients.size();
    if (segments == 0 || segments > S || curve.temperature_ranges.size() != segments + 1) {
        throw ParameterError(std::format(
          "{} needs between 1 and {} segments with one more temperature range, got {} and {}",
          role,
          S,
          segments,
          curve.temperature_ranges.size()));
    }
    TransportCurveInput<S> out {};
    out.segments = static_cast<int>(segments);
    for (std::size_t k = 0; k <= segments; ++k) { out.breaks[k] = curve.temperature_ranges[k]; }
    for (std::size_t k = 0; k < segments; ++k) { out.coefficients[k] = curve.coefficients[k]; }
    return out;
}

}  // namespace detail

// Transport inputs for the named species of file, in the given order, for
// make_wilke_transport, with at most S segments per curve. Each name points
// to the species' key in file, so the result is valid while file exists and its records remain
// unmodified. Throws DatabaseError for a name that file lacks, and ParameterError (with the species
// name) when a species has no transport record, lacks its viscosity, thermal conductivity or
// Lennard-Jones record, or has a curve that detail::transport_curve rejects.
template <std::size_t S, std::size_t N>
std::array<SpeciesTransportInput<S>, N>
  species_transport(const SpeciesFile& file, const std::array<std::string_view, N>& names)
{
    std::array<SpeciesTransportInput<S>, N> out {};
    for (std::size_t i = 0; i < N; ++i) {
        const auto entry = file.species.find(std::string(names[i]));
        if (entry == file.species.end()) {
            throw DatabaseError(std::format("species \"{}\" is not in the file", names[i]));
        }
        const std::optional<TransportRecord>& transport = entry->second.transport;
        try {
            if (!transport || !transport->viscosity || !transport->thermal_conductivity ||
                !transport->lennard_jones) {
                throw ParameterError("needs viscosity, thermal_conductivity and lennard_jones "
                                     "transport records");
            }
            out[i] = SpeciesTransportInput<S> {
              .name = entry->first.c_str(),
              .viscosity = detail::transport_curve<S>(*transport->viscosity, "viscosity"),
              .conductivity = detail::transport_curve<S>(*transport->thermal_conductivity,
                                                         "thermal_conductivity"),
              .lennard_jones =
                {
                  .sigma = transport->lennard_jones->sigma,
                  .epsilon = transport->lennard_jones->epsilon,
                },
            };
        } catch (const ParameterError& error) {
            throw ParameterError(std::format("species {}: {}", entry->first, error.what()));
        }
    }
    return out;
}

// Electronic levels of the named species of file, in the given order, for
// make_multi_temperature_mixture, with at most L levels each; a species
// without an electronic_levels record has zero levels. Each name points to
// the species' key in file, so the result is valid while file exists and its records remain
// unmodified. Throws DatabaseError for a name that file lacks, and ParameterError (with the species
// name) when a species has more than L levels or different numbers of degeneracies and
// temperatures. The values themselves are checked by make_multi_temperature_mixture.
template <std::size_t L, std::size_t N>
std::array<ElectronicLevelsInput<L>, N>
  electronic_levels(const SpeciesFile& file, const std::array<std::string_view, N>& names)
{
    std::array<ElectronicLevelsInput<L>, N> out {};
    for (std::size_t i = 0; i < N; ++i) {
        const auto entry = file.species.find(std::string(names[i]));
        if (entry == file.species.end()) {
            throw DatabaseError(std::format("species \"{}\" is not in the file", names[i]));
        }
        out[i].name = entry->first.c_str();
        if (!entry->second.electronic_levels) { continue; }
        const ElectronicLevelsRecord& levels = *entry->second.electronic_levels;
        const std::size_t count = levels.degeneracies.size();
        if (levels.temperatures.size() != count || count > L) {
            throw ParameterError(
              std::format("species {}: needs at most {} levels with as many degeneracies as "
                          "temperatures, got {} and {}",
                          entry->first,
                          L,
                          count,
                          levels.temperatures.size()));
        }
        out[i].levels = static_cast<int>(count);
        for (std::size_t k = 0; k < count; ++k) {
            out[i].degeneracies[k] = levels.degeneracies[k];
            out[i].temperatures[k] = levels.temperatures[k];
        }
    }
    return out;
}

// Model inputs of the K pairs of file, in file order, for
// make_gupta_yos_transport. Names point into file, so the result is valid while file exists and its
// records remain unmodified. Throws ParameterError when file does not have K pairs.
template <std::size_t K>
std::array<CollisionPairInput, K> collision_integral_pairs(const CollisionIntegralsFile& file)
{
    if (file.pairs.size() != K) {
        throw ParameterError(
          std::format("expected {} pairs, the file has {}", K, file.pairs.size()));
    }
    const auto fit = [](const CollisionFitRecord& r) {
        return CollisionFitInput {.A = r.A, .B = r.B, .C = r.C, .D = r.D};
    };
    std::array<CollisionPairInput, K> out {};
    for (std::size_t k = 0; k < K; ++k) {
        const CollisionPairRecord& pair = file.pairs[k];
        out[k] = CollisionPairInput {
          .first = pair.species[0].c_str(),
          .second = pair.species[1].c_str(),
          .pi_omega_11 = fit(pair.pi_omega_11),
          .pi_omega_22 = fit(pair.pi_omega_22),
        };
    }
    return out;
}

}  // namespace canary::database

#endif  // CANARY_DATABASE_HPP
