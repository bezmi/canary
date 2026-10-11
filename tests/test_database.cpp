// Tests of the species database library: reading and writing species,
// reactions and energy-exchange files, errors, the generated schemas,
// conversion to model inputs, cp/R of the shipped five-species air file
// against reference values, and the compiled Kim-Jo mechanism, Park
// mechanism, Park energy exchange, species transport data, collision
// integrals and electronic levels against their files.

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <gaslib/build.hpp>
#include <gaslib/data/air11_kim_jo.hpp>
#include <gaslib/data/air11_kim_jo_exchange.hpp>
#include <gaslib/data/air11_species.hpp>
#include <gaslib/data/air5.hpp>
#include <gaslib/data/air5_electronic.hpp>
#include <gaslib/data/air5_kim_jo.hpp>
#include <gaslib/data/air5_park.hpp>
#include <gaslib/data/air5_park_exchange.hpp>
#include <gaslib/data/air5_transport.hpp>
#include <gaslib/database.hpp>
#include <gaslib/energy_exchange.hpp>
#include <gaslib/nasa9.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/reactions.hpp>
#include <gaslib/transport/mixture.hpp>
#include <gaslib/units.hpp>
#include <ios>
#include <iterator>
#include <map>
#include <optional>
#include <string>
#include <string_view>

namespace
{

namespace db = gaslib::database;

std::filesystem::path source_path(std::string_view relative)
{ return std::filesystem::path(GASLIB_SOURCE_DIR) / relative; }

std::filesystem::path air5_path() { return source_path("database/data/air5.json"); }

std::filesystem::path air11_path() { return source_path("database/data/air11.json"); }

std::filesystem::path air11_kim_jo_path() { return source_path("database/data/air11_kim_jo.json"); }

std::filesystem::path air11_exchange_path()
{ return source_path("database/data/air11_kim_jo_exchange.json"); }

std::filesystem::path kim_jo_path() { return source_path("database/data/air5_kim_jo.json"); }

std::filesystem::path park_path() { return source_path("database/data/air5_park.json"); }

std::filesystem::path exchange_path()
{ return source_path("database/data/air5_park_exchange.json"); }

std::filesystem::path collision_path()
{ return source_path("database/data/air5_collision_integrals.json"); }

std::string read_text(const std::filesystem::path& path)
{
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

// air5.json with one textual replacement.
std::string air5_with(std::string_view from, std::string_view to)
{
    std::string text = read_text(air5_path());
    const auto at = text.find(from);
    if (at != std::string::npos) { text.replace(at, from.size(), to); }
    return text;
}

// air5_kim_jo.json with one textual replacement.
std::string kim_jo_with(std::string_view from, std::string_view to)
{
    std::string text = read_text(kim_jo_path());
    const auto at = text.find(from);
    if (at != std::string::npos) { text.replace(at, from.size(), to); }
    return text;
}

// The same species coefficients in any order.
template <std::size_t K>
bool same_coefficients(const gaslib::SpeciesCoefficients<K>& a,
                       const gaslib::SpeciesCoefficients<K>& b)
{
    if (a.count != b.count) { return false; }
    for (std::size_t i = 0; i < static_cast<std::size_t>(a.count); ++i) {
        bool found = false;
        for (std::size_t j = 0; j < static_cast<std::size_t>(b.count); ++j) {
            found = found || (std::string_view(a.entries[i].name) == b.entries[j].name &&
                              a.entries[i].value == b.entries[j].value);
        }
        if (!found) { return false; }
    }
    return true;
}

// The same other temperature (both absent, or the same name).
bool same_other(const gaslib::ArrheniusInput& a, const gaslib::ArrheniusInput& b)
{
    if (a.other_temperature == nullptr || b.other_temperature == nullptr) {
        return a.other_temperature == b.other_temperature;
    }
    return std::string_view(a.other_temperature) == b.other_temperature;
}

bool same_correction(const gaslib::RelaxationCorrectionInput& a,
                     const gaslib::RelaxationCorrectionInput& b)
{
    return a.model == b.model && a.cross_section == b.cross_section &&
           a.reference_temperature == b.reference_temperature && a.exponent == b.exponent;
}

bool same_rate(const gaslib::ArrheniusInput& a, const gaslib::ArrheniusInput& b)
{
    return a.A == b.A && a.n == b.n && a.activation_temperature == b.activation_temperature &&
           a.translational_exponent == b.translational_exponent;
}

// air5_collision_integrals.json with one textual replacement.
std::string collision_with(std::string_view from, std::string_view to)
{
    std::string text = read_text(collision_path());
    const auto at = text.find(from);
    if (at != std::string::npos) { text.replace(at, from.size(), to); }
    return text;
}

// air5_park_exchange.json with one textual replacement.
std::string exchange_with(std::string_view from, std::string_view to)
{
    std::string text = read_text(exchange_path());
    const auto at = text.find(from);
    if (at != std::string::npos) { text.replace(at, from.size(), to); }
    return text;
}

template <typename F> bool throws_database_error(const F& f)
{
    try {
        f();
    } catch (const db::DatabaseError&) {
        return true;
    }
    return false;
}

template <typename F> bool throws_parameter_error(const F& f)
{
    try {
        f();
    } catch (const gaslib::ParameterError&) {
        return true;
    }
    return false;
}

constexpr std::array<std::string_view, 5> air5_names {"N2", "O2", "N", "O", "NO"};
constexpr std::array<std::string_view, 11> air11_names {
  "N2",
  "O2",
  "N",
  "O",
  "NO",
  "N+",
  "O+",
  "N2+",
  "O2+",
  "NO+",
  "e-",
};

double cp_over_r(const gaslib::Nasa9Curve<double, 4>& curve, double t)
{ return gaslib::nasa9_cp_h(curve, t, std::log(t)).cp_over_r; }

}  // namespace

int main()
{
    using namespace boost::ut;

    "the shipped five-species air file reads with its five species"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        expect(file.format == db::species_format);
        expect(file.units == db::Units {});
        expect(file.species.size() == 5_u);
        for (const auto name : air5_names) {
            expect(file.species.contains(std::string(name))) << name;
        }
        expect(file.species.at("N2").thermo.model == "nasa9");
    };

    "read, write and read again gives the same data and the same text"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        const std::string text = db::write_species(file);
        const db::SpeciesFile again = db::parse_species(text);
        expect(again == file);
        expect(db::write_species(again) == text);
    };

    "writing to a file and reading it back gives the same data"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        const auto path = std::filesystem::temp_directory_path() / "gaslib_test_species.json";
        db::write_species_file(path, file);
        expect(db::read_species_file(path) == file);
        std::filesystem::remove(path);
    };

    "unknown keys are errors at every level"_test = [] {
        expect(throws_database_error(
          [] { db::parse_species(air5_with("\"format\"", R"("extra": 1, "format")")); }));
        expect(throws_database_error([] {
            db::parse_species(air5_with("\"molar_mass\": 0.0280134", "\"mass\": 0.0280134"));
        }));
        expect(throws_database_error(
          [] { db::parse_species(air5_with("\"model\"", R"("curve": 1, "model")")); }));
    };

    "a missing required key is an error and missing notes are not"_test = [] {
        expect(throws_database_error(
          [] { db::parse_species(air5_with("\"molar_mass\": 0.0280134,", "")); }));
        expect(
          throws_database_error([] { db::parse_species(air5_with(R"("model": "nasa9",)", "")); }));
        const db::SpeciesFile file = db::read_species_file(air5_path());
        db::SpeciesFile no_notes = file;
        for (auto& [name, record] : no_notes.species) { record.notes.reset(); }
        expect(db::parse_species(db::write_species(no_notes)) == no_notes);
    };

    "another format, other units, invalid JSON or a missing file are errors"_test = [] {
        expect(throws_database_error(
          [] { db::parse_species(air5_with("gaslib-species-1", "gaslib-species-2")); }));
        expect(
          throws_database_error([] { db::parse_species(air5_with("\"kg/mol\"", "\"g/mol\"")); }));
        expect(throws_database_error([] { db::parse_species(air5_with("\"K\"", "\"degC\"")); }));
        expect(throws_database_error([] { db::parse_species("{\"format\": "); }));
        expect(
          throws_database_error([] { db::read_species_file(source_path("no_such_file.json")); }));
    };

    "the shipped schema is the schema of the current file structs"_test = [] {
        expect(db::species_schema() == read_text(source_path("database/species.schema.json")));
    };

    "model inputs from the file equal the compiled five-species air data"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        const auto from_file = db::nasa9_species<4>(file, air5_names);
        for (std::size_t i = 0; i < 5; ++i) {
            const auto& a = from_file[i];
            const auto& b = gaslib::air5_species[i];
            expect(std::string_view(a.name) == std::string_view(b.name));
            expect(a.molar_mass == b.molar_mass && a.charge == b.charge) << b.name;
            expect(a.fit.segments == b.fit.segments) << b.name;
            const auto n = static_cast<std::size_t>(b.fit.segments);
            for (std::size_t k = 0; k <= n; ++k) {
                expect(a.fit.breaks[k] == b.fit.breaks[k]) << b.name;
            }
            for (std::size_t k = 0; k + 1 < n; ++k) {
                expect(a.fit.blends[k] == b.fit.blends[k]) << b.name;
            }
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t k = 0; k < 9; ++k) {
                    expect(a.fit.a[j][k] == b.fit.a[j][k]) << b.name;
                }
            }
        }
    };

    "conversion errors name the problem"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        expect(throws_database_error(
          [&] { (void)db::nasa9_species<4>(file, std::array<std::string_view, 1> {"Ar"}); }));
        expect(throws_parameter_error(
          [&] { (void)db::nasa9_species<3>(file, std::array<std::string_view, 1> {"N"}); }));
        db::ThermoRecord thermo = file.species.at("N2").thermo;
        thermo.model = "nasa7";
        expect(throws_parameter_error([&] { (void)db::nasa9_data<4>(thermo); }));
        thermo = file.species.at("N2").thermo;
        thermo.blend_widths.pop_back();
        expect(throws_parameter_error([&] { (void)db::nasa9_data<4>(thermo); }));
        thermo = file.species.at("N2").thermo;
        thermo.temperature_ranges.pop_back();
        expect(throws_parameter_error([&] { (void)db::nasa9_data<4>(thermo); }));
    };

    // cp/R computed by GDTk (libgas.so, gdtk-uq/gdtk commit bf559c3f) from
    // the same N2 and O2 data, with GDTk's own gas constant: an independent
    // evaluation of the same fits, including the blending regions at 1000 K.
    "cp/R of N2 and O2 equals GDTk's evaluation of the same data"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        const auto species =
          db::nasa9_species<4>(file, std::array<std::string_view, 2> {"N2", "O2"});
        constexpr std::array temps {300.0, 1000.0, 5000.0, 15000.0};
        constexpr std::array<std::array<double, 4>, 2> gdtk {
          {
            {3.5029350227462328, 3.9324555662029987, 4.5621214841400004, 7.9038660786805508},
            {3.534485021788023, 4.1953819290074996, 5.1713062458949999, 3.9642060171736162},
          },
        };
        for (std::size_t s = 0; s < 2; ++s) {
            const auto curve = gaslib::make_nasa9_curve(species[s].fit, 1.0);
            for (std::size_t k = 0; k < temps.size(); ++k) {
                const double cp = cp_over_r(curve, temps[k]);
                expect(std::abs(cp - gdtk[s][k]) <= 1e-12 * gdtk[s][k])
                  << species[s].name << temps[k] << cp;
            }
        }
    };

    // NIST-JANAF Cp in J/(mol K). The fits differ from JANAF by up to 0.75 %
    // here (O2 at 5000 K); the tolerance matches the curve tests'.
    "Cp of N2 and O2 agrees with JANAF up to its 6000 K limit"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        const auto species =
          db::nasa9_species<4>(file, std::array<std::string_view, 2> {"N2", "O2"});
        constexpr std::array temps {300.0, 1000.0, 5000.0, 6000.0};
        constexpr std::array<std::array<double, 4>, 2> janaf {
          {
            {29.125, 32.697, 37.912, 38.276},
            {29.385, 34.870, 42.675, 44.387},
          },
        };
        const double r_u = gaslib::constants::universal_gas_constant<double>;
        for (std::size_t s = 0; s < 2; ++s) {
            const auto curve = gaslib::make_nasa9_curve(species[s].fit, 1.0);
            for (std::size_t k = 0; k < temps.size(); ++k) {
                const double cp = r_u * cp_over_r(curve, temps[k]);
                expect(std::abs(cp - janaf[s][k]) <= 1.5e-2 * janaf[s][k])
                  << species[s].name << temps[k] << cp;
            }
        }
    };

    "the shipped Kim-Jo file reads with its 13 reactions"_test = [] {
        const db::ReactionsFile file = db::read_reactions_file(kim_jo_path());
        expect(file.format == db::reactions_format);
        expect(file.units == db::ReactionUnits {});
        expect(file.temperature_limits == db::TemperatureLimits {.lower = 300.0, .upper = 30000.0});
        expect(file.reactions.size() == 13_u);
        expect(file.reactions[0].reactants == std::map<std::string, double> {{"N2", 2.0}});
        expect(file.reactions[10].third_body.has_value());
        db::RateRecord equilibrium {};
        equilibrium.model = "equilibrium";
        expect(file.reactions[12].backward == equilibrium);
    };

    "a reactions file reads, writes and reads again to the same data and text"_test = [] {
        const db::ReactionsFile file = db::read_reactions_file(kim_jo_path());
        const std::string text = db::write_reactions(file);
        const db::ReactionsFile again = db::parse_reactions(text);
        expect(again == file);
        expect(db::write_reactions(again) == text);
        const auto path = std::filesystem::temp_directory_path() / "gaslib_test_reactions.json";
        db::write_reactions_file(path, file);
        expect(db::read_reactions_file(path) == file);
        std::filesystem::remove(path);
    };

    "unknown and missing keys of a reactions file are errors"_test = [] {
        expect(throws_database_error(
          [] { db::parse_reactions(kim_jo_with("\"format\"", R"("extra": 1, "format")")); }));
        expect(throws_database_error(
          [] { db::parse_reactions(kim_jo_with("\"equation\"", R"("note": "", "equation")")); }));
        expect(throws_database_error(
          [] { db::parse_reactions(kim_jo_with("\"activation_temperature\"", "\"theta\"")); }));
        expect(throws_database_error(
          [] { db::parse_reactions(kim_jo_with("\"efficiencies\"", "\"efficiency\"")); }));
        expect(throws_database_error(
          [] { db::parse_reactions(kim_jo_with(R"("lower": 300.0,)", "")); }));
        expect(throws_database_error([] {
            db::parse_reactions(kim_jo_with(R"("equation": "N2 + N2 <=> N + N + N2",)", ""));
        }));
        expect(throws_database_error(
          [] { db::parse_reactions(kim_jo_with(R"("model": "equilibrium")", "")); }));
        db::ReactionsFile file = db::read_reactions_file(kim_jo_path());
        for (auto& reaction : file.reactions) { reaction.label.reset(); }
        file.reference.reset();
        expect(db::parse_reactions(db::write_reactions(file)) == file);
    };

    "another reactions format, other units, or a missing file are errors"_test = [] {
        expect(throws_database_error(
          [] { db::parse_reactions(kim_jo_with("gaslib-reactions-1", "gaslib-reactions-2")); }));
        expect(throws_database_error([] {
            db::parse_reactions(kim_jo_with(R"("quantity": "mol")", R"("quantity": "kmol")"));
        }));
        expect(throws_database_error(
          [] { db::parse_reactions(kim_jo_with(R"("length": "m")", R"("length": "cm")")); }));
        expect(throws_database_error(
          [] { db::parse_reactions(kim_jo_with(R"("time": "s")", R"("time": "ms")")); }));
        expect(throws_database_error([] {
            db::parse_reactions(kim_jo_with(R"("temperature": "K")", R"("temperature": "degC")"));
        }));
        expect(
          throws_database_error([] { db::read_reactions_file(source_path("no_such_file.json")); }));
        expect(throws_database_error([] { db::parse_reactions(read_text(air5_path())); }));
    };

    "the shipped reactions schema is the schema of the current file structs"_test = [] {
        expect(db::reactions_schema() == read_text(source_path("database/reactions.schema.json")));
    };

    "model inputs from the Kim-Jo file equal the compiled data"_test = [] {
        const db::ReactionsFile file = db::read_reactions_file(kim_jo_path());
        const auto input = db::mechanism_input<13, 5>(file);
        const auto& compiled = gaslib::air5_kim_jo;
        expect(input.temperature_lower == compiled.temperature_lower);
        expect(input.temperature_upper == compiled.temperature_upper);
        for (std::size_t k = 0; k < 13; ++k) {
            const auto& a = input.reactions[k];
            const auto& b = compiled.reactions[k];
            expect(std::string_view(a.equation) == b.equation) << k;
            expect(same_coefficients(a.reactants, b.reactants)) << k;
            expect(same_coefficients(a.products, b.products)) << k;
            expect(same_rate(a.forward, b.forward)) << k;
            expect(a.backward_model == b.backward_model) << k;
            expect(same_rate(a.backward, b.backward)) << k;
            expect(a.third_body == b.third_body) << k;
            expect(same_coefficients(a.efficiencies, b.efficiencies)) << k;
        }
        // The mechanism built from the file is the compiled one.
        const gaslib::UnitScales<double, gaslib::si> scales {};
        const auto built = gaslib::make_air5<double>(scales);
        const auto from_file =
          gaslib::make_reaction_mechanism(input, air5_names, built.species, scales);
        const auto from_header = gaslib::make_air5_kim_jo(built.species, scales);
        for (std::size_t k = 0; k < 13; ++k) {
            const auto& a = from_file.reactions[k];
            const auto& b = from_header.reactions[k];
            expect(a.reactants == b.reactants && a.products == b.products &&
                   a.coefficient_change == b.coefficient_change)
              << k;
            expect(a.forward.A == b.forward.A && a.forward.n == b.forward.n &&
                   a.forward.activation_temperature == b.forward.activation_temperature)
              << k;
            expect(a.backward_model == b.backward_model && a.third_body == b.third_body) << k;
            expect(a.efficiencies == b.efficiencies) << k;
        }
    };

    "reactions conversion errors name the reaction and the problem"_test = [] {
        const db::ReactionsFile file = db::read_reactions_file(kim_jo_path());
        expect(throws_parameter_error([&] { (void)db::mechanism_input<12, 5>(file); }));
        // NO + M lists five efficiencies.
        expect(throws_parameter_error([&] { (void)db::mechanism_input<13, 4>(file); }));
        const auto message = [](const db::ReactionsFile& f) -> std::string {
            try {
                (void)db::mechanism_input<13, 5>(f);
            } catch (const gaslib::ParameterError& error) {
                return error.what();
            }
            return "no exception";
        };
        db::ReactionsFile changed = file;
        changed.reactions[2].forward = file.reactions[2].backward;
        expect(message(changed).starts_with(
          R"(reaction 2 (N2 + NO <=> N + N + NO): forward rate model must be "arrhenius")"))
          << message(changed);
        changed = file;
        changed.reactions[4].forward.A.reset();
        expect(message(changed).starts_with("reaction 4")) << message(changed);
        changed = file;
        changed.reactions[6].backward.A = 1.0;
        expect(message(changed).starts_with("reaction 6")) << message(changed);
        changed = file;
        changed.reactions[7].backward.model = "power";
        expect(message(changed).starts_with("reaction 7")) << message(changed);
        // Values are checked by the builder.
        changed = file;
        changed.reactions[8].reactants = {{"N2", 1.0}, {"Ar", 1.0}};
        const auto input = db::mechanism_input<13, 5>(changed);
        const gaslib::UnitScales<double, gaslib::si> scales {};
        const auto built = gaslib::make_air5<double>(scales);
        expect(throws_parameter_error([&] {
            (void)gaslib::make_reaction_mechanism(input, air5_names, built.species, scales);
        }));
        changed = file;
        changed.reactions[8].forward.A = -1.0;
        const auto negative = db::mechanism_input<13, 5>(changed);
        expect(throws_parameter_error([&] {
            (void)gaslib::make_reaction_mechanism(negative, air5_names, built.species, scales);
        }));
    };

    "model inputs from the Park file equal the compiled data, with translational exponents"_test =
      [] {
          const db::ReactionsFile file = db::read_reactions_file(park_path());
          expect(file.reactions[0].forward.translational_exponent == 0.5);
          expect(!file.reactions[3].forward.translational_exponent.has_value());
          const auto input = db::mechanism_input<5, 5>(file);
          const auto& compiled = gaslib::air5_park;
          expect(input.temperature_lower == compiled.temperature_lower);
          expect(input.temperature_upper == compiled.temperature_upper);
          for (std::size_t k = 0; k < 5; ++k) {
              const auto& a = input.reactions[k];
              const auto& b = compiled.reactions[k];
              expect(std::string_view(a.equation) == b.equation) << k;
              expect(same_coefficients(a.reactants, b.reactants)) << k;
              expect(same_coefficients(a.products, b.products)) << k;
              expect(same_rate(a.forward, b.forward)) << k;
              expect(a.backward_model == b.backward_model) << k;
              expect(a.third_body == b.third_body) << k;
              expect(same_coefficients(a.efficiencies, b.efficiencies)) << k;
          }
          expect(input.reactions[0].forward.translational_exponent == 0.5);
          expect(input.reactions[3].forward.translational_exponent == 1.0);
          // The exponent survives writing and reading.
          expect(db::parse_reactions(db::write_reactions(file)) == file);
      };

    "an equilibrium backward rate with a translational exponent is an error"_test = [] {
        db::ReactionsFile file = db::read_reactions_file(park_path());
        file.reactions[1].backward.translational_exponent = 0.5;
        expect(throws_parameter_error([&] { (void)db::mechanism_input<5, 5>(file); }));
    };

    "the shipped energy-exchange file reads with its pairs and couplings"_test = [] {
        const db::EnergyExchangeFile file = db::read_energy_exchange_file(exchange_path());
        expect(file.format == db::energy_exchange_format);
        expect(file.units == db::ExchangeUnits {});
        expect(file.vibration_translation.size() == 15_u);
        expect(file.chemistry_vibration.size() == 3_u);
        const db::VibrationTranslationRecord& first = file.vibration_translation[0];
        expect(first.species == "N2" && first.partner == "N" && first.model == "landau_teller");
        expect(first.relaxation_time.model == "millikan_white");
        expect(first.relaxation_time.a == 180.0_d && first.relaxation_time.b == 0.0262_d);
        expect(first.relaxation_time.high_temperature_correction.has_value());
        expect(first.relaxation_time.high_temperature_correction->cross_section == 1.0e-21);
        const db::ChemistryVibrationRecord& coupling = file.chemistry_vibration[1];
        expect(coupling.reaction == 1_i && coupling.species == "O2");
        expect(coupling.model == "impartial_dissociation");
        expect(coupling.vibrational_temperature == 2273.53_d);
    };

    "an energy-exchange file reads, writes and reads again to the same data and text"_test = [] {
        const db::EnergyExchangeFile file = db::read_energy_exchange_file(exchange_path());
        const std::string text = db::write_energy_exchange(file);
        const db::EnergyExchangeFile again = db::parse_energy_exchange(text);
        expect(again == file);
        expect(db::write_energy_exchange(again) == text);
        const auto path = std::filesystem::temp_directory_path() / "gaslib_test_exchange.json";
        db::write_energy_exchange_file(path, file);
        expect(db::read_energy_exchange_file(path) == file);
        std::filesystem::remove(path);
    };

    "unknown and missing keys, another format or other units of an energy-exchange file are errors"_test =
      [] {
          expect(throws_database_error([] {
              db::parse_energy_exchange(exchange_with("\"format\"", R"("extra": 1, "format")"));
          }));
          expect(throws_database_error(
            [] { db::parse_energy_exchange(exchange_with("\"partner\"", "\"collider\"")); }));
          expect(throws_database_error(
            [] { db::parse_energy_exchange(exchange_with(R"("model": "landau_teller",)", "")); }));
          expect(throws_database_error(
            [] { db::parse_energy_exchange(exchange_with(R"("reaction": 0,)", "")); }));
          expect(throws_database_error([] {
              db::parse_energy_exchange(
                exchange_with("gaslib-energy-exchange-1", "gaslib-energy-exchange-2"));
          }));
          expect(throws_database_error([] {
              db::parse_energy_exchange(
                exchange_with(R"("energy": "J/mol")", R"("energy": "J/kmol")"));
          }));
          expect(throws_database_error([] {
              db::parse_energy_exchange(exchange_with(R"("area": "m^2")", R"("area": "cm^2")"));
          }));
          expect(throws_database_error(
            [] { db::read_energy_exchange_file(source_path("no_such_file.json")); }));
          expect(throws_database_error([] { db::parse_energy_exchange(read_text(park_path())); }));
          // The correction is optional.
          db::EnergyExchangeFile file = db::read_energy_exchange_file(exchange_path());
          file.vibration_translation[2].relaxation_time.high_temperature_correction.reset();
          expect(db::parse_energy_exchange(db::write_energy_exchange(file)) == file);
      };

    "the shipped energy-exchange schema is the schema of the current file structs"_test = [] {
        expect(db::energy_exchange_schema() ==
               read_text(source_path("database/energy_exchange.schema.json")));
    };

    "model inputs from the energy-exchange file equal the compiled data"_test = [] {
        const db::EnergyExchangeFile file = db::read_energy_exchange_file(exchange_path());
        const auto input = db::energy_exchange_input<15, 3>(file);
        const auto& compiled = gaslib::air5_park_exchange;
        for (std::size_t k = 0; k < 15; ++k) {
            const auto& a = input.vibration_translation[k];
            const auto& b = compiled.vibration_translation[k];
            expect(std::string_view(a.species) == b.species &&
                   std::string_view(a.partner) == b.partner)
              << k;
            expect(a.millikan_white.a == b.millikan_white.a &&
                   a.millikan_white.b == b.millikan_white.b)
              << k;
            expect(same_correction(a.correction, b.correction)) << k;
        }
        for (std::size_t k = 0; k < 3; ++k) {
            const auto& a = input.chemistry_vibration[k];
            const auto& b = compiled.chemistry_vibration[k];
            expect(a.reaction == b.reaction && std::string_view(a.species) == b.species) << k;
            expect(a.dissociation_energy == b.dissociation_energy &&
                   a.vibrational_temperature == b.vibrational_temperature)
              << k;
        }
        // The energy exchange built from the file builds like the compiled one.
        const gaslib::UnitScales<double, gaslib::si> scales {};
        const auto built = gaslib::make_air5<double>(scales);
        const auto mechanism = gaslib::make_air5_park(built.species, scales);
        (void)gaslib::make_energy_exchange(input, air5_names, built.species, mechanism, scales);
        (void)gaslib::make_air5_park_exchange(built.species, mechanism, scales);
    };

    "energy-exchange conversion errors name the entry and the problem"_test = [] {
        const db::EnergyExchangeFile file = db::read_energy_exchange_file(exchange_path());
        expect(throws_parameter_error([&] { (void)db::energy_exchange_input<14, 3>(file); }));
        expect(throws_parameter_error([&] { (void)db::energy_exchange_input<15, 2>(file); }));
        const auto message = [](const db::EnergyExchangeFile& f) -> std::string {
            try {
                (void)db::energy_exchange_input<15, 3>(f);
            } catch (const gaslib::ParameterError& error) {
                return error.what();
            }
            return "no exception";
        };
        db::EnergyExchangeFile changed = file;
        changed.vibration_translation[1].model = "ssh";
        expect(message(changed).starts_with(
          R"(vibration_translation 1 (N2-O): model must be "landau_teller")"))
          << message(changed);
        changed = file;
        changed.vibration_translation[4].relaxation_time.model = "schwartz";
        expect(message(changed).starts_with("vibration_translation 4 (N2-NO)")) << message(changed);
        changed = file;
        changed.vibration_translation[5].relaxation_time.high_temperature_correction->model =
          "park_2000";
        expect(message(changed).starts_with("vibration_translation 5 (O2-N)")) << message(changed);
        changed = file;
        changed.chemistry_vibration[2].model = "marrone_treanor";
        expect(message(changed).starts_with(
          R"(chemistry_vibration 2 (NO): model must be "impartial_dissociation")"))
          << message(changed);
    };

    "the shipped air file holds transport data for every species"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        for (const auto name : air5_names) {
            const auto& transport = file.species.at(std::string(name)).transport;
            expect(transport.has_value() && transport->viscosity.has_value() &&
                   transport->thermal_conductivity.has_value() &&
                   transport->lennard_jones.has_value())
              << name;
        }
        const db::TransportRecord& n2 = *file.species.at("N2").transport;
        expect(n2.viscosity->model == "cea");
        expect(n2.viscosity->temperature_ranges.size() == 4_u);
        expect(n2.lennard_jones->sigma == 3.621e-10);
        expect(n2.lennard_jones->epsilon == 97.53_d);
        expect(file.species.at("N").transport->viscosity->temperature_ranges.front() == 1000.0_d);
    };

    "unknown keys, a missing length unit and other lengths are errors of a species file"_test = [] {
        expect(throws_database_error(
          [] { db::parse_species(air5_with("\"lennard_jones\"", "\"lj\"")); }));
        expect(throws_database_error(
          [] { db::parse_species(air5_with(R"("sigma")", R"("diameter": 1.0, "sigma")")); }));
        expect(throws_database_error(
          [] { db::parse_species(air5_with(R"("length": "m")", R"("length": "cm")")); }));
        expect(throws_database_error([] {
            db::parse_species(air5_with(R"(,
      "length": "m")",
                                        ""));
        }));
        // Transport data are optional, as is each of their records.
        db::SpeciesFile file = db::read_species_file(air5_path());
        file.species.at("O").transport.reset();
        file.species.at("N").transport->lennard_jones.reset();
        expect(db::parse_species(db::write_species(file)) == file);
    };

    "transport inputs from the air file equal the compiled data"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        const auto input = db::species_transport<3>(file, air5_names);
        for (std::size_t i = 0; i < 5; ++i) {
            const auto& a = input[i];
            const auto& b = gaslib::air5_species_transport[i];
            expect(std::string_view(a.name) == b.name) << i;
            expect(a.viscosity.segments == b.viscosity.segments &&
                   a.viscosity.breaks == b.viscosity.breaks &&
                   a.viscosity.coefficients == b.viscosity.coefficients)
              << i;
            expect(a.conductivity.segments == b.conductivity.segments &&
                   a.conductivity.breaks == b.conductivity.breaks &&
                   a.conductivity.coefficients == b.conductivity.coefficients)
              << i;
            expect(a.lennard_jones.sigma == b.lennard_jones.sigma &&
                   a.lennard_jones.epsilon == b.lennard_jones.epsilon)
              << i;
        }
    };

    "transport conversion errors name the species and the problem"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        const auto message = [](const db::SpeciesFile& f) -> std::string {
            try {
                (void)db::species_transport<3>(f, air5_names);
            } catch (const gaslib::ParameterError& error) {
                return error.what();
            }
            return "no exception";
        };
        // N2 has three segments.
        expect(throws_parameter_error([&] { (void)db::species_transport<2>(file, air5_names); }));
        db::SpeciesFile changed = file;
        changed.species.at("O2").transport.reset();
        expect(message(changed).starts_with("species O2: needs viscosity")) << message(changed);
        changed = file;
        changed.species.at("NO").transport->thermal_conductivity->model = "chemkin";
        expect(
          message(changed).starts_with(R"(species NO: thermal_conductivity model must be "cea")"))
          << message(changed);
        changed = file;
        changed.species.at("N").transport->viscosity->temperature_ranges.pop_back();
        expect(message(changed).starts_with("species N: viscosity needs")) << message(changed);
        changed = file;
        changed.species.erase("O");
        expect(throws_database_error([&] { (void)db::species_transport<3>(changed, air5_names); }));
    };

    "the shipped collision-integral file reads with its 15 pairs"_test = [] {
        const db::CollisionIntegralsFile file = db::read_collision_integrals_file(collision_path());
        expect(file.format == db::collision_integrals_format);
        expect(file.units == db::CollisionUnits {});
        expect(file.pairs.size() == 15_u);
        expect(file.pairs[1].species == std::array<std::string, 2> {"O2", "N2"});
        expect(file.pairs[1].pi_omega_11.C == 0.5729_d);
    };

    "a collision-integral file reads, writes and reads again to the same data and text"_test = [] {
        const db::CollisionIntegralsFile file = db::read_collision_integrals_file(collision_path());
        const std::string text = db::write_collision_integrals(file);
        const db::CollisionIntegralsFile again = db::parse_collision_integrals(text);
        expect(again == file);
        expect(db::write_collision_integrals(again) == text);
        const auto path = std::filesystem::temp_directory_path() / "gaslib_test_collisions.json";
        db::write_collision_integrals_file(path, file);
        expect(db::read_collision_integrals_file(path) == file);
        std::filesystem::remove(path);
    };

    "unknown and missing keys, another format or other units of a collision-integral file are errors"_test =
      [] {
          expect(throws_database_error([] {
              db::parse_collision_integrals(
                collision_with("\"format\"", R"("extra": 1, "format")"));
          }));
          expect(throws_database_error([] {
              db::parse_collision_integrals(collision_with("\"pi_omega_22\"", "\"omega_22\""));
          }));
          expect(throws_database_error(
            [] { db::parse_collision_integrals(collision_with(R"("A": 0.0,)", "")); }));
          expect(throws_database_error([] {
              db::parse_collision_integrals(
                collision_with("gaslib-collision-integrals-1", "gaslib-collision-integrals-2"));
          }));
          expect(throws_database_error([] {
              db::parse_collision_integrals(collision_with(R"("area": "m^2")", R"("area": "A^2")"));
          }));
          expect(throws_database_error(
            [] { db::read_collision_integrals_file(source_path("no_such_file.json")); }));
          expect(
            throws_database_error([] { db::parse_collision_integrals(read_text(air5_path())); }));
      };

    "the shipped collision-integral schema is the schema of the current file structs"_test = [] {
        expect(db::collision_integrals_schema() ==
               read_text(source_path("database/collision_integrals.schema.json")));
    };

    "collision-pair inputs from the file equal the compiled data"_test = [] {
        const db::CollisionIntegralsFile file = db::read_collision_integrals_file(collision_path());
        const auto input = db::collision_integral_pairs<15>(file);
        const auto same = [](const gaslib::CollisionFitInput& a,
                             const gaslib::CollisionFitInput& b) {
            return a.A == b.A && a.B == b.B && a.C == b.C && a.D == b.D;
        };
        for (std::size_t k = 0; k < 15; ++k) {
            const auto& a = input[k];
            const auto& b = gaslib::air5_collision_integrals[k];
            expect(std::string_view(a.first) == b.first && std::string_view(a.second) == b.second)
              << k;
            expect(same(a.pi_omega_11, b.pi_omega_11) && same(a.pi_omega_22, b.pi_omega_22)) << k;
        }
        expect(throws_parameter_error([&] { (void)db::collision_integral_pairs<14>(file); }));
    };

    "the shipped air files hold electronic levels"_test = [] {
        const db::SpeciesFile air5 = db::read_species_file(air5_path());
        const auto& n2 = air5.species.at("N2").electronic_levels;
        expect(n2.has_value() && n2->degeneracies.size() == 7_u && n2->temperatures.size() == 7_u);
        expect(n2->temperatures[0] == 0.0_d);
        // 50203.66 cm^-1 times h c / k_B = 1.438776877... cm K.
        expect(std::abs(n2->temperatures[1] - (50203.66 * 1.4387768775039338)) <= 1.0e-9 * 7.2e4);
        const db::SpeciesFile air11 = db::read_species_file(air11_path());
        expect(air11.species.size() == 11_u);
        expect(air11.species.at("e-").charge == -1.0_d);
        expect(!air11.species.at("e-").electronic_levels.has_value());
        expect(air11.species.at("NO+").electronic_levels->notes.has_value());
        expect(air11.species.at("N2").thermo == air5.species.at("N2").thermo);
    };

    "electronic-level inputs from the air files equal the compiled data"_test = [] {
        const auto same = [](const auto& a, const auto& b) {
            return std::string_view(a.name) == b.name && a.levels == b.levels &&
                   a.degeneracies == b.degeneracies && a.temperatures == b.temperatures;
        };
        const db::SpeciesFile air5 = db::read_species_file(air5_path());
        const auto input5 = db::electronic_levels<8>(air5, air5_names);
        for (std::size_t i = 0; i < 5; ++i) {
            expect(same(input5[i], gaslib::air5_electronic_levels[i])) << i;
        }
        const db::SpeciesFile air11 = db::read_species_file(air11_path());
        const auto input11 = db::electronic_levels<8>(air11, air11_names);
        for (std::size_t i = 0; i < 11; ++i) {
            expect(same(input11[i], gaslib::air11_electronic_levels[i])) << i;
        }
        const auto species11 = db::nasa9_species<4>(air11, air11_names);
        for (std::size_t i = 0; i < 11; ++i) {
            const auto& a = species11[i];
            const auto& b = gaslib::air11_species[i];
            expect(std::string_view(a.name) == b.name && a.molar_mass == b.molar_mass &&
                   a.charge == b.charge && a.fit.segments == b.fit.segments &&
                   a.fit.breaks == b.fit.breaks && a.fit.blends == b.fit.blends &&
                   a.fit.a == b.fit.a)
              << i;
        }
        // N2 has 7 levels.
        expect(throws_parameter_error([&] { (void)db::electronic_levels<6>(air5, air5_names); }));
        db::SpeciesFile changed = air5;
        changed.species.at("O").electronic_levels->temperatures.pop_back();
        expect(
          throws_parameter_error([&] { (void)db::electronic_levels<8>(changed, air5_names); }));
        changed = air5;
        changed.species.erase("NO");
        expect(throws_database_error([&] { (void)db::electronic_levels<8>(changed, air5_names); }));
    };

    "the eleven-species Kim-Jo mechanism reads with electron temperatures"_test = [] {
        const db::ReactionsFile file = db::read_reactions_file(air11_kim_jo_path());
        expect(file.reactions.size() == 40_u);
        const db::ReactionRecord& n2e = file.reactions[21];
        expect(n2e.equation == "N2 + e- <=> N + N + e-");
        expect(n2e.forward.other_temperature == std::optional<std::string> {"electron"});
        expect(n2e.backward.other_temperature == std::optional<std::string> {"electron"});
        expect(!file.reactions[0].forward.other_temperature.has_value());
        const auto input = db::mechanism_input<40, 11>(file);
        for (std::size_t k = 0; k < 40; ++k) {
            const auto& a = input.reactions[k];
            const auto& b = gaslib::air11_kim_jo.reactions[k];
            expect(std::string_view(a.equation) == b.equation) << k;
            expect(same_rate(a.forward, b.forward)) << k;
            expect(same_other(a.forward, b.forward) && same_other(a.backward, b.backward)) << k;
            expect(same_coefficients(a.efficiencies, b.efficiencies)) << k;
        }
        expect(std::string_view(input.reactions[38].backward.other_temperature) == "electron");
    };

    "the eleven-species exchange file reads and maps to the compiled electron exchanges"_test = [] {
        const db::EnergyExchangeFile file = db::read_energy_exchange_file(air11_exchange_path());
        expect(file.vibration_translation.size() == 30_u &&
               file.chemistry_vibration.size() == 21_u);
        expect(file.electron_translation.size() == 10_u && file.vibration_electron.size() == 1_u &&
               file.electron_impact_ionisation.size() == 2_u);
        expect(db::parse_energy_exchange(db::write_energy_exchange(file)) == file);
        const db::SpeciesFile species = db::read_species_file(air11_path());
        expect(species.species.at("N").ionisation_energy.has_value());
        const auto electrons = db::electron_exchange_input<10, 1, 2>(file, species);
        const auto& compiled = gaslib::air11_kim_jo_exchange_electrons;
        for (std::size_t k = 0; k < 10; ++k) {
            const auto& a = electrons.electron_translation[k];
            const auto& b = compiled.electron_translation[k];
            expect(std::string_view(a.partner) == b.partner && a.cross_section == b.cross_section &&
                   a.a == b.a && a.b == b.b && a.c == b.c &&
                   a.maximum_temperature == b.maximum_temperature)
              << k;
        }
        const auto& ve = electrons.vibration_electron[0];
        expect(std::string_view(ve.species) == "N2" && ve.relaxation_time.a_low == 7.5_d &&
               ve.relaxation_time.blend_upper == 7050.0_d);
        for (std::size_t k = 0; k < 2; ++k) {
            const auto& a = electrons.electron_impact_ionisation[k];
            const auto& b = compiled.electron_impact_ionisation[k];
            expect(a.reaction == b.reaction && std::string_view(a.species) == b.species &&
                   a.ionisation_energy == b.ionisation_energy)
              << k;
        }
        const auto input = db::energy_exchange_input<30, 21>(file);
        for (std::size_t k = 0; k < 30; ++k) {
            const auto& a = input.vibration_translation[k];
            const auto& b = gaslib::air11_kim_jo_exchange.vibration_translation[k];
            expect(same_correction(a.correction, b.correction) &&
                   a.millikan_white.a == b.millikan_white.a)
              << k;
        }
        // Errors: counts, an unknown cross-section, a missing ionisation energy.
        expect(throws_parameter_error(
          [&] { (void)db::electron_exchange_input<9, 1, 2>(file, species); }));
        db::EnergyExchangeFile changed = file;
        changed.electron_translation[3].cross_section.model = "maxwell";
        expect(throws_parameter_error(
          [&] { (void)db::electron_exchange_input<10, 1, 2>(changed, species); }));
        changed = file;
        changed.vibration_electron[0].relaxation_time.model = "millikan_white";
        expect(throws_parameter_error(
          [&] { (void)db::electron_exchange_input<10, 1, 2>(changed, species); }));
        db::SpeciesFile no_energy = species;
        no_energy.species.at("O").ionisation_energy.reset();
        expect(throws_parameter_error(
          [&] { (void)db::electron_exchange_input<10, 1, 2>(file, no_energy); }));
    };
}
