// Tests of the species database library: reading and writing species files,
// errors, the generated schema, conversion to model inputs, and cp/R of the
// shipped five-species air file against reference values.

#include <array>
#include <boost/ut.hpp>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <gaslib/build.hpp>
#include <gaslib/data/air5.hpp>
#include <gaslib/database.hpp>
#include <gaslib/nasa9.hpp>
#include <gaslib/units.hpp>
#include <ios>
#include <iterator>
#include <string>
#include <string_view>

namespace
{

namespace db = gaslib::database;

std::filesystem::path source_path(std::string_view relative)
{ return std::filesystem::path(GASLIB_SOURCE_DIR) / relative; }

std::filesystem::path air5_path() { return source_path("database/data/air5.json"); }

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
}
