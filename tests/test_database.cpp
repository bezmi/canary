// Tests of the species database library: reading and writing species,
// reactions and energy-exchange files, errors, the generated schemas,
// conversion to model inputs, cp/R of the shipped five-species air file
// against reference values, and the compiled Kim-Jo mechanism, Park
// mechanism, Park energy exchange, species transport data, collision
// integrals and electronic levels against their files.

#include "../database/file_io.hpp"

#include <algorithm>
#include <array>
#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/data/air11_kim_jo.hpp>
#include <canary/data/air11_kim_jo_exchange.hpp>
#include <canary/data/air11_species.hpp>
#include <canary/data/air5.hpp>
#include <canary/data/air5_electronic.hpp>
#include <canary/data/air5_kim_jo.hpp>
#include <canary/data/air5_park.hpp>
#include <canary/data/air5_park_exchange.hpp>
#include <canary/data/air5_transport.hpp>
#include <canary/database.hpp>
#include <canary/energy_exchange.hpp>
#include <canary/nasa9.hpp>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/transport/mixture.hpp>
#include <canary/units.hpp>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <streambuf>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>

namespace
{

namespace db = canary::database;

// Each fixture owns a directory created atomically, including across processes.
class TemporaryDirectory
{
  public:
    TemporaryDirectory()
    {
        const auto parent = std::filesystem::temp_directory_path();
        std::random_device random;
        for (int attempt = 0; attempt < 128; ++attempt) {
            auto candidate =
              parent / ("canary-test-" + std::to_string(random()) + "-" + std::to_string(random()));
            std::error_code error;
            if (std::filesystem::create_directory(candidate, error)) {
                path_ = std::move(candidate);
                return;
            }
            if (error && error != std::errc::file_exists) {
                throw std::filesystem::filesystem_error("cannot create test directory",
                                                        candidate,
                                                        error);
            }
        }
        throw std::runtime_error("cannot create a unique test directory");
    }

    TemporaryDirectory(const TemporaryDirectory&) = delete;
    TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

    ~TemporaryDirectory()
    {
        std::error_code error;
        std::filesystem::remove_all(path_, error);
    }

    const std::filesystem::path& path() const { return path_; }

  private:
    std::filesystem::path path_;
};

std::filesystem::path source_path(std::string_view relative)
{ return std::filesystem::path(CANARY_SOURCE_DIR) / relative; }

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
    boost::ut::expect(boost::ut::fatal(at != std::string::npos))
      << "fixture text to replace is missing: " << from;
    text.replace(at, from.size(), to);
    return text;
}

// air5_kim_jo.json with one textual replacement.
std::string kim_jo_with(std::string_view from, std::string_view to)
{
    std::string text = read_text(kim_jo_path());
    const auto at = text.find(from);
    boost::ut::expect(boost::ut::fatal(at != std::string::npos))
      << "fixture text to replace is missing: " << from;
    text.replace(at, from.size(), to);
    return text;
}

// The same species coefficients in any order.
template <std::size_t K>
bool same_coefficients(const canary::SpeciesCoefficients<K>& a,
                       const canary::SpeciesCoefficients<K>& b)
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

// Compare physical identities and species text independently of ownership.
bool same_temperature(const std::optional<canary::PhysicalTemperatureInput>& a,
                      const std::optional<canary::PhysicalTemperatureInput>& b)
{
    if (!a || !b) { return a.has_value() == b.has_value(); }
    if (a->index() != b->index()) { return false; }
    if (const auto* v = std::get_if<canary::VibrationTemperatureInput>(&*a)) {
        const auto* w = std::get_if<canary::VibrationTemperatureInput>(&*b);
        if (!v->species || !w->species) { return v->species == w->species; }
        return std::string_view(v->species) == w->species;
    }
    if (const auto* v = std::get_if<canary::BoundElectronicTemperatureInput>(&*a)) {
        const auto* w = std::get_if<canary::BoundElectronicTemperatureInput>(&*b);
        if (!v->species || !w->species) { return v->species == w->species; }
        return std::string_view(v->species) == w->species;
    }
    return true;
}

bool same_correction(const canary::RelaxationCorrectionInput& a,
                     const canary::RelaxationCorrectionInput& b)
{
    return a.model == b.model && a.cross_section == b.cross_section &&
           a.reference_temperature == b.reference_temperature && a.exponent == b.exponent;
}

bool same_rate(const canary::ArrheniusInput& a, const canary::ArrheniusInput& b)
{
    return a.A == b.A && a.n == b.n && a.activation_temperature == b.activation_temperature &&
           a.translational_exponent == b.translational_exponent &&
           same_temperature(a.other_temperature, b.other_temperature);
}

// air5_collision_integrals.json with one textual replacement.
std::string collision_with(std::string_view from, std::string_view to)
{
    std::string text = read_text(collision_path());
    const auto at = text.find(from);
    boost::ut::expect(boost::ut::fatal(at != std::string::npos))
      << "fixture text to replace is missing: " << from;
    text.replace(at, from.size(), to);
    return text;
}

// air5_park_exchange.json with one textual replacement.
std::string exchange_with(std::string_view from, std::string_view to)
{
    std::string text = read_text(exchange_path());
    const auto at = text.find(from);
    boost::ut::expect(boost::ut::fatal(at != std::string::npos))
      << "fixture text to replace is missing: " << from;
    text.replace(at, from.size(), to);
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
    } catch (const canary::ParameterError&) {
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

double cp_over_r(const canary::Nasa9Curve<double, 4>& curve, double t)
{ return canary::nasa9_cp_h(curve, t, std::log(t)).cp_over_r; }

}  // namespace

int main()
{
    using namespace boost::ut;

    "converted owners survive source destruction and value operations"_test = [] {
        const canary::UnitScales<double, canary::si> scales {};
        const auto gas5 = canary::make_air5<double>(scales);
        const auto mechanism5 = canary::make_air5_park(scales);
        const auto gas11 = canary::make_air11_three_temperature(scales);
        const auto mechanism11 =
          canary::make_air11_kim_jo(scales,
                                    canary::Pressure<double, canary::si> {1e5},
                                    gas11.model);
        const auto exercise = [](const auto& factory, const auto& consume) {
            using Owner = decltype(factory());
            static_assert(std::is_copy_constructible_v<Owner> && std::is_copy_assignable_v<Owner>);
            static_assert(std::is_move_constructible_v<Owner> && std::is_move_assignable_v<Owner>);
            consume(factory());
            const auto copied = [&] {
                const auto source = factory();
                return Owner(source);
            }();
            consume(copied);
            auto copy_assigned = factory();
            {
                const auto source = factory();
                copy_assigned = source;
            }
            consume(copy_assigned);
            const auto moved = [&] {
                auto source = factory();
                return Owner(std::move(source));
            }();
            consume(moved);
            auto move_assigned = factory();
            {
                auto source = factory();
                move_assigned = std::move(source);
            }
            consume(move_assigned);
        };
        exercise(
          [] {
              auto file = db::read_species_file(air5_path());
              auto owner = db::nasa9_species<4>(file, air5_names);
              file.species.clear();
              return owner;
          },
          [&](const auto& owner) {
              const auto built =
                canary::make_thermally_perfect_mixture<double, canary::si>(owner, scales);
              expect(built.species.mass == gas5.species.mass);
              const auto two = canary::make_two_temperature_mixture(owner,
                                                                    canary::air5_rotational_degrees,
                                                                    scales);
              expect(two.species.mass == gas5.species.mass);
          });
        exercise(
          [] { return db::mechanism_input<40, 11>(db::read_reactions_file(air11_kim_jo_path())); },
          [&](const auto& owner) {
              const auto built =
                canary::make_reaction_mechanism(owner,
                                                canary::reaction_species(canary::air11_species),
                                                scales,
                                                canary::Pressure<double, canary::si> {1e5},
                                                gas11.model);
              for (const auto i : canary::species_indices(0, 11)) {
                  expect(built.reactions()[canary::ReactionIndex {21}].reactants()[i] ==
                         mechanism11.reactions()[canary::ReactionIndex {21}].reactants()[i]);
              }
              expect(built.reactions()[canary::ReactionIndex {38}].backward.A ==
                     mechanism11.reactions()[canary::ReactionIndex {38}].backward.A);
          });
        exercise(
          [] {
              return db::energy_exchange_input<15, 3>(
                db::read_energy_exchange_file(exchange_path()));
          },
          [&](const auto& owner) {
              const auto built =
                canary::make_energy_exchange(owner, air5_names, gas5.species, mechanism5, scales);
              const auto reference =
                canary::make_air5_park_exchange(gas5.species, mechanism5, scales);
              expect(built.chemistry_vibration()[0].species ==
                     reference.chemistry_vibration()[0].species);
              for (const auto i : canary::species_indices(0, 5)) {
                  expect(built.mass()[i] == reference.mass()[i]) << i.value();
              }
          });
        exercise(
          [] { return db::species_transport<3>(db::read_species_file(air5_path()), air5_names); },
          [&](const auto& owner) {
              const auto built =
                canary::make_wilke_transport(owner, air5_names, gas5.species, scales);
              expect(built.mass == gas5.species.mass);
          });
        exercise(
          [] {
              return db::collision_integral_pairs<15>(
                db::read_collision_integrals_file(collision_path()));
          },
          [&](const auto& owner) {
              const auto built = canary::make_gupta_yos_transport(owner,
                                                                  air5_names,
                                                                  gas5.species,
                                                                  canary::air5_rotational_degrees,
                                                                  scales);
              expect(built.mass == gas5.species.mass);
          });
        exercise(
          [] { return db::electronic_levels<8>(db::read_species_file(air11_path()), air11_names); },
          [&](const auto& owner) {
              const auto species =
                db::nasa9_species<4>(db::read_species_file(air11_path()), air11_names);
              const auto built =
                canary::make_multi_temperature_mixture(species,
                                                       canary::air11_rotational_degrees,
                                                       owner,
                                                       canary::air11_three_temperature_modes,
                                                       scales);
              const auto mixed =
                canary::make_multi_temperature_mixture(canary::air11_species,
                                                       canary::air11_rotational_degrees,
                                                       owner,
                                                       canary::air11_three_temperature_modes,
                                                       scales);
              expect(built.species.mass == gas11.species.mass &&
                     mixed.species.mass == built.species.mass);
          });
        exercise(
          [] {
              return db::electron_exchange_input<10, 1, 2>(
                db::read_energy_exchange_file(air11_exchange_path()),
                db::read_species_file(air11_path()));
          },
          [&](const auto& owner) {
              const auto heavy = db::energy_exchange_input<30, 21>(
                db::read_energy_exchange_file(air11_exchange_path()));
              const auto built = canary::make_multi_temperature_exchange(heavy,
                                                                         owner,
                                                                         air11_names,
                                                                         gas11.species,
                                                                         mechanism11,
                                                                         gas11.model,
                                                                         scales);
              const auto mixed =
                canary::make_multi_temperature_exchange(canary::air11_kim_jo_exchange,
                                                        owner,
                                                        air11_names,
                                                        gas11.species,
                                                        mechanism11,
                                                        gas11.model,
                                                        scales);
              for (const auto i : canary::species_indices(0, 11)) {
                  expect(built.mass()[i] == gas11.species.mass[i.value()] &&
                         mixed.mass()[i] == built.mass()[i])
                    << i.value();
              }
          });
    };

    "owning species retain names and optional element metadata independently"_test = [] {
        using Input = std::remove_cv_t<decltype(canary::air5_species)>;
        using Owner = canary::OwnedInput<Input>;
        static_assert(!std::is_convertible_v<Owner, Input>);
        static_assert(!std::is_constructible_v<canary::detail::OwnedSpeciesName,
                                               canary::detail::OwnedReactionEquation>);
        static_assert(!std::is_constructible_v<canary::detail::OwnedSpeciesName,
                                               canary::detail::OwnedElementName>);
        const auto make_owner = [] {
            std::string name = "nitrogen_species_with_a_long_owned_name";
            std::string element = "element_identifier_with_a_long_owned_name";
            std::string short_element = "N";
            std::array<canary::ElementCount, 2> elements {{{element, 2}, {short_element, 1}}};
            auto input = canary::air5_species;
            input[0].name = name.c_str();
            input[0].metadata.elements = std::span<const canary::ElementCount> {elements};
            input[1].metadata.elements = std::span<const canary::ElementCount> {};
            input[2].metadata.elements = std::nullopt;
            Owner owner {input};
            name.assign("changed");
            element.assign("changed");
            short_element.assign("changed");
            elements[0].count = 99;
            input[0].molar_mass = 99;
            return owner;
        };
        const auto copied = [&] {
            const auto source = make_owner();
            return Owner(source);
        }();
        const auto moved = [&] {
            auto source = make_owner();
            return Owner(std::move(source));
        }();
        auto copy_assigned = make_owner();
        {
            const auto source = make_owner();
            copy_assigned = source;
        }
        auto move_assigned = make_owner();
        {
            auto source = make_owner();
            move_assigned = std::move(source);
        }
        for (const auto* owner : {&copied,
                                  &moved,
                                  static_cast<const Owner*>(&copy_assigned),
                                  static_cast<const Owner*>(&move_assigned)}) {
            canary::detail::with_input(*owner, [&](const auto& input) {
                expect(std::string_view(input[0].name) ==
                       "nitrogen_species_with_a_long_owned_name");
                expect(input[0].molar_mass == canary::air5_species[0].molar_mass);
                const auto elements = *input[0].metadata.elements;
                expect(elements.size() == 2_u);
                expect(elements[0].element == "element_identifier_with_a_long_owned_name" &&
                       elements[0].count == 2_u);
                expect(elements[1].element == "N" && elements[1].count == 1_u);
                expect(input[1].metadata.elements.has_value() &&
                       input[1].metadata.elements->empty());
                expect(!input[2].metadata.elements.has_value());
            });
            const auto built =
              canary::make_thermally_perfect_mixture(*owner,
                                                     canary::UnitScales<double, canary::si> {});
            expect(built.species.n == 5_u);
        }
        auto input = canary::air5_species;
        input[0].name = nullptr;
        input[1].name = "";
        const Owner owner {input};
        canary::detail::with_input(owner, [&](const auto& borrowed) {
            expect(borrowed[0].name == nullptr);
            expect(borrowed[1].name != nullptr && std::string_view(borrowed[1].name).empty());
        });
        expect(throws_parameter_error([&] {
            (void)canary::make_thermally_perfect_mixture(owner,
                                                         canary::UnitScales<double, canary::si> {});
        }));
    };

    "temporary directories are distinct and removed after scope exit"_test = [] {
        std::filesystem::path first_path;
        std::filesystem::path second_path;
        {
            const TemporaryDirectory first;
            const TemporaryDirectory second;
            first_path = first.path();
            second_path = second.path();
            expect(first_path != second_path);
            expect(std::filesystem::is_directory(first_path));
            expect(std::filesystem::is_directory(second_path));
        }
        expect(!std::filesystem::exists(first_path));
        expect(!std::filesystem::exists(second_path));
    };

    "temporary directories remove files during exception unwinding"_test = [] {
        std::filesystem::path path;
        expect(throws_database_error([&] {
            const TemporaryDirectory directory;
            path = directory.path();
            db::write_species_file(path / "species.json", db::read_species_file(air5_path()));
            throw db::DatabaseError("test exception");
        }));
        expect(!path.empty());
        expect(!std::filesystem::exists(path));
    };

    "all database writers reject invalid destinations without changing existing files"_test = [] {
        const TemporaryDirectory directory;
        const auto occupied = directory.path() / "occupied";
        expect(fatal(std::filesystem::create_directory(occupied)));
        const auto marker = occupied / "marker.txt";
        {
            std::ofstream out(marker);
            out << "retained content";
            out.close();
            expect(fatal(bool(out)));
        }
        const auto check = [&](const auto& file, auto write) {
            expect(throws_database_error([&] { write(occupied, file); }));
            expect(throws_database_error(
              [&] { write(directory.path() / "missing" / "output.json", file); }));
            expect(throws_database_error([&] { write(marker / "output.json", file); }));
            expect(read_text(marker) == "retained content");
            expect(std::filesystem::is_directory(occupied));
            expect(!std::filesystem::exists(directory.path() / "missing"));
        };
        check(db::read_species_file(air5_path()), db::write_species_file);
        check(db::read_reactions_file(kim_jo_path()), db::write_reactions_file);
        check(db::read_energy_exchange_file(exchange_path()), db::write_energy_exchange_file);
        check(db::read_collision_integrals_file(collision_path()),
              db::write_collision_integrals_file);
    };

    "database reads require completion and preserve binary contents"_test = [] {
        for (const std::size_t size : {0U, 1U, 16383U, 16384U, 16385U, 32768U}) {
            std::string contents(size, 'x');
            if (size > 3) { contents[3] = '\0'; }
            std::istringstream input(contents);
            expect(db::detail::read_contents(input, "sample.json", "species") == contents);
        }
        class FailedRead : public std::streambuf {
            std::streamsize xsgetn(char* output, std::streamsize count) override
            {
                if (read_once_) { throw std::ios_base::failure("injected read failure"); }
                read_once_ = true;
                std::fill_n(output, count, 'x');
                return count;
            }
            bool read_once_ = false;
        } buffer;
        std::istream input(&buffer);
        try {
            static_cast<void>(db::detail::read_contents(input, "sample.json", "species"));
            expect(false) << "partial contents were accepted after a read failure";
        } catch (const db::DatabaseError& error) {
            expect(std::string_view(error.what()) == "cannot read species file sample.json");
        }
    };

    "database writers replace complete files and preserve permissions"_test = [] {
        const TemporaryDirectory directory;
        const auto destination = directory.path() / "output.json";
        const auto retained = directory.path() / "retained.json";
        const auto permissions = std::filesystem::perms::owner_read |
                                 std::filesystem::perms::owner_write |
                                 std::filesystem::perms::group_read;
        const auto check = [&](const auto& file, auto write, auto serialize) {
            const auto expected = serialize(file);
            write(destination, file);
            expect(read_text(destination) == expected);
            std::filesystem::remove(destination);
            {
                std::ofstream old(destination);
                old << "previous contents";
            }
            std::filesystem::permissions(destination, permissions);
            std::filesystem::create_hard_link(destination, retained);
            write(destination, file);
            expect(read_text(destination) == expected);
            expect(read_text(retained) == "previous contents");
            expect(!std::filesystem::equivalent(destination, retained));
            expect(std::filesystem::status(destination).permissions() == permissions);
            std::filesystem::remove(destination);
            std::filesystem::remove(retained);
            expect(std::filesystem::is_empty(directory.path()));
        };
        check(db::read_species_file(air5_path()), db::write_species_file, db::write_species);
        check(db::read_reactions_file(kim_jo_path()), db::write_reactions_file, db::write_reactions);
        check(db::read_energy_exchange_file(exchange_path()), db::write_energy_exchange_file,
              db::write_energy_exchange);
        check(db::read_collision_integrals_file(collision_path()), db::write_collision_integrals_file,
              db::write_collision_integrals);
    };

    "database writers reject existing and dangling destination symlinks"_test = [] {
        const TemporaryDirectory directory;
        const auto target = directory.path() / "target.json";
        const auto link = directory.path() / "output.json";
        const auto check = [&](const auto& file, auto write) {
            for (const bool exists : {true, false}) {
                if (exists) {
                    std::ofstream output(target);
                    output << "retained contents";
                }
                std::filesystem::create_symlink(target, link);
                expect(throws_database_error([&] { write(link, file); }));
                expect(std::filesystem::is_symlink(link));
                expect(std::filesystem::exists(target) == exists);
                if (exists) { expect(read_text(target) == "retained contents"); }
                std::filesystem::remove(link);
                std::filesystem::remove(target);
                expect(std::filesystem::is_empty(directory.path()));
            }
        };
        check(db::read_species_file(air5_path()), db::write_species_file);
        check(db::read_reactions_file(kim_jo_path()), db::write_reactions_file);
        check(db::read_energy_exchange_file(exchange_path()), db::write_energy_exchange_file);
        check(db::read_collision_integrals_file(collision_path()), db::write_collision_integrals_file);
    };

    "database serialization failure preserves the destination"_test = [] {
        const TemporaryDirectory directory;
        const auto path = directory.path() / "output.json";
        const auto file = db::read_species_file(air5_path());
        db::write_species_file(path, file);
        const auto original = read_text(path);
        auto invalid = file;
        invalid.format = "unsupported";
        expect(throws_database_error([&] { db::write_species_file(path, invalid); }));
        expect(read_text(path) == original);
        std::filesystem::remove(path);
        expect(std::filesystem::is_empty(directory.path()));
    };

    "database parse diagnostics retain file kind and source path"_test = [] {
        const TemporaryDirectory directory;
        const auto path = directory.path() / "malformed.json";
        {
            std::ofstream out(path);
            out << "{";
            out.close();
            expect(fatal(bool(out)));
        }
        const auto error_message = [](auto action) {
            try {
                action();
            } catch (const db::DatabaseError& error) {
                return std::string(error.what());
            }
            return std::string {};
        };
        const auto check = [&](std::string_view kind, auto parse, auto read) {
            const auto diagnostic = error_message([&] { (void)parse("{"); });
            expect(diagnostic.starts_with("invalid " + std::string(kind) + " file: "));
            expect(error_message([&] { (void)read(path); }) == path.string() + ": " + diagnostic);
            const auto missing = directory.path() / "missing.json";
            expect(error_message([&] { (void)read(missing); }) ==
                   "cannot open " + std::string(kind) + " file " + missing.string());
        };
        check("species", db::parse_species, db::read_species_file);
        check("reactions", db::parse_reactions, db::read_reactions_file);
        check("energy-exchange", db::parse_energy_exchange, db::read_energy_exchange_file);
        check("collision-integral",
              db::parse_collision_integrals,
              db::read_collision_integrals_file);
    };

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
        const TemporaryDirectory directory;
        const auto path = directory.path() / "species.json";
        db::write_species_file(path, file);
        expect(db::read_species_file(path) == file);
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
          [] { db::parse_species(air5_with("canary-species-2", "canary-species-1")); }));
        expect(
          throws_database_error([] { db::parse_species(air5_with("\"kg/mol\"", "\"g/mol\"")); }));
        expect(throws_database_error([] { db::parse_species(air5_with("\"K\"", "\"degC\"")); }));
        expect(throws_database_error([] { db::parse_species("{\"format\": "); }));
        expect(
          throws_database_error([] { db::read_species_file(source_path("no_such_file.json")); }));
    };

    "species records require explicit roles and diagnose old formats"_test = [] {
        for (const std::string_view replacement :
             {"", "\"role\": \"unspecified\",", "\"role\": \"unknown\","}) {
            expect(throws_database_error(
              [&] { db::parse_species(air5_with("\"role\": \"heavy\",", replacement)); }));
        }
        expect(throws_database_error(
          [] { db::parse_species(air5_with("\"role\": \"heavy\"", "\"role\": \"electron\"")); }));
        try {
            db::parse_species(R"({"format":"canary-species-1","species":{}})");
            expect(false) << "old species format accepted";
        } catch (const db::DatabaseError& error) {
            expect(std::string_view(error.what())
                     .find("regenerate species files with an explicit heavy/electron role") !=
                   std::string_view::npos)
              << error.what();
        }
        const auto file = db::read_species_file(air11_path());
        expect(file.species.at("e-").role == canary::SpeciesRole::electron);
        expect(file.species.at("N+").role == canary::SpeciesRole::heavy);
        const auto input_owner = db::nasa9_species<4>(file, air11_names);
        canary::detail::with_input(input_owner, [&](const auto& input) {
            for (std::size_t i = 0; i < input.size(); ++i) {
                expect(input[i].metadata.role == canary::air11_species[i].metadata.role);
            }
        });
    };

    "element objects preserve missing empty and populated composition"_test = [] {
        static_assert(!std::is_convertible_v<std::string, db::ElementName>);
        auto file = db::read_species_file(air11_path());
        expect(file.species.at("N2").elements->entries ==
               db::ElementComposition {{{db::ElementName {"N"}, 2}}}.entries);
        expect(file.species.at("e-").elements.has_value());
        expect(file.species.at("e-").elements->entries.empty());
        file.species.at("N").elements.reset();
        file.species.at("O").elements = db::ElementComposition {};
        const auto text = db::write_species(file);
        expect(db::parse_species(text) == file);
        const auto owner = [&] {
            auto source = db::parse_species(text);
            const auto converted = db::nasa9_species<4>(source, air11_names);
            source.species.clear();
            return converted;
        }();
        canary::detail::with_input(owner, [&](const auto& input) {
            expect((*input[0].metadata.elements)[0].element == "N");
            expect((*input[0].metadata.elements)[0].count == 2_u);
            expect(!input[2].metadata.elements.has_value());
            expect(input[3].metadata.elements.has_value() && input[3].metadata.elements->empty());
            expect(input[10].metadata.elements.has_value() && input[10].metadata.elements->empty());
        });
        const auto built =
          canary::make_thermally_perfect_mixture(owner, canary::UnitScales<double, canary::si> {});
        expect(built.species.n == 11_u);
    };

    "element objects reject invalid counts identifiers and repeated decoded keys"_test = [] {
        auto file = db::read_species_file(air5_path());
        for (auto& [name, record] : file.species) { record.elements.reset(); }
        const auto base = db::write_species(file);
        const auto with_elements = [&](std::string_view elements) {
            auto text = base;
            const auto position = text.find("\"thermo\"");
            text.insert(position, "\"elements\": " + std::string(elements) + ",");
            return text;
        };
        for (const std::string_view invalid : {R"({"N":0})",
                                               R"({"N":-1})",
                                               R"({"N":1.5})",
                                               R"({"N":4294967296})",
                                               R"({"N":true})",
                                               R"({"N":"2"})",
                                               R"({"":1})",
                                               R"({"N":2,"N":3})",
                                               R"({"N":2,"\u004e":3})",
                                               R"([{"N":2}])"}) {
            expect(throws_database_error([&] { (void)db::parse_species(with_elements(invalid)); }))
              << invalid;
        }
        expect(nothrow([&] { (void)db::parse_species(with_elements(R"({"N":4294967295})")); }));
        for (const std::string_view duplicate : {R"({"N":2,"N":3})", R"({"N":2,"\u004e":3})"}) {
            try {
                (void)db::parse_species(with_elements(duplicate));
                expect(false) << "duplicate element accepted";
            } catch (const db::DatabaseError& error) {
                expect(std::string_view(error.what()).find("element N: given twice") !=
                       std::string_view::npos);
            }
        }
        auto& entries = file.species.at("N2").elements.emplace().entries;
        entries = {{db::ElementName {"N"}, 2}, {db::ElementName {"N"}, 3}};
        expect(throws_database_error([&] { (void)db::write_species(file); }));
        const auto owner = db::nasa9_species<4>(file, air5_names);
        expect(throws_parameter_error([&] {
            (void)canary::make_thermally_perfect_mixture(owner,
                                                         canary::UnitScales<double, canary::si> {});
        }));
        auto electrons = db::read_species_file(air11_path());
        electrons.species.at("e-").elements = db::ElementComposition {{{db::ElementName {"N"}, 1}}};
        expect(throws_database_error([&] { (void)db::write_species(electrons); }));
    };

    "the shipped schema is the schema of the current file structs"_test = [] {
        expect(db::species_schema() == read_text(source_path("database/species.schema.json")));
    };

    "model inputs from the file equal the compiled five-species air data"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        const auto from_file_owner = db::nasa9_species<4>(file, air5_names);
        canary::detail::with_input(from_file_owner, [&](const auto& from_file) {
            for (std::size_t i = 0; i < 5; ++i) {
                const auto& a = from_file[i];
                const auto& b = canary::air5_species[i];
                expect(std::string_view(a.name) == std::string_view(b.name));
                expect(a.molar_mass == b.molar_mass && a.charge == b.charge) << b.name;
                expect(a.metadata.role == b.metadata.role) << b.name;
                expect(a.metadata.elements.has_value() == b.metadata.elements.has_value());
                if (a.metadata.elements && b.metadata.elements) {
                    expect(a.metadata.elements->size() == b.metadata.elements->size());
                    for (std::size_t e = 0; e < a.metadata.elements->size(); ++e) {
                        expect((*a.metadata.elements)[e].element ==
                               (*b.metadata.elements)[e].element);
                        expect((*a.metadata.elements)[e].count == (*b.metadata.elements)[e].count);
                    }
                }
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
        });
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
        const auto species_owner =
          db::nasa9_species<4>(file, std::array<std::string_view, 2> {"N2", "O2"});
        canary::detail::with_input(species_owner, [&](const auto& species) {
            constexpr std::array temps {300.0, 1000.0, 5000.0, 15000.0};
            constexpr std::array<std::array<double, 4>, 2> gdtk {
              {
                {3.5029350227462328, 3.9324555662029987, 4.5621214841400004, 7.9038660786805508},
                {3.534485021788023, 4.1953819290074996, 5.1713062458949999, 3.9642060171736162},
              },
            };
            for (std::size_t s = 0; s < 2; ++s) {
                const auto curve = canary::make_nasa9_curve(species[s].fit, 1.0);
                for (std::size_t k = 0; k < temps.size(); ++k) {
                    const double cp = cp_over_r(curve, temps[k]);
                    expect(std::abs(cp - gdtk[s][k]) <= 1e-12 * gdtk[s][k])
                      << species[s].name << temps[k] << cp;
                }
            }
        });
    };

    // NIST-JANAF Cp in J/(mol K). The fits differ from JANAF by up to 0.75 %
    // here (O2 at 5000 K); the tolerance matches the curve tests'.
    "Cp of N2 and O2 agrees with JANAF up to its 6000 K limit"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        const auto species_owner =
          db::nasa9_species<4>(file, std::array<std::string_view, 2> {"N2", "O2"});
        canary::detail::with_input(species_owner, [&](const auto& species) {
            constexpr std::array temps {300.0, 1000.0, 5000.0, 6000.0};
            constexpr std::array<std::array<double, 4>, 2> janaf {
              {
                {29.125, 32.697, 37.912, 38.276},
                {29.385, 34.870, 42.675, 44.387},
              },
            };
            const double r_u = canary::constants::universal_gas_constant<double>;
            for (std::size_t s = 0; s < 2; ++s) {
                const auto curve = canary::make_nasa9_curve(species[s].fit, 1.0);
                for (std::size_t k = 0; k < temps.size(); ++k) {
                    const double cp = r_u * cp_over_r(curve, temps[k]);
                    expect(std::abs(cp - janaf[s][k]) <= 1.5e-2 * janaf[s][k])
                      << species[s].name << temps[k] << cp;
                }
            }
        });
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
        equilibrium.temperature =
          db::PhysicalTemperatureRecord {.kind = "translation", .species = {}};
        expect(file.reactions[12].backward == equilibrium);
    };

    "a reactions file reads, writes and reads again to the same data and text"_test = [] {
        const db::ReactionsFile file = db::read_reactions_file(kim_jo_path());
        const std::string text = db::write_reactions(file);
        const db::ReactionsFile again = db::parse_reactions(text);
        expect(again == file);
        expect(db::write_reactions(again) == text);
        const TemporaryDirectory directory;
        const auto path = directory.path() / "reactions.json";
        db::write_reactions_file(path, file);
        expect(db::read_reactions_file(path) == file);
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
          [] { db::parse_reactions(kim_jo_with("canary-reactions-2", "canary-reactions-1")); }));
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
        const auto input_owner = db::mechanism_input<13, 5>(file);
        canary::detail::with_input(input_owner, [&](const auto& input) {
            const auto& compiled = canary::air5_kim_jo;
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
                expect(same_temperature(a.equilibrium_temperature, b.equilibrium_temperature)) << k;
                expect(a.third_body == b.third_body) << k;
                expect(same_coefficients(a.efficiencies, b.efficiencies)) << k;
            }
            // The mechanism built from the file is the compiled one.
            const canary::UnitScales<double, canary::si> scales {};
            const auto from_file =
              canary::make_reaction_mechanism(input,
                                              canary::reaction_species(canary::air5_species),
                                              scales);
            const auto from_header = canary::make_air5_kim_jo(scales);
            for (const auto k : canary::reaction_indices(0, 13)) {
                const auto a_records = from_file.reactions();
                const auto& a = a_records[k];
                const auto b_records = from_header.reactions();
                const auto& b = b_records[k];
                for (const auto i : canary::species_indices(0, 5)) {
                    expect(a.reactants()[i] == b.reactants()[i] &&
                           a.products()[i] == b.products()[i])
                      << k.value() << i.value();
                    expect(a.efficiencies()[i] == b.efficiencies()[i]) << k.value() << i.value();
                }
                expect(a.coefficient_change == b.coefficient_change) << k.value();
                expect(a.forward.A == b.forward.A && a.forward.n == b.forward.n &&
                       a.forward.activation_temperature == b.forward.activation_temperature)
                  << k.value();
                expect(a.backward_model == b.backward_model && a.third_body == b.third_body)
                  << k.value();
            }
        });
    };

    "reactions conversion errors name the reaction and the problem"_test = [] {
        const db::ReactionsFile file = db::read_reactions_file(kim_jo_path());
        expect(throws_parameter_error([&] { (void)db::mechanism_input<12, 5>(file); }));
        // NO + M lists five efficiencies.
        expect(throws_parameter_error([&] { (void)db::mechanism_input<13, 4>(file); }));
        const auto message = [](const db::ReactionsFile& f) -> std::string {
            try {
                (void)db::mechanism_input<13, 5>(f);
            } catch (const canary::ParameterError& error) {
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
        const auto input_owner = db::mechanism_input<13, 5>(changed);
        canary::detail::with_input(input_owner, [&](const auto& input) {
            const canary::UnitScales<double, canary::si> scales {};
            expect(throws_parameter_error([&] {
                (void)canary::make_reaction_mechanism(
                  input,
                  canary::reaction_species(canary::air5_species),
                  scales);
            }));
            changed = file;
            changed.reactions[8].forward.A = -1.0;
            const auto negative_owner = db::mechanism_input<13, 5>(changed);
            canary::detail::with_input(negative_owner, [&](const auto& negative) {
                expect(throws_parameter_error([&] {
                    (void)canary::make_reaction_mechanism(
                      negative,
                      canary::reaction_species(canary::air5_species),
                      scales);
                }));
            });
        });
    };

    "model inputs from the Park file equal the compiled data, with translational exponents"_test =
      [] {
          const db::ReactionsFile file = db::read_reactions_file(park_path());
          expect(file.reactions[0].forward.translational_exponent == 0.5);
          expect(!file.reactions[3].forward.translational_exponent.has_value());
          const auto input_owner = db::mechanism_input<5, 5>(file);
          canary::detail::with_input(input_owner, [&](const auto& input) {
              const auto& compiled = canary::air5_park;
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
                  expect(same_rate(a.backward, b.backward)) << k;
                  expect(same_temperature(a.equilibrium_temperature, b.equilibrium_temperature))
                    << k;
                  expect(a.third_body == b.third_body) << k;
                  expect(same_coefficients(a.efficiencies, b.efficiencies)) << k;
              }
              expect(input.reactions[0].forward.translational_exponent == 0.5);
              expect(input.reactions[3].forward.translational_exponent == 1.0);
              // The exponent survives writing and reading.
              expect(db::parse_reactions(db::write_reactions(file)) == file);
          });
      };

    "explicit backward Arrhenius inputs preserve every rate field"_test = [] {
        auto file = db::read_reactions_file(air11_kim_jo_path());
        file.reactions[38].backward = {
          .model = "arrhenius",
          .A = 2.5e8,
          .n = -0.75,
          .activation_temperature = 12500.0,
          .translational_exponent = 0.25,
          .other_temperature =
            db::PhysicalTemperatureRecord {.kind = "free_electron", .species = {}},
          .temperature = {},
        };
        const auto owner = db::mechanism_input<40, 11>(
          db::parse_reactions(db::write_reactions(file)));
        canary::detail::with_input(owner, [&](const auto& input) {
            const canary::ArrheniusInput expected {
              .A = 2.5e8,
              .n = -0.75,
              .activation_temperature = 12500.0,
              .translational_exponent = 0.25,
              .other_temperature = canary::FreeElectronTemperatureInput {},
            };
            expect(input.reactions[38].backward_model == canary::BackwardRateModel::arrhenius);
            expect(same_rate(input.reactions[38].backward, expected));
        });
    };

    "physical temperature conversion validates kind and applicable species fields"_test = [] {
        const auto original = db::read_reactions_file(air11_kim_jo_path());
        for (const auto& bad : std::array {
               db::PhysicalTemperatureRecord {.kind = "unknown", .species = {}},
               db::PhysicalTemperatureRecord {.kind = "vibration", .species = {}},
               db::PhysicalTemperatureRecord {.kind = "bound_electronic", .species = ""},
               db::PhysicalTemperatureRecord {.kind = "translation", .species = "N2"},
               db::PhysicalTemperatureRecord {.kind = "rotation", .species = "N2"},
               db::PhysicalTemperatureRecord {.kind = "free_electron", .species = "e-"}}) {
            for (const bool equilibrium : {false, true}) {
                auto file = original;
                if (equilibrium) {
                    file.reactions[0].backward.temperature = bad;
                } else {
                    file.reactions[0].forward.other_temperature = bad;
                }
                expect(throws_parameter_error([&] { (void)db::mechanism_input<40, 11>(file); }));
            }
        }
        auto file = original;
        file.reactions[0].forward.temperature = file.reactions[0].backward.temperature;
        expect(throws_parameter_error([&] { (void)db::mechanism_input<40, 11>(file); }));
        file = original;
        file.reactions[0].backward.other_temperature = file.reactions[0].forward.other_temperature;
        expect(throws_parameter_error([&] { (void)db::mechanism_input<40, 11>(file); }));
    };

    "physical temperature owners retain short and long species names"_test = [] {
        using namespace canary;
        const UnitScales<double, si> scales {};
        for (const bool electronic : {false, true}) {
            for (const std::string& name :
                 {std::string {"X"}, std::string {"nitrogen_species_with_a_long_owned_name"}}) {
                const auto factory = [&] {
                    auto species = air5_species;
                    std::string transient = name;
                    species[0].name = transient.c_str();
                    auto r = air5_park.reactions[0];
                    r.reactants.entries[0].name = transient.c_str();
                    for (auto& entry : r.efficiencies.entries) {
                        if (entry.name && std::string_view(entry.name) == "N2") {
                            entry.name = transient.c_str();
                        }
                    }
                    const PhysicalTemperatureInput requirement =
                      electronic
                        ? PhysicalTemperatureInput {BoundElectronicTemperatureInput {
                            transient.c_str()}}
                        : PhysicalTemperatureInput {VibrationTemperatureInput {transient.c_str()}};
                    r.forward.other_temperature = requirement;
                    r.equilibrium_temperature = requirement;
                    const MechanismInput<1, 5> input {.temperature_lower = 300.0,
                                                      .temperature_upper = 30000.0,
                                                      .reactions = {r}};
                    return std::pair {OwnedInput {input}, OwnedInput {species}};
                };
                const auto copied = [&] {
                    const auto source = factory();
                    return decltype(source)(source);
                }();
                const auto moved = [&] {
                    auto source = copied;
                    return decltype(source)(std::move(source));
                }();
                const auto gas =
                  make_two_temperature_mixture(moved.second, air5_rotational_degrees, scales);
                const auto mechanism = make_reaction_mechanism(moved.first,
                                                               reaction_species(moved.second),
                                                               scales,
                                                               Pressure<double, si> {101325.0},
                                                               gas.model);
                expect(mechanism.reactions()[ReactionIndex {0}].forward.other_mode ==
                       OptionalModeIndex {ModeIndex {0}});
                expect(mechanism.reactions()[ReactionIndex {0}].backward.other_mode ==
                       OptionalModeIndex {ModeIndex {0}});
            }
        }
    };

    "a forward rate requires an explicit model"_test = [] {
        expect(throws_database_error([] {
            (void)db::parse_reactions(kim_jo_with(R"("model": "arrhenius",)", ""));
        }));
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
        const TemporaryDirectory directory;
        const auto path = directory.path() / "exchange.json";
        db::write_energy_exchange_file(path, file);
        expect(db::read_energy_exchange_file(path) == file);
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
                exchange_with("canary-energy-exchange-1", "canary-energy-exchange-2"));
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

    "every exchange model field is required"_test = [] {
        for (const auto& path : {exchange_path(), air11_exchange_path()}) {
            const auto text = read_text(path);
            expect(db::parse_energy_exchange(text) == db::read_energy_exchange_file(path));
            std::size_t position = text.find("\"model\"");
            expect(fatal(position != std::string::npos));
            while (position != std::string::npos) {
                const auto value = text.find('"', text.find(':', position) + 1);
                const auto end = text.find('"', value + 1) + 1;
                const auto next = text.find_first_not_of(" \t\r\n", end);
                expect(fatal(next != std::string::npos));
                expect(fatal(text[next] == ',' || text[next] == '}'));
                auto changed = text;
                changed.erase(position, (text[next] == ',' ? next + 1 : end) - position);
                expect(throws_database_error([&] { (void)db::parse_energy_exchange(changed); }))
                  << path << " model at byte " << position;
                position = text.find("\"model\"", end);
            }
        }
    };

    "exchange units require every field in the declared units"_test = [] {
        for (const std::string_view field : {"energy", "temperature", "area"}) {
            auto file = db::read_energy_exchange_file(air11_exchange_path());
            if (field == "energy") { file.units.energy = "J/kmol"; }
            if (field == "temperature") { file.units.temperature = "degC"; }
            if (field == "area") { file.units.area = "cm^2"; }
            expect(throws_database_error([&] {
                (void)db::parse_energy_exchange(db::write_energy_exchange(file));
            })) << field;
        }
        for (const std::string_view token : {R"("energy": "J/mol",)",
                                             R"("temperature": "K",)",
                                             R"(,
      "area": "m^2")"}) {
            expect(throws_database_error([&] {
                (void)db::parse_energy_exchange(exchange_with(token, ""));
            })) << token;
        }
    };

    "the shipped energy-exchange schema is the schema of the current file structs"_test = [] {
        expect(db::energy_exchange_schema() ==
               read_text(source_path("database/energy_exchange.schema.json")));
    };

    "model inputs from the energy-exchange file equal the compiled data"_test = [] {
        const db::EnergyExchangeFile file = db::read_energy_exchange_file(exchange_path());
        const auto input_owner = db::energy_exchange_input<15, 3>(file);
        canary::detail::with_input(input_owner, [&](const auto& input) {
            const auto& compiled = canary::air5_park_exchange;
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
            const canary::UnitScales<double, canary::si> scales {};
            const auto built = canary::make_air5<double>(scales);
            const auto mechanism = canary::make_air5_park(scales);
            (void)canary::make_energy_exchange(input, air5_names, built.species, mechanism, scales);
            (void)canary::make_air5_park_exchange(built.species, mechanism, scales);
        });
    };

    "energy-exchange corrections require only their model-specific fields"_test = [] {
        const auto original = db::read_energy_exchange_file(exchange_path());
        for (const std::string_view model : {"park_1985", "park_1993", "kim_2021"}) {
            for (unsigned fields = 0; fields < 4; ++fields) {
                auto file = original;
                auto& correction =
                  *file.vibration_translation[0].relaxation_time.high_temperature_correction;
                correction.model = model;
                correction.reference_temperature =
                  (fields & 1U) ? std::optional {3000.0} : std::nullopt;
                correction.exponent = (fields & 2U) ? std::optional {0.0} : std::nullopt;
                const bool valid = (model == "park_1985" && fields == 0) ||
                                   (model == "park_1993" && fields == 1) ||
                                   (model == "kim_2021" && fields == 2);
                if (valid) {
                    const auto parsed = db::parse_energy_exchange(db::write_energy_exchange(file));
                    const auto input_owner = db::energy_exchange_input<15, 3>(parsed);
                    canary::detail::with_input(input_owner, [&](const auto& input) {
                        expect(input.vibration_translation[0].correction.reference_temperature ==
                               correction.reference_temperature.value_or(0.0));
                        expect(input.vibration_translation[0].correction.exponent == 0.0);
                    });
                } else {
                    std::string message;
                    try {
                        (void)db::energy_exchange_input<15, 3>(file);
                    } catch (const canary::ParameterError& error) {
                        message = error.what();
                    }
                    expect(message.starts_with("vibration_translation 0 (N2-N): ")) << message;
                    expect(message.find(model) != std::string::npos) << message;
                    expect(message.find("requires") != std::string::npos ||
                           message.find("forbids") != std::string::npos)
                      << message;
                }
            }
        }
    };

    "energy-exchange conversion errors name the entry and the problem"_test = [] {
        const db::EnergyExchangeFile file = db::read_energy_exchange_file(exchange_path());
        expect(throws_parameter_error([&] { (void)db::energy_exchange_input<14, 3>(file); }));
        expect(throws_parameter_error([&] { (void)db::energy_exchange_input<15, 2>(file); }));
        const auto message = [](const db::EnergyExchangeFile& f) -> std::string {
            try {
                (void)db::energy_exchange_input<15, 3>(f);
            } catch (const canary::ParameterError& error) {
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
        const auto input_owner = db::species_transport<3>(file, air5_names);
        canary::detail::with_input(input_owner, [&](const auto& input) {
            for (std::size_t i = 0; i < 5; ++i) {
                const auto& a = input[i];
                const auto& b = canary::air5_species_transport[i];
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
        });
    };

    "transport conversion errors name the species and the problem"_test = [] {
        const db::SpeciesFile file = db::read_species_file(air5_path());
        const auto message = [](const db::SpeciesFile& f) -> std::string {
            try {
                (void)db::species_transport<3>(f, air5_names);
            } catch (const canary::ParameterError& error) {
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
        const TemporaryDirectory directory;
        const auto path = directory.path() / "collisions.json";
        db::write_collision_integrals_file(path, file);
        expect(db::read_collision_integrals_file(path) == file);
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
                collision_with("canary-collision-integrals-1", "canary-collision-integrals-2"));
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
        const auto input_owner = db::collision_integral_pairs<15>(file);
        canary::detail::with_input(input_owner, [&](const auto& input) {
            const auto same = [](const canary::CollisionFitInput& a,
                                 const canary::CollisionFitInput& b) {
                return a.A == b.A && a.B == b.B && a.C == b.C && a.D == b.D;
            };
            for (std::size_t k = 0; k < 15; ++k) {
                const auto& a = input[k];
                const auto& b = canary::air5_collision_integrals[k];
                expect(std::string_view(a.first) == b.first &&
                       std::string_view(a.second) == b.second)
                  << k;
                expect(same(a.pi_omega_11, b.pi_omega_11) && same(a.pi_omega_22, b.pi_omega_22))
                  << k;
            }
            expect(throws_parameter_error([&] { (void)db::collision_integral_pairs<14>(file); }));
        });
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
        const auto input5_owner = db::electronic_levels<8>(air5, air5_names);
        canary::detail::with_input(input5_owner, [&](const auto& input5) {
            for (std::size_t i = 0; i < 5; ++i) {
                expect(same(input5[i], canary::air5_electronic_levels[i])) << i;
            }
            const db::SpeciesFile air11 = db::read_species_file(air11_path());
            const auto input11_owner = db::electronic_levels<8>(air11, air11_names);
            canary::detail::with_input(input11_owner, [&](const auto& input11) {
                for (std::size_t i = 0; i < 11; ++i) {
                    expect(same(input11[i], canary::air11_electronic_levels[i])) << i;
                }
                const auto species11_owner = db::nasa9_species<4>(air11, air11_names);
                canary::detail::with_input(species11_owner, [&](const auto& species11) {
                    for (std::size_t i = 0; i < 11; ++i) {
                        const auto& a = species11[i];
                        const auto& b = canary::air11_species[i];
                        expect(std::string_view(a.name) == b.name && a.molar_mass == b.molar_mass &&
                               a.charge == b.charge && a.fit.segments == b.fit.segments &&
                               a.fit.breaks == b.fit.breaks && a.fit.blends == b.fit.blends &&
                               a.fit.a == b.fit.a)
                          << i;
                    }
                    // N2 has 7 levels.
                    expect(throws_parameter_error(
                      [&] { (void)db::electronic_levels<6>(air5, air5_names); }));
                    db::SpeciesFile changed = air5;
                    changed.species.at("O").electronic_levels->temperatures.pop_back();
                    expect(throws_parameter_error(
                      [&] { (void)db::electronic_levels<8>(changed, air5_names); }));
                    changed = air5;
                    changed.species.erase("NO");
                    expect(throws_database_error(
                      [&] { (void)db::electronic_levels<8>(changed, air5_names); }));
                });
            });
        });
    };

    "the eleven-species Kim-Jo mechanism reads with electron temperatures"_test = [] {
        const db::ReactionsFile file = db::read_reactions_file(air11_kim_jo_path());
        expect(file.reactions.size() == 40_u);
        const db::ReactionRecord& n2e = file.reactions[21];
        expect(n2e.equation == "N2 + e- <=> N + N + e-");
        expect(
          n2e.forward.other_temperature ==
          std::optional {db::PhysicalTemperatureRecord {.kind = "free_electron", .species = {}}});
        expect(
          n2e.backward.temperature ==
          std::optional {db::PhysicalTemperatureRecord {.kind = "free_electron", .species = {}}});
        expect(
          file.reactions[0].forward.other_temperature ==
          std::optional {db::PhysicalTemperatureRecord {.kind = "vibration", .species = "N2"}});
        const auto input_owner = db::mechanism_input<40, 11>(file);
        canary::detail::with_input(input_owner, [&](const auto& input) {
            expect(input.temperature_lower == canary::air11_kim_jo.temperature_lower);
            expect(input.temperature_upper == canary::air11_kim_jo.temperature_upper);
            for (std::size_t k = 0; k < 40; ++k) {
                const auto& a = input.reactions[k];
                const auto& b = canary::air11_kim_jo.reactions[k];
                expect(std::string_view(a.equation) == b.equation) << k;
                expect(same_rate(a.forward, b.forward)) << k;
                expect(same_coefficients(a.reactants, b.reactants)) << k;
                expect(same_coefficients(a.products, b.products)) << k;
                expect(a.backward_model == b.backward_model) << k;
                expect(same_rate(a.backward, b.backward)) << k;
                expect(same_temperature(a.equilibrium_temperature, b.equilibrium_temperature)) << k;
                expect(a.third_body == b.third_body) << k;
                expect(same_coefficients(a.efficiencies, b.efficiencies)) << k;
            }
            expect(std::holds_alternative<canary::FreeElectronTemperatureInput>(
              *input.reactions[38].equilibrium_temperature));
        });
    };

    "reaction species adapters retain owned definitions through source destruction"_test = [] {
        const canary::UnitScales<double, canary::si> scales {};
        const auto adapted = [] {
            auto file = db::read_species_file(air11_path());
            auto owner = db::nasa9_species<4>(file, air11_names);
            return canary::reaction_species(std::move(owner));
        }();
        const auto copied = adapted;
        const auto mechanism =
          db::mechanism_input<40, 11>(db::read_reactions_file(air11_kim_jo_path()));
        const auto gas = canary::make_air11_three_temperature(scales);
        const auto reference = canary::make_air11_kim_jo(scales);
        const auto result = canary::make_reaction_mechanism(mechanism, copied, scales);
        for (const auto i : canary::species_indices(0, 11)) {
            expect(result.mass()[i] == gas.species.mass[i.value()]);
            expect(result.reactions()[canary::ReactionIndex {38}].products()[i] ==
                   reference.reactions()[canary::ReactionIndex {38}].products()[i]);
        }
        auto file = db::read_species_file(air11_path());
        file.species.at("O+").charge = 2;
        const auto invalid = canary::reaction_species(db::nasa9_species<4>(file, air11_names));
        expect(throws_parameter_error(
          [&] { (void)canary::make_reaction_mechanism(mechanism, invalid, scales); }));
        file = db::read_species_file(air11_path());
        file.species.at("O").elements.reset();
        expect(throws_parameter_error([&] {
            (void)canary::make_reaction_mechanism(
              mechanism,
              canary::reaction_species(db::nasa9_species<4>(file, air11_names)),
              scales);
        }));
    };

    "exchange JSON rejects noninteger and unrepresentable ionisation reaction indices"_test = [] {
        const auto json = db::write_energy_exchange(
          db::read_energy_exchange_file(air11_exchange_path()));
        const auto section = json.find("\"electron_impact_ionisation\"");
        expect(fatal(section != std::string::npos));
        const auto key = json.find("\"reaction\"", section);
        expect(fatal(key != std::string::npos));
        const auto start = json.find_first_not_of(" \t\r\n", json.find(':', key) + 1);
        const auto end = json.find_first_not_of("0123456789", start);
        expect(fatal(start != std::string::npos && end != std::string::npos && end > start));
        for (const auto token : {"0.5", "true", "\"38\"", "2147483648", "-2147483649"}) {
            auto changed = json;
            changed.replace(start, end - start, token);
            expect(throws_database_error([&] { (void)db::parse_energy_exchange(changed); })) << token;
        }
    };

    "file-backed ionisation construction checks reaction bounds and charge transfer"_test = [] {
        const auto file = db::read_energy_exchange_file(air11_exchange_path());
        const auto species_file = db::read_species_file(air11_path());
        const auto heavy = db::energy_exchange_input<30, 21>(file);
        const canary::UnitScales<double, canary::si> scales {};
        const auto gas = canary::make_air11_three_temperature(scales);
        const auto mechanism =
          canary::make_air11_kim_jo(scales, canary::Pressure<double, canary::si> {1e5}, gas.model);
        const auto electron = *canary::electron_species_index(gas.model);
        const auto ion = *canary::detail::species_index(air11_names, "O+");
        const auto spectator = *canary::detail::species_index(air11_names, "N2");
        const auto check = [&](const auto& electrons,
                               const auto& species,
                               const auto& reactions,
                               std::string_view reason) {
            std::string message;
            try {
                (void)canary::make_multi_temperature_exchange(heavy,
                                                              electrons,
                                                              air11_names,
                                                              species,
                                                              reactions,
                                                              gas.model,
                                                              scales);
            } catch (const canary::ParameterError& error) {
                message = error.what();
            }
            if (reason.empty()) {
                expect(message.empty()) << message;
            } else {
                expect(message.starts_with("electron-impact ionisation 0 (O): ")) << message;
                expect(message.find(reason) != std::string::npos) << message;
            }
        };
        const auto electrons = [&](const auto& input_file) {
            return db::electron_exchange_input<10, 1, 2>(input_file, species_file);
        };
        check(electrons(file), gas.species, mechanism, "");
        for (const int index : {-1, 40}) {
            auto changed = file;
            changed.electron_impact_ionisation[0].reaction = index;
            check(electrons(changed),
                  gas.species,
                  mechanism,
                  "reaction index must be between 0 and 39");
        }
        const auto source = *canary::detail::species_index(air11_names, "O");
        const auto check_record = [&](const auto& record, std::string_view reason) {
            std::string message;
            try {
                canary::detail::check_ionisation_participants(record,
                                                              air11_names,
                                                              gas.species,
                                                              source,
                                                              electron,
                                                              canary::ReactionIndex {38});
            } catch (const canary::ParameterError& error) {
                message = error.what();
            }
            expect(message.find(reason) != std::string::npos) << message;
        };
        auto changed = mechanism.reactions()[canary::ReactionIndex {38}];
        changed.products()[electron] = 1;
        check_record(changed, "net electron production +1");
        auto species = gas.species;
        species.charge[ion.value()] =
          canary::Charge<double, canary::si> {2.0 * species.charge[ion.value()].value()};
        // O+ is also a coulomb partner, whose charge is checked before the
        // ionisations; without electron-translation pairs the ionisation
        // check is reached.
        auto ionisation_only = file;
        ionisation_only.electron_translation.clear();
        check(db::electron_exchange_input<0, 1, 2>(ionisation_only, species_file),
              species,
              mechanism,
              "must have charge opposite to the electron charge");
        changed = mechanism.reactions()[canary::ReactionIndex {38}];
        changed.products()[spectator] += 1;
        check_record(changed, "no other net species changes");
    };

    "electron exchange conversion rejects incomplete parameters and unsupported models"_test = [] {
        const auto file = db::read_energy_exchange_file(air11_exchange_path());
        const auto species = db::read_species_file(air11_path());
        for (const auto member : {&db::ElectronCrossSectionRecord::a,
                                  &db::ElectronCrossSectionRecord::b,
                                  &db::ElectronCrossSectionRecord::c,
                                  &db::ElectronCrossSectionRecord::maximum_temperature}) {
            auto changed = file;
            (changed.electron_translation[0].cross_section.*member).reset();
            expect(throws_parameter_error([&] {
                (void)db::electron_exchange_input<10, 1, 2>(changed, species);
            }));
        }
        auto changed = file;
        changed.electron_translation[0].model = "inelastic";
        expect(throws_parameter_error([&] {
            (void)db::electron_exchange_input<10, 1, 2>(changed, species);
        }));
        changed = file;
        changed.vibration_electron[0].model = "ssh";
        expect(throws_parameter_error([&] {
            (void)db::electron_exchange_input<10, 1, 2>(changed, species);
        }));
        expect(throws_parameter_error([&] {
            (void)db::electron_exchange_input<10, 0, 2>(file, species);
        }));
        expect(throws_parameter_error([&] {
            (void)db::electron_exchange_input<10, 1, 1>(file, species);
        }));
        auto absent = species;
        absent.species.erase("O");
        expect(throws_database_error([&] {
            (void)db::electron_exchange_input<10, 1, 2>(file, absent);
        }));
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
        const auto electrons_owner = db::electron_exchange_input<10, 1, 2>(file, species);
        canary::detail::with_input(electrons_owner, [&](const auto& electrons) {
            const auto& compiled = canary::air11_kim_jo_exchange_electrons;
            for (std::size_t k = 0; k < 10; ++k) {
                const auto& a = electrons.electron_translation[k];
                const auto& b = compiled.electron_translation[k];
                expect(std::string_view(a.partner) == b.partner &&
                       a.cross_section == b.cross_section && a.a == b.a && a.b == b.b &&
                       a.c == b.c && a.maximum_temperature == b.maximum_temperature)
                  << k;
            }
            for (std::size_t k = 0; k < electrons.vibration_electron.size(); ++k) {
                const auto& a = electrons.vibration_electron[k];
                const auto& b = compiled.vibration_electron[k];
                expect(std::string_view(a.species) == b.species) << k;
                expect(a.relaxation_time.a_low == b.relaxation_time.a_low &&
                       a.relaxation_time.b_low == b.relaxation_time.b_low &&
                       a.relaxation_time.c_low == b.relaxation_time.c_low &&
                       a.relaxation_time.a_high == b.relaxation_time.a_high &&
                       a.relaxation_time.b_high == b.relaxation_time.b_high &&
                       a.relaxation_time.c_high == b.relaxation_time.c_high &&
                       a.relaxation_time.blend_lower == b.relaxation_time.blend_lower &&
                       a.relaxation_time.blend_upper == b.relaxation_time.blend_upper)
                  << k;
            }
            for (std::size_t k = 0; k < 2; ++k) {
                const auto& a = electrons.electron_impact_ionisation[k];
                const auto& b = compiled.electron_impact_ionisation[k];
                expect(a.reaction == b.reaction && std::string_view(a.species) == b.species &&
                       a.ionisation_energy == b.ionisation_energy)
                  << k;
            }
            const auto input_owner = db::energy_exchange_input<30, 21>(file);
            canary::detail::with_input(input_owner, [&](const auto& input) {
                for (std::size_t k = 0; k < 30; ++k) {
                    const auto& a = input.vibration_translation[k];
                    const auto& b = canary::air11_kim_jo_exchange.vibration_translation[k];
                    expect(std::string_view(a.species) == b.species &&
                           std::string_view(a.partner) == b.partner)
                      << k;
                    expect(same_correction(a.correction, b.correction) &&
                           a.millikan_white.a == b.millikan_white.a &&
                           a.millikan_white.b == b.millikan_white.b)
                      << k;
                }
                for (std::size_t k = 0; k < input.chemistry_vibration.size(); ++k) {
                    const auto& a = input.chemistry_vibration[k];
                    const auto& b = canary::air11_kim_jo_exchange.chemistry_vibration[k];
                    expect(a.reaction == b.reaction && std::string_view(a.species) == b.species &&
                           a.dissociation_energy == b.dissociation_energy &&
                           a.vibrational_temperature == b.vibrational_temperature)
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
            });
        });
    };
}
