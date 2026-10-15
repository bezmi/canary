// Construct models from headers emitted by the Python database tools.
#include "generated_exchange.hpp"
#include "generated_reactions.hpp"
#include "generated_species.hpp"
#include "generated_transport.hpp"
#include "negative_capacity.hpp"
#include "small_mass.hpp"

#include <boost/ut.hpp>
#include <canary/build.hpp>
#include <canary/data/air5.hpp>

namespace
{
using namespace boost::ut;

template <typename P, canary::UnitSystem U>
void construct_generated(const canary::UnitScales<P, U>& scales)
{
    const auto built =
      canary::make_thermally_perfect_mixture(canary::generated_species_species, scales);
    const auto mechanism = canary::make_generated_reactions(built.species, scales);
    expect(
      nothrow([&] { (void)canary::make_generated_exchange(built.species, mechanism, scales); }));
    expect(nothrow([&] {
        (void)canary::make_wilke_transport(canary::generated_species_transport,
                                           canary::detail::generated_reactions_species_names,
                                           built.species,
                                           scales);
    }));
    expect(nothrow([&] {
        (void)canary::make_gupta_yos_transport(canary::generated_collision_integrals,
                                               canary::detail::generated_reactions_species_names,
                                               built.species,
                                               canary::air5_rotational_degrees,
                                               scales);
    }));
    expect(nothrow([&] {
        (void)canary::make_multi_temperature_mixture(canary::generated_species_species,
                                                     canary::air5_rotational_degrees,
                                                     canary::generated_species_electronic_levels,
                                                     canary::air5_four_temperature_modes,
                                                     scales);
    }));
}

suite generated_database = [] {
    "generated inputs construct with float and double in SI and scaled units"_test = [] {
        construct_generated(canary::UnitScales<float, canary::si> {});
        construct_generated(canary::UnitScales<double, canary::si> {});
        construct_generated(canary::UnitScales<float, canary::scaled> {
          {.mass = 4.65e-26, .density = 1.0e-3, .temperature = 1000.0, .length = 1.0}});
        construct_generated(canary::UnitScales<double, canary::scaled> {
          {.mass = 4.65e-26, .density = 1.0e-3, .temperature = 1000.0, .length = 1.0}});
    };

    "construction rejects nonphysical source-valid generated coefficients"_test = [] {
        expect(throws<canary::ParameterError>([] {
            (void)canary::make_thermally_perfect_mixture(canary::negative_capacity_species,
                                                         canary::UnitScales<double, canary::si> {});
        }));
    };

    "generated parameters are checked for the selected storage precision"_test = [] {
        expect(nothrow([] {
            (void)canary::make_thermally_perfect_mixture(canary::small_mass_species,
                                                         canary::UnitScales<double, canary::si> {});
        }));
        expect(throws<canary::ParameterError>([] {
            (void)canary::make_thermally_perfect_mixture(canary::small_mass_species,
                                                         canary::UnitScales<float, canary::si> {});
        }));
    };
};
}  // namespace

int main() {}
