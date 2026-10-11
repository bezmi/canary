#ifndef GASLIB_DATA_AIR11_HPP
#define GASLIB_DATA_AIR11_HPP

// Eleven-species air: N2, O2, N, O, NO, N+, O+, N2+, O2+, NO+ and e-, in
// this order (air11_species and air11_electronic_levels).

#include <array>
#include <concepts>
#include <gaslib/build.hpp>
#include <gaslib/data/air11_species.hpp>
#include <gaslib/models/multi_temperature.hpp>
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/species.hpp>
#include <gaslib/units.hpp>

namespace gaslib
{

// Rotational degrees of freedom of air11_species: 2 for the linear
// molecules N2, O2, NO, N2+, O2+ and NO+, 0 for the atoms, the atomic ions
// and e-.
inline constexpr std::array<int, 11> air11_rotational_degrees {2, 2, 0, 0, 2, 0, 0, 2, 2, 2, 0};

// Thermally perfect eleven-species air.
template <std::floating_point P, UnitSystem U> using Air11 = ThermallyPerfectMixture<P, U, 11, 4>;

// Air11 model of air11_species in the unit system of the scales, with its
// species data. See make_thermally_perfect_mixture.
template <std::floating_point P, UnitSystem U>
Built<Air11<P, U>, SpeciesData<P, U, 11>> make_air11(const UnitScales<P, U>& scales)
{ return make_thermally_perfect_mixture(air11_species, scales); }

// Energy modes of three-temperature eleven-species air, in this order:
// - vibration: the vibration of N2, O2, NO, N2+, O2+ and NO+;
// - electron: the electronic energy of the ten heavy species and the
//   electron energy of e-.
inline constexpr std::array<EnergyModeInput<11>, 2> air11_three_temperature_modes {
  {
    {
      .count = 6,
      .components =
        {
          {
            {.species = "N2", .component = ModeComponent::vibration},
            {.species = "O2", .component = ModeComponent::vibration},
            {.species = "NO", .component = ModeComponent::vibration},
            {.species = "N2+", .component = ModeComponent::vibration},
            {.species = "O2+", .component = ModeComponent::vibration},
            {.species = "NO+", .component = ModeComponent::vibration},
            {.species = nullptr, .component = ModeComponent::vibration},
            {.species = nullptr, .component = ModeComponent::vibration},
            {.species = nullptr, .component = ModeComponent::vibration},
            {.species = nullptr, .component = ModeComponent::vibration},
            {.species = nullptr, .component = ModeComponent::vibration},
          },
        },
    },
    {
      .count = 11,
      .components =
        {
          {
            {.species = "N2", .component = ModeComponent::electronic},
            {.species = "O2", .component = ModeComponent::electronic},
            {.species = "N", .component = ModeComponent::electronic},
            {.species = "O", .component = ModeComponent::electronic},
            {.species = "NO", .component = ModeComponent::electronic},
            {.species = "N+", .component = ModeComponent::electronic},
            {.species = "O+", .component = ModeComponent::electronic},
            {.species = "N2+", .component = ModeComponent::electronic},
            {.species = "O2+", .component = ModeComponent::electronic},
            {.species = "NO+", .component = ModeComponent::electronic},
            {.species = "e-", .component = ModeComponent::electron},
          },
        },
    },
  },
};

// Three-temperature eleven-species air.
template <std::floating_point P, UnitSystem U>
using Air11ThreeTemperature = MultiTemperatureMixture<P, U, 11, 2, 4, 8>;

// Air11ThreeTemperature model of air11_species, air11_rotational_degrees,
// air11_electronic_levels and air11_three_temperature_modes in the unit
// system of the scales, with its species data. See
// make_multi_temperature_mixture.
template <std::floating_point P, UnitSystem U>
Built<Air11ThreeTemperature<P, U>, SpeciesData<P, U, 11>>
  make_air11_three_temperature(const UnitScales<P, U>& scales)
{
    return make_multi_temperature_mixture(air11_species,
                                          air11_rotational_degrees,
                                          air11_electronic_levels,
                                          air11_three_temperature_modes,
                                          scales);
}

}  // namespace gaslib

#endif  // GASLIB_DATA_AIR11_HPP
