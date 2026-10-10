#ifndef GASLIB_DATA_AIR5_HPP
#define GASLIB_DATA_AIR5_HPP

// Five-species air: N2, O2, N, O and NO, in this order.
// Each species cites the source of its fit and of its molar mass. The
// values were taken from the species database distributed with GDTk
// (gdtk-uq/gdtk, commit bf559c3f, src/gas/species-database), which cites
// the same sources; they are data only.

#include <array>
#include <concepts>
#include <gaslib/build.hpp>
#include <gaslib/models/thermally_perfect.hpp>
#include <gaslib/models/two_temperature.hpp>
#include <gaslib/quantity.hpp>
#include <gaslib/species.hpp>
#include <gaslib/units.hpp>

namespace gaslib
{

// Molar masses in kg/mol, charges in elementary charges, fits in SI
// (kelvin). Entries past a species' segment count are zero.
inline constexpr std::array<Nasa9Species<4>, 5> air5_species {
  {
    // N2. Fit: NASA CEA2 thermo.inp (McBride, Zehe and Gordon, NASA/TP-2002-211556).
    // Molar mass: NASA CEA2 thermo.inp.
    {
      .name = "N2",
      .molar_mass = 28.0134e-3,
      .charge = 0,
      .fit =
        {
          .segments = 3,
          .breaks =
            {
              200.0,
              1000.0,
              6000.0,
              20000.0,
            },
          .blends =
            {
              400.0,
              1000.0,
            },
          .a =
            {
              {
                {
                  2.210371497e+04,
                  -3.818461820e+02,
                  6.082738360e+00,
                  -8.530914410e-03,
                  1.384646189e-05,
                  -9.625793620e-09,
                  2.519705809e-12,
                  7.108460860e+02,
                  -1.076003744e+01,
                },
                {
                  5.877124060e+05,
                  -2.239249073e+03,
                  6.066949220e+00,
                  -6.139685500e-04,
                  1.491806679e-07,
                  -1.923105485e-11,
                  1.061954386e-15,
                  1.283210415e+04,
                  -1.586640027e+01,
                },
                {
                  8.310139160e+08,
                  -6.420733540e+05,
                  2.020264635e+02,
                  -3.065092046e-02,
                  2.486903333e-06,
                  -9.705954110e-11,
                  1.437538881e-15,
                  4.938707040e+06,
                  -1.672099740e+03,
                },
              },
            },
        },
    },
    // O2. Fit: NASA CEA2 thermo.inp (McBride, Zehe and Gordon, NASA/TP-2002-211556).
    // Molar mass: NASA CEA2 thermo.inp.
    {
      .name = "O2",
      .molar_mass = 0.03199880,
      .charge = 0,
      .fit =
        {
          .segments = 3,
          .breaks =
            {
              200.0,
              1000.0,
              6000.0,
              20000.0,
            },
          .blends =
            {
              400.0,
              1000.0,
            },
          .a =
            {
              {
                {
                  -3.425563420e+04,
                  4.847000970e+02,
                  1.119010961e+00,
                  4.293889240e-03,
                  -6.836300520e-07,
                  -2.023372700e-09,
                  1.039040018e-12,
                  -3.391454870e+03,
                  1.849699470e+01,
                },
                {
                  -1.037939022e+06,
                  2.344830282e+03,
                  1.819732036e+00,
                  1.267847582e-03,
                  -2.188067988e-07,
                  2.053719572e-11,
                  -8.193467050e-16,
                  -1.689010929e+04,
                  1.738716506e+01,
                },
                {
                  4.975294300e+08,
                  -2.866106874e+05,
                  6.690352250e+01,
                  -6.169959020e-03,
                  3.016396027e-07,
                  -7.421416600e-12,
                  7.278175770e-17,
                  2.293554027e+06,
                  -5.530621610e+02,
                },
              },
            },
        },
    },
    // N. Fit: Johnston et al., AIAA Paper 2015-3110.
    // Molar mass: NASA CEA2 thermo.inp.
    {
      .name = "N",
      .molar_mass = 14.0067e-3,
      .charge = 0,
      .fit =
        {
          .segments = 4,
          .breaks =
            {
              200.0,
              1000.0,
              6000.0,
              20000.0,
              50000.0,
            },
          .blends =
            {
              400.0,
              1000.0,
              1000.0,
            },
          .a =
            {
              {
                {
                  0.000000000e+00,
                  0.000000000e+00,
                  2.500000000e+00,
                  0.000000000e+00,
                  0.000000000e+00,
                  0.000000000e+00,
                  0.000000000e+00,
                  5.61061063e+04,
                  4.19425139e+00,
                },
                {
                  -2.27073277e+05,
                  8.14052944e+02,
                  1.32705137e+00,
                  8.62721731e-04,
                  -3.35747089e-07,
                  6.29010687e-11,
                  -3.90674587e-15,
                  5.10943141e+04,
                  1.22823713e+01,
                },
                {
                  -2.04738994e+09,
                  1.45842847e+06,
                  -4.18833824e+02,
                  6.25994407e-02,
                  -4.96542822e-06,
                  1.98247052e-10,
                  -3.05470194e-15,
                  -1.12727730e+07,
                  3.58487417e+03,
                },
                {
                  5.74291902e+11,
                  -1.29039294e+08,
                  1.15381467e+04,
                  -5.25078568e-01,
                  1.29219090e-05,
                  -1.63974231e-10,
                  8.41878585e-16,
                  1.15261836e+09,
                  -1.11649232e+05,
                },
              },
            },
        },
    },
    // O. Fit: Johnston et al., AIAA Paper 2015-3110.
    // Molar mass: NASA CEA2 thermo.inp.
    {
      .name = "O",
      .molar_mass = 0.01599940,
      .charge = 0,
      .fit =
        {
          .segments = 4,
          .breaks =
            {
              200.0,
              1000.0,
              6000.0,
              20000.0,
              50000.0,
            },
          .blends =
            {
              400.0,
              1000.0,
              1000.0,
            },
          .a =
            {
              {
                {
                  0.00000000e+00,
                  0.00000000e+00,
                  2.50000000e+00,
                  0.00000000e+00,
                  0.00000000e+00,
                  0.00000000e+00,
                  0.00000000e+00,
                  2.92727578e+04,
                  5.20440739e+00,
                },
                {
                  -4.59493000e+05,
                  1.36004299e+03,
                  9.30779889e-01,
                  9.08157296e-04,
                  -2.81746654e-07,
                  4.49451675e-11,
                  -2.68564490e-15,
                  2.06167786e+04,
                  1.63928984e+01,
                },
                {
                  -1.36180488e+09,
                  9.17683684e+05,
                  -2.46884942e+02,
                  3.48520165e-02,
                  -2.61713430e-06,
                  9.91544159e-11,
                  -1.44192000e-15,
                  -7.15328977e+06,
                  2.14015002e+03,
                },
                {
                  4.38121309e+11,
                  -8.98310623e+07,
                  7.33626546e+03,
                  -3.05369185e-01,
                  6.90318054e-06,
                  -8.09257822e-11,
                  3.86322459e-16,
                  8.10576891e+08,
                  -7.16401144e+04,
                },
              },
            },
        },
    },
    // NO. Fit: NASA CEA2 thermo.inp (McBride, Zehe and Gordon, NASA/TP-2002-211556).
    // Molar mass: NASA CEA2 thermo.inp.
    {
      .name = "NO",
      .molar_mass = 0.03000610,
      .charge = 0,
      .fit =
        {
          .segments = 3,
          .breaks =
            {
              200.0,
              1000.0,
              6000.0,
              20000.0,
            },
          .blends =
            {
              400.0,
              1000.0,
            },
          .a =
            {
              {
                {
                  -1.143916503e+04,
                  1.536467592e+02,
                  3.431468730e+00,
                  -2.668592368e-03,
                  8.481399120e-06,
                  -7.685111050e-09,
                  2.386797655e-12,
                  9.098214410e+03,
                  6.728725490e+00,
                },
                {
                  2.239018716e+05,
                  -1.289651623e+03,
                  5.433936030e+00,
                  -3.656034900e-04,
                  9.880966450e-08,
                  -1.416076856e-11,
                  9.380184620e-16,
                  1.750317656e+04,
                  -8.501669090e+00,
                },
                {
                  -9.575303540e+08,
                  5.912434480e+05,
                  -1.384566826e+02,
                  1.694339403e-02,
                  -1.007351096e-06,
                  2.912584076e-11,
                  -3.295109350e-16,
                  -4.677501240e+06,
                  1.242081216e+03,
                },
              },
            },
        },
    },
  },
};

// Thermally perfect five-species air.
template <std::floating_point P, UnitSystem U> using Air5 = ThermallyPerfectMixture<P, U, 5, 4>;

// Air5 model of air5_species in the unit system of the scales, with its
// species data. See make_thermally_perfect_mixture.
template <std::floating_point P, UnitSystem U>
Built<Air5<P, U>, SpeciesData<P, U, 5>> make_air5(const UnitScales<P, U>& scales)
{ return make_thermally_perfect_mixture(air5_species, scales); }

// Rotational degrees of freedom of air5_species: 2 for the linear molecules
// N2, O2 and NO, 0 for the atoms N and O.
inline constexpr std::array<int, 5> air5_rotational_degrees {2, 2, 0, 0, 2};

// Two-temperature five-species air.
template <std::floating_point P, UnitSystem U>
using Air5TwoTemperature = TwoTemperatureMixture<P, U, 5, 4>;

// Air5TwoTemperature model of air5_species and air5_rotational_degrees in
// the unit system of the scales, with its species data. See
// make_two_temperature_mixture.
template <std::floating_point P, UnitSystem U>
Built<Air5TwoTemperature<P, U>, SpeciesData<P, U, 5>>
  make_air5_two_temperature(const UnitScales<P, U>& scales)
{ return make_two_temperature_mixture(air5_species, air5_rotational_degrees, scales); }

}  // namespace gaslib

#endif  // GASLIB_DATA_AIR5_HPP
