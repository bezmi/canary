#ifndef GASLIB_DATA_AIR5_ELECTRONIC_HPP
#define GASLIB_DATA_AIR5_ELECTRONIC_HPP

// Generated from database/data/air5.json by
// tools/species_database/species_to_header.py; do not edit.
// Reference: GDTk species database (gdtk-uq/gdtk commit bf559c3f,
// src/gas/species-database); fits: NASA CEA2 thermo.inp (McBride, Zehe and
// Gordon, NASA/TP-2002-211556) for N2, O2, NO and Johnston et al., AIAA Paper
// 2015-3110 for N and O; transport: NASA CEA trans.inp fits as in the GDTk
// database, Lennard-Jones parameters from GRI-Mech 3.0 (Svehla 1962 for N);
// electronic levels: NIST, as in the GDTk database

#include <array>
#include <gaslib/models/multi_temperature.hpp>

namespace gaslib
{

// Electronic levels of N2, O2, N, O, NO, in this order: degeneracies
// and level temperatures E_k / k_B in K; entries past a species' level count
// are zero.
inline constexpr std::array<ElectronicLevelsInput<8>, 5> air5_electronic_levels {
  {
    {
      .name = "N2",
      .levels = 7,
      .degeneracies =
        {
          {
            1,
            3,
            6,
            6,
            3,
            1,
            2,
            0,
          },
        },
      .temperatures =
        {
          {
            0.0,
            72231.86517406914,
            85778.56814982601,
            86050.36748975527,
            95351.34061437946,
            98056.47134838725,
            99682.8647307177,
            0.0,
          },
        },
    },
    {
      .name = "O2",
      .levels = 6,
      .degeneracies =
        {
          {
            3,
            2,
            1,
            1,
            6,
            3,
            0,
            0,
          },
        },
      .temperatures =
        {
          {
            0.0,
            11392.292867151247,
            18984.804776352157,
            47562.078872710794,
            49911.601513674716,
            50927.80962225574,
            0.0,
            0.0,
          },
        },
    },
    {
      .name = "N",
      .levels = 7,
      .degeneracies =
        {
          {
            4,
            10,
            6,
            12,
            6,
            12,
            2,
            0,
          },
        },
      .temperatures =
        {
          {
            0.0,
            27664.729861801763,
            41493.145350173894,
            119901.33435291682,
            124012.1932595523,
            126802.93121777158,
            134642.97030097825,
            0.0,
          },
        },
    },
    {
      .name = "O",
      .levels = 6,
      .degeneracies =
        {
          {
            9,
            5,
            1,
            5,
            3,
            15,
            0,
            0,
          },
        },
      .temperatures =
        {
          {
            0.0,
            22830.31006346957,
            48619.982735201884,
            106135.9804550857,
            110490.84153237705,
            124639.93161120726,
            0.0,
            0.0,
          },
        },
    },
    {
      .name = "NO",
      .levels = 6,
      .degeneracies =
        {
          {
            4,
            8,
            2,
            4,
            4,
            4,
            0,
            0,
          },
        },
      .temperatures =
        {
          {
            0.0,
            55306.583171251215,
            63256.8325632747,
            66086.33117057393,
            70039.6583968915,
            75069.190727582,
            0.0,
            0.0,
          },
        },
    },
  },
};

}  // namespace gaslib

#endif  // GASLIB_DATA_AIR5_ELECTRONIC_HPP
