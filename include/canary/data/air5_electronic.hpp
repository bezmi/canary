#ifndef CANARY_DATA_AIR5_ELECTRONIC_HPP
#define CANARY_DATA_AIR5_ELECTRONIC_HPP

// Generated from database/data/air5.json by
// tools/species_database/species_to_header.py; do not edit.
// Reference: GDTk species database (gdtk-uq/gdtk commit bf559c3f,
// src/gas/species-database); fits: NASA CEA2 thermo.inp (McBride, Zehe and
// Gordon, NASA/TP-2002-211556) for N2, O2, NO and Johnston et al., AIAA Paper
// 2015-3110 for N and O; transport: NASA CEA trans.inp fits as in the GDTk
// database, Lennard-Jones parameters from GRI-Mech 3.0 (Svehla 1962 for N);
// electronic levels: NIST, as in the GDTk database

#include <array>
#include <canary/models/multi_temperature.hpp>

namespace canary
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
            72231.82201076281,
            85778.94223181416,
            86050.36748975527,
            95351.29745107313,
            98056.47134838725,
            99682.86473071769,
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
            11392.379193763898,
            18984.804776352157,
            47562.078872710794,
            49911.169880611465,
            50929.53615450875,
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
      .levels = 8,
      .degeneracies =
        {
          {
            5,
            3,
            1,
            5,
            1,
            5,
            3,
            15,
          },
        },
      .temperatures =
        {
          {
            0.0,
            227.70802251816005,
            326.5692593252104,
            22830.312941023323,
            48619.98705153252,
            106135.98045508568,
            110490.8386548233,
            124639.93055610423,
          },
        },
    },
    {
      .name = "NO",
      .levels = 7,
      .degeneracies =
        {
          {
            2,
            2,
            8,
            2,
            4,
            4,
            4,
            0,
          },
        },
      .temperatures =
        {
          {
            0.0,
            172.39424546252133,
            55306.583171251215,
            63256.8325632747,
            66101.15057241223,
            70039.6583968915,
            74997.68351677005,
            0.0,
          },
        },
    },
  },
};

}  // namespace canary

#endif  // CANARY_DATA_AIR5_ELECTRONIC_HPP
