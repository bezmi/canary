#ifndef CANARY_DATA_AIR5_TRANSPORT_HPP
#define CANARY_DATA_AIR5_TRANSPORT_HPP

// Generated from database/data/air5.json and
// database/data/air5_collision_integrals.json by
// tools/species_database/transport_to_header.py; do not edit.
// References: Molar masses and NASA9 fits of N2, O2 and NO: NASA CEA
// thermo.inp (McBride, Zehe and Gordon, NASA/TP-2002-211556); NASA9 fits of N
// and O: C. O. Johnston, J. Samareh and A. M. Brandis, AIAA Paper 2015-3110,
// Tables 1 and 4; viscosity and conductivity fits: NASA CEA trans.inp;
// Lennard-Jones parameters: GRI-Mech 3.0 transport file; electronic levels as
// cited per species; Gupta, Yos, Thompson and Lee, A review of reaction rates
// and thermodynamic and transport properties for an 11-species air model for
// chemical and thermal nonequilibrium calculations to 30000 K, NASA RP-1232,
// 1990, Tables VIII and IX; ln(1e-20) is added to D for cross-sections in m^2

#include <array>
#include <canary/transport/mixture.hpp>

namespace canary
{

// Transport data of N2, O2, N, O, NO, in this order, in SI units:
// CEA viscosity and conductivity fits (ln(x / SI unit) = A ln T + B / T +
// C / T^2 + D) and Lennard-Jones sigma in m and epsilon / k_B in K.
inline constexpr std::array<SpeciesTransportInput<3>, 5> air5_species_transport {
  {
    {
      .name = "N2",
      .viscosity =
        {
          .segments = 3,
          .breaks =
            {
              {
                200.0,
                1000.0,
                5000.0,
                15000.0,
              },
            },
          .coefficients =
            {
              {
                {
                  {
                    0.62526577,
                    -31.779652,
                    -1640.7983,
                    -14.372596450958321,
                  },
                },
                {
                  {
                    0.87395209,
                    561.52222,
                    -173948.09,
                    -16.51145523095832,
                  },
                },
                {
                  {
                    0.88503551,
                    909.02171,
                    -731290.61,
                    -16.65313403095832,
                  },
                },
              },
            },
        },
      .conductivity =
        {
          .segments = 3,
          .breaks =
            {
              {
                200.0,
                1000.0,
                5000.0,
                15000.0,
              },
            },
          .coefficients =
            {
              {
                {
                  {
                    0.85439436,
                    105.73224,
                    -12347.848,
                    -8.732409091976182,
                  },
                },
                {
                  {
                    0.88407146,
                    133.57293,
                    -11429.64,
                    -8.966170181976182,
                  },
                },
                {
                  {
                    2.4176185,
                    8047.7749,
                    3105580.2,
                    -23.72810137197618,
                  },
                },
              },
            },
        },
      .lennard_jones =
        {
          .sigma = 3.621e-10,
          .epsilon = 97.53,
        },
    },
    {
      .name = "O2",
      .viscosity =
        {
          .segments = 3,
          .breaks =
            {
              {
                200.0,
                1000.0,
                5000.0,
                15000.0,
              },
            },
          .coefficients =
            {
              {
                {
                  {
                    0.6091618,
                    -52.244847,
                    -599.74009,
                    -14.07701555095832,
                  },
                },
                {
                  {
                    0.72216486,
                    175.50839,
                    -57974.816,
                    -15.027991250958321,
                  },
                },
                {
                  {
                    0.73981127,
                    391.94906,
                    -378331.68,
                    -15.20877785095832,
                  },
                },
              },
            },
        },
      .conductivity =
        {
          .segments = 3,
          .breaks =
            {
              {
                200.0,
                1000.0,
                5000.0,
                15000.0,
              },
            },
          .coefficients =
            {
              {
                {
                  {
                    0.77229167,
                    6.846321,
                    -5893.3377,
                    -7.9893038719761815,
                  },
                },
                {
                  {
                    0.90917351,
                    291.24182,
                    -79650.171,
                    -9.145488740976182,
                  },
                },
                {
                  {
                    -1.1218262,
                    -19286.378,
                    23295011.0,
                    11.131702628023818,
                  },
                },
              },
            },
        },
      .lennard_jones =
        {
          .sigma = 3.458e-10,
          .epsilon = 107.4,
        },
    },
    {
      .name = "N",
      .viscosity =
        {
          .segments = 2,
          .breaks =
            {
              {
                1000.0,
                5000.0,
                15000.0,
                0.0,
              },
            },
          .coefficients =
            {
              {
                {
                  {
                    0.83724737,
                    439.9715,
                    -174507.53,
                    -16.01443876095832,
                  },
                },
                {
                  {
                    0.89986588,
                    1411.2801,
                    -1820047.8,
                    -16.67621281095832,
                  },
                },
                {
                  {
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                  },
                },
              },
            },
        },
      .conductivity =
        {
          .segments = 2,
          .breaks =
            {
              {
                1000.0,
                5000.0,
                15000.0,
                0.0,
              },
            },
          .coefficients =
            {
              {
                {
                  {
                    0.83771661,
                    442.4327,
                    -175784.46,
                    -8.310911221976182,
                  },
                },
                {
                  {
                    0.9000171,
                    1414.1175,
                    -1826240.3,
                    -8.969855241976182,
                  },
                },
                {
                  {
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                  },
                },
              },
            },
        },
      .lennard_jones =
        {
          .sigma = 3.298e-10,
          .epsilon = 71.4,
        },
    },
    {
      .name = "O",
      .viscosity =
        {
          .segments = 2,
          .breaks =
            {
              {
                1000.0,
                5000.0,
                15000.0,
                0.0,
              },
            },
          .coefficients =
            {
              {
                {
                  {
                    0.77269241,
                    83.842977,
                    -58502.098,
                    -15.267087380958321,
                  },
                },
                {
                  {
                    0.87669586,
                    1015.842,
                    -1088456.6,
                    -16.29810642095832,
                  },
                },
                {
                  {
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                  },
                },
              },
            },
        },
      .conductivity =
        {
          .segments = 2,
          .breaks =
            {
              {
                1000.0,
                5000.0,
                15000.0,
                0.0,
              },
            },
          .coefficients =
            {
              {
                {
                  {
                    0.77271664,
                    83.9891,
                    -58580.966,
                    -7.692350371976182,
                  },
                },
                {
                  {
                    0.87676666,
                    1017.0744,
                    -1090669.0,
                    -8.723898051976182,
                  },
                },
                {
                  {
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                  },
                },
              },
            },
        },
      .lennard_jones =
        {
          .sigma = 2.75e-10,
          .epsilon = 80.0,
        },
    },
    {
      .name = "NO",
      .viscosity =
        {
          .segments = 3,
          .breaks =
            {
              {
                200.0,
                1000.0,
                5000.0,
                15000.0,
              },
            },
          .coefficients =
            {
              {
                {
                  {
                    0.60262029,
                    -62.017783,
                    -139.54524,
                    -14.09126245095832,
                  },
                },
                {
                  {
                    0.7800905,
                    304.86891,
                    -94847.722,
                    -15.58936184095832,
                  },
                },
                {
                  {
                    0.80580582,
                    624.27878,
                    -578792.1,
                    -15.85293115095832,
                  },
                },
              },
            },
        },
      .conductivity =
        {
          .segments = 3,
          .breaks =
            {
              {
                200.0,
                1000.0,
                5000.0,
                15000.0,
              },
            },
          .coefficients =
            {
              {
                {
                  {
                    0.95028758,
                    76.667058,
                    -9989.4764,
                    -9.216618043676181,
                  },
                },
                {
                  {
                    0.86215238,
                    445.68223,
                    -238564.66,
                    -8.748241611976182,
                  },
                },
                {
                  {
                    -1.0377865,
                    -34486.864,
                    67451187.0,
                    11.718408628023818,
                  },
                },
              },
            },
        },
      .lennard_jones =
        {
          .sigma = 3.621e-10,
          .epsilon = 97.53,
        },
    },
  },
};

// Collision integrals of the 15 pairs of these species, in the order of
// air5_collision_integrals.json: ln(pi Omega / m^2) = A (ln T)^3 + B (ln T)^2 + C ln T + D.
inline constexpr std::array<CollisionPairInput, 15> air5_collision_integrals {
  {
    {
      .first = "N2",
      .second = "N2",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0112,
          .C = -0.1182,
          .D = -41.20530185988091,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0203,
          .C = 0.0683,
          .D = -41.96170185988092,
        },
    },
    {
      .first = "O2",
      .second = "N2",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0465,
          .C = 0.5729,
          .D = -44.43320185988092,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0558,
          .C = 0.759,
          .D = -45.156201859880916,
        },
    },
    {
      .first = "O2",
      .second = "O2",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.041,
          .C = 0.4977,
          .D = -44.22150185988092,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0485,
          .C = 0.6475,
          .D = -44.791001859880915,
        },
    },
    {
      .first = "N",
      .second = "N2",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0194,
          .C = 0.0119,
          .D = -41.946201859880915,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.019,
          .C = 0.0239,
          .D = -41.87350185988092,
        },
    },
    {
      .first = "N",
      .second = "O2",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0179,
          .C = 0.0152,
          .D = -42.05210185988091,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0203,
          .C = 0.073,
          .D = -42.169901859880916,
        },
    },
    {
      .first = "N",
      .second = "N",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0033,
          .C = -0.0572,
          .D = -41.00650185988091,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0118,
          .C = -0.096,
          .D = -41.72650185988091,
        },
    },
    {
      .first = "O",
      .second = "N2",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0139,
          .C = -0.0825,
          .D = -41.473201859880916,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0169,
          .C = -0.0143,
          .D = -41.632201859880915,
        },
    },
    {
      .first = "O",
      .second = "O2",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0226,
          .C = 0.13,
          .D = -42.71540185988091,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0247,
          .C = 0.1783,
          .D = -42.800001859880915,
        },
    },
    {
      .first = "O",
      .second = "N",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = 0.0048,
          .C = -0.4195,
          .D = -40.274301859880914,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = 0.0065,
          .C = -0.4467,
          .D = -40.009101859880914,
        },
    },
    {
      .first = "O",
      .second = "O",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0034,
          .C = -0.0572,
          .D = -41.061601859880916,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0207,
          .C = 0.078,
          .D = -42.48590185988091,
        },
    },
    {
      .first = "NO",
      .second = "N2",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0291,
          .C = 0.2324,
          .D = -42.84350185988092,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0385,
          .C = 0.4226,
          .D = -43.60100185988092,
        },
    },
    {
      .first = "NO",
      .second = "O2",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0438,
          .C = 0.5352,
          .D = -44.32650185988091,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0522,
          .C = 0.7045,
          .D = -44.977901859880916,
        },
    },
    {
      .first = "NO",
      .second = "N",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0185,
          .C = 0.0118,
          .D = -41.99270185988092,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0196,
          .C = 0.0478,
          .D = -42.019601859880915,
        },
    },
    {
      .first = "NO",
      .second = "O",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0179,
          .C = 0.0152,
          .D = -42.05210185988091,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0203,
          .C = 0.073,
          .D = -42.169901859880916,
        },
    },
    {
      .first = "NO",
      .second = "NO",
      .pi_omega_11 =
        {
          .A = 0.0,
          .B = -0.0364,
          .C = 0.3825,
          .D = -43.57990185988091,
        },
      .pi_omega_22 =
        {
          .A = 0.0,
          .B = -0.0453,
          .C = 0.5624,
          .D = -44.284801859880915,
        },
    },
  },
};

}  // namespace canary

#endif  // CANARY_DATA_AIR5_TRANSPORT_HPP
