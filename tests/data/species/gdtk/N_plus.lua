-- Test input in the layout of GDTk species-database files (bracket form of
-- the species name), written for gaslib's converter test. Values: NASA CEA
-- thermo.inp (Apache License 2.0, see ../NOTICE.md); electronic levels from
-- NIST (cm^-1).
db['N+'] = {}
db['N+'].charge = 1
db['N+'].M = {
   value = 0.0140061514,
   units = 'kg/mol',
}
db['N+'].electronic_levels = {
  Te = {
    value = {0.0, 70.07, 188.19, 22036.47, 47031.63, 67312.23},
    units = 'cm^(-1)',
  },
  g = {
    value = {1, 3, 5, 5, 1, 5},
  }
}
db['N+'].ceaThermoCoeffs = {
   nsegments = 4,
   T_break_points = {298.15, 1000.0, 6000.0, 20000.0, 50000.0},
   T_blend_ranges = {400.0, 1000.0, 1000.0},
   segment0 = {
      5237.07921, 2.29995832, 2.48237856, 2.73749076e-05, -3.13444758e-08, 1.85011133e-11, -4.44735098e-15, 225626.95, 5.0477146
   },
   segment1 = {
      351764.922, -1044.68568, 3.68498994, -0.000641320819, 1.66764737e-07, -1.82509101e-11, 7.37893962e-16, 232355.85, -3.4854045
   },
   segment2 = {
      20893183.1, -13621.3868, 5.54565345, -0.000268374643, 1.46878484e-08, -4.22829994e-13, 5.70569781e-18, 333626.095, -22.0294391
   },
   segment3 = {
      1483513370.0, -269100.02, 19.1340133, -0.00021338085, -1.51589694e-08, 5.96512853e-13, -5.46224064e-18, 2700086.12, -164.949269
   },
   reference = 'Johnston et al. AIAA Paper 2015-3110'
}
