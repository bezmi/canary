# Third-party data in this directory

- `thermo_inp_excerpt.inp`, and the coefficients in `gdtk/N2.lua` and
  `gdtk/N_plus.lua`: excerpts
  of `data/thermo.inp` (and, for the viscosity and conductivity fits in
  `gdtk/N2.lua`, `data/trans.inp`) from NASA CEA, https://github.com/nasa/cea, licensed
  under the Apache License, Version 2.0
  (https://www.apache.org/licenses/LICENSE-2.0). Data: McBride, Zehe and
  Gordon, "NASA Glenn Coefficients for Calculating Thermodynamic Properties
  of Individual Species", NASA/TP-2002-211556 (2002).
- `cantera_excerpt.yaml`: excerpt of `data/airNASA9.yaml` from Cantera,
  https://github.com/Cantera/cantera, licensed under the BSD 3-Clause
  License. Copyright (c) 2001-2026, Cantera Developers. Redistribution of
  this excerpt retains this notice; the license text is at
  https://github.com/Cantera/cantera/blob/main/License.txt.
- The Lennard-Jones parameters in `gdtk/N2.lua` are those of N2 in the
  GRI-Mech 3.0 transport data (Smith et al., http://combustion.berkeley.edu/gri-mech/).
- `gdtk_collision_integrals_sample.lua`: collision-integral curve fits of
  three pairs from Gupta, Yos, Thompson and Lee, NASA RP-1232 (1990), a
  publication of the United States government.
- The electronic levels in `gdtk/N2.lua` and `gdtk/N_plus.lua` are NIST
  atomic and molecular data (level energies in cm^-1).
