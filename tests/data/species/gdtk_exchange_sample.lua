-- Test input in the layout of GDTk energy-exchange files as written by
-- prep-kinetics, written for gaslib's converter test. Millikan-White a, b:
-- Park (1993), Table 1; dissociation energy and vibrational temperature of
-- N2 as GDTk's prep-kinetics derives them.

mechanism = {}

mechanism[1] = {
  type = 'V-T',
  p = 'N2', q = 'N',
  mode_p = 0,
  rate = 'Landau-Teller',
  relaxation_time = {model='ParkHTC2', sigma=1.0000e-21, submodel={model='Millikan-White', a=180.000, b=0.026200}}
}

mechanism[2] = {
  type = 'V-T',
  p = 'O2', q = 'O2',
  mode_p = 0,
  rate = 'Landau-Teller',
  relaxation_time = {model='Millikan-White', a=138.000, b=0.030000}
}

mechanism[3] = {
  type = 'C-V',
  p = 'N2',
  rate = 'Marrone-Treanor',
  mode_p = 0,
  reaction_index = 0,
  coupling_model = {model='ImpartialDissociation', D=945360.000000, Thetav=3393.440000}
}

mechanism[4] = {
  type = 'V-T',
  p = 'N2', q = 'N2',
  mode_p = 0,
  rate = 'Landau-Teller',
  relaxation_time = {model='KimHTC', sigma=1.7000e-23, exponent=0.5500, submodel={model='Millikan-White', a=221.530, b=0.029000}}
}

mechanism[5] = {
  type = 'V-T',
  p = 'NO', q = 'N',
  mode_p = 0,
  rate = 'Landau-Teller',
  relaxation_time = {model='ParkHTC', sigma=5.8000e-20, submodel={model='Millikan-White', a=49.500, b=0.042000}}
}

mechanism[6] = {
  type = 'E-T',
  p = 'e-', q = 'N2',
  rate = 'ElectronExchange',
  mode_p = 1,
  exchange_cross_section = {c = -1e-28, b = 5.5e-24, a = 7.5e-20, type = 'GnoffoNeutral'}
}

mechanism[7] = {
  type = 'E-T',
  p = 'e-', q = 'N+',
  rate = 'ElectronExchange',
  mode_p = 1,
  exchange_cross_section = {type = 'Coulomb'}
}

mechanism[8] = {
  type = 'E-V',
  p = 'N2', q = 'e-',
  mode_p = 0,
  rate = 'Landau-Teller',
  relaxation_time = {model='Candler', a_low=7.5000e+00, b_low=-5.7000e+01, c_low=9.8700e+01, a_high=2.3600e+00, b_high=-1.7900e+01, c_high=2.4350e+01 }
}

mechanism[9] = {
  type = 'C-V',
  p = 'N2+',
  rate = 'Marrone-Treanor',
  mode_p = 0,
  reaction_index = 0,
  coupling_model = {model='ImpartialDissociation', D=-564148.420000, Thetav=3393.440000}
}
