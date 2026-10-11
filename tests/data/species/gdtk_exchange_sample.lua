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
