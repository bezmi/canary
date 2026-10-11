-- Test input in the layout of GDTk collision-integral data files, written
-- for gaslib's converter test. Values: Gupta, Yos, Thompson and Lee (1990),
-- NASA RP-1232.
cis = {}
cis['N2:N2'] = {
  pi_Omega_11 = {A=  0.0000, B= -0.0112, C= -0.1182, D=  4.8464},
  pi_Omega_22 = {A=  0.0000, B= -0.0203, C=  0.0683, D=  4.0900}
}
cis['O2:N2'] = {
  pi_Omega_11 = {A=  0.0000, B= -0.0465, C=  0.5729, D=  1.6185},
  pi_Omega_22 = {A=  0.0000, B= -0.0558, C=  0.7590, D=  0.8955}
}
cis['O2:O2'] = {
  pi_Omega_11 = {A=  0.0000, B= -0.0410, C=  0.4977, D=  1.8302},
  pi_Omega_22 = {A=  0.0000, B= -0.0485, C=  0.6475, D=  1.2607}
}
