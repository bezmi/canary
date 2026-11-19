# Electron-vibration relaxation time of N2

`fit.py` computes the polynomial that `database/data/air11_kim_jo_exchange.json`
stores as the `relaxation_time` of the `vibration_electron` record of N2.

## Source

J.-H. Lee, "Electron-impact vibrational excitation rates in the flow field of
aeroassisted orbital transfer vehicles", AIAA Paper 85-1035, 1985.

- `lee_1985_table2.txt` is Table 2 of the paper: the excitation rate
  coefficients k_0,v of N2 by electron impact in cm^3/s for v = 1 to 10, at 50
  electron temperatures from 0.1 eV to 5.0 eV in steps of 0.1 eV. The numbers
  were transcribed from the printed table.
- Eq. (19) of the paper gives the relaxation time of the Landau-Teller rate of
  its Eq. (16):

  ```text
  p_e tau = k_B T_e / ((1 - exp(-theta_v / T_e))^2 (1/2) sum_v k_0,v v^2)
  ```

  with theta_v = 3395 K and the sum over the ten tabulated levels.

## Fit

```text
ln(p_e tau / (1 atm s)) = a + n ln(T_e / T_ref) + B u + C u^2 + D u^3
u = T_ref / T_e,  T_ref = 10000 K
```

The values of a, n, B, C and D are the least-squares solution over the 50
temperatures, rounded to 7 significant digits. The rounded fit differs from
Eq. (19) by at most 0.76 % of tau at the tabulated temperatures, 1160.45 K
(0.1 eV) to 58022.59 K (5.0 eV).

Outside the table the form has these limits:

- high temperature: p_e tau proportional to T_e^n with n = 4.475. For rate
  coefficients proportional to T_e^(-3/2) far above the resonance energy,
  Eq. (19) gives T_e^4.5;
- low temperature: B + 2 C u + 3 D u^2 is positive for every u > 0, so
  p_e tau rises monotonically as T_e falls. The script stops when a fit does
  not have this property.

The paper gives no values outside the table to compare with.

## Run

Python 3.9 or later; only the standard library is used. From the repository
root:

```sh
python scripts/electron_vibration_relaxation/fit.py
python scripts/electron_vibration_relaxation/fit.py --check
```

The first command prints the reference temperature, the constant a, the
exponent n, the coefficients B, C, D and the largest error. With `--check` the script also
compares them with the values in `database/data/air11_kim_jo_exchange.json`
and exits with status 1 when they differ.

The script reads the table and checks its layout: 50 rows, temperatures 0.1
to 5.0 eV, and k_0,v decreasing with v in every row. It writes no file.
