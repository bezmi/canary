# SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
# Copyright (C) 2026 Shahzeb Imran and the canary contributors
"""Fit the electron-vibration relaxation time of N2 to Lee's rate coefficients.

Input: lee_1985_table2.txt, the excitation rate coefficients k_0,v (cm^3/s,
v = 1 to 10) of N2 by electron impact at 50 electron temperatures from 0.1 eV
to 5.0 eV (J.-H. Lee, AIAA Paper 85-1035, Table 2).

Relaxation time (Lee, Eq. (19)), with the integral over v taken as the sum
over the tabulated levels:

    p_e tau = k_B T_e / ((1 - exp(-theta_v / T_e))^2 (1/2) sum_v k_0,v v^2)

The script fits

    ln(p_e tau / (1 atm s)) = a + n ln(T_e / T_ref) + B u + C u^2 + D u^3,
    u = T_ref / T_e,

by least squares over the 50 temperatures, rounds the coefficients to 7
significant digits and prints them with the errors of the rounded fit. It
stops when B + 2 C u + 3 D u^2 has a zero for u > 0, because the fitted
relaxation time would then not rise monotonically towards low temperature.

Usage: python scripts/electron_vibration_relaxation/fit.py [--check]

--check exits with status 1 unless the printed values equal those of the
"relaxation_time" record of N2 in database/data/air11_kim_jo_exchange.json.
The script needs only the Python standard library.
"""

import json
import math
import sys
from pathlib import Path

# Exact SI values (2019 redefinition).
BOLTZMANN = 1.380649e-23  # J/K
ELEMENTARY_CHARGE = 1.602176634e-19  # C
ATMOSPHERE = 101325.0  # Pa

# Characteristic vibrational temperature of N2 in K used in Eq. (19).
THETA_V = 3395.0
REFERENCE_TEMPERATURE = 1.0e4  # K
SIGNIFICANT_DIGITS = 7

HERE = Path(__file__).resolve().parent


def read_table(path):
    """Rows (T_e in eV, k_0,1 ... k_0,10 in cm^3/s)."""
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        values = [float(word) for word in line.split()]
        if len(values) != 11:
            raise SystemExit(f"expected 11 numbers per row, got {len(values)}: {line}")
        rows.append(values)
    if len(rows) != 50:
        raise SystemExit(f"expected 50 rows, got {len(rows)}")
    for index, row in enumerate(rows):
        if abs(row[0] - 0.1 * (index + 1)) > 1.0e-9:
            raise SystemExit(f"row {index + 1}: unexpected temperature {row[0]} eV")
        if any(row[v + 1] >= row[v] for v in range(1, 10)):
            raise SystemExit(f"row {index + 1}: k_0,v does not decrease with v")
    return rows


def relaxation_points(rows):
    """Pairs (T_e in K, ln(p_e tau / (1 atm s))) of Eq. (19)."""
    points = []
    for row in rows:
        temperature = row[0] * ELEMENTARY_CHARGE / BOLTZMANN
        # k in cm^3/s to m^3/s.
        moment = sum(row[v] * 1.0e-6 * v * v for v in range(1, 11))
        factor = (1.0 - math.exp(-THETA_V / temperature)) ** 2
        p_tau = BOLTZMANN * temperature / (factor * 0.5 * moment)  # Pa s
        points.append((temperature, math.log(p_tau / ATMOSPHERE)))
    return points


def basis(temperature):
    """The five functions multiplied by a, n, B, C, D."""
    u = REFERENCE_TEMPERATURE / temperature
    return [1.0, -math.log(u), u, u * u, u * u * u]


def solve(matrix, rhs):
    """Gaussian elimination with partial pivoting."""
    n = len(rhs)
    a = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda r: abs(a[r][column]))
        a[column], a[pivot] = a[pivot], a[column]
        for r in range(column + 1, n):
            ratio = a[r][column] / a[column][column]
            for c in range(column, n + 1):
                a[r][c] -= ratio * a[column][c]
    x = [0.0] * n
    for r in reversed(range(n)):
        x[r] = (a[r][n] - sum(a[r][c] * x[c] for c in range(r + 1, n))) / a[r][r]
    return x


def fit(points):
    """Least-squares values of a, n, B, C, D."""
    rows = [basis(t) for t, _ in points]
    values = [value for _, value in points]
    n = len(rows[0])
    normal = [[math.fsum(r[i] * r[j] for r in rows) for j in range(n)] for i in range(n)]
    rhs = [math.fsum(r[i] * v for r, v in zip(rows, values)) for i in range(n)]
    return solve(normal, rhs)


def evaluate(coefficients, temperature):
    return math.fsum(c * b for c, b in zip(coefficients, basis(temperature)))


def main():
    rows = read_table(HERE / "lee_1985_table2.txt")
    points = relaxation_points(rows)
    coefficients = [float(f"{c:.{SIGNIFICANT_DIGITS}g}") for c in fit(points)]
    constant, exponent, b, c, d = coefficients
    # B + 2 C u + 3 D u^2 > 0 for all u > 0.
    if not (b > 0.0 and d > 0.0 and (c >= 0.0 or c * c < 3.0 * b * d)):
        raise SystemExit("the fitted relaxation time is not monotone in 1 / T_e")
    errors = [evaluate(coefficients, t) - value for t, value in points]
    worst = max(range(len(errors)), key=lambda i: abs(errors[i]))
    print(f"reference_temperature {REFERENCE_TEMPERATURE!r}")
    print(f"constant {constant!r}")
    print(f"exponent {exponent!r}")
    print(f"inverse_coefficients {[b, c, d]!r}")
    print(
        f"largest error of ln(p_e tau): {errors[worst]:+.5f} at {rows[worst][0]} eV, "
        f"{100.0 * (math.exp(abs(errors[worst])) - 1.0):.2f} % of tau"
    )
    if "--check" in sys.argv[1:]:
        path = HERE.parents[1] / "database/data/air11_kim_jo_exchange.json"
        records = json.loads(path.read_text())["vibration_electron"]
        stored = next(r for r in records if r["species"] == "N2")["relaxation_time"]
        expected = {
            "reference_temperature": REFERENCE_TEMPERATURE,
            "constant": constant,
            "exponent": exponent,
            "inverse_coefficients": [b, c, d],
        }
        for key, value in expected.items():
            if stored[key] != value:
                raise SystemExit(f"{path.name}: {key} is {stored[key]!r}, the fit gives {value!r}")
        print(f"{path.name} holds this fit")


if __name__ == "__main__":
    main()
