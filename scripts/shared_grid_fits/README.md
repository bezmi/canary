# Atomic thermodynamic fit generation

Generate shared-grid polynomial fits for N, O, N+ and O+ from the electronic
levels in `levels.json`. The calculation uses ideal translation, Boltzmann
electronic populations and a fixed finite level list. It does not impose
chemical equilibrium or charge neutrality.

## Run

Python 3.9 or later is required. Install the dependencies in a virtual
environment, then run from the repository root:

```sh
python -m venv /tmp/canary-fit-venv
/tmp/canary-fit-venv/bin/python -m pip install -r scripts/shared_grid_fits/requirements.txt
/tmp/canary-fit-venv/bin/python scripts/shared_grid_fits/verify_levels.py --output-dir /tmp/canary-atomic-fits
/tmp/canary-fit-venv/bin/python scripts/shared_grid_fits/fit.py --output-dir /tmp/canary-atomic-fits
/tmp/canary-fit-venv/bin/python scripts/shared_grid_fits/verify_fits.py --output-dir /tmp/canary-atomic-fits
```

The scripts resolve their input data relative to their own directory. They
need no network access after dependency installation and no files from a
developer's plan directory. The output directory is created when necessary;
rerunning replaces generated files with the same names. Run the checks with
ordinary Python, without `-O`, because verification uses assertions.

The generator uses mpmath; the independent level verification also uses NumPy
and SciPy. The scripts run sequentially and do not invoke a CFD simulation.

## Representation and outputs

All four species share 44 intervals over 1000–50000 K. The grid has four
intervals from 1000 to 1200 K and eight intervals in each of 1200–5500,
5500–6500, 6500–19500, 19500–20500 and 20500–50000 K.

At each endpoint the generator stores internal energy and heat capacity:

```text
energy = (u(T) - u(1000 K)) / R
capacity = cv(T) / R
```

Within each interval, degree-seven energy consists of cubic Hermite
interpolation plus `x^2 (1-x)^2 q(2x-1)`, where `x` is the local interval
coordinate and `q` is cubic. Four interior Chebyshev samples determine `q`.
The energy and its first derivative share endpoint data. Degree-six capacity
is the derivative of this energy. Positive Bernstein coefficients verify
positive fitted `cv` throughout each interval in the high-precision fit.

| Output | Contents |
|---|---|
| `fits.json` | Grid, endpoint values, interior coefficients and entropy increments |
| `fitted-data.hpp` | C++ coefficient arrays in species order N, O, N+, O+ |
| `fit-summary.json` | Sampled approximation errors, positivity and endpoint residuals |
| `fit-comparison.csv` | 65 high-precision reference comparisons per interval and species |
| `runtime-reference.txt` | Reference values for float/double temperatures and temperature scales 1 and 1000 |
| `level-verification.json` | Independent level-sum quadrature and precision checks |
| `fit-verification.json` | Independent entropy quadrature, joins and derivative checks |

Calculations use 70 decimal digits; entropy verification uses 80. The entropy
integral is evaluated offline and is not a device evaluator. Sampled errors
are not continuous worst-case error bounds. Positivity before coefficient
rounding does not establish positivity for every runtime representation.

The output does not supply absolute formation enthalpies or standard entropy
constants, a model below 1000 K, or an extrapolation policy. These require
consistent treatment when integrating the fits into a reacting-gas model.

## Input provenance

The level energies and degeneracies are retained exactly as decimal inputs.
The universal gas constant is 8.31446261815324 J/(mol K); conversions from eV
and inverse centimetres use exact SI constants.

| Species | Entries | Source |
|---|---:|---|
| N | 35 | Johnston's 2006 dissertation, Table 3.1 |
| O | 32 | Johnston's 2006 dissertation, Table 3.2 |
| N+ | 26 | Johnston, Samareh and Brandis (2015), Table 7 |
| O+ | 275 | NIST ASD O II numerical levels, based on Martin, Kaufman and Musgrove (1993) |

- [Johnston's dissertation](https://vtechworks.lib.vt.edu/handle/10919/29769).
- [Johnston, Samareh and Brandis, AIAA 2015-3110](https://doi.org/10.2514/6.2015-3110).
- [NIST Atomic Spectra Database](https://physics.nist.gov/asd).
- [Martin, Kaufman and Musgrove's O II compilation](https://www.nist.gov/system/files/documents/srd/jpcrd459.pdf).

N and O use the published grouped levels, including N entry 32's degeneracy
of 822. Exact identity with the unpublished inputs used for the later
Johnston polynomial fits is unverified. In particular, grouped O levels do
not resolve the ground-term fine structure at low temperature.

For O+, `nist-current-OII.tsv` retains the source export and `levels.json`
retains per-level flags. The calculation includes 16 tentative levels and
10 levels with an unspecified common offset, evaluated at the nominal offset
zero. It includes 28 levels above the first ionization limit. Five continuum
limit rows and seven states without numerical energies are excluded. No
density-dependent cutoff, ionization-potential correction or completion of
unlisted Rydberg levels is applied.

These input choices define the fitted model. Small polynomial approximation
errors do not quantify uncertainties in that physical model.
