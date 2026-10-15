"""Evaluate published electronic-level sets independently of polynomial fits.

Returned molar cp, h and standard-pressure s omit additive reference
constants. The fit generator uses h and s increments from 1000 K. Fixed
composition, ideal translation, Boltzmann electronic populations and fixed
level lists apply.
No density-dependent ionization correction is included.
"""
from pathlib import Path
import json
import mpmath as mp

ROOT = Path(__file__).resolve().parent
mp.mp.dps = 70
R = mp.mpf('8.31446261815324')
CM_OVER_K = mp.mpf('6.62607015e-34') * 299792458 * 100 / mp.mpf('1.380649e-23')
EV_OVER_K = mp.mpf('1.602176634e-19') / mp.mpf('1.380649e-23')
LEVELS = json.loads((ROOT / 'levels.json').read_text())


def electronic(name, temperature):
    t = mp.mpf(temperature)
    assert t > 0
    rows = LEVELS[name]['levels']
    theta = [mp.mpf(row['energy_cm_inverse']) * CM_OVER_K if 'energy_cm_inverse' in row
             else mp.mpf(row['energy_eV']) * EV_OVER_K for row in rows]
    weights = [row['degeneracy'] * mp.exp(-x / t) for row, x in zip(rows, theta)]
    z = mp.fsum(weights)
    mean = mp.fsum(w*x for w, x in zip(weights, theta)) / z
    variance = mp.fsum(w*(x-mean)**2 for w, x in zip(weights, theta)) / z
    return (R * (mp.mpf('2.5') + variance / t**2),
            R * (mp.mpf('2.5') * t + mean),
            R * (mp.mpf('2.5') * mp.log(t) + mp.log(z) + mean / t))
