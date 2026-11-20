# SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
# Copyright (C) 2026 Shahzeb Imran and the canary contributors
"""Independently integrate fitted heat capacity and check entropy joins."""
import argparse
from pathlib import Path
import json
import mpmath as mp
import fit

mp.mp.dps = 80
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', type=Path, required=True,
                    help='Directory containing the generated fits.json.')
output_dir = parser.parse_args().output_dir
data = json.loads((output_dir/'fits.json').read_text())
knots = list(map(mp.mpf, data['knots_K']))
out = {}
for name, table in data['species'].items():
    max_integral = mp.mpf(0)
    max_join = mp.mpf(0)
    max_derivative = mp.mpf(0)
    for i, p in enumerate(table['pieces']):
        lo, hi = knots[i:i+2]
        width = hi-lo
        a = list(map(mp.mpf, p['energy_power_u']))
        cp = [(k+1)*a[k+1]/width for k in range(7)]
        cp[0] += 1
        for u in (mp.mpf('.25'), mp.mpf('.5'), mp.mpf(1)):
            exact = mp.quad(lambda x: fit.horner(cp, x)/(x+lo/width), [0,u])
            analytic = fit.entropy_integral(cp,lo,width,u)
            max_integral = max(max_integral, abs(exact-analytic)*fit.ref.R)
            derivative = mp.diff(lambda x: fit.entropy_integral(cp,lo,width,x),u)
            max_derivative = max(max_derivative,abs(derivative*(u+lo/width)/fit.horner(cp,u)-1))
        if i+1<len(table['pieces']):
            right = mp.mpf(p['entropy_left_over_R'])+fit.entropy_integral(cp,lo,width,mp.mpf(1))
            next_left = mp.mpf(table['pieces'][i+1]['entropy_left_over_R'])
            max_join = max(max_join,abs(right-next_left)*fit.ref.R)
    assert max_integral < mp.mpf('1e-60')
    assert max_join < mp.mpf('1e-55')
    assert max_derivative < mp.mpf('1e-60')
    out[name] = dict(max_entropy_quadrature_error_J_mol_K=str(max_integral),
                     max_entropy_join_J_mol_K=str(max_join),
                     max_entropy_identity_relative_residual=str(max_derivative))
(output_dir/'fit-verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
