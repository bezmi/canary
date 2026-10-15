"""Check partition sums against independent binary64 quadrature and higher precision."""
import argparse
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
import calculate as ref

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir', type=Path, required=True,
                    help='Directory for the independent level verification report.')
output_dir = parser.parse_args().output_dir
output_dir.mkdir(parents=True, exist_ok=True)
results = {}
for name, data in ref.LEVELS.items():
    theta = np.array([float(x['energy_cm_inverse'])*6.62607015e-34*299792458*100/1.380649e-23
                      if 'energy_cm_inverse' in x else float(x['energy_eV'])*1.602176634e-19/1.380649e-23
                      for x in data['levels']])
    weights = np.array([x['degeneracy'] for x in data['levels']])

    def capacity(t):
        p = weights * np.exp(-theta/t)
        p /= p.sum()
        return 8.31446261815324 * (2.5 + (p @ (theta*theta) - (p @ theta)**2)/t**2)

    checks = []
    for t in (6000, 20000, 47400, 50000):
        cp, h, s = map(float, ref.electronic(name,t))
        anchor = list(map(float,ref.electronic(name,1000)))
        ih, eh = quad(capacity,1000,t,epsabs=1e-7,epsrel=1e-11)
        iss, es = quad(lambda x: capacity(x)/x,1000,t,epsabs=1e-10,epsrel=1e-11)
        checks.append(dict(temperature_K=t, cp_difference_J_mol_K=capacity(t)-cp,
                           integrated_h_difference_J_mol=ih-(h-anchor[1]),
                           integrated_s_difference_J_mol_K=iss-(s-anchor[2]),
                           quadrature_h_error_estimate_J_mol=eh,
                           quadrature_s_error_estimate_J_mol_K=es))
        assert abs(capacity(t)-cp) < 1e-10
        assert abs(ih-(h-anchor[1])) < 1e-7
        assert abs(iss-(s-anchor[2])) < 1e-10
    precision = []
    for t in (1000,6000,20000,47400,50000):
        with ref.mp.workdps(40):
            a = ref.electronic(name,t)
        with ref.mp.workdps(90):
            b = ref.electronic(name,t)
            residual = max(abs(x/y-1) for x,y in zip(a,b))
            assert residual < ref.mp.mpf('1e-38')
            precision.append(str(residual))
    results[name] = dict(quadrature=checks, precision_40_vs_90_relative=precision)

(output_dir/'level-verification.json').write_text(json.dumps(results, indent=2)+'\n')
print(f'Passed: {len(ref.LEVELS)*4} quadrature cases and {len(ref.LEVELS)*5} precision comparisons.')
