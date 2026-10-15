"""Generate atomic thermodynamic fits with shared energy and capacity endpoints.

The fixed 44-interval grid covers 1000--50000 K for N, O, N+ and O+. Energy is (u(T)-u(1000 K))/R. Degree-seven
energy gives degree-six capacity. Entropy is the analytic integral of fitted
cp/T, evaluated offline at 70 digits; this is not a device entropy evaluator.
"""
import argparse
from pathlib import Path
import json
import csv
import struct
from math import comb
from functools import lru_cache
import mpmath as mp
import calculate as ref

mp.mp.dps = 70
NAMES = tuple(ref.LEVELS)


def convolve(a,b):
    c = [mp.mpf(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j] += x*y
    return c


def horner(a,x):
    return mp.polyval(list(reversed(a)),x)


def power_energy(left,right,width,q):
    e0,c0 = left
    e1,c1 = right
    # Cubic Hermite endpoint constraints plus u^2(1-u)^2 q(2u-1).
    a = [e0,width*c0,3*(e1-e0)-width*(2*c0+c1),-2*(e1-e0)+width*(c0+c1)]
    qu = [mp.mpf(0)]*4
    for k in range(4):
        for j in range(k+1): qu[j] += q[k]*comb(k,j)*2**j*(-1)**(k-j)
    bubble = convolve([0,0,1,-2,1],qu)
    return [bubble[i]+(a[i] if i<len(a) else 0) for i in range(8)]


def entropy_integral(cp,lo,width,u):
    # Polynomial division by u+lo/width leaves a log1p term.
    ratio = lo/width
    quotient = [mp.mpf(0)]*6
    remainder = cp[:]
    for k in range(6,0,-1):
        quotient[k-1] = remainder[k]
        remainder[k-1] -= ratio*remainder[k]
    return mp.fsum(quotient[k]*u**(k+1)/(k+1) for k in range(6))+remainder[0]*mp.log1p(u/ratio)


def rounded(x,bits):
    f = float(x)
    return mp.mpf(struct.unpack('f',struct.pack('f',f))[0] if bits==32 else f)


def main(output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    regions = [(mp.mpf(1000),mp.mpf(1200),4)]
    regions += [(mp.mpf(a),mp.mpf(b),8) for a,b in zip((1200,5500,6500,19500,20500),(5500,6500,19500,20500,50000))]
    knots = [regions[0][0]]
    for lo,hi,n in regions:
        knots.extend(lo+(hi-lo)*k/n for k in range(1,n+1))
    tables, summary, sampled, runtime = {}, {}, [], []
    for name in NAMES:
        print('Fitting', name, flush=True)
        electronic = lru_cache(maxsize=None)(lambda t: ref.electronic(name,t))
        anchor = electronic(1000)
        def energy(t):
            cp,h,_ = electronic(t)
            return ((h-anchor[1])/ref.R-(t-1000), cp/ref.R-1)
        nodes = [energy(t) for t in knots]
        pieces, s_anchor = [], mp.mpf(0)
        max_cp,max_h,max_s,min_cv = mp.mpf(0),mp.mpf(0),mp.mpf(0),mp.inf
        endpoint_e,endpoint_cv = mp.mpf(0),mp.mpf(0)
        for i,(lo,hi) in enumerate(zip(knots[:-1],knots[1:])):
            width = hi-lo
            base = power_energy(nodes[i],nodes[i+1],width,[0]*4)
            xs = [mp.cos(mp.pi*(j+mp.mpf('.5'))/4) for j in range(4)]
            values=[]
            for x in xs:
                u=(x+1)/2
                values.append((energy(lo+width*u)[0]-horner(base,u))/(u*(1-u))**2)
            q = list(mp.lu_solve(mp.matrix([[x**k for k in range(4)] for x in xs]),mp.matrix(values)))
            a = power_energy(nodes[i],nodes[i+1],width,q)
            cv = [(k+1)*a[k+1]/width for k in range(7)]
            cp = cv[:];cp[0]+=1
            bernstein = [mp.fsum(cv[k]*mp.mpf(comb(j,k))/comb(6,k) for k in range(j+1)) for j in range(7)]
            min_cv = min(min_cv,min(bernstein))
            assert min(bernstein)>0
            for u,node in ((mp.mpf(0),nodes[i]),(mp.mpf(1),nodes[i+1])):
                endpoint_e=max(endpoint_e,abs(horner(a,u)-node[0]))
                endpoint_cv=max(endpoint_cv,abs(horner(cv,u)-node[1]))
            for j in range(65):
                u=mp.mpf(j)/64;t=lo+width*u
                rcp,rh,rs=electronic(t)
                e=horner(a,u);cpfit=horner(cp,u)*ref.R
                hfit=(e+t-1000)*ref.R
                sfit=(s_anchor+entropy_integral(cp,lo,width,u))*ref.R
                cperr=abs(cpfit/rcp-1)
                herr=abs(hfit-(rh-anchor[1]))
                serr=abs(sfit-(rs-anchor[2]))
                max_cp,max_h,max_s=max(max_cp,cperr),max(max_h,herr),max(max_s,serr)
                sampled.append(dict(species=name,temperature_K=str(t),cp_fit_J_mol_K=str(cpfit),
                                    h_increment_fit_J_mol=str(hfit),s_increment_fit_J_mol_K=str(sfit),
                                    cp_relative_error=str(cperr),h_absolute_error_J_mol=str(herr),
                                    s_absolute_error_J_mol_K=str(serr)))
                if j%4==0:
                    for bits in (32,64):
                        for scale in (1,1000):
                            native=rounded(t/scale,bits)
                            actual=native*scale
                            rcp,rh,_=electronic(actual)
                            runtime.append([bits,scale,NAMES.index(name),str(native),str(rcp/ref.R-1),str(((rh-anchor[1])/ref.R-(actual-1000))/scale)])
            pieces.append(dict(q=list(map(str,q)),entropy_left_over_R=str(s_anchor),
                               energy_power_u=list(map(str,a))))
            s_anchor+=entropy_integral(cp,lo,width,mp.mpf(1))
        tables[name]=dict(nodes=[[str(x) for x in n] for n in nodes],pieces=pieces)
        summary[name]=dict(max_relative_cp_error=str(max_cp),max_h_error_J_mol=str(max_h),
                           max_s_error_J_mol_K=str(max_s),min_bernstein_cv_over_R=str(min_cv),
                           max_endpoint_energy_residual_K=str(endpoint_e),max_endpoint_cv_residual=str(endpoint_cv))
    data=dict(precision_digits=70,energy_degree=7,capacity_degree=6,reference_temperature_K=1000,
              knots_K=list(map(str,knots)),regions=[[str(a),str(b),n] for a,b,n in regions],species=tables)
    (output_dir/'fits.json').write_text(json.dumps(data,indent=2)+'\n')
    (output_dir/'fit-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    with (output_dir/'fit-comparison.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(sampled[0]));w.writeheader();w.writerows(sampled)
    with (output_dir/'runtime-reference.txt').open('w') as f:
        for row in runtime:f.write(' '.join(map(str,row))+'\n')
    def literal(x):return mp.nstr(mp.mpf(x),28)+'L'
    def initializer(rows):return ',\n'.join('{'+','.join(literal(x) for x in row)+'}' for row in rows)
    header=['#pragma once','namespace fitted_levels {',f'inline constexpr int count={len(knots)-1};',
            f'inline constexpr int species_count={len(NAMES)};',
            'inline constexpr long double knots[]={'+','.join(literal(x) for x in knots)+'};',
            'inline constexpr long double nodes[species_count][count+1][2]={'+',\n'.join('{'+initializer(tables[n]['nodes'])+'}' for n in NAMES)+'};',
            'inline constexpr long double coefficients[species_count][count][4]={'+',\n'.join('{'+initializer([p['q'] for p in tables[n]['pieces']])+'}' for n in NAMES)+'};','}']
    (output_dir/'fitted-data.hpp').write_text('\n'.join(header)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True,
                        help='Directory for generated coefficients and verification data.')
    main(parser.parse_args().output_dir)
