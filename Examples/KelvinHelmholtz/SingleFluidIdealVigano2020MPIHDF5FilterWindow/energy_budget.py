#!/usr/bin/env python
"""
Budget of the resolved magnetic energy from frames: dE/dt = P - X - D (lab frame, flat, ideal MHD, mu0 = 1).

For frames of one run, optionally block-averaged by r (the filtered fine run; r = 1 for a coarse run itself):
    E = int |<B>|^2/2 dV
    P = int J.(<v> x <B>) dV      production by the resolved flow            (J = curl <B>, centred differences on the coarse grid)
    X = int J.( -(<v x B> - <v> x <B>) ) dV   drain by the sub-filter EMF   (zero for r = 1: nothing is sub-filter)
    D = P - X - dE/dt             what is left: numerical dissipation and discretisation error, in the same units
dE/dt is the centred difference of E between the neighbouring frames, so the frames should be close in time compared with 1/(growth rate).
Everything is printed divided by E (a rate per unit time), against time. The closure of SubfilterGradient models X.

    python energy_budget.py --block 3 --dirs RUN384 [--tmin 2] [--out budget_b3.txt]        # the filtered fine run
    python energy_budget.py --block 1 --dirs SCAN/c0_d0_l0 [--tmin 2] [--out budget_coarse.txt]  # the coarse baseline

For the reading of the numbers: if the coarse run's D/E is as large as the filtered fine run's X/E (or larger), the numerics already remove what
the closure would, and adding the closure double counts.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import viz_common as vc
from apriori_emf import cross, dot, bm


def one(path, r):
    S = vc.load(path)
    P, dx = S['P'], S['dx'] * r
    v = [P['vx'], P['vy'], P['vz']]
    B = [P['Bx'], P['By'], P['Bz']]
    if r > 1:
        vxB = [bm(c, r) for c in cross(v, B)]
        vm = [bm(c, r) for c in v]
        Bm = [bm(c, r) for c in B]
    else:
        vm, Bm = v, B
        vxB = cross(v, B)
    vmBm = cross(vm, Bm)
    J = [vc.ddx(Bm[2], 1, dx) - vc.ddx(Bm[1], 2, dx), vc.ddx(Bm[0], 2, dx) - vc.ddx(Bm[2], 0, dx), vc.ddx(Bm[1], 0, dx) - vc.ddx(Bm[0], 1, dx)]
    dV = dx ** 3
    E = 0.5 * dV * np.sum(dot(Bm, Bm))
    prod = dV * np.sum(dot(J, vmBm))
    drain = dV * np.sum(dot(J, [-(a - b) for a, b in zip(vxB, vmBm)])) if r > 1 else 0.0
    return S['t'], E, prod, drain


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('--block', type=int, required=True)
    p.add_argument('--dirs', nargs='+', required=True)
    p.add_argument('--tmin', type=float, default=0.0)
    p.add_argument('--tmax', type=float, default=1e9)
    p.add_argument('--out', default='budget.txt')
    a = p.parse_args()
    files = []
    for d in a.dirs:
        files += vc.frames(d)
    rows = []
    for f in files:
        t, E, Pr, X = one(f, a.block)
        rows.append((t, E, Pr, X))
        print("read t = %.3f" % t, file=sys.stderr, flush=True)
    rows.sort()
    rows = np.array(rows)
    lines = ["%-8s %-12s %-9s %-9s %-9s %-9s %-9s" % ("t", "E", "dlnE/dt", "P/E", "X/E", "D/E", "D/(P-X)")]
    for i in range(1, len(rows) - 1):
        t, E, Pr, X = rows[i]
        if not (a.tmin <= t <= a.tmax):
            continue
        dEdt = (rows[i + 1, 1] - rows[i - 1, 1]) / (rows[i + 1, 0] - rows[i - 1, 0])
        D = Pr - X - dEdt
        lines.append("%-8.3f %-12.4e %-9.3f %-9.3f %-9.3f %-9.3f %-9.3f" % (t, E, dEdt / E, Pr / E, X / E, D / E, D / (Pr - X) if Pr != X else np.nan))
    print("\n".join(lines))
    open(a.out, 'w').write("\n".join(lines) + "\n")


if __name__ == '__main__':
    main()
