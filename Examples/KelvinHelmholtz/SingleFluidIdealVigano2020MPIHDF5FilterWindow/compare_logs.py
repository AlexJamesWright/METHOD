#!/usr/bin/env python
"""
Compare a closure run with the baseline from their run logs, at matched times (works on a run in progress).

    python compare_logs.py BASELINE_LOG CLOSURE_LOG [--dt 0.5] [--seed 5e-7]

At every --dt time units up to the end of the shorter log, prints: E_mag and E_kin of the closure run relative to the baseline, the same for the
perturbation energy E_mag - seed (the uniform B_x seed carries --seed, 5e-7 for the Vigano initial data, and dominates E_mag early, so this ratio is the
sensitive one), the largest Lorentz factor of both, the number of failed conserved-to-primitive recoveries so far in each, and the closure's energy exchange
`net/abs` (positive: the closure drains the resolved field) and the exchange per unit E_mag, int J.EMF dV / E_mag (a rate), from the `closure` lines.
Values are linearly interpolated between the log lines (every 20 iterations). A ratio different from 1 by less than the run-to-run spread of
nominally identical runs (c = 1e-6 run of the scan) is not an effect.
"""
import argparse
import re

import numpy as np


def parse(log):
    diag, clo, fails = [], [], []
    t_now = 0.0
    n_fail = 0
    for line in open(log):
        m = re.match(r'diag t=(\S+) iter=\d+ Emag=(\S+) Ekin=(\S+) mass=\S+ Wmax=(\S+)', line)
        if m:
            diag.append([float(x) for x in m.groups()])
            fails.append(n_fail)
            continue
        m = re.match(r'closure t=(\S+) exchange=(\S+) abs=(\S+) net/abs=(\S+)', line)
        if m:
            clo.append([float(x) for x in m.groups()])
            continue
        if line.startswith('t = '):
            t_now = float(line.split()[2])
        elif 'failed' in line:
            n_fail += 1
    return np.array(diag), np.array(clo), np.array(fails)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('baseline_log')
    p.add_argument('closure_log')
    p.add_argument('--dt', type=float, default=0.5)
    p.add_argument('--seed', type=float, default=5e-7)
    a = p.parse_args()
    Db, Cb, Fb = parse(a.baseline_log)
    Dc, Cc, Fc = parse(a.closure_log)
    tend = min(Db[-1, 0], Dc[-1, 0])
    print("compared up to t = %.3f (baseline log ends %.3f, closure log %.3f)" % (tend, Db[-1, 0], Dc[-1, 0]))
    print("%-6s %-9s %-9s %-10s %-9s %-11s %-11s %-9s %-9s %-10s" % ("t", "Emag", "Ekin", "pert.Emag", "Wmax b/c", "fails b/c", "net/abs", "X/Emag", "", ""))
    for t in np.arange(a.dt, tend + 1e-9, a.dt):
        ib = lambda col: np.interp(t, Db[:, 0], Db[:, col])
        ic = lambda col: np.interp(t, Dc[:, 0], Dc[:, col])
        pert_b, pert_c = ib(1) - a.seed, ic(1) - a.seed
        fb = int(np.interp(t, Db[:, 0], Fb)); fc = int(np.interp(t, Dc[:, 0], Fc))
        na = np.interp(t, Cc[:, 0], Cc[:, 3]) if len(Cc) else float('nan')
        X = np.interp(t, Cc[:, 0], Cc[:, 1]) / ic(1) if len(Cc) else float('nan')
        print("%-6.2f %-9.4f %-9.4f %-10.4f %-4.2f/%-4.2f %-5d/%-5d %-11.3f %-9.3f" % (t, ic(1) / ib(1), ic(2) / ib(2), pert_c / pert_b if pert_b > 0 else float('nan'),
              ib(3), ic(3), fb, fc, na, X))


if __name__ == '__main__':
    main()
