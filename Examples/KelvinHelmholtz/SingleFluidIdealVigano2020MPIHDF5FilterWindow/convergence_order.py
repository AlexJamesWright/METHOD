#!/usr/bin/env python
"""
Convergence of the energy histories with resolution, from the `diag` lines of runs of the same problem at several grids.

    python convergence_order.py 128:SCAN/c0_d0_l0/run.log 192:CONV/N192_c0/run.log 256:CONV/N256_c0/run.log 384:RUN384/khi384-run-4534288.out [--tmin 0]

For each level (E_kin falling through 0.7 ... 0.2, E_mag rising through 1e-6 ... 2e-5) the time t_N at which the run of N cells per side reaches it
is fitted by  t_N = t_inf + A N^(-p)  (p searched on a grid, t_inf and A by least squares; needs at least three resolutions, and with
three the fit is exact). Printed: the crossing times, p, the extrapolated t_inf (the continuum value), and how late each grid is
relative to it. The shift between two grids at a level is the quantity of time_shift.py. Not a proof of asymptotic convergence:
the fit assumes the runs are in the regime where one power law holds, and the run-to-run spread of nominally identical runs
(0.1-0.6 time units, design doc section 11) is the uncertainty of every t_N.
"""
import argparse

import numpy as np

from time_shift import crossing, read


def fit(Ns, ts):
    """Best p, t_inf, A of t = t_inf + A N^-p over a grid of p (least squares for the other two); exact for three points."""
    Ns, ts = np.asarray(Ns, float), np.asarray(ts, float)
    best = None
    for p in np.arange(0.2, 6.0, 0.01):
        X = np.vstack([np.ones_like(Ns), Ns ** -p]).T
        coef, res, *_ = np.linalg.lstsq(X, ts, rcond=None)
        r = np.sum((X @ coef - ts) ** 2)
        if best is None or r < best[0]:
            best = (r, p, coef[0], coef[1])
    return best[1], best[2], best[3], best[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument('runs', nargs='+', help='N:LOG pairs')
    ap.add_argument('--tmin', type=float, default=0.0)
    a = ap.parse_args()
    runs = []
    for item in a.runs:
        n, log = item.split(':', 1)
        D = read(log)
        runs.append((int(n), D[D[:, 0] >= a.tmin]))
    runs.sort(key=lambda r: r[0])
    Ns = [r[0] for r in runs]
    print("%-6s %-9s " % ("qty", "level") + " ".join("t(%d)" % n for n in Ns).ljust(8 * len(Ns)) + "   p      t_inf    late: " + " ".join("N=%d" % n for n in Ns))
    ps = []
    for name, col, levels, falling in (("Ekin", 2, [0.7, 0.6, 0.5, 0.4, 0.3, 0.2], True), ("Emag", 1, [1e-6, 2e-6, 5e-6, 1e-5, 2e-5], False)):
        for lev in levels:
            ts = [crossing(D[:, 0], D[:, col], lev, falling) for _, D in runs]
            if any(t is None for t in ts):
                continue
            if len(ts) < 3:
                print("%-6s %-9.3g " % (name, lev) + " ".join("%.3f" % t for t in ts))
                continue
            p, tinf, A, r = fit(Ns, ts)
            ps.append(p)
            print("%-6s %-9.3g " % (name, lev) + " ".join("%-7.3f" % t for t in ts) + "   %-5.2f  %-7.3f     " % (p, tinf) + " ".join("%-6.2f" % (t - tinf) for t in ts))
    if ps:
        print("\nfitted order p: median %.2f, range %.2f-%.2f over %d levels (p = 1: first order, 2: second order)" % (np.median(ps), min(ps), max(ps), len(ps)))


if __name__ == '__main__':
    main()
