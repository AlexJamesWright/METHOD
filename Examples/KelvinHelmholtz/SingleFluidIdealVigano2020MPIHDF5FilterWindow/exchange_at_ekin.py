#!/usr/bin/env python
"""
The closure's own energy exchange in a coarse run, at equal E_kin, against the exact sub-filter exchange of the filtered fine run.

    python exchange_at_ekin.py --fine-log FINE_LOG --exact BUDGET_FILE [BUDGET_FILE ...] LABEL=C:RUN_LOG [LABEL=C:RUN_LOG ...] [--window 0.5] [--tmin 5]

An amplitude of a flux does not depend on when the run's field happens to take off, so it is not limited by the run-to-run noise that the levels of E_mag are (design doc
sections 15-18). The exact exchange X/E_mag of the filtered fine run comes from the budget tables of `energy_budget.py --block R` (columns t, E, dlnE/dt, P/E, X/E, D/E, D/(P-X); X
positive: the resolved field loses energy to the sub-filter scales; give the tables of the same R, and of the run whose log is --fine-log, for example R = 3 for the 384^3 run).
For each coarse run (closure parameter C, the `closure t=... exchange=...` and `diag` lines of its log): at each row of the table the fine E_kin is read from --fine-log, the time at which
the coarse run's E_kin first falls through that value is found, and the closure's exchange per unit energy there, `exchange / E_mag` averaged over +- window/2 time units
(the lines are 20 iterations apart and the exchange instantaneous), is compared with the exact one. The exchange of the closure is proportional to C, so
`C_eff = C x exact / closure` is the C that would reproduce the exact exchange at that stage; if the closure is consistent with the data it is the same for every run
(a-priori best C: 0.10-0.125). Printed per run and, in the last lines, the median of C_eff over the rows with t_fine >= --tmin.
"""
import argparse
import re

import numpy as np

from time_shift import crossing, read


def read_closure(log):
    rows = []
    with open(log, errors='replace') as fh:
        for line in fh:
            if not line.startswith('closure'):
                continue
            m = re.match(r'closure t=(\S+) exchange=(\S+) abs=(\S+) net/abs=(\S+)', line)
            if m:
                rows.append([float(x) for x in m.groups()])
    return np.array(rows)


def read_table(path):
    rows = []
    for line in open(path):
        parts = line.split()
        try:
            vals = [float(x) for x in parts]
        except ValueError:
            continue
        if len(vals) >= 5:
            rows.append(vals)
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('runs', nargs='+', help='LABEL=C:LOG')
    p.add_argument('--fine-log', required=True)
    p.add_argument('--exact', nargs='+', required=True, help='budget tables of energy_budget.py (t, E, dlnE/dt, P/E, X/E, ...)')
    p.add_argument('--window', type=float, default=0.5, help='averaging window of the closure exchange, time units')
    p.add_argument('--tmin', type=float, default=5.0, help='first fine time of the summary median')
    a = p.parse_args()
    F = read(a.fine_log)
    exact = []
    for path in a.exact:
        exact += read_table(path)
    exact = sorted({row[0]: row for row in exact}.values())
    summary = []
    for item in a.runs:
        label, rest = item.split('=', 1)
        c, log = rest.split(':', 1)
        c = float(c)
        D, C = read(log), read_closure(log)
        if len(C) == 0:
            print("== %s: no closure lines in %s" % (label, log))
            continue
        xe = C[:, 1] / np.interp(C[:, 0], D[:, 0], D[:, 1])            # exchange per unit E_mag of the coarse run
        print("== %s (C = %g)" % (label, c))
        print("%-8s %-8s %-8s %-10s %-10s %-8s %-8s" % ("t_fine", "Ekin_f", "t_coarse", "X/E exact", "X/E coarse", "coarse/exact", "C_eff"))
        ceff = []
        for row in exact:
            tf, x_exact = row[0], row[4]
            lev = np.interp(tf, F[:, 0], F[:, 2])
            tc = crossing(D[:, 0], D[:, 2], lev, True)
            if tc is None or tc > C[-1, 0] - a.window / 2 or tc < C[0, 0] + a.window / 2:
                continue
            sel = np.abs(C[:, 0] - tc) <= a.window / 2
            if sel.sum() < 2:
                continue
            x_c = float(np.mean(xe[sel]))
            ratio = x_c / x_exact if x_exact > 0 else float('nan')
            ce = c / ratio if ratio > 0 else float('nan')
            print("%-8.3f %-8.4f %-8.3f %-10.4f %-10.4f %-8.3f %-8.4f" % (tf, lev, tc, x_exact, x_c, ratio, ce))
            if tf >= a.tmin and np.isfinite(ce):
                ceff.append(ce)
        if ceff:
            summary.append((label, c, np.median(ceff), np.min(ceff), np.max(ceff), len(ceff)))
    if summary:
        print("\nC_eff (the C that reproduces the exact exchange at equal E_kin), median [min, max] over the rows with t_fine >= %g:" % a.tmin)
        for label, c, med, lo, hi, n in summary:
            print("  %-10s C = %-10g C_eff = %.3f [%.3f, %.3f]  (%d rows)" % (label, c, med, lo, hi, n))


if __name__ == '__main__':
    main()
