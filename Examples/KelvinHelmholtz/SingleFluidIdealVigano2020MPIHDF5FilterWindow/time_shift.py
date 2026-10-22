#!/usr/bin/env python
"""
How far one run's history is delayed relative to a reference run, from the `diag t=... Emag=... Ekin=...` lines of the two run logs.

    python time_shift.py REFERENCE_LOG OTHER_LOG [--tmin 0]

For a set of levels, prints the time at which each run reaches the level, and the difference: Ekin falling through 0.7 ... 0.1 and Emag rising through
1e-6 ... 1e-4 (levels outside the range of either run are skipped). A constant difference means the same evolution, delayed; a difference that
grows means a different rate. Used on the 384^3 run and the 128^3 baseline: the kinetic-energy history is delayed by 2.7-3.0 time units and the
magnetic one by about 4.7.
"""
import argparse
import re

import numpy as np


def read(log):
    rows = []
    with open(log, errors='replace') as fh:
        for line in fh:
            if not line.startswith('diag'):          # skip the (possibly millions of) `failed` lines without a regex
                continue
            m = re.match(r'diag t=(\S+) iter=\d+ Emag=(\S+) Ekin=(\S+)', line)
            if m:
                rows.append([float(x) for x in m.groups()])
    return np.array(rows)


def crossing(t, y, level, falling):
    """First time y passes level (interpolated), or None."""
    s = (y - level) * (-1.0 if falling else 1.0)
    idx = np.where((s[:-1] < 0) & (s[1:] >= 0))[0]
    if not len(idx):
        return None
    i = idx[0]
    return t[i] + (t[i + 1] - t[i]) * (0.0 - s[i]) / (s[i + 1] - s[i])


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('reference_log')
    p.add_argument('other_log')
    p.add_argument('--tmin', type=float, default=0.0)
    a = p.parse_args()
    R, O = read(a.reference_log), read(a.other_log)
    R, O = R[R[:, 0] >= a.tmin], O[O[:, 0] >= a.tmin]
    print("%-12s %-10s %-10s %-10s %-10s" % ("quantity", "level", "t_ref", "t_other", "shift"))
    for name, col, levels, falling in (("Ekin", 2, [0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1], True),
                                       ("Emag", 1, [1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4], False)):
        for lev in levels:
            tr, to = crossing(R[:, 0], R[:, col], lev, falling), crossing(O[:, 0], O[:, col], lev, falling)
            if tr is None or to is None:
                continue
            print("%-12s %-10.3g %-10.3f %-10.3f %-10.3f" % (name, lev, tr, to, to - tr))


if __name__ == '__main__':
    main()
