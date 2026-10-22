#!/usr/bin/env python
"""
One line per run of a scan: the key numbers from the run logs.

    python scan_summary.py REFERENCE_LOG BASELINE_LOG LABEL:LOG [LABEL:LOG ...]

REFERENCE_LOG is the fine run (384^3), BASELINE_LOG the coarse run without the closure, then the closure runs. For each run: the times at which E_mag first reaches
1e-6 and 2e-5, their delay against the reference, the mean growth rate between them (ln 20 / time), the time at which E_kin falls to 0.2 and its delay,
E_mag at the end relative to the baseline, the failed recoveries (cumulative, `failed` lines) and the peak of the closure drain X/E_mag. Missing levels are shown as `-`.
"""
import argparse
import re

import numpy as np

from time_shift import crossing, read


def failures(log):
    return sum(1 for line in open(log) if 'failed' in line)


def peak_drain(log):
    best = None
    last_E = None
    for line in open(log):
        m = re.match(r'diag t=\S+ iter=\d+ Emag=(\S+)', line)
        if m:
            last_E = float(m.group(1))
            continue
        m = re.match(r'closure t=\S+ exchange=(\S+)', line)
        if m and last_E:
            x = float(m.group(1)) / last_E
            best = x if best is None else max(best, x)
    return best


def row(D, ref, base_end):
    e1, e2 = crossing(D[:, 0], D[:, 1], 1e-6, False), crossing(D[:, 0], D[:, 1], 2e-5, False)
    k2 = crossing(D[:, 0], D[:, 2], 0.2, True)
    f = lambda x: "%7.2f" % x if x is not None else "%7s" % "-"
    d = lambda x, r: "%6.2f" % (x - r) if x is not None and r is not None else "%6s" % "-"
    rate = "%6.2f" % (np.log(20.0) / (e2 - e1)) if e1 is not None and e2 is not None else "%6s" % "-"
    return e1, e2, k2, f(e1) + " " + d(e1, ref[0]) + "   " + f(e2) + " " + d(e2, ref[1]) + "  " + rate + "   " + f(k2) + " " + d(k2, ref[2]), D[-1, 1] / base_end if base_end else None


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('reference_log')
    p.add_argument('baseline_log')
    p.add_argument('runs', nargs='*')
    a = p.parse_args()
    R = read(a.reference_log)
    ref = (crossing(R[:, 0], R[:, 1], 1e-6, False), crossing(R[:, 0], R[:, 1], 2e-5, False), crossing(R[:, 0], R[:, 2], 0.2, True))
    B = read(a.baseline_log)
    print("%-14s %-15s %-15s %-6s %-15s %-9s %-7s %-8s" % ("run", "t(Emag=1e-6) dly", "t(Emag=2e-5) dly", "rate", "t(Ekin=.2) dly", "Emag/base", "fails", "X/E pk"))
    print("%-14s %7.2f %6s   %7.2f %6s  %6.2f   %7.2f %6s" % ("reference", ref[0], "", ref[1], "", np.log(20.0) / (ref[1] - ref[0]), ref[2], ""))
    items = [("baseline", a.baseline_log, B)] + [(r.split(':', 1)[0], r.split(':', 1)[1], read(r.split(':', 1)[1])) for r in a.runs]
    for label, log, D in items:
        _, _, _, s, ratio = row(D, ref, B[-1, 1])
        pk = peak_drain(log)
        print("%-14s %s  %-9s %-7d %-8s" % (label, s, "%.3f" % ratio if ratio else "-", failures(log), "%.3f" % pk if pk is not None else "-"))


if __name__ == '__main__':
    main()
