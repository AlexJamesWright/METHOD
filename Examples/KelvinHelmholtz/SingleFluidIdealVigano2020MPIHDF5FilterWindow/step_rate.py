#!/usr/bin/env python
"""
Seconds per step of runs, from the `timing t=... iter=... wall=...` lines that the ClosureScan `main` prints with each `diag` line (since commit "timing lines").

    python step_rate.py LOG [LOG ...]

On a shared node the wall time of a whole run is not a measure of the cost of a scheme (identical WENO5 runs took 38 to 104 minutes): the fastest segments are. Printed per log: the iterations,
the failed-recovery lines per step (they cost time too), and the seconds per step over each segment between timing lines: the minimum, the 10th percentile and the median; the minimum
(or the 10th percentile) is the best estimate of the uncontended rate. Compare schemes by that, not by the elapsed time of the job.
"""
import re
import sys

import numpy as np


def main():
    print("%-40s %-8s %-10s %-9s %-9s %-9s %-9s" % ("log", "iters", "fail/step", "min s/st", "p10 s/st", "median", "whole run"))
    for log in sys.argv[1:]:
        it, wall, fails = [], [], 0
        with open(log, errors='replace') as fh:
            for line in fh:
                if line.startswith('timing'):
                    m = re.match(r'timing t=\S+ iter=(\d+) wall=(\S+)', line)
                    if m:
                        it.append(int(m.group(1)))
                        wall.append(float(m.group(2)))
                elif 'failed' in line:
                    fails += 1
        if len(it) < 3:
            print("%-40s no timing lines (older run)" % log[-40:])
            continue
        it, wall = np.array(it), np.array(wall)
        rate = np.diff(wall) / np.diff(it)
        print("%-40s %-8d %-10.3f %-9.3f %-9.3f %-9.3f %-9.3f" % (log[-40:], it[-1], fails / it[-1], rate.min(), np.percentile(rate, 10), np.median(rate), wall[-1] / it[-1]))


if __name__ == '__main__':
    main()
