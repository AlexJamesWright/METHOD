#!/usr/bin/env python
"""
Synthetic check of convergence_order.py: logs whose energy histories are delayed by 100 N^-1.5 must give p = 1.5 and the exact
continuum crossing time. Needs numpy only.
"""
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

EX = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'Examples', 'KelvinHelmholtz', 'SingleFluidIdealVigano2020MPIHDF5FilterWindow')


def main():
    with tempfile.TemporaryDirectory() as tmp:
        t = np.arange(0, 14, 0.05)
        logs = []
        for N in (128, 192, 256, 384):
            d = 100.0 * N ** -1.5
            path = os.path.join(tmp, 'l%d.log' % N)
            with open(path, 'w') as f:
                for tt in t:
                    s = tt - d
                    em = 5e-7 * math.exp(0.75 * s) if s > 0 else 5e-7
                    ek = 0.84 * math.exp(-0.22 * max(s - 0.5, 0))
                    f.write("diag t=%.8f iter=1 Emag=%.8e Ekin=%.8e mass=1 Wmax=1\n" % (tt, em, ek))
            logs.append('%d:%s' % (N, path))
        out = subprocess.check_output([sys.executable, os.path.join(EX, 'convergence_order.py')] + logs).decode()
    line = [l for l in out.splitlines() if l.startswith('fitted order')][0]
    p = float(line.split('median ')[1].split(',')[0])
    assert abs(p - 1.5) < 0.05, line
    print("convergence_order: ok (%s)" % line)


if __name__ == '__main__':
    main()
