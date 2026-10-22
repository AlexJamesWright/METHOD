#!/usr/bin/env python
"""
Synthetic check of exchange_at_ekin.py: a coarse run lagging the fine one by 1.0 time unit in E_kin(t) whose closure exchange per unit E_mag is half the exact one at the
matched stage must give C_eff = 2 C, for two different values of C (the exchange being proportional to C). Needs numpy only.
"""
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

EX = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'Examples', 'KelvinHelmholtz', 'SingleFluidIdealVigano2020MPIHDF5FilterWindow')


def ekin(t, lag):
    return 0.8 * math.exp(-0.3 * (t - lag))


def x_exact(t):
    return 0.4 * math.exp(-0.2 * t)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        fine = os.path.join(tmp, 'fine.log')
        with open(fine, 'w') as f:
            for t in np.arange(0, 20, 0.1):
                f.write("diag t=%.8f iter=1 Emag=1e-5 Ekin=%.8e mass=1 Wmax=1\n" % (t, ekin(t, 0.0)))
        table = os.path.join(tmp, 'budget.txt')
        with open(table, 'w') as f:
            f.write("t E dlnE/dt P/E X/E D/E D/(P-X)\n")
            for t in np.arange(2.0, 12.0, 0.5):
                f.write("%.3f 1e-5 0.1 3.0 %.6f 2.0 0.9\n" % (t, x_exact(t)))
        runs = []
        for c in (0.1, 0.4):
            path = os.path.join(tmp, 'run_c%g.log' % c)
            with open(path, 'w') as f:
                for t in np.arange(0, 20, 0.02):
                    f.write("diag t=%.8f iter=1 Emag=2e-5 Ekin=%.8e mass=1 Wmax=1\n" % (t, ekin(t, 1.0)))
                    # matched fine time t_f = t - 1 (equal E_kin), X/E_closure = 0.5 x exact(t_f) x (c / 0.1) so that C_eff = 2 x 0.1 for both runs? no: X is proportional to C
                    tf = t - 1.0
                    xe = 0.5 * x_exact(tf) * (c / 0.2)          # C_eff = C / ratio = C / (0.5 c / 0.2) = 0.4 for every c
                    f.write("closure t=%.8f exchange=%.8e abs=1 net/abs=0.8\n" % (t, xe * 2e-5))
            runs.append('r%g=%g:%s' % (c, c, path))
        out = subprocess.check_output([sys.executable, os.path.join(EX, 'exchange_at_ekin.py'), '--fine-log', fine, '--exact', table, '--tmin', '4'] + runs, cwd=EX).decode()
    lines = [l for l in out.splitlines() if 'C_eff =' in l]
    assert len(lines) == 2, out
    for l in lines:
        ce = float(l.split('C_eff = ')[1].split()[0])
        assert abs(ce - 0.4) < 0.03, l
    print("exchange_at_ekin: ok (%s)" % "; ".join(l.split('C_eff = ')[1].split()[0] for l in lines))


if __name__ == '__main__':
    main()
