#!/usr/bin/env python
"""
Mean and spread of the key numbers over groups of nominally identical runs, from the run logs (the noise ensemble and the ensemble at one closure setting).

    python ensemble_stats.py NAME=LOG1,LOG2,... [NAME=LOG1,LOG2,...]  [--reference LOG]

Per group (n runs): mean, standard deviation (n - 1) and range of: the time at which E_mag reaches 1e-6 and 2e-5, the growth rate ln 20 / (difference), the time at
which E_kin falls to 0.2, the largest E_mag over the run and E_mag at the end (and, with --reference, the delays of the first two against the reference log). With the
standard deviation of the noise runs one can say whether the difference of the means of two groups exceeds its standard error (sd_a^2/n_a + sd_b^2/n_b)^(1/2).
"""
import argparse

import numpy as np

from time_shift import crossing, read


def metrics(D):
    e1, e2 = crossing(D[:, 0], D[:, 1], 1e-6, False), crossing(D[:, 0], D[:, 1], 2e-5, False)
    return {'t(Emag=1e-6)': e1, 't(Emag=2e-5)': e2, 'rate 1e-6..2e-5': (np.log(20.0) / (e2 - e1)) if e1 is not None and e2 is not None else None,
            't(Ekin=0.2)': crossing(D[:, 0], D[:, 2], 0.2, True), 'Emag max': float(D[:, 1].max()), 'Emag end': float(D[-1, 1])}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('groups', nargs='+', help='NAME=LOG1,LOG2,...')
    p.add_argument('--reference')
    a = p.parse_args()
    ref = metrics(read(a.reference)) if a.reference else None
    results = {}
    for g in a.groups:
        name, logs = g.split('=', 1)
        rows = [metrics(read(l)) for l in logs.split(',')]
        results[name] = rows
        print("== %s (n = %d)" % (name, len(rows)))
        for key in rows[0]:
            v = [r[key] for r in rows if r[key] is not None]
            if not v:
                print("  %-18s -" % key)
                continue
            sd = np.std(v, ddof=1) if len(v) > 1 else float('nan')
            extra = ""
            if ref and key in ('t(Emag=1e-6)', 't(Emag=2e-5)') and ref[key] is not None:
                extra = "   delay vs reference: %.2f" % (np.mean(v) - ref[key])
            fmt = "%.3e" if key.startswith('Emag') else "%.3f"
            print(("  %-18s mean " + fmt + "  sd " + fmt + "  range " + fmt + " - " + fmt + "  (%d/%d runs)%s") % (key, np.mean(v), sd, min(v), max(v), len(v), len(rows), extra))
    names = list(results)
    if len(names) == 2:
        print("\n== difference of the means, in units of the standard error, group %s minus %s" % (names[1], names[0]))
        for key in results[names[0]][0]:
            va = [r[key] for r in results[names[0]] if r[key] is not None]
            vb = [r[key] for r in results[names[1]] if r[key] is not None]
            if len(va) > 1 and len(vb) > 1 and (np.var(va, ddof=1) + np.var(vb, ddof=1)) > 1e-12 * (abs(np.mean(va)) + abs(np.mean(vb)) + 1e-30) ** 2:
                se = np.sqrt(np.var(va, ddof=1) / len(va) + np.var(vb, ddof=1) / len(vb))
                fmt = "%.3e" if key.startswith('Emag') else "%.3f"
                print(("  %-18s " + fmt + "  (%.1f standard errors)") % (key, np.mean(vb) - np.mean(va), (np.mean(vb) - np.mean(va)) / se if se > 0 else float('nan')))


if __name__ == '__main__':
    main()
