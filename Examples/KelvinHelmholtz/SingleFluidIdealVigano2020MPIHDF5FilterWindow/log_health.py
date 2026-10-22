#!/usr/bin/env python
"""
Health of run logs: failed recoveries, non-finite diagnostics and Lorentz-factor spikes.

    python log_health.py LOG [LOG ...] [--wmax 5]

Per log: the number of `diag` lines, the failed-recovery lines, the `diag` lines with a non-finite value (a NaN cell existed at that moment; the times of the first few), the largest W_max with
its time, and the number of `diag` lines with W_max above --wmax. Isolated events (one non-finite line, a spike of a few lines) are repaired by the fallback within a few steps and are of the
kind every run has; a burst of non-finite lines or a W_max of tens (the WENO9 run at 128^3: ten lines, W_max 61, design doc section 23) means a corrupted solution. Reads only the log.
"""
import argparse
import re


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('logs', nargs='+')
    p.add_argument('--wmax', type=float, default=5.0)
    a = p.parse_args()
    print("%-52s %6s %8s %5s %9s %7s %s" % ("log", "diag", "failed", "nan", "Wmax max", "at t", "non-finite diag at t (first 4) | n(Wmax>%g)" % a.wmax))
    for log in a.logs:
        nd = nf = 0
        nan_t, big = [], 0
        wmax, wt = 0.0, 0.0
        with open(log, errors='replace') as fh:
            for line in fh:
                if line.startswith('diag'):
                    nd += 1
                    m = re.match(r'diag t=(\S+) iter=\d+ Emag=(\S+) Ekin=(\S+) mass=(\S+) Wmax=(\S+)', line)
                    if not m:
                        continue
                    t = float(m.group(1))
                    vals = m.groups()[1:4]
                    if any('nan' in v.lower() or 'inf' in v.lower() for v in vals):
                        nan_t.append(t)
                    w = float(m.group(5))
                    if w > wmax:
                        wmax, wt = w, t
                    if w > a.wmax:
                        big += 1
                elif 'failed' in line:
                    nf += 1
        print("%-52s %6d %8d %5d %9.2f %7.2f %s | %d" % (log[-52:], nd, nf, len(nan_t), wmax, wt, " ".join("%.2f" % t for t in nan_t[:4]), big))


if __name__ == '__main__':
    main()
