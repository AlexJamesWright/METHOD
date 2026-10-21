#!/usr/bin/env python
"""
Experiment A driver: compare a scan of coarse closure runs with the 384^3 reference at every frame time they share.

    python compare_scan.py --fine-dir RUN384 --coarse c=0:SCAN/c0_d0_l0 c=1/12:SCAN/c0.0833333333_d0_l0 ... --out cmp [--tmin 4]

For each frame.checkpoint.*.hdf5 of the fine run that has a frame of the same time (within 1e-6) in EVERY coarse directory, runs
compare_coarse.py (the first coarse entry is the baseline) into OUT/t_<t>/, then prints and writes OUT/summary.txt: against time and
closure setting, the resolved magnetic energy relative to the filtered fine run, the high-k tail of the B spectrum relative to it, and
the relative error of B_z. Coarse frames come from the ClosureScan example, which writes them at the fine frame times (frame_times.py).
The magnetic-energy history figure is made once, from the first time.
"""
import argparse
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import viz_common as vc


def frame_times(d):
    import glob
    import h5py
    out = {}
    for f in glob.glob(os.path.join(d, 'frame.checkpoint.*.hdf5')):
        try:
            with h5py.File(f, 'r') as h:
                out[float(np.ravel(h.attrs['t'])[0])] = f
        except OSError:
            pass
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('--fine-dir', nargs='+', required=True, help='directories of the fine run (frames)')
    p.add_argument('--coarse', nargs='+', required=True, help='LABEL:DIR, first is the baseline')
    p.add_argument('--out', default='cmp')
    p.add_argument('--tmin', type=float, default=0.0)
    p.add_argument('--tol', type=float, default=1e-6)
    a = p.parse_args()
    fine = {}
    for d in a.fine_dir:
        fine.update(frame_times(d))
    coarse = []
    for item in a.coarse:
        label, d = item.split(':', 1)
        coarse.append((label, frame_times(d)))
    rows, first = [], True
    for tf in sorted(fine):
        if tf < a.tmin:
            continue
        match = []
        for label, fr in coarse:
            hit = [t for t in fr if abs(t - tf) < a.tol]
            match.append(fr[hit[0]] if hit else None)
        if any(m is None for m in match):
            print("t = %.6f: not in every coarse run, skipped" % tf)
            continue
        out = os.path.join(a.out, 't_%.4f' % tf)
        cmd = [sys.executable, os.path.join(HERE, 'compare_coarse.py'), '--fine', fine[tf], '--out', out, '--time-tol', str(a.tol),
               '--coarse'] + ['%s:%s' % (c[0], m) for c, m in zip(coarse, match)]
        if first:
            cmd += ['--fine-frames'] + a.fine_dir
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL)
        z = np.load(os.path.join(out, 'comparison.npz'))
        rows.append((tf, z['Emag'] / float(z['Emag_ref']), z['tail'], z['err_Bz']))
        first = False
    if not rows:
        sys.exit("no common frame times")
    labels = [c[0] for c in coarse]
    lines = []
    for title, idx in (("RESOLVED MAGNETIC ENERGY / filtered fine (1 = agreement)", 1), ("B-SPECTRUM TAIL (k > k_Nyquist/2) / filtered fine", 2),
                       ("RELATIVE ERROR OF B_z against the filtered fine field (pointwise; dominated by eddy positions)", 3)):
        lines += ["", title, "%-9s" % "t" + "".join("%12s" % l for l in labels)]
        for r in rows:
            lines.append("%-9.4f" % r[0] + "".join("%12.4f" % x for x in r[idx]))
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(a.out, 'summary.txt'), 'w') as fh:
        fh.write(text + "\n")


if __name__ == '__main__':
    main()
