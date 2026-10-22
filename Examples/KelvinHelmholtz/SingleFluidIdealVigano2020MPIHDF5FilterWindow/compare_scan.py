#!/usr/bin/env python
"""
Experiment A driver: compare a scan of coarse closure runs with the 384^3 reference at every frame time they share.

    python compare_scan.py --fine-dir RUN384 --coarse c=0:SCAN/c0_d0_l0 c=1/12:SCAN/c0.0833333333_d0_l0 ... --out cmp [--tmin 4]

For each frame.checkpoint.*.hdf5 of the fine run that has a frame of the same time (within 1e-6) in EVERY coarse directory, runs
compare_coarse.py (the first coarse entry is the baseline) into OUT/t_<t>/, then prints and writes OUT/summary.txt: against time and
closure setting, the resolved magnetic energy relative to the filtered fine run, the high-k tail of the B spectrum relative to it, and
the relative error of B_z. Coarse frames come from the ClosureScan example, which writes them at the fine frame times (frame_times.py).
The magnetic-energy history figure is made once, from the first time.

Equal-E_kin mode (--match ekin --fine-log LOG): the runs of a chaotic flow are at different stages at the same time (the coarse runs lag the fine one by
2-5 time units, and by differing amounts from run to run), but E_kin, the integral of rho W (W - 1) printed in the `diag` lines, is a clock that the hydrodynamics
fixes (it converges with resolution, the magnetic field does not). For each fine frame, the coarse frame of each run whose run-log E_kin is closest to the fine
E_kin at that frame is compared with it (the first time each coarse run's E_kin falls through the fine value, then the nearest frame; frames more than --ekin-tol
(fractional, default 0.1) off in E_kin are not used and the fine frame is skipped). The summary then lists the matched coarse times and the E_kin mismatch too.
Needs the `run.log` of each coarse run next to its frames and --fine-log for the fine run. Spacing of the coarse frames sets the accuracy of the match.
"""
import argparse
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import viz_common as vc
from time_shift import crossing, read


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
    p.add_argument('--match', choices=['time', 'ekin'], default='time', help='pair frames at equal time (default) or at equal E_kin from the run logs')
    p.add_argument('--fine-log', help='run log of the fine run (--match ekin)')
    p.add_argument('--ekin-tol', type=float, default=0.1, help='largest fractional E_kin mismatch of a matched frame (--match ekin)')
    a = p.parse_args()
    fine = {}
    for d in a.fine_dir:
        fine.update(frame_times(d))
    coarse = []
    for item in a.coarse:
        label, d = item.split(':', 1)
        coarse.append((label, frame_times(d), d))
    if a.match == 'ekin':
        if not a.fine_log:
            sys.exit("--match ekin needs --fine-log")
        fineD = read(a.fine_log)
        coarseD = []
        for label, fr, d in coarse:
            lg = os.path.join(d, 'run.log')
            if not os.path.exists(lg):
                sys.exit("--match ekin needs %s" % lg)
            coarseD.append(read(lg))
    rows, first = [], True
    for tf in sorted(fine):
        if tf < a.tmin:
            continue
        match, tc_used, mis = [], [], []
        if a.match == 'time':
            for label, fr, d in coarse:
                hit = [t for t in fr if abs(t - tf) < a.tol]
                match.append(fr[hit[0]] if hit else None)
        else:
            level = np.interp(tf, fineD[:, 0], fineD[:, 2])
            for (label, fr, d), D in zip(coarse, coarseD):
                tc = crossing(D[:, 0], D[:, 2], level, True)
                best = min(fr, key=lambda t: abs(t - tc)) if (tc is not None and fr) else None
                if best is None or abs(np.interp(best, D[:, 0], D[:, 2]) / level - 1.0) > a.ekin_tol:
                    match.append(None)
                    continue
                match.append(fr[best])
                tc_used.append(best)
                mis.append(np.interp(best, D[:, 0], D[:, 2]) / level - 1.0)
        if any(m is None for m in match):
            print("t = %.6f: not in every coarse run%s, skipped" % (tf, "" if a.match == 'time' else " (no frame with E_kin within %.0f %%)" % (100 * a.ekin_tol)))
            continue
        out = os.path.join(a.out, 't_%.4f' % tf)
        cmd = [sys.executable, os.path.join(HERE, 'compare_coarse.py'), '--fine', fine[tf], '--out', out, '--time-tol', str(a.tol if a.match == 'time' else 1e9),
               '--coarse'] + ['%s:%s' % (c[0], m) for c, m in zip(coarse, match)]
        if first and a.match == 'time':
            cmd += ['--fine-frames'] + a.fine_dir
        subprocess.check_call(cmd, stdout=subprocess.DEVNULL)
        z = np.load(os.path.join(out, 'comparison.npz'))
        rows.append((tf, z['Emag'] / float(z['Emag_ref']), z['tail'], z['err_Bz'], tc_used, mis))
        first = False
    if not rows:
        sys.exit("no common frame times")
    labels = [c[0] for c in coarse]
    if a.match == 'ekin':
        titles = (("RESOLVED MAGNETIC ENERGY / filtered fine, at equal E_kin (1 = agreement)", 1), ("B-SPECTRUM TAIL (k > k_Nyquist/2) / filtered fine, at equal E_kin", 2),
                  ("RELATIVE ERROR OF B_z against the filtered fine field (pointwise; dominated by eddy positions)", 3))
    else:
        titles = (("RESOLVED MAGNETIC ENERGY / filtered fine (1 = agreement)", 1), ("B-SPECTRUM TAIL (k > k_Nyquist/2) / filtered fine", 2),
                  ("RELATIVE ERROR OF B_z against the filtered fine field (pointwise; dominated by eddy positions)", 3))
    lines = []
    for title, idx in titles:
        lines += ["", title, "%-9s" % "t" + "".join("%12s" % l for l in labels)]
        for r in rows:
            lines.append("%-9.4f" % r[0] + "".join("%12.4f" % x for x in r[idx]))
    if a.match == 'ekin':
        lines += ["", "MATCHED COARSE TIMES (and fractional E_kin mismatch of the frame)", "%-9s" % "t_fine" + "".join("%16s" % l for l in labels)]
        for r in rows:
            lines.append("%-9.4f" % r[0] + "".join("%9.3f(%+4.0f%%)" % (tc, 100 * m) for tc, m in zip(r[4], r[5])))
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(a.out, 'summary.txt'), 'w') as fh:
        fh.write(text + "\n")


if __name__ == '__main__':
    main()
