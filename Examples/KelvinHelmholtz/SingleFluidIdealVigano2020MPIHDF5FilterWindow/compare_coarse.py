#!/usr/bin/env python
"""
Compare coarse (sub-filter closure) runs of the Vigano2020 KHI with the filtered fine run.

The fine run is coarse-grained by averaging its cell-averaged conserved variables over r^3 blocks of cells
(r = N_fine / N_coarse); this is the exact finite-volume restriction and a top-hat filter of width one coarse cell
(variance Dx^2/12 per direction, i.e. c = 1/12 in sigma^2 = c Dx^2). A coarse run with the right closure should
reproduce it better than the same run without the closure. Lab-frame, no Lagrangian filter: the real-data version
(filtered 384^3 pickles) is the same comparison with the reference supplied by the filter pipeline.

    python compare_coarse.py --fine FINE_FILE --coarse c=0:DIR0/data_final.hdf5 c=1/3:DIR1/data_final.hdf5 ... \
           [--fine-frames DIR [DIR ...]] [--out OUTDIR]

  --fine         a fine-run file (frame or restart frame) at the SAME time as the coarse data_final files (primitives
                 only is fine: conserved variables are rebuilt from them; a stored Conserved group is used when present)
  --coarse       LABEL:FILE pairs; the first is the baseline (normally the run without closure). The run.log beside each
                 file, if present, supplies the Emag(t) diagnostics for the time-series figure
  --fine-frames  directories of frame.checkpoint.*.hdf5 of the fine run, for the fine curve of the time series

Writes comparison.txt (the tables), comparison.npz, spectra.png, slices.png, timeseries.png.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import viz_common as vc


def read_diag(log):
    """(t, Emag, Ekin) from the `diag` lines of a METHOD run log."""
    rows = []
    for line in open(log):
        m = re.match(r'diag t=(\S+) iter=\d+ Emag=(\S+) Ekin=(\S+)', line)
        if m:
            rows.append([float(x) for x in m.groups()])
    return np.array(rows)


def reference(fine_path, r):
    """Block-averaged conserved variables of the fine file, and its time."""
    S = vc.load(fine_path)
    C = vc.load_cons(fine_path) or vc.conserved_from_prims(S['P'])
    return S['t'], S['dx'], {k: vc.block_mean(v, r) for k, v in C.items()}, S


def fluct_norm(a):
    return np.sqrt(np.sum((a - a.mean()) ** 2))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('--fine', required=True)
    p.add_argument('--coarse', nargs='+', required=True, help='LABEL:FILE, first is the baseline')
    p.add_argument('--fine-frames', nargs='*', default=[])
    p.add_argument('--out', default='.')
    p.add_argument('--time-tol', type=float, default=1e-6)
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)

    runs = []
    for item in a.coarse:
        label, path = item.split(':', 1)
        S = vc.load(path)
        C = vc.load_cons(path)
        if C is None:
            C = vc.conserved_from_prims(S['P'])
        runs.append({'label': label, 'path': path, 'S': S, 'C': C})
    n_c = runs[0]['S']['P']['rho'].shape[0]
    with_fine = vc.load(a.fine)
    n_f = with_fine['P']['rho'].shape[0]
    if n_f % n_c:
        sys.exit("fine grid %d is not a multiple of the coarse grid %d" % (n_f, n_c))
    r = n_f // n_c
    t, dx_f, ref, Sf = reference(a.fine, r)
    dx_c = runs[0]['S']['dx']
    for R in runs:
        if abs(R['S']['t'] - t) > a.time_tol:
            sys.exit("time mismatch: fine t=%.9f, %s t=%.9f" % (t, R['label'], R['S']['t']))
    out = []
    pr = lambda s='': (print(s), out.append(s))
    pr("fine %s (N=%d) at t = %.6f block-averaged by r = %d onto N = %d (coarse dx = %.5g)" % (os.path.basename(a.fine), n_f, t, r, n_c, dx_c))
    pr("baseline: %s\n" % runs[0]['label'])

    # 1. relative error of each conserved variable: |c - ref| / |ref - mean(ref)|
    pr("RELATIVE ERROR  |c - ref| / |ref - <ref>|   (lower is better; % change relative to the baseline in brackets)")
    pr("%-10s" % "variable" + "".join("%18s" % R['label'] for R in runs))
    err = {}
    for k in vc.CONS:
        nrm = fluct_norm(ref[k])
        err[k] = [np.sqrt(np.sum((R['C'][k] - ref[k]) ** 2)) / nrm for R in runs]
        row = "%-10s" % k
        for i, e in enumerate(err[k]):
            row += "%10.4f" % e + ("" if i == 0 or err[k][0] < 1e-12 else "(%+5.1f%%)" % (100 * (e / err[k][0] - 1)))
            row = row.ljust(10 + 18 * (i + 1))
        pr(row)

    # 1b. size of the change the closure makes, |c_i - c_baseline| / |c_baseline - ref|: pointwise errors of a turbulent flow are
    # dominated by the resolution difference itself (the eddies are not in the same place), so this says how much of that error
    # the closure could possibly have addressed
    pr("\nCLOSURE SHIFT  |c - c_baseline| / |c_baseline - ref|   (how large the change made by the closure is compared with the baseline error)")
    pr("%-10s" % "variable" + "".join("%12s" % R['label'] for R in runs[1:]))
    for k in vc.CONS:
        den = np.sqrt(np.sum((runs[0]['C'][k] - ref[k]) ** 2))
        pr("%-10s" % k + "".join("%12.4f" % (np.sqrt(np.sum((R['C'][k] - runs[0]['C'][k]) ** 2)) / max(den, 1e-300)) for R in runs[1:]))

    # 2. amplitude of the fluctuations: rms(c - <c>) / rms(ref - <ref>)
    pr("\nFLUCTUATION AMPLITUDE  rms(c - <c>) / rms(ref - <ref>)   (1 = the coarse run keeps the filtered fine amplitude)")
    pr("%-10s" % "variable" + "".join("%12s" % R['label'] for R in runs))
    amp = {}
    for k in vc.CONS:
        amp[k] = [fluct_norm(R['C'][k]) / fluct_norm(ref[k]) for R in runs]
        pr("%-10s" % k + "".join("%12.4f" % x for x in amp[k]))

    # 3. integrated resolved magnetic energy  int |B|^2/2 dV  (the quantity the closure is meant to feed)
    vol = dx_c ** 3
    Emag_ref = 0.5 * vol * sum(np.sum(ref[k] ** 2) for k in ('Bx', 'By', 'Bz'))
    pr("\nRESOLVED MAGNETIC ENERGY  int |B|^2/2 dV;  reference = filtered fine field: %.6e" % Emag_ref)
    Em = []
    for R in runs:
        e = 0.5 * vol * sum(np.sum(R['S']['P'][k] ** 2) for k in ('Bx', 'By', 'Bz'))
        Em.append(e)
        pr("  %-12s %.6e   ratio to reference %.4f" % (R['label'], e, e / Emag_ref))

    # 4. spectra of B and of S on the coarse grid
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    spec = {}
    for ax_, name, keys in ((ax[0], 'B', ('Bx', 'By', 'Bz')), (ax[1], 'S', ('Sx', 'Sy', 'Sz'))):
        kk, Pr = vc.shell_spectrum([ref[k] - ref[k].mean() for k in keys], dx_c)
        ax_.loglog(kk[1:], Pr[1:], 'k-', lw=2, label='filtered fine (ref)')
        spec[name + '_ref'] = Pr
        for R in runs:
            kk, Pc = vc.shell_spectrum([R['C'][k] - R['C'][k].mean() for k in keys], dx_c)
            ax_.loglog(kk[1:], Pc[1:], label=R['label'])
            spec[name + '_' + R['label']] = Pc
        ax_.set_xlabel('k'); ax_.set_ylabel('P(k)  of ' + name); ax_.legend(fontsize=8); ax_.set_title('shell spectrum of ' + name + ', t = %.3f' % t)
    spec['k'] = kk
    fig.tight_layout(); fig.savefig(os.path.join(a.out, 'spectra.png'), dpi=130); plt.close(fig)
    # high-k part of the spectrum: k above half the coarse Nyquist
    pr("\nSPECTRAL TAIL  sum of P_B(k) for k > k_Nyquist/2, relative to the reference")
    hi = kk > 0.5 * kk[-1]
    for R in runs:
        pr("  %-12s %.4f" % (R['label'], spec['B_' + R['label']][hi].sum() / spec['B_ref'][hi].sum()))

    # 5. mid-plane slices of |B|
    kmid = n_c // 2
    mag = lambda B: np.sqrt(B['Bx'] ** 2 + B['By'] ** 2 + B['Bz'] ** 2)
    fields = [('filtered fine', mag(ref)[:, :, kmid])] + [(R['label'], mag(R['C'])[:, :, kmid]) for R in runs]
    vmax = np.percentile(fields[0][1], 99.5)
    fig, axs = plt.subplots(2, len(fields), figsize=(3.1 * len(fields), 6.2))
    for i, (name, f) in enumerate(fields):
        axs[0, i].imshow(f.T, origin='lower', vmin=0, vmax=vmax, cmap='viridis'); axs[0, i].set_title(name, fontsize=9)
        d = f - fields[0][1]
        dm = max(np.abs(fields[1][1] - fields[0][1]).max(), 1e-30)
        axs[1, i].imshow(d.T, origin='lower', vmin=-dm, vmax=dm, cmap='seismic')
        axs[1, i].set_title('minus ref' if i else '', fontsize=9)
        for r_ in (0, 1): axs[r_, i].set_xticks([]); axs[r_, i].set_yticks([])
    fig.suptitle('|B| in the z mid-plane, t = %.3f' % t)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, 'slices.png'), dpi=130); plt.close(fig)

    # 6. time series of the magnetic energy: coarse diagnostics against the block-averaged fine frames
    fig, axt = plt.subplots(figsize=(6.5, 4.4))
    series = {}
    pts = []
    for d in a.fine_frames:
        for f in vc.frames(d):
            S = vc.load(f)
            Bm = [vc.block_mean(S['P'][k], r) for k in ('Bx', 'By', 'Bz')]
            full = 0.5 * S['dx'] ** 3 * sum(np.sum(S['P'][k] ** 2) for k in ('Bx', 'By', 'Bz'))
            pts.append((S['t'], 0.5 * vol * sum(np.sum(b ** 2) for b in Bm), full))
    if pts:
        pts = np.array(sorted(pts))
        axt.semilogy(pts[:, 0], pts[:, 2], 'k:', label='fine, unfiltered')
        axt.semilogy(pts[:, 0], pts[:, 1], 'ko-', ms=3, label='fine, block-averaged (reference)')
        series['fine'] = pts
    for R in runs:
        log = os.path.join(os.path.dirname(R['path']), 'run.log')
        if os.path.exists(log):
            D = read_diag(log)
            if len(D):
                axt.semilogy(D[:, 0], D[:, 1], label=R['label']); series[R['label']] = D
    axt.set_xlabel('t'); axt.set_ylabel('resolved magnetic energy'); axt.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(a.out, 'timeseries.png'), dpi=130); plt.close(fig)
    if 'fine' in series:
        pr("\nMAGNETIC ENERGY HISTORY, coarse / block-averaged fine at the fine frame times (log-interpolated)")
        pr("%-12s" % "t" + "".join("%12s" % R['label'] for R in runs))
        for tt, ef, _ in series['fine'][2:]:
            row = "%-12.4f" % tt
            for R in runs:
                D = series.get(R['label'])
                row += "%12.4f" % (np.exp(np.interp(tt, D[:, 0], np.log(D[:, 1]))) / ef) if D is not None and tt <= D[-1, 0] else "%12s" % "-"
            pr(row)

    with open(os.path.join(a.out, 'comparison.txt'), 'w') as fh:
        fh.write("\n".join(out) + "\n")
    np.savez(os.path.join(a.out, 'comparison.npz'), t=t, r=r, labels=np.array([R['label'] for R in runs]),
             **{'err_' + k: np.array(v) for k, v in err.items()}, **{'amp_' + k: np.array(v) for k, v in amp.items()},
             Emag_ref=Emag_ref, Emag=np.array(Em), **{'spec_' + k: v for k, v in spec.items()})
    print("\nWrote comparison.txt, comparison.npz and the figures to %s" % a.out)


if __name__ == '__main__':
    main()
