#!/usr/bin/env python
"""
Slices of coarse runs side by side at equal development: the baseline against closure runs (and, optionally, the filtered fine run), only resolved quantities.

    python compare_slices.py --ekin 0.4 0.2 0.1 --outdir slices \
        baseline=DIR0 c1/8=DIR1 c0.6=DIR2 [w7=DIR3 ...] \
        [--reference DIR384 --reference-log LOG384 --block 3 --reference-label "384^3, filtered"] \
        [--vars rho,wz,Bz,Bmag,Jmag] [--plane xy|xz|yz] [--index N] [--crop x0 x1 y0 y1]

One figure per stage, `slices_ekin_<level>.png`: one row per variable, one column per run (the reference first), panels touching (hairline gaps), one thin colour bar per row, the SAME
colour scale in every panel of a row (the range is set from all the panels of that stage together). The run name and frame time are above each column, the variable name beside its colour bar on the right. Standard colour map viridis, symmetric seismic for signed quantities (log10 values with a linear colour bar for the strictly positive quantities that span decades).

Stages are defined by E_kin, not by time: the runs are at different stages at equal time (delays of 2-5 time units, differing from run to run), and E_kin is a clock that the
hydrodynamics fixes (design doc sections 14-18). For each level the frame of each run whose run-log E_kin is closest is used (the first time the run's E_kin falls through the level, then
the nearest frame); a run with no frame within --ekin-tol (fractional, default 0.1) gets an empty panel. Needs `run.log` (the `diag` lines) next to the frames of each run and, for the
reference, --reference-log. The reference frames are block-averaged (primitives) by --block to the coarse grid. Pointwise differences are deliberately not drawn: the eddies of runs that
differ by noise are not in the same place. Writes `slices_match.txt` with the frames used and their E_kin mismatch.

Variables (--vars, comma separated): rho, p, vmag, W, wz (vorticity z), wy, Bx, By, Bz, Bmag, beta (log), Jmag (log), Emag (log, |v x B|), emag (magnetic energy density), ekin (kinetic energy density).
Defaults: rho,wz,Bz,Bmag,Jmag. Everything is computed from the primitives on the grid of the run.
"""
import argparse
import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import viz_common as vc
from compare_scan import frame_times
from time_shift import crossing, read

# name -> (label, getter(P, D), colour map, scale: lin | sym | log)
VARIABLES = {
    'rho': ('$\\rho$', lambda P, D: P['rho'], 'viridis', 'lin'),
    'p': ('$p$', lambda P, D: P['p'], 'viridis', 'lin'),
    'vmag': ('$|v|$', lambda P, D: D['vmag'], 'viridis', 'lin'),
    'W': ('$W$', lambda P, D: D['W'], 'viridis', 'lin'),
    'wz': ('$\\omega_z$', lambda P, D: D['wz'], 'seismic', 'sym'),
    'wy': ('$\\omega_y$', lambda P, D: D['wy'], 'seismic', 'sym'),
    'Bx': ('$B_x$', lambda P, D: P['Bx'], 'seismic', 'sym'),
    'By': ('$B_y$', lambda P, D: P['By'], 'seismic', 'sym'),
    'Bz': ('$B_z$', lambda P, D: P['Bz'], 'seismic', 'sym'),
    'Bmag': ('$|B|$', lambda P, D: D['Bmag'], 'viridis', 'lin'),
    'beta': ('$\\beta$', lambda P, D: D['beta'], 'viridis', 'log'),
    'Jmag': ('$|J|$', lambda P, D: D['Jmag'], 'viridis', 'log'),
    'Emag': ('$|v\\times B|$', lambda P, D: D['Emag'], 'viridis', 'log'),
    'emag': ('$B^2/2$', lambda P, D: D['emag'], 'viridis', 'lin'),
    'ekin': ('$\\rho W(W-1)$', lambda P, D: D['ekin'], 'viridis', 'lin'),
}


def find_frame(frames, D, level, tol):
    """(time, path, mismatch) of the frame whose run-log E_kin is closest to `level` (first crossing, then nearest frame), or None."""
    tc = crossing(D[:, 0], D[:, 2], level, True)
    if tc is None or not frames:
        return None
    t = min(frames, key=lambda x: abs(x - tc))
    mis = np.interp(t, D[:, 0], D[:, 2]) / level - 1.0
    return None if abs(mis) > tol else (t, frames[t], mis)


def slice_fields(path, names, plane, index, block):
    """2D slices of the requested variables from one frame; the frame is block-averaged by `block` first if block > 1."""
    S = vc.load(path)
    if block > 1:
        S['P'] = {k: vc.block_mean(v, block) for k, v in S['P'].items()}
        S['dx'] = S['dx'] * block
    D = vc.derived(S, subfilter=False)
    n = S['P']['rho'].shape
    ax = {'xy': 2, 'xz': 1, 'yz': 0}[plane]
    i = n[ax] // 2 if index is None else index
    take = {'xy': lambda f: f[:, :, i].T, 'xz': lambda f: f[:, i, :].T, 'yz': lambda f: f[i, :, :].T}[plane]
    out = {nm: take(VARIABLES[nm][1](S['P'], D)).copy() for nm in names}
    coords = vc.axes_of(S)
    ext = {'xy': (0, 1), 'xz': (0, 2), 'yz': (1, 2)}[plane]
    extent = [S['lo'][ext[0]], S['hi'][ext[0]], S['lo'][ext[1]], S['hi'][ext[1]]]
    return out, extent, S['t'], coords[ax][i]


def norm_for(mode, panels):
    allv = np.concatenate([p.ravel() for p in panels])
    if mode == 'sym':
        lim = float(np.percentile(np.abs(allv), 99.5)) or 1.0
        return dict(vmin=-lim, vmax=lim)
    if mode == 'log':
        pos = np.log10(allv[allv > 0])
        lo, hi = (np.percentile(pos, 1), np.percentile(pos, 99.8)) if pos.size else (-12.0, 0.0)
        return dict(vmin=float(lo), vmax=float(hi))
    return dict(vmin=float(np.percentile(allv, 0.2)), vmax=float(np.percentile(allv, 99.8)))


def draw_stage(level, columns, names, plane, outdir, crop, dpi, size):
    """columns: list of (label, None or (slices, extent, t, mis)); writes slices_ekin_<level>.png.
    Run titles above the top row, the variable name beside its colour bar on the right, hairline gaps between the panels."""
    nr, nc = len(names), len(columns)
    head, foot, cb_gap, cb_w, right = 0.32, 0.2, 0.06, 0.12, 1.05           # inches: title strip, footer strip, gap to the colour bar, its width, tick labels + variable name
    width, height = size * nc + cb_gap + cb_w + right, size * nr + head + foot
    fig = plt.figure(figsize=(width, height))
    gs = fig.add_gridspec(nr, nc, wspace=0.012, hspace=0.012, left=0.0, right=size * nc / width, top=1.0 - head / height, bottom=foot / height)
    axes = [[None] * nc for _ in range(nr)]
    for r, nm in enumerate(names):
        label, cmap, mode = VARIABLES[nm][0], VARIABLES[nm][2], VARIABLES[nm][3]
        if mode == 'log':
            label = '$\\log_{10}$ ' + label                                   # log quantities are drawn as log10 values with a linear colour bar
        panels = [c[1][0][nm] for c in columns if c[1] is not None]
        kw = norm_for(mode, panels) if panels else {}
        im = None
        for c, (clabel, data) in enumerate(columns):
            ax = fig.add_subplot(gs[r, c])
            axes[r][c] = ax
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_linewidth(0.4)
            if data is None:
                ax.set_facecolor('0.85')
                ax.text(0.5, 0.5, 'no frame', transform=ax.transAxes, ha='center', va='center', fontsize=8, color='0.3')
            else:
                sl, extent, t, mis = data
                field = np.log10(np.maximum(sl[nm], 1e-300)) if mode == 'log' else sl[nm]
                im = ax.imshow(field, origin='lower', extent=extent, cmap=cmap, interpolation='nearest', aspect='equal', **kw)
                if crop:
                    ax.set_xlim(crop[0], crop[1]); ax.set_ylim(crop[2], crop[3])
            if r == 0:
                ax.set_title(clabel + ('' if data is None else ',  $t$ = %.2f' % data[2]), fontsize=10, pad=3)
        pos = axes[r][nc - 1].get_position()
        cax = fig.add_axes([pos.x1 + cb_gap / width, pos.y0, cb_w / width, pos.height])
        if im is not None:
            cb = fig.colorbar(im, cax=cax)
            cb.ax.tick_params(labelsize=7, length=2, pad=1)
            cb.ax.yaxis.set_label_position('right')
            cb.set_label(label, rotation=270, labelpad=14, fontsize=11)
        else:
            cax.axis('off')
    fig.text(0.002, 0.003, '$E_{kin}$ = %g,  %s plane' % (level, plane), ha='left', va='bottom', fontsize=7, color='0.35')
    out = os.path.join(outdir, 'slices_ekin_%g.png' % level)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('runs', nargs='+', help='LABEL=DIR (frames and run.log); the first is the baseline')
    p.add_argument('--ekin', nargs='+', type=float, required=True, help='E_kin levels (one figure each)')
    p.add_argument('--outdir', default='slices')
    p.add_argument('--vars', default='rho,wz,Bz,Bmag,Jmag')
    p.add_argument('--plane', choices=['xy', 'xz', 'yz'], default='xy')
    p.add_argument('--index', type=int, help='slice index along the third axis of the run grid (default: the middle); the reference uses the same position')
    p.add_argument('--crop', nargs=4, type=float, metavar=('X0', 'X1', 'Y0', 'Y1'), help='zoom to this part of the plane (domain coordinates)')
    p.add_argument('--reference', help='directory of frames of the fine run')
    p.add_argument('--reference-log', help='run log of the fine run')
    p.add_argument('--block', type=int, default=1, help='block size that brings the fine frames to the grid of the runs')
    p.add_argument('--reference-label', default='filtered fine')
    p.add_argument('--ekin-tol', type=float, default=0.1)
    p.add_argument('--dpi', type=int, default=130)
    p.add_argument('--panel-size', type=float, default=2.2, help='inches per panel')
    a = p.parse_args()
    names = a.vars.split(',')
    for nm in names:
        if nm not in VARIABLES:
            sys.exit("unknown variable %s (known: %s)" % (nm, ", ".join(VARIABLES)))
    os.makedirs(a.outdir, exist_ok=True)

    runs = []
    for item in a.runs:
        label, d = item.split('=', 1)
        lg = os.path.join(d, 'run.log')
        if not os.path.exists(lg):
            sys.exit("%s needs %s" % (label, lg))
        runs.append((label, frame_times(d), read(lg), 1))
    if a.reference:
        if not a.reference_log:
            sys.exit("--reference needs --reference-log")
        runs.insert(0, (a.reference_label, frame_times(a.reference), read(a.reference_log), a.block))

    report = []
    for level in a.ekin:
        columns = []
        report.append("E_kin = %g" % level)
        index_ref = None
        for label, frames, D, block in runs:
            hit = find_frame(frames, D, level, a.ekin_tol)
            if hit is None:
                columns.append((label, None))
                report.append("  %-18s no frame within %.0f %%" % (label, 100 * a.ekin_tol))
                continue
            t, path, mis = hit
            # the same physical position for every run: the index of the first run's grid, scaled for a reference on a finer grid
            idx = a.index
            sl, extent, tt, pos = slice_fields(path, names, a.plane, idx, block)
            columns.append((label, (sl, extent, t, mis)))
            report.append("  %-18s t = %8.4f  (%s)  E_kin mismatch %+.1f %%" % (label, t, os.path.basename(path), 100 * mis))
        out = draw_stage(level, columns, names, a.plane, a.outdir, a.crop, a.dpi, a.panel_size)
        print("wrote", out)
    with open(os.path.join(a.outdir, 'slices_match.txt'), 'w') as fh:
        fh.write("\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == '__main__':
    main()
