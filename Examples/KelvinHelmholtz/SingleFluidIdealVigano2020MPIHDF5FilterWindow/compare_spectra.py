#!/usr/bin/env python
"""
Kinetic and magnetic spectra of several runs at equal development, in the layout of Figs. 5 and 6 of Vigano et al., arXiv:2004.00870, with the conventions of plot_quantities.py.

    python compare_spectra.py --ekin 0.4 0.2 0.1 0.05 --out spectra_stages.png \
        baseline=DIR0 c0.6=DIR1 w7=DIR2 [...] \
        [--reference DIR384 --reference-log LOG384 --block 3 --reference-label "384^3"] [--reference-mode native|filtered|both]
    python compare_spectra.py --layout groups --ekin 0.1 --group "WENO5=baseline,c0.6" --group "WENO7=w7,w7c0.6" --out spectra_groups.png  baseline=... c0.6=... w7=... w7c0.6=...

Layout `stages` (default, as Fig. 5): one figure, one panel per E_kin level (left to right, in the order given), every run in every panel; layout `groups` (as Fig. 6): one figure for the first
--ekin level, one panel per `--group NAME=label1,label2` (repeat the flag), the reference in every panel. Kinetic spectra dash-dot, magnetic dashed (the dash pattern of each curve is shifted by a fraction of its period, so curves that lie on top of each
other show their colours in turn), one colour per run (viridis; the reference black and grey); thin solid reference slopes k^-5/3 (black) and k^3/2 (dark grey), anchored to the baseline's spectra. All panels of a row share the y-axis (--no-sharey for independent axes).

Convention of plot_quantities.py (viz_common.spectrum), not the one of the paper: E(k) = sum over the integer shell k of |f_k|^2 / 2 with f the FFT divided by N^3, k in units of 2 pi / L, for the
velocity v (`--kinetic sqrtrho` uses sqrt(rho) v as the paper does) and for B, from the primitives. The paper averages 4 pi k^2 |f|^2 over the shell instead of summing (a power density per unit wave
number): the two differ by the number of modes in the shell, so slopes are not comparable, only the shape within one convention.

Stages are defined by E_kin, not time (the runs are at different stages at equal time): each run's frame is the one whose run-log E_kin is closest to the level (first crossing, then the nearest
frame; no frame within --ekin-tol, default 0.1, leaves the run out of that panel), as compare_slices.py. The reference (filtered fine run): `native` its spectrum on its own grid, out to its own
Nyquist; `filtered` after block-averaging its primitives by --block onto the grid of the runs (the top-hat is a poor filter near the Nyquist wavenumber of the coarse grid); `both` (default).
`--compensate` makes two rows, the kinetic spectrum times k^(5/3) above and the magnetic one times k^(-3/2) below (a plateau is the slope; no guide lines are drawn); the `.npz` always holds the raw spectra.
Writes the figure and, next to it, `<out>.npz` with the spectra and `<out>_match.txt` with the frames used.
"""
import argparse
import os
import sys

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import viz_common as vc
from compare_scan import frame_times
from compare_slices import find_frame
from time_shift import read


def spectra_of(path, block, kinetic):
    """(k, E_kinetic, E_magnetic) of one frame, the primitives block-averaged by `block` first if block > 1."""
    S = vc.load(path)
    P = S['P']
    if block > 1:
        P = {k: vc.block_mean(v, block) for k, v in P.items()}
    w = np.sqrt(P['rho']) if kinetic == 'sqrtrho' else 1.0
    k, Ek = vc.spectrum(w * P['vx'], w * P['vy'], w * P['vz'], S['dx'])
    k, Em = vc.spectrum(P['Bx'], P['By'], P['Bz'], S['dx'])
    return k, Ek, Em, S['t']


KIN_PATTERN = (6.0, 2.0, 1.0, 2.0)          # kinetic: dash-dot
MAG_PATTERN = (5.0, 3.0)                    # magnetic: dashed


def dashes(pattern, i, n):
    """Line style with the dash pattern shifted by i/n of its period: curves that lie on top of each other then show their colours in turn."""
    return (i * sum(pattern) / max(n, 1), pattern)


def draw_panel(ax, curves, title, slopes, first, which='both', compensate=False):
    """which: 'both' (kinetic dash-dot, magnetic dashed on one axis), 'kin' or 'mag' (one row of the compensated layout). compensate: multiply by k^(5/3) (kinetic) or k^(-3/2) (magnetic).
    Reference slopes are thin solid black lines."""
    n = len(curves)
    for i, c in enumerate(curves):
        k = c['k'][1:]
        if which in ('both', 'kin'):
            ax.loglog(k, c['Ek'][1:] * (k ** (5.0 / 3.0) if compensate else 1.0), color=c['color'], ls=dashes(KIN_PATTERN, i, n), lw=c['lw'], label=c['label'])
        if which in ('both', 'mag'):
            ax.loglog(k, c['Em'][1:] * (k ** -1.5 if compensate else 1.0), color=c['color'], ls=dashes(MAG_PATTERN, i, n), lw=c['lw'], label=c['label'] if which == 'mag' else None)
    if slopes is not None and not compensate:
        k0, ek0, em0 = slopes
        kk = np.array([k0, 6.0 * k0])
        ax.loglog(kk, ek0 * (kk / k0) ** (-5.0 / 3.0), 'k-', lw=0.8, label='$k^{-5/3}$' if first else None)
        ax.loglog(kk, em0 * (kk / k0) ** 1.5 * 0.35, '-', color='0.25', lw=0.8, label='$k^{3/2}$' if first else None)
    if title:
        ax.set_title(title, fontsize=10)
    ax.set_xlabel('k (units of $2\\pi/L$)')
    ax.grid(True, which='major', lw=0.3, alpha=0.4)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('runs', nargs='+', help='LABEL=DIR (frames and run.log); the first is the baseline')
    p.add_argument('--ekin', nargs='+', type=float, required=True)
    p.add_argument('--out', default='spectra.png')
    p.add_argument('--layout', choices=['stages', 'groups'], default='stages')
    p.add_argument('--group', action='append', default=[], help='NAME=label1,label2: one panel of the groups layout; repeat the flag for each panel (put it anywhere)')
    p.add_argument('--groups', nargs='*', default=[], help='as --group, several in a row; must come AFTER the runs (it takes every following argument)')
    p.add_argument('--reference')
    p.add_argument('--reference-log')
    p.add_argument('--block', type=int, default=1)
    p.add_argument('--reference-label', default='filtered fine')
    p.add_argument('--reference-mode', choices=['native', 'filtered', 'both'], default='both')
    p.add_argument('--kinetic', choices=['v', 'sqrtrho'], default='v')
    p.add_argument('--ekin-tol', type=float, default=0.1)
    p.add_argument('--no-sharey', dest='sharey', action='store_false', help='independent y-axes (default: one y-axis for all panels of a row)')
    p.add_argument('--compensate', action='store_true', help='two rows: kinetic x k^(5/3) above, magnetic x k^(-3/2) below (a plateau is the slope)')
    p.add_argument('--panel-size', type=float, default=3.4, help='inches per panel')
    p.add_argument('--dpi', type=int, default=130)
    a = p.parse_args()

    runs = []
    for item in a.runs:
        label, d = item.split('=', 1)
        lg = os.path.join(d, 'run.log')
        if not os.path.exists(lg):
            sys.exit("%s needs %s" % (label, lg))
        runs.append((label, frame_times(d), read(lg)))
    ref = None
    if a.reference:
        if not a.reference_log:
            sys.exit("--reference needs --reference-log")
        ref = (a.reference_label, frame_times(a.reference), read(a.reference_log))
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, max(len(runs), 2)))
    cache, report = {}, []

    def curves_at(level, labels):
        curves = []
        report.append("E_kin = %g" % level)
        if ref is not None:
            hit = find_frame(ref[1], ref[2], level, a.ekin_tol)
            if hit is None:
                report.append("  %-18s no frame within %.0f %%" % (ref[0], 100 * a.ekin_tol))
            else:
                t, path, mis = hit
                report.append("  %-18s t = %8.4f  (%s)  E_kin mismatch %+.1f %%" % (ref[0], t, os.path.basename(path), 100 * mis))
                if a.reference_mode in ('native', 'both'):
                    key = (path, 1)
                    cache[key] = cache.get(key) or spectra_of(path, 1, a.kinetic)
                    k, Ek, Em, _ = cache[key]
                    curves.append(dict(k=k, Ek=Ek, Em=Em, color='k', lw=1.6, label='%s%s' % (ref[0], ' (native)' if a.reference_mode == 'both' else '')))
                if a.reference_mode in ('filtered', 'both') and a.block > 1:
                    key = (path, a.block)
                    cache[key] = cache.get(key) or spectra_of(path, a.block, a.kinetic)
                    k, Ek, Em, _ = cache[key]
                    curves.append(dict(k=k, Ek=Ek, Em=Em, color='0.45', lw=1.2, label='%s, filtered' % ref[0]))
        for i, (label, frames, D) in enumerate(runs):
            if labels is not None and label not in labels:
                continue
            hit = find_frame(frames, D, level, a.ekin_tol)
            if hit is None:
                report.append("  %-18s no frame within %.0f %%" % (label, 100 * a.ekin_tol))
                continue
            t, path, mis = hit
            report.append("  %-18s t = %8.4f  (%s)  E_kin mismatch %+.1f %%" % (label, t, os.path.basename(path), 100 * mis))
            key = (path, 1)
            cache[key] = cache.get(key) or spectra_of(path, 1, a.kinetic)
            k, Ek, Em, _ = cache[key]
            curves.append(dict(k=k, Ek=Ek, Em=Em, color=colors[i], lw=1.3, label='%s' % label if a.layout == 'stages' else '%s  $t$=%.1f' % (label, t)))
        return curves

    def slope_anchor(curves):
        for c in curves:
            if not isinstance(c['color'], str):                       # the first run (a numpy colour), not the reference
                k0 = 4.0
                return k0, float(c['Ek'][4]) * 1.8, float(c['Em'][4])
        return None

    if a.layout == 'stages':
        panels = [("$E_{kin}$ = %g" % lv, curves_at(lv, None)) for lv in a.ekin]
    else:
        allgroups = a.group + a.groups
        if not allgroups:
            sys.exit("--layout groups needs --group NAME=label1,label2 (repeat for each panel)")
        panels = []
        for g in allgroups:
            name, labs = g.split('=', 1)
            panels.append((name, curves_at(a.ekin[0], labs.split(','))))
    n = len(panels)
    from matplotlib.lines import Line2D
    styles = [Line2D([], [], color='0.3', ls=(0, KIN_PATTERN), lw=1.3), Line2D([], [], color='0.3', ls=(0, MAG_PATTERN), lw=1.3)]
    if a.compensate:
        fig, axes = plt.subplots(2, n, figsize=(a.panel_size * n, a.panel_size * 1.9), sharey='row' if a.sharey else False, sharex='col', squeeze=False,
                                 gridspec_kw=dict(wspace=0.04 if a.sharey else 0.3, hspace=0.06))
        for i, (title, curves) in enumerate(panels):
            draw_panel(axes[0][i], curves, title, None, i == 0, 'kin', True)
            draw_panel(axes[1][i], curves, None, None, i == 0, 'mag', True)
            axes[0][i].set_xlabel('')
            if a.sharey and i > 0:
                for r in (0, 1):
                    axes[r][i].tick_params(labelleft=False)
        axes[0][0].set_ylabel('$k^{5/3}E_K(k)$ (kinetic)')
        axes[1][0].set_ylabel('$k^{-3/2}E_M(k)$ (magnetic)')
        if a.layout == 'stages':
            h, l = axes[0][0].get_legend_handles_labels()
            fig.legend(h, l, fontsize=8, loc='upper center', ncol=min(len(h), 6), bbox_to_anchor=(0.5, 0.0), frameon=False)
        else:
            for i in range(n):
                axes[0][i].legend(fontsize=7, loc='lower left')
    else:
        fig, axes = plt.subplots(1, n, figsize=(a.panel_size * n, a.panel_size * 1.05), sharey=a.sharey, squeeze=False, gridspec_kw=dict(wspace=0.04 if a.sharey else 0.22))
        for i, (ax, (title, curves)) in enumerate(zip(axes[0], panels)):
            draw_panel(ax, curves, title, slope_anchor(curves), i == 0)
            if i == 0:
                ax.set_ylabel('$E(k)$')
            if a.sharey and i > 0:
                ax.tick_params(labelleft=False)
        if a.layout == 'stages':                       # one legend below the panels (the labels carry the frame times of the first panel; see the match file for the others)
            h, l = axes[0][0].get_legend_handles_labels()
            fig.legend(h + styles, l + ['kinetic', 'magnetic'], fontsize=8, loc='upper center', ncol=min(len(h) + 2, 6), bbox_to_anchor=(0.5, -0.04), frameon=False)
        else:
            for ax in axes[0]:
                h, l = ax.get_legend_handles_labels()
                ax.legend(h + styles, l + ['kinetic', 'magnetic'], fontsize=7, loc='lower left')
    fig.savefig(a.out, dpi=a.dpi, bbox_inches='tight')
    plt.close(fig)
    base = os.path.splitext(a.out)[0]
    save = {}
    for pi, (title, curves) in enumerate(panels):
        for c in curves:
            tag = '%d_%s' % (pi, c['label'].replace(' ', '_'))
            save[tag + '_k'], save[tag + '_Ek'], save[tag + '_Em'] = c['k'], c['Ek'], c['Em']
    np.savez(base + '.npz', **save)
    with open(base + '_match.txt', 'w') as fh:
        fh.write("\n".join(report) + "\n")
    print("wrote", a.out)
    print("\n".join(report))


if __name__ == '__main__':
    main()
