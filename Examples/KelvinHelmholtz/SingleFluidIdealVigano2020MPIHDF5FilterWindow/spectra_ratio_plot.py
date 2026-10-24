#!/usr/bin/env python
"""
Ratio plot of the spectra of runs against a reference, with the noise band of nominally identical runs, from the .npz written by compare_spectra.py (the table of the same numbers: spectra_ratio_table.py).

    python spectra_ratio_plot.py spectra_w7_noise.npz --out spectra_ratio.png \
        --noise baseline,n1,n2,n3 --group "c=0.2:c0.2a,c0.2b,c0.2c,c0.2d" --single c0.6 n256 [--ref native] [--kinetic] [--levels 0 1 2]

One column per E_kin stage (panel of the .npz; --levels selects panels by index), one row per ratio:
  row 1: magnetic spectrum over the reference spectrum (default: the native fine run), log-log; the noise group as a band (mean +- one standard deviation over its runs), each --group as
         a line (mean) with a lighter band (+- sd), each --single run as a line;
  row 2: the same spectra over the MEAN OF THE NOISE GROUP: the band is then the noise around 1, and a group or run outside it is an effect of the closure (or of the resolution) at that
         wavenumber; semi-log;
  row 3 (--kinetic): the kinetic spectra over the mean of the noise group.
Labels are matched exactly to the curve labels of the .npz (the labels given to compare_spectra.py; a curve missing from a panel is left out of it). A group or the noise band needs at least two
runs in the panel. The runs of a panel are at the frames of the same E_kin level, not necessarily at the same time (the match file next to the .npz lists the frames): use the levels at which
all the runs have the same frame for a clean comparison of the noise band.
"""
import argparse
import os

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from spectra_ratio_table import load_npz, read_match


def curve(panel, label):
    """(k, E_kinetic, E_magnetic) of the curve whose label is `label` (underscores stand for spaces in the keys of the .npz), or None."""
    for key, val in panel.items():
        if key.replace('_', ' ') == label.replace('_', ' '):
            return val
    return None


def find_ref(panel, ref):
    """The reference curve: 'native' -> the label containing '(native)', 'filtered' -> containing 'filtered', else a substring of the label; the first curve if there is no hit."""
    pat = {'native': '(native)', 'filtered': 'filtered'}.get(ref, ref)
    hits = [l for l in panel if pat in l.replace('_', ' ')]
    return panel[hits[0] if hits else list(panel)[0]]


def on_grid(k_from, E, k_to):
    """E at the integer wavenumbers k_to that are present in k_from (nan elsewhere)."""
    out = np.full(len(k_to), np.nan)
    pos = {int(k): i for i, k in enumerate(k_from)}
    for j, k in enumerate(k_to):
        if int(k) in pos:
            out[j] = E[pos[int(k)]]
    return out


def stack(panel, labels, idx):
    """(k, array [n_runs, n_k]) of the spectra idx (1 kinetic, 2 magnetic) of the labels present in the panel, on the integer wavenumbers of the first of them."""
    cs = [curve(panel, l) for l in labels]
    cs = [c for c in cs if c is not None]
    if not cs:
        return None, None
    k = cs[0][0]
    return k, np.array([on_grid(c[0], c[idx], k) for c in cs])


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('npz')
    p.add_argument('--out', default='spectra_ratio.png')
    p.add_argument('--noise', required=True, help='comma separated labels of the nominally identical runs (the noise group)')
    p.add_argument('--group', action='append', default=[], help='NAME:label1,label2,...  a group of runs shown by its mean and sd; repeat the flag')
    p.add_argument('--single', nargs='*', default=[], help='labels of single runs shown as lines')
    p.add_argument('--ref', default='native', help="reference curve: 'native' (default), 'filtered' or a substring of its label")
    p.add_argument('--kinetic', action='store_true', help='add a row with the kinetic spectra over the noise mean')
    p.add_argument('--levels', nargs='*', type=int, default=None, help='indices of the panels (E_kin stages) to draw; default all')
    p.add_argument('--panel-size', type=float, default=3.4)
    a = p.parse_args()

    panels = load_npz(a.npz)
    match = read_match(os.path.splitext(a.npz)[0] + '_match.txt')
    idxs = a.levels if a.levels else sorted(panels)
    noise = [s for s in a.noise.split(',') if s]
    groups = []
    for g in a.group:
        name, labs = g.split(':', 1)
        groups.append((name, [s for s in labs.split(',') if s]))
    nrow = 3 if a.kinetic else 2
    fig, axes = plt.subplots(nrow, len(idxs), figsize=(a.panel_size * len(idxs) + 1.0, 2.7 * nrow), squeeze=False, sharex='col')
    cmap = plt.get_cmap('viridis')
    colours = {}
    names = [g[0] for g in groups] + list(a.single)
    for i, n in enumerate(names):
        colours[n] = cmap(0.15 + 0.7 * i / max(1, len(names) - 1))
    for ci, pi in enumerate(idxs):
        panel = panels[pi]
        title = "$E_{kin}$ = %g" % match[pi][0] if pi < len(match) else "panel %d" % pi
        ref = find_ref(panel, a.ref)
        kr, Er = ref[0], ref[2]
        kn, Mn = stack(panel, noise, 2)
        if Mn is None or len(Mn) < 2:
            axes[0][ci].set_title(title + " (noise group missing)")
            continue
        mean_n, sd_n = np.nanmean(Mn, axis=0), np.nanstd(Mn, axis=0, ddof=1)
        Rr = on_grid(kr, Er, kn)
        ax = axes[0][ci]
        ax.fill_between(kn, (mean_n - sd_n) / Rr, (mean_n + sd_n) / Rr, color='0.7', alpha=0.8, label='noise, $\\pm$ sd (n = %d)' % len(Mn), lw=0)
        ax.plot(kn, mean_n / Rr, color='0.35', lw=1.0)
        ax2 = axes[1][ci]
        ax2.axhline(1.0, color='0.35', lw=0.8)
        ax2.fill_between(kn, 1 - sd_n / mean_n, 1 + sd_n / mean_n, color='0.7', alpha=0.8, lw=0)
        if a.kinetic:
            _, Kn = stack(panel, noise, 1)
            mean_k, sd_k = np.nanmean(Kn, axis=0), np.nanstd(Kn, axis=0, ddof=1)
            ax3 = axes[2][ci]
            ax3.axhline(1.0, color='0.35', lw=0.8)
            ax3.fill_between(kn, 1 - sd_k / mean_k, 1 + sd_k / mean_k, color='0.7', alpha=0.8, lw=0)
        for name, labs in groups:
            k, M = stack(panel, labs, 2)
            if M is None:
                continue
            m = np.nanmean(M, axis=0)
            c = colours[name]
            if len(M) > 1:
                s = np.nanstd(M, axis=0, ddof=1)
                ax.fill_between(k, (m - s) / on_grid(kr, Er, k), (m + s) / on_grid(kr, Er, k), color=c, alpha=0.25, lw=0)
                ax2.fill_between(k, (m - s) / mean_n, (m + s) / mean_n, color=c, alpha=0.25, lw=0)
            ax.plot(k, m / on_grid(kr, Er, k), color=c, lw=1.6, label="%s (n = %d)" % (name, len(M)))
            ax2.plot(k, m / mean_n, color=c, lw=1.6)
            if a.kinetic:
                _, K = stack(panel, labs, 1)
                axes[2][ci].plot(k, np.nanmean(K, axis=0) / mean_k, color=c, lw=1.6)
        for name in a.single:
            c = curve(panel, name)
            if c is None:
                continue
            col = colours[name]
            ax.plot(c[0], c[2] / on_grid(kr, Er, c[0]), color=col, lw=1.6, ls='--', label=name)
            ax2.plot(c[0], on_grid(c[0], c[2], kn) / mean_n, color=col, lw=1.6, ls='--')
            if a.kinetic:
                axes[2][ci].plot(kn, on_grid(c[0], c[1], kn) / mean_k, color=col, lw=1.6, ls='--')
        ax.set_title(title)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax2.set_xscale('log')
        for row in range(nrow):
            axes[row][ci].grid(alpha=0.25)
        axes[-1][ci].set_xlabel('k (units of $2\\pi/L$)')
    axes[0][0].set_ylabel('$E_m$ / reference')
    axes[1][0].set_ylabel('$E_m$ / noise mean')
    if a.kinetic:
        axes[2][0].set_ylabel('$E_k$ / noise mean')
    axes[0][0].legend(fontsize=7, loc='lower left')
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print('wrote', a.out)


if __name__ == '__main__':
    main()
