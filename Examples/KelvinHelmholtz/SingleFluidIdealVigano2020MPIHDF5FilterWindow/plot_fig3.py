#!/usr/bin/env python
"""
Recreation of Fig. 3 of Vigano et al. (arXiv:2004.00870): density (top row) and the z-component of the magnetic field
B^z (bottom row) in the x-y plane, one column per snapshot, for the 3D KHI of the Vigano et al. (2020) initial data.
Standard colour map viridis, diverging colour map seismic (B^z, symmetric about zero). The paper does not say which z
plane it shows; the default is the middle one.

    python plot_fig3.py data_window_00.hdf5 [more snapshots ...] --labels "N=128" --out fig3.png

Colour limits: --rho-lim and --bz-lim (the paper's are 0.8-2.7 for rho with rho_0=1.5, rho_1=0.5, and +-0.02 for B^z); with
--bz-lim auto the B^z limits are the 99.5th percentile of |B^z| over all panels. Writes <out> and <out stem>_autoscale.png
when the limits were fixed, so both can be compared.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from viz_common import axes_of, load


def draw(snaps, labels, zi, rho_lim, bz_lim, out):
    n = len(snaps)
    fig, axes = plt.subplots(2, n, figsize=(4.3 * n + 0.8, 8.4), squeeze=False, constrained_layout=True)
    S0 = snaps[0]
    x, y, z = axes_of(S0)
    ext = [x[0] - S0['dx'] / 2, x[-1] + S0['dx'] / 2, y[0] - S0['dx'] / 2, y[-1] + S0['dx'] / 2]
    for c, (S, lab) in enumerate(zip(snaps, labels)):
        for r, (name, cmap, lim) in enumerate((('rho', 'viridis', rho_lim), ('Bz', 'seismic', bz_lim))):
            f = S['P'][name][:, :, zi].T                               # rows = y, columns = x
            if r == 1:
                vmin, vmax = -lim, lim
            else:
                vmin, vmax = lim
            im = axes[r, c].imshow(f, origin='lower', extent=ext, cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
            axes[r, c].set_xlabel('x')
            if c == 0:
                axes[r, c].set_ylabel('y')
            if r == 0:
                axes[r, c].set_title('%s, t = %.3f, z = %.3f' % (lab, S['t'], z[zi]), fontsize=10)
            if c == n - 1:
                fig.colorbar(im, ax=axes[r, :].tolist(), shrink=0.85, label=('density $\\rho$' if r == 0 else '$B^z$'))
    fig.savefig(out, dpi=140)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("files", nargs="+")
    p.add_argument("--labels", nargs="*", default=None)
    p.add_argument("--z-index", type=int, default=None)
    p.add_argument("--rho-lim", type=float, nargs=2, default=None, metavar=("MIN", "MAX"))
    p.add_argument("--bz-lim", default="0.02", help="symmetric limit for B^z, or 'auto'")
    p.add_argument("--out", default="fig3.png")
    a = p.parse_args()
    snaps = [load(f) for f in a.files]
    labels = a.labels if a.labels and len(a.labels) == len(snaps) else [os.path.basename(f) for f in a.files]
    zi = a.z_index if a.z_index is not None else snaps[0]['P']['rho'].shape[2] // 2
    rho_all = np.concatenate([s['P']['rho'].ravel() for s in snaps])
    rho_lim = tuple(a.rho_lim) if a.rho_lim else (float(rho_all.min()), float(rho_all.max()))
    bz_auto = float(np.percentile(np.abs(np.concatenate([s['P']['Bz'][:, :, zi].ravel() for s in snaps])), 99.5))
    print("rho range %.3f - %.3f; max |Bz| on the plane %.3g (99.5th percentile %.3g)" % (
        rho_all.min(), rho_all.max(), max(np.abs(s['P']['Bz'][:, :, zi]).max() for s in snaps), bz_auto))
    if a.bz_lim == 'auto':
        draw(snaps, labels, zi, rho_lim, bz_auto, a.out)
    else:
        draw(snaps, labels, zi, rho_lim, float(a.bz_lim), a.out)
        stem, ext = os.path.splitext(a.out)
        draw(snaps, labels, zi, rho_lim, bz_auto, stem + "_autoscale" + ext)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
