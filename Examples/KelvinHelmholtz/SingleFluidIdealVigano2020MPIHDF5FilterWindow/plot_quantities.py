#!/usr/bin/env python
"""
Snapshot figures of the quantities used in the filtering work, for one run of the Vigano et al. (2020) 3D KHI.

    python plot_quantities.py RUN_DIR [CONTINUATION_DIR ...] --snapshot data_window_00.hdf5 --outdir figures --filter-cells 4

RUN_DIR holds the checkpoint frames (frame.checkpoint.<t>.hdf5, written every `safety` iterations) and the window
snapshots. Writes, to --outdir:
  maps_xy.png, maps_xz.png   12 maps each of the final snapshot in the x-y plane (z middle) and the x-z plane (at the
                             shear layer y = 1/4): rho, p, |v|, Lorentz factor W, vorticity w_z (w_y), B_x, B_z,
                             |B|, plasma beta (log), |J| = |curl B| (log), |E| = |v x B| (log), and the sub-filter
                             magnetic energy <B^2> - |<B>|^2 for a top-hat of --filter-cells cells (log)
  energy_history.png         total magnetic energy (1/2 int B^2) and kinetic energy (int rho W (W-1)) against time, from
                             the frames (compare Fig. 4 of Vigano et al., arXiv:2004.00870)
  spectra.png                kinetic (v) and magnetic (B) shell spectra of the final snapshot, and of three earlier frames
  profiles.png               x-z averages of rho, v_x, |B|^2 and |J| against y at several times
  jcurrent_pdf.png           distribution of |J| and of the sub-filter magnetic energy, and where the 5 % largest |J|
                             values are in y
Standard colour map viridis; diverging colour map seismic for signed quantities.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from viz_common import axes_of, derived, frames, load, spectrum


def load_time(path):
    import h5py
    with h5py.File(path, 'r') as f:
        return float(np.ravel(f.attrs['t'])[0]) if 't' in f.attrs else float(np.ravel(f['Domain'].attrs['endTime'])[0])


def field_list(S, D):
    P = S['P']
    return [('density $\\rho$', P['rho'], 'viridis', 'lin'), ('pressure $p$', P['p'], 'viridis', 'lin'),
            ('speed $|v|$', D['vmag'], 'viridis', 'lin'), ('Lorentz factor $W$', D['W'], 'viridis', 'lin'),
            ('vorticity $\\omega_{n}$', None, 'seismic', 'sym'), ('$B_x$', P['Bx'], 'seismic', 'sym'),
            ('$B_z$', P['Bz'], 'seismic', 'sym'), ('$|B|$', D['Bmag'], 'viridis', 'lin'),
            ('plasma $\\beta = p/(B^2/2)$', D['beta'], 'viridis', 'log'), ('$|J| = |\\nabla\\times B|$', D['Jmag'], 'viridis', 'log'),
            ('$|E| = |v\\times B|$', D['Emag'], 'viridis', 'log'), ('sub-filter $\\langle B^2\\rangle-|\\langle B\\rangle|^2$', D['Bsub'], 'viridis', 'log')]


def maps(S, D, plane, index, out, filter_cells):
    x, y, z = axes_of(S)
    dx = S['dx']
    if plane == 'xy':
        take, a1, a2, labs, ort = (lambda f: f[:, :, index].T), x, y, ('x', 'y'), 'wz'
        title = "x-y plane, z = %.3f" % z[index]
    else:
        take, a1, a2, labs, ort = (lambda f: f[:, index, :].T), x, z, ('x', 'z'), 'wy'
        title = "x-z plane, y = %.3f" % y[index]
    ext = [a1[0] - dx / 2, a1[-1] + dx / 2, a2[0] - dx / 2, a2[-1] + dx / 2]
    fig, axes = plt.subplots(3, 4, figsize=(19, 14.5 if plane == 'xy' else 10.5), constrained_layout=True)
    for ax, (name, f, cmap, mode) in zip(axes.ravel(), field_list(S, D)):
        if f is None:
            f = D[ort]
            name = name.replace('{n}', '{%s}' % ort[1])
        im2 = take(f)
        if mode == 'sym':
            lim = float(np.percentile(np.abs(im2), 99.5)) or 1.0
            kw = dict(vmin=-lim, vmax=lim)
        elif mode == 'log':
            pos = im2[im2 > 0]
            lo, hi = (np.percentile(pos, 1), np.percentile(pos, 99.8)) if pos.size else (1e-12, 1.0)
            kw = dict(norm=LogNorm(vmin=lo, vmax=hi))
        else:
            kw = dict(vmin=float(np.percentile(im2, 0.2)), vmax=float(np.percentile(im2, 99.8)))
        im = ax.imshow(im2, origin='lower', extent=ext, cmap=cmap, interpolation='nearest', aspect='equal', **kw)
        ax.set_title(name, fontsize=10)
        ax.set_xlabel(labs[0]); ax.set_ylabel(labs[1])
        fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle("%s, t = %.3f, %s (sub-filter energy: top-hat of %d cells)" % (os.path.basename(S['path']), S['t'], title, filter_cells))
    fig.savefig(out, dpi=110)
    plt.close(fig)


def energy_history(files, out):
    t, em, ek = [], [], []
    for f in files:
        S = load(f)
        D = derived(S, 2)
        vol = S['dx'] ** 3
        t.append(S['t']); em.append(float(D['emag'].sum() * vol)); ek.append(float(D['ekin'].sum() * vol))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    ax[0].semilogy(t, em, 'o-'); ax[0].set_xlabel('time'); ax[0].set_ylabel('magnetic energy $\\frac{1}{2}\\int B^2$')
    ax[0].set_title('magnetic energy (initial $B_0^2/2 = 5\\times10^{-7}$)')
    ax[1].semilogy(t, ek, 'o-'); ax[1].set_xlabel('time'); ax[1].set_ylabel('kinetic energy $\\int \\rho W(W-1)$')
    ax[1].set_title('kinetic energy')
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return np.array(t), np.array(em), np.array(ek)


def spectra(files, final, out):
    sel = [files[len(files) // 4], files[len(files) // 2], files[-1]] if len(files) >= 3 else files
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), constrained_layout=True)
    cols = plt.cm.viridis(np.linspace(0.1, 0.9, len(sel) + 1))
    for c, f in zip(cols, sel + [None]):
        S = final if f is None else load(f)
        k, Ek = spectrum(S['P']['vx'], S['P']['vy'], S['P']['vz'], S['dx'])
        k, Em = spectrum(S['P']['Bx'], S['P']['By'], S['P']['Bz'], S['dx'])
        ax[0].loglog(k[1:], Ek[1:], color=c, label='t = %.2f' % S['t'])
        ax[1].loglog(k[1:], Em[1:], color=c, label='t = %.2f' % S['t'])
    kk = np.array([3.0, 15.0])
    ax[0].loglog(kk, 3e-3 * (kk / 3.0) ** (-5.0 / 3.0), 'k--', lw=0.8, label='$k^{-5/3}$ (guide)')
    for a_, ttl in zip(ax, ('velocity spectrum $\\frac{1}{2}|\\hat v|^2$', 'magnetic spectrum $\\frac{1}{2}|\\hat B|^2$')):
        a_.set_xlabel('k (units of $2\\pi/L$)'); a_.set_title(ttl); a_.legend(fontsize=8)
    fig.savefig(out, dpi=130)
    plt.close(fig)


def profiles(files, out):
    sel = [files[len(files) // 4], files[len(files) // 2], files[-1]] if len(files) >= 3 else files
    fig, ax = plt.subplots(1, 4, figsize=(19, 4.2), constrained_layout=True)
    cols = plt.cm.viridis(np.linspace(0.1, 0.9, len(sel)))
    for c, f in zip(cols, sel):
        S = load(f); D = derived(S, 2); y = axes_of(S)[1]
        for a_, q in zip(ax, (S['P']['rho'], S['P']['vx'], D['B2'], D['Jmag'])):
            a_.plot(y, q.mean(axis=(0, 2)), color=c, label='t = %.2f' % S['t'])
    for a_, ttl, lg in zip(ax, ('$\\langle\\rho\\rangle_{xz}$', '$\\langle v_x\\rangle_{xz}$', '$\\langle |B|^2\\rangle_{xz}$', '$\\langle |J|\\rangle_{xz}$'), (False, False, True, True)):
        a_.set_xlabel('y'); a_.set_title(ttl)
        if lg:
            a_.set_yscale('log')
    ax[0].legend(fontsize=8)
    fig.savefig(out, dpi=130)
    plt.close(fig)


def jcurrent(S, D, out):
    J, Bs = D['Jmag'].ravel(), D['Bsub'].ravel()
    y = axes_of(S)[1]
    thr = np.percentile(J, 95)
    big = (D['Jmag'] >= thr)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2), constrained_layout=True)
    ax[0].hist(np.log10(J[J > 0]), bins=100); ax[0].set_xlabel('$\\log_{10}|J|$'); ax[0].set_title('distribution of $|J|$')
    ax[1].hist(np.log10(Bs[Bs > 0]), bins=100); ax[1].set_xlabel('$\\log_{10}$ sub-filter $B^2$'); ax[1].set_title('distribution of the sub-filter $B^2$')
    ax[2].plot(y, big.mean(axis=(0, 2))); ax[2].set_xlabel('y'); ax[2].set_title('fraction of the cells in the top 5 % of $|J|$')
    fig.savefig(out, dpi=130)
    plt.close(fig)
    rank = np.argsort(np.argsort(J)).astype(float); rs = np.argsort(np.argsort(Bs)).astype(float)
    corr = np.corrcoef(rank, rs)[0, 1]
    return thr, corr


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("run_dir", nargs="+", help="one or more run directories (e.g. a run and its continuations); their frames are merged in time")
    p.add_argument("--snapshot", default="data_window_00.hdf5", help="snapshot for the maps: a path, or a file name in the first RUN_DIR")
    p.add_argument("--outdir", default="figures")
    p.add_argument("--filter-cells", type=int, default=4)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    S = load(a.snapshot if os.path.exists(a.snapshot) else os.path.join(a.run_dir[0], a.snapshot))
    D = derived(S, a.filter_cells)
    N = S['P']['rho'].shape
    maps(S, D, 'xy', N[2] // 2, os.path.join(a.outdir, 'maps_xy.png'), a.filter_cells)
    y = axes_of(S)[1]
    maps(S, D, 'xz', int(np.argmin(np.abs(y - 0.25))), os.path.join(a.outdir, 'maps_xz.png'), a.filter_cells)
    files = sorted((f for d in a.run_dir for f in frames(d)), key=lambda f: load_time(f))
    if files:
        t, em, ek = energy_history(files, os.path.join(a.outdir, 'energy_history.png'))
        print("magnetic energy: %.3g at t = %.3f -> %.3g at t = %.3f; kinetic %.3g -> %.3g" % (em[0], t[0], em[-1], t[-1], ek[0], ek[-1]))
        np.savez(os.path.join(a.outdir, 'energy_history.npz'), t=t, magnetic=em, kinetic=ek)
        spectra(files + [], S, os.path.join(a.outdir, 'spectra.png'))
        profiles(files, os.path.join(a.outdir, 'profiles.png'))
    thr, corr = jcurrent(S, D, os.path.join(a.outdir, 'jcurrent_pdf.png'))
    print("final snapshot t = %.3f: |J| 95th percentile %.3g; rank correlation of |J| with the sub-filter B^2 (top-hat %d cells): %.2f" % (S['t'], thr, a.filter_cells, corr))
    print("max |v| %.3f (W %.2f), rho [%.3f, %.3f], max |B| %.3g, min beta %.3g" % (D['vmag'].max(), D['W'].max(), S['P']['rho'].min(), S['P']['rho'].max(), D['Bmag'].max(), D['beta'].min()))


if __name__ == "__main__":
    main()
