"""Shared helpers for the Vigano2020 KHI plots: reading METHOD snapshots and derived quantities."""
import glob
import os
import re

import h5py
import numpy as np

PRIMS = ('rho', 'p', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz')


def load(path):
    """Primitive variables, grid spacing and time of one METHOD HDF5 file (snapshot or checkpoint frame)."""
    with h5py.File(path, 'r') as f:
        dom = f['Domain'].attrs
        dx = float(np.ravel(dom['dx'])[0])
        lo = [float(np.ravel(dom[k + 'min'])[0]) for k in 'xyz']
        hi = [float(np.ravel(dom[k + 'max'])[0]) for k in 'xyz']
        if 't' in f.attrs:                                   # checkpoint frames: t at the file root
            t = float(np.ravel(f.attrs['t'])[0])
        else:                                                # window snapshots: endTime was set to t
            t = float(np.ravel(dom['endTime'])[0])
        P = {k: f['Primitive'][k][()].astype(np.float64) for k in PRIMS}
    return {'t': t, 'dx': dx, 'lo': lo, 'hi': hi, 'P': P, 'path': path}


def frames(directory):
    """Checkpoint frames in a run directory, sorted by time (file name frame.checkpoint.<t>.hdf5)."""
    out = []
    for p in glob.glob(os.path.join(directory, 'frame.checkpoint.*.hdf5')):
        m = re.search(r'frame\.checkpoint\.([0-9.]+)\.hdf5$', p)
        out.append((float(m.group(1).rstrip('.')), p))
    good = []
    for _, p in sorted(out):
        try:
            with h5py.File(p, 'r'):
                good.append(p)
        except OSError:      # the last frame is never closed by METHOD's saveAll(true) and may be unreadable
            print("skipping unreadable frame", os.path.basename(p))
    return good


def ddx(f, axis, dx):
    """Centred difference on the periodic box."""
    return (np.roll(f, -1, axis) - np.roll(f, 1, axis)) / (2.0 * dx)


def curl(vx, vy, vz, dx):
    return (ddx(vz, 1, dx) - ddx(vy, 2, dx), ddx(vx, 2, dx) - ddx(vz, 0, dx), ddx(vy, 0, dx) - ddx(vx, 1, dx))


def tophat(f, cells):
    """Box average over cells^3 grid cells (periodic): a lab-frame stand-in for the filter of the pipeline."""
    from scipy.ndimage import uniform_filter
    return uniform_filter(f, size=cells, mode='wrap')


def derived(S, filter_cells=4):
    """Derived fields used in the figures. Lab-frame, flat metric. J = curl B (displacement current neglected)."""
    P, dx = S['P'], S['dx']
    v2 = P['vx'] ** 2 + P['vy'] ** 2 + P['vz'] ** 2
    W = 1.0 / np.sqrt(1.0 - v2)
    B2 = P['Bx'] ** 2 + P['By'] ** 2 + P['Bz'] ** 2
    Ex = -(P['vy'] * P['Bz'] - P['vz'] * P['By'])
    Ey = -(P['vz'] * P['Bx'] - P['vx'] * P['Bz'])
    Ez = -(P['vx'] * P['By'] - P['vy'] * P['Bx'])
    Jx, Jy, Jz = curl(P['Bx'], P['By'], P['Bz'], dx)
    wx, wy, wz = curl(P['vx'], P['vy'], P['vz'], dx)
    # sub-filter magnetic energy <B^2> - |<B>|^2 for a top-hat of `filter_cells` cells (lab frame)
    Bsub = tophat(B2, filter_cells) - sum(tophat(P[k], filter_cells) ** 2 for k in ('Bx', 'By', 'Bz'))
    D = {'W': W, 'vmag': np.sqrt(v2), 'B2': B2, 'Bmag': np.sqrt(B2), 'beta': P['p'] / np.maximum(0.5 * B2, 1e-300),
         'E2': Ex ** 2 + Ey ** 2 + Ez ** 2, 'Emag': np.sqrt(Ex ** 2 + Ey ** 2 + Ez ** 2),
         'Jmag': np.sqrt(Jx ** 2 + Jy ** 2 + Jz ** 2), 'Jz': Jz, 'wx': wx, 'wy': wy, 'wz': wz,
         'Bsub': np.maximum(Bsub, 0.0), 'ekin': P['rho'] * W * (W - 1.0), 'emag': 0.5 * B2}
    return D


def axes_of(S):
    N = S['P']['rho'].shape
    return [S['lo'][a] + (np.arange(N[a]) + 0.5) * S['dx'] for a in range(3)]


def spectrum(fx, fy, fz, dx):
    """Shell-summed power spectrum of a vector field: E(k) = sum over the shell of |f_k|^2 / 2 (k in units of 2 pi / L)."""
    N = fx.shape[0]
    power = sum(np.abs(np.fft.fftn(f)) ** 2 for f in (fx, fy, fz)) / (2.0 * N ** 6)
    kk = np.fft.fftfreq(N, d=1.0 / N)
    k = np.sqrt(kk[:, None, None] ** 2 + kk[None, :, None] ** 2 + kk[None, None, :] ** 2)
    shell = np.rint(k).astype(int)
    kmax = N // 2
    E = np.bincount(shell.ravel(), weights=power.ravel(), minlength=kmax + 1)[:kmax + 1]
    return np.arange(kmax + 1), E


CONS = ('D', 'Sx', 'Sy', 'Sz', 'tau', 'Bx', 'By', 'Bz')


def conserved_from_prims(P, gamma=4.0 / 3.0):
    """Conserved variables (D, S, tau, B) of ideal SRMHD from primitives, as SRMHD::primsToAll (mu0 = 1, B^2 without 1/2 in h)."""
    v2 = P['vx'] ** 2 + P['vy'] ** 2 + P['vz'] ** 2
    W = 1.0 / np.sqrt(1.0 - v2)
    vB = P['vx'] * P['Bx'] + P['vy'] * P['By'] + P['vz'] * P['Bz']
    B2 = P['Bx'] ** 2 + P['By'] ** 2 + P['Bz'] ** 2
    b0 = W * vB
    bsq = (B2 + b0 ** 2) / W ** 2
    h = 1.0 + P['p'] / P['rho'] * (gamma / (gamma - 1.0))
    wtot = (P['rho'] * h + bsq) * W ** 2
    C = {'D': P['rho'] * W}
    for a in 'xyz':
        b_a = P['B' + a] / W + b0 * P['v' + a]
        C['S' + a] = wtot * P['v' + a] - b0 * b_a
    C['tau'] = wtot - (P['p'] + bsq / 2.0) - b0 ** 2 - C['D']
    for a in 'xyz':
        C['B' + a] = P['B' + a]
    return C


def block_mean(f, r):
    """Average over r^3 blocks of cells: the conservative restriction of a cell-averaged field to a grid r times coarser."""
    n = f.shape
    return f.reshape(n[0] // r, r, n[1] // r, r, n[2] // r, r).mean(axis=(1, 3, 5))


def load_cons(path):
    """Conserved variables stored in a METHOD file (restart frames, data_final), or None."""
    with h5py.File(path, 'r') as f:
        if 'Conserved' not in f:
            return None
        return {k: f['Conserved'][k][()].astype(np.float64) for k in CONS}


def shell_spectrum(fields, dx):
    """Shell-summed power spectrum sum_i |f_i(k)|^2 of vector/scalar fields on the periodic cube; returns (k, P(k))."""
    n = fields[0].shape[0]
    kx = np.fft.fftfreq(n, d=dx) * 2.0 * np.pi
    K = np.sqrt(kx[:, None, None] ** 2 + kx[None, :, None] ** 2 + kx[None, None, :] ** 2)
    dk = 2.0 * np.pi / (n * dx)
    power = sum(np.abs(np.fft.fftn(f)) ** 2 for f in fields) / float(n ** 6)
    idx = np.rint(K / dk).astype(int)
    P = np.bincount(idx.ravel(), weights=power.ravel())[:n // 2 + 1]
    return dk * np.arange(n // 2 + 1), P
