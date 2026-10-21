#!/usr/bin/env python
"""
Synthetic check of compare_coarse.py (no METHOD build, no real data). Needs numpy, h5py, scipy, matplotlib.

A smooth random fine state is written as a METHOD-like file; two "coarse runs" are made from it:
  exact : conserved variables = the block average of the fine ones (the error must be ~0, energy ratio 1)
  noisy : the same plus noise                                         (the error must be clearly positive)
The script is run on them and its npz output is checked. Also checks conserved_from_prims against hand values.
"""
import os
import subprocess
import sys
import tempfile

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.join(HERE, '..', '..', 'Examples', 'KelvinHelmholtz', 'SingleFluidIdealVigano2020MPIHDF5FilterWindow')
sys.path.insert(0, EX)
import viz_common as vc

PRIMS = ('rho', 'p', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz')


def smooth(rng, n, amp):
    f = np.fft.fftn(rng.standard_normal((n, n, n)))
    k = np.fft.fftfreq(n) * n
    K = np.sqrt(k[:, None, None] ** 2 + k[None, :, None] ** 2 + k[None, None, :] ** 2)
    f *= np.exp(-(K / 4.0) ** 2)
    g = np.fft.ifftn(f).real
    return amp * g / g.std()


def write(path, P, t, C=None):
    n = P['rho'].shape[0]
    with h5py.File(path, 'w') as f:
        f.attrs['t'] = np.array([t])
        d = f.create_group('Domain')
        d.attrs['dx'] = np.array([1.0 / n])
        for k in 'xyz':
            d.attrs[k + 'min'] = np.array([-0.5])
            d.attrs[k + 'max'] = np.array([0.5])
        g = f.create_group('Primitive')
        for k in PRIMS:
            g[k] = P[k]
        if C is not None:
            g = f.create_group('Conserved')
            for k in vc.CONS:
                g[k] = C[k]


def check_conversion():
    P = {k: np.full((2, 2, 2), v) for k, v in dict(rho=1.0, p=0.5, vx=0.0, vy=0.0, vz=0.0, Bx=0.0, By=0.0, Bz=0.0).items()}
    C = vc.conserved_from_prims(P, gamma=4.0 / 3.0)
    h = 1.0 + 0.5 / 1.0 * 4.0       # gamma/(gamma-1) = 4
    assert np.allclose(C['D'], 1.0) and np.allclose(C['Sx'], 0.0)
    assert np.allclose(C['tau'], 1.0 * h - 0.5 - 1.0), C['tau']      # rho h - p - D at rest, no field
    P['vx'][:] = 0.6; P['Bx'][:] = 0.3                               # B parallel to v: S = (rho h + B^2/W^2 ...) W^2 v - b0 b_x
    C = vc.conserved_from_prims(P, gamma=4.0 / 3.0)
    W = 1.0 / np.sqrt(1 - 0.36)
    b0 = W * 0.6 * 0.3
    bx = 0.3 / W + b0 * 0.6
    bsq = (0.09 + b0 ** 2) / W ** 2
    assert np.allclose(C['Sx'], (h + bsq) * W ** 2 * 0.6 - b0 * bx)
    print("conserved_from_prims: ok")


def main():
    check_conversion()
    assert np.allclose(vc.block_mean(np.arange(64.0).reshape(4, 4, 4), 2)[0, 0, 0], np.mean([0, 1, 4, 5, 16, 17, 20, 21]))
    rng = np.random.default_rng(1)
    n, r = 32, 2
    P = {'rho': 1.0 + 0.1 * smooth(rng, n, 1.0), 'p': 1.0 + 0.0 * smooth(rng, n, 1.0)}
    P['p'] = 0.5 + 0.05 * smooth(rng, n, 1.0)
    for k in ('vx', 'vy', 'vz'):
        P[k] = 0.2 * smooth(rng, n, 1.0)
    for k in ('Bx', 'By', 'Bz'):
        P[k] = 0.01 * smooth(rng, n, 1.0)
    Cf = vc.conserved_from_prims(P)
    ref = {k: vc.block_mean(v, r) for k, v in Cf.items()}
    Pc = {k: vc.block_mean(v, r) for k, v in P.items()}      # prims of the coarse run: block means (only B is used by the script)
    with tempfile.TemporaryDirectory() as tmp:
        t = 1.25
        write(os.path.join(tmp, 'fine.hdf5'), P, t)
        write(os.path.join(tmp, 'exact.hdf5'), Pc, t, C=ref)
        noisy = {k: v + 0.02 * np.std(v) * rng.standard_normal(v.shape) for k, v in ref.items()}
        noisy.update({k: ref[k] for k in ('Bx', 'By', 'Bz')})
        write(os.path.join(tmp, 'noisy.hdf5'), Pc, t, C=noisy)
        out = os.path.join(tmp, 'out')
        subprocess.check_call([sys.executable, os.path.join(EX, 'compare_coarse.py'), '--fine', os.path.join(tmp, 'fine.hdf5'),
                               '--coarse', 'exact:' + os.path.join(tmp, 'exact.hdf5'), 'noisy:' + os.path.join(tmp, 'noisy.hdf5'),
                               '--out', out])
        z = np.load(os.path.join(out, 'comparison.npz'))
        assert int(z['r']) == r
        for k in vc.CONS:
            e = z['err_' + k]
            assert e[0] < 1e-10, (k, e)
            if k not in ('Bx', 'By', 'Bz'):
                assert 0.005 < e[1] < 0.2, (k, e)
            assert abs(z['amp_' + k][0] - 1.0) < 1e-10
        assert abs(z['Emag'][0] / float(z['Emag_ref']) - 1.0) < 1e-12
        for f in ('comparison.txt', 'spectra.png', 'slices.png', 'timeseries.png'):
            assert os.path.exists(os.path.join(out, f)), f
        # mismatched times must be refused
        write(os.path.join(tmp, 'late.hdf5'), Pc, t + 0.1, C=ref)
        rc = subprocess.call([sys.executable, os.path.join(EX, 'compare_coarse.py'), '--fine', os.path.join(tmp, 'fine.hdf5'),
                              '--coarse', 'late:' + os.path.join(tmp, 'late.hdf5'), '--out', out],
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        assert rc != 0
        # the scan driver: two frame times, a baseline and a noisy run; the second noisy frame is missing in one run -> skipped
        for sub in ('fine', 'base', 'run1'):
            os.makedirs(os.path.join(tmp, sub))
        for tt in (0.5, 1.0):
            write(os.path.join(tmp, 'fine', 'frame.checkpoint.%.6f.hdf5' % tt), P, tt)
            write(os.path.join(tmp, 'base', 'frame.checkpoint.%.6f.hdf5' % tt), Pc, tt, C=ref)
        write(os.path.join(tmp, 'run1', 'frame.checkpoint.0.500000.hdf5'), Pc, 0.5, C=noisy)
        out2 = os.path.join(tmp, 'scan')
        txt = subprocess.check_output([sys.executable, os.path.join(EX, 'compare_scan.py'), '--fine-dir', os.path.join(tmp, 'fine'),
                                       '--coarse', 'base:' + os.path.join(tmp, 'base'), 'run1:' + os.path.join(tmp, 'run1'), '--out', out2]).decode()
        assert os.path.exists(os.path.join(out2, 'summary.txt')) and os.path.exists(os.path.join(out2, 't_0.5000', 'comparison.npz'))
        assert not os.path.exists(os.path.join(out2, 't_1.0000')) and 'skipped' in txt
    print("compare_coarse: all checks passed")


if __name__ == '__main__':
    main()
