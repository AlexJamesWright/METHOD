#!/usr/bin/env python
"""
Synthetic check of compare_coarse.py (no METHOD build, no real data). Needs numpy, h5py, scipy, matplotlib.

A smooth random fine state is written as a METHOD-like file; two "coarse runs" are made from it:
  exact : conserved variables = the block average of the fine ones (the error must be ~0, energy ratio 1)
  noisy : the same plus noise                                         (the error must be clearly positive)
The script is run on them and its npz output is checked. Also checks conserved_from_prims against hand values.
"""
import os
import shutil
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


def smooth(rng, n, amp, kc=4.0):
    f = np.fft.fftn(rng.standard_normal((n, n, n)))
    k = np.fft.fftfreq(n) * n
    K = np.sqrt(k[:, None, None] ** 2 + k[None, :, None] ** 2 + k[None, None, :] ** 2)
    f *= np.exp(-(K / kc) ** 2)
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


def check_apriori():
    """On smooth fields the second-order EMF closure must predict the exact block-average EMF of the discrete data up to the within-cell variance."""
    sys.path.insert(0, EX)
    import apriori_emf
    rng = np.random.default_rng(3)
    n = 64
    P = {'rho': 1.0 + 0.05 * smooth(rng, n, 1.0, kc=1.5), 'p': 0.5 + 0.0 * smooth(rng, n, 1.0, kc=1.5)}
    for k in ('vx', 'vy', 'vz'):
        P[k] = 0.2 * smooth(rng, n, 1.0, kc=1.5)
    for k in ('Bx', 'By', 'Bz'):
        P[k] = 0.01 * smooth(rng, n, 1.0, kc=1.5)
    with tempfile.TemporaryDirectory() as tmp:
        write(os.path.join(tmp, 'f.hdf5'), P, 1.0)
        res = apriori_emf.analyse(os.path.join(tmp, 'f.hdf5'), 2)
    assert res['lab']['R2best'] > 0.999, res['lab']                       # the direction is exact on smooth fields
    assert abs(res['lab']['scale'] - 0.75) < 0.03, res['lab']            # (r^2 - 1)/r^2 for point values and r = 2
    assert abs(res['lab']['scale_c'] - 1.0) < 0.03, res['lab']           # with the within-cell variance added back
    assert res['favre']['R2best'] > 0.99, res['favre']
    print("apriori_emf: ok (R2 %.3f, scale %.3f)" % (res['lab']['R2'], res['lab']['scale']))


def main():
    check_conversion()
    check_apriori()
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
        # equal-E_kin mode: the coarse runs lag the fine one by 0.2 in E_kin(t); coarse frames at 0.7 and 1.2 must be matched to the fine frames at 0.5 and 1.0
        import math
        for sub in ('finelog', 'cb', 'cr'):
            os.makedirs(os.path.join(tmp, sub))
        def mklog(path, lag):
            with open(path, 'w') as f:
                for tt in np.arange(0.0, 2.0, 0.05):
                    f.write("diag t=%.8f iter=1 Emag=5e-07 Ekin=%.8e mass=1 Wmax=1\n" % (tt, 0.8 * math.exp(-0.3 * (tt - lag))))
        mklog(os.path.join(tmp, 'finelog', 'fine.out'), 0.0)
        mklog(os.path.join(tmp, 'cb', 'run.log'), 0.2)
        mklog(os.path.join(tmp, 'cr', 'run.log'), 0.2)
        for tt in (0.5, 1.0):
            write(os.path.join(tmp, 'finelog', 'frame.checkpoint.%.6f.hdf5' % tt), P, tt)
        for sub, C_ in (('cb', ref), ('cr', noisy)):
            for tt in (0.7, 1.2, 1.6):
                write(os.path.join(tmp, sub, 'frame.checkpoint.%.6f.hdf5' % tt), Pc, tt, C=C_)
        out3 = os.path.join(tmp, 'scan_ekin')
        txt = subprocess.check_output([sys.executable, os.path.join(EX, 'compare_scan.py'), '--fine-dir', os.path.join(tmp, 'finelog'), '--coarse',
                                       'base:' + os.path.join(tmp, 'cb'), 'run1:' + os.path.join(tmp, 'cr'), '--match', 'ekin',
                                       '--fine-log', os.path.join(tmp, 'finelog', 'fine.out'), '--out', out3]).decode()
        assert 'MATCHED COARSE TIMES' in txt, txt
        lines = [l for l in txt.splitlines() if l.startswith('0.5000') or l.startswith('1.0000')]
        assert any('0.700' in l for l in lines) and any('1.200' in l for l in lines), txt
        assert os.path.exists(os.path.join(out3, 't_0.5000', 'comparison.npz'))
        # frames that are 6 % off in E_kin (coarse frames at 0.9 and 1.4 instead of 0.7 and 1.2) are refused at a 1 % tolerance and accepted at 10 %
        os.makedirs(os.path.join(tmp, 'coff'))
        shutil.copy(os.path.join(tmp, 'cb', 'run.log'), os.path.join(tmp, 'coff', 'run.log'))
        for tt in (0.9, 1.4):
            write(os.path.join(tmp, 'coff', 'frame.checkpoint.%.6f.hdf5' % tt), Pc, tt, C=ref)
        cmd = [sys.executable, os.path.join(EX, 'compare_scan.py'), '--fine-dir', os.path.join(tmp, 'finelog'), '--coarse', 'off:' + os.path.join(tmp, 'coff'),
               '--match', 'ekin', '--fine-log', os.path.join(tmp, 'finelog', 'fine.out')]
        rc = subprocess.run(cmd + ['--ekin-tol', '0.01', '--out', os.path.join(tmp, 'scan_tight')], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert rc.returncode != 0 and b'no common frame times' in rc.stderr + rc.stdout
        rc = subprocess.run(cmd + ['--ekin-tol', '0.1', '--out', os.path.join(tmp, 'scan_loose')], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert rc.returncode == 0, rc.stderr
        # compare_slices.py: one figure per E_kin level, a column per run (the second run lags by 0.2 in E_kin(t)), an empty panel when no frame is close enough
        os.makedirs(os.path.join(tmp, 'sa'))
        os.makedirs(os.path.join(tmp, 'sb'))
        mklog(os.path.join(tmp, 'sa', 'run.log'), 0.0)
        mklog(os.path.join(tmp, 'sb', 'run.log'), 0.2)
        for tt in (0.5, 0.7, 1.0):
            write(os.path.join(tmp, 'sa', 'frame.checkpoint.%.6f.hdf5' % tt), P, tt)
            write(os.path.join(tmp, 'sb', 'frame.checkpoint.%.6f.hdf5' % tt), P, tt)
        out4 = os.path.join(tmp, 'slices')
        txt = subprocess.check_output([sys.executable, os.path.join(EX, 'compare_slices.py'), '--ekin', '0.6', '0.01', '--outdir', out4, '--vars', 'rho,wz,Bz,Jmag',
                                       'base=' + os.path.join(tmp, 'sa'), 'lag=' + os.path.join(tmp, 'sb')]).decode()
        assert os.path.exists(os.path.join(out4, 'slices_ekin_0.6.png')) and os.path.exists(os.path.join(out4, 'slices_match.txt')), txt
        assert 'no frame within' in txt, txt                       # E_kin = 0.01 is never reached by these logs within the frames
        # compare_spectra.py: stages layout and groups layout, with a reference, one figure and the npz / match files each
        out5 = os.path.join(tmp, 'spec_stages.png')
        txt = subprocess.check_output([sys.executable, os.path.join(EX, 'compare_spectra.py'), '--ekin', '0.6', '0.4', '--ekin-tol', '0.5', '--out', out5, '--reference', os.path.join(tmp, 'sa'),
                                       '--reference-log', os.path.join(tmp, 'sa', 'run.log'), '--block', '2', 'base=' + os.path.join(tmp, 'sa'), 'lag=' + os.path.join(tmp, 'sb')]).decode()
        assert os.path.exists(out5) and os.path.exists(out5[:-4] + '.npz') and os.path.exists(out5[:-4] + '_match.txt'), txt
        out5c = os.path.join(tmp, 'spec_comp.png')                  # compensated: two rows
        subprocess.check_output([sys.executable, os.path.join(EX, 'compare_spectra.py'), '--compensate', '--ekin', '0.6', '0.4', '--ekin-tol', '0.5', '--out', out5c,
                                 'base=' + os.path.join(tmp, 'sa'), 'lag=' + os.path.join(tmp, 'sb')])
        assert os.path.exists(out5c)
        out6 = os.path.join(tmp, 'spec_groups.png')
        subprocess.check_output([sys.executable, os.path.join(EX, 'compare_spectra.py'), '--layout', 'groups', '--groups', 'A=base', 'B=lag,base', '--ekin', '0.6', '--ekin-tol', '0.5', '--out', out6,
                                 'base=' + os.path.join(tmp, 'sa'), 'lag=' + os.path.join(tmp, 'sb')])
        assert os.path.exists(out6)
    print("compare_coarse: all checks passed")


if __name__ == '__main__':
    main()
