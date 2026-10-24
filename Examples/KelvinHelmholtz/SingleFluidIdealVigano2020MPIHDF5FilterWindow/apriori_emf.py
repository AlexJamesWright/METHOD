#!/usr/bin/env python
"""
A-priori test of the closure's electromotive force on fine-run frames, as a function of time (lab frame, flat, no Lagrangian filter).

For every frame the fine primitives are block-averaged over r^3 cells (a top-hat of one coarse cell, Delta = r dx, sigma^2 = Delta^2/12,
i.e. c = 1/12). The EXACT sub-filter EMF of the ideal field E = -v x B is

    E_true = <E> - E_res = -( <v x B> - <v> x <B> )                        (lab-frame velocity average)
    E_true,F = -( <v x B> - v~ x <B> ),  v~ = <D v>/<D>, D = rho W         (Favre average, as the observer of the covariant filter)

and the closure of SubfilterGradient (c = 1/12 here; `scale` rescales it) predicts

    E_model   = -sigma^2 sum_l d_l<v> x d_l<B>
    E_model,F = E_model + sigma^2 sum_l (d_l<D> d_l<v> / <D>) x <B>        (the density term, densityTerm = 1)

with centred differences on the coarse grid. The fine data are cell values, so the block average of r of them has variance (r^2 - 1)/12 dx_f^2 and
misses the variance dx_f^2/12 inside one fine cell (1/r^2 of Delta^2/12): on smooth fields the best scale is (r^2 - 1)/r^2 (0.75 for r = 2, 0.89 for
r = 3), and `scale_c` = scale + 1/r^2 adds the missing part back (1 on smooth fields; test: Tests/CompareCoarse). Printed per frame: R^2 = 1 - |E_true - E_model|^2/|E_true|^2 of the parameter-free
prediction, the best scale s = <E_true.E_model>/<E_model.E_model> (the c that fits best is s/12 in units of the coarse dx^2), the
amplitude ratio, the rms sub-filter EMF relative to the resolved E_res = -<v> x <B>, and the energy exchange with the resolved field,
X = sum J.E_true (J = curl <B>; dE_mag/dt = -X from the sub-filter EMF, so X/E_mag is the rate at which it drains the resolved field, to be compared with the growth rate of E_mag): net/abs and the backscatter fraction (negative part / abs) for the true and the modelled EMF.

    python apriori_emf.py --block 3 --dirs RUN384 [RUN384_restart ...] [--tmin 1.5] [--every 1] [--out apriori_emf.txt]
    python apriori_emf.py --block 3 --files data_window_15.hdf5 --crop 144 240 144 240 168 216    # one snapshot, only a box (the crop of a covariant filtering run, for a like-for-like comparison)

--block 3 for a 384^3 run (onto 128^3), 2 for a 128^3 run (onto 64^3). Needs about 8 fine-grid arrays of memory (about 4 GB at 384^3): run it
on a compute node. Frames of the first --dirs entry are read first; unreadable frames (the newest one of a running job) are skipped.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import viz_common as vc


def bm(f, r):
    return vc.block_mean(f, r)


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def grad_terms(f, dx):
    return [vc.ddx(f, ax, dx) for ax in range(3)]


def analyse(path, r, crop=None):  # noqa: C901
    S = vc.load(path)
    P, dx = S['P'], S['dx'] * r
    sl = (slice(None),) * 3
    if crop is not None:                                    # fine-cell box i0 i1 j0 j1 k0 k1, multiples of r; the outermost coarse cell layer is left out of the sums (periodic differences wrap)
        i0, i1, j0, j1, k0, k1 = crop
        if any(c % r for c in crop):
            raise SystemExit('the crop %s must be a multiple of the block %d' % (crop, r))
        P = {k: v[i0:i1, j0:j1, k0:k1] for k, v in P.items()}
        sl = (slice(1, -1),) * 3
    v = [P['vx'], P['vy'], P['vz']]
    B = [P['Bx'], P['By'], P['Bz']]
    W = 1.0 / np.sqrt(1.0 - (v[0] ** 2 + v[1] ** 2 + v[2] ** 2))
    D = P['rho'] * W
    vxB = cross(v, B)
    vxB_m = [bm(c, r) for c in vxB]
    del vxB
    vm = [bm(c, r) for c in v]
    Bm = [bm(c, r) for c in B]
    Dm = bm(D, r)
    Dvm = [bm(D * c, r) for c in v]
    del v, B, D, W, P
    vF = [c / Dm for c in Dvm]
    vmBm = cross(vm, Bm)
    vFBm = cross(vF, Bm)
    E_true = [-(a - b) for a, b in zip(vxB_m, vmBm)]
    E_trueF = [-(a - b) for a, b in zip(vxB_m, vFBm)]
    sigma2 = dx ** 2 / 12.0
    # model: -sigma^2 sum_l d_l v x d_l B
    dv = [grad_terms(c, dx) for c in vm]       # dv[m][l]
    dB = [grad_terms(c, dx) for c in Bm]
    E_model = [np.zeros_like(Dm) for _ in range(3)]
    for l in range(3):
        c = cross([dv[m][l] for m in range(3)], [dB[m][l] for m in range(3)])
        for m in range(3):
            E_model[m] -= sigma2 * c[m]
    dD = grad_terms(Dm, dx)
    zeta = [np.zeros_like(Dm) for _ in range(3)]
    for l in range(3):
        for m in range(3):
            zeta[m] += sigma2 * dD[l] * dv[m][l] / Dm
    zB = cross(zeta, Bm)
    E_modelF = [a + b for a, b in zip(E_model, zB)]
    Jm = [vc.ddx(Bm[2], 1, dx) - vc.ddx(Bm[1], 2, dx), vc.ddx(Bm[0], 2, dx) - vc.ddx(Bm[2], 0, dx), vc.ddx(Bm[1], 0, dx) - vc.ddx(Bm[0], 1, dx)]

    def score(Et, Em):  # r is the block size
        Et, Em = [a[sl] for a in Et], [a[sl] for a in Em]
        n2, m2, c = np.sum(dot(Et, Et)), np.sum(dot(Em, Em)), np.sum(dot(Et, Em))
        err = [a - b for a, b in zip(Et, Em)]
        s = c / m2 if m2 > 0 else np.nan
        return {'R2': 1 - np.sum(dot(err, err)) / n2, 'scale': s, 'scale_c': s + 1.0 / r ** 2, 'R2best': c * c / (n2 * m2) if m2 > 0 else np.nan, 'amp': np.sqrt(m2 / n2)}

    def exchange(E):
        x = dot([a[sl] for a in Jm], [a[sl] for a in E])
        return np.sum(x) / np.sum(np.abs(x)), np.sum(np.abs(x) - x) / (2 * np.sum(np.abs(x)))

    res = {'t': S['t'], 'lab': score(E_true, E_model), 'favre': score(E_trueF, E_modelF)}
    Eres = [-a[sl] for a in vmBm]
    res['rel_size'] = np.sqrt(np.sum(dot([a[sl] for a in E_true], [a[sl] for a in E_true])) / np.sum(dot(Eres, Eres)))
    res['X_true'] = exchange(E_true)
    res['X_model'] = exchange(E_model)
    Js = [a[sl] for a in Jm]
    res['X_ratio'] = np.sum(dot(Js, [a[sl] for a in E_model])) / np.sum(dot(Js, [a[sl] for a in E_true]))
    # d/dt int |B|^2/2 dV = -int J.E_total dV, so the sub-filter EMF changes the resolved energy at the rate -X; X / E_mag is a rate per unit time
    res['rate'] = np.sum(dot(Js, [a[sl] for a in E_true])) / (0.5 * np.sum(dot([a[sl] for a in Bm], [a[sl] for a in Bm])))
    return res


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('--block', type=int, required=True)
    p.add_argument('--dirs', nargs='*', default=[])
    p.add_argument('--files', nargs='*', default=[], help='HDF5 files (snapshots or frames) instead of the frame.checkpoint.* of --dirs, e.g. data_window_15.hdf5')
    p.add_argument('--crop', nargs=6, type=int, default=None, metavar=('I0', 'I1', 'J0', 'J1', 'K0', 'K1'), help='analyse only this box of fine cells (multiples of --block); the outermost coarse layer is dropped from the sums')
    p.add_argument('--tmin', type=float, default=0.0)
    p.add_argument('--tmax', type=float, default=1e9)
    p.add_argument('--every', type=int, default=1, help='use every n-th frame')
    p.add_argument('--out', default='apriori_emf.txt')
    a = p.parse_args()
    files = list(a.files)
    for d in a.dirs:
        files += vc.frames(d)
    if not files:
        raise SystemExit('give --dirs or --files')
    rows = []
    hdr = "%-8s %-8s %-8s %-8s %-8s %-8s | %-8s %-8s %-9s | %-9s | %-16s | %-16s | %-8s %-8s" % (
        "t", "R2", "scale", "scale_c", "R2best", "amp", "R2(F)", "scale(F)", "scale_c(F)", "|E|/|Eres|", "net/abs,bs(true)", "net/abs,bs(mod)", "X mod/true", "X/Emag")
    lines = [hdr]
    for i, f in enumerate(files):
        if i % a.every:
            continue
        res = analyse(f, a.block, a.crop)
        if not (a.tmin <= res['t'] <= a.tmax):
            continue
        L, F = res['lab'], res['favre']
        line = "%-8.3f %-8.3f %-8.3f %-8.3f %-8.3f %-8.3f | %-8.3f %-8.3f %-9.3f | %-9.4f | %6.2f %6.2f     | %6.2f %6.2f     | %-8.3f %-8.3f" % (
            res['t'], L['R2'], L['scale'], L['scale_c'], L['R2best'], L['amp'], F['R2'], F['scale'], F['scale_c'], res['rel_size'],
            res['X_true'][0], res['X_true'][1], res['X_model'][0], res['X_model'][1], res['X_ratio'], res['rate'])
        print(line, flush=True)
        lines.append(line)
    with open(a.out, 'w') as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == '__main__':
    main()
