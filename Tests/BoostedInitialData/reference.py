#!/usr/bin/env python3
"""
Independent numpy check of KHInstabilitySingleFluidVigano2020Boosted: the rest-frame data are evaluated from the formulas, the state is boosted with the Lorentz matrix
acting on u^mu and on F^{mu nu} (F^{0i} = E_i, F^{ij} = eps_ijk B_k, E = -v x B), and compared with the primitives written by test.cc; also checks F u = 0 in the boosted
frame, the scalars rho and p, and the continuity of the data at y = 0 and across the periodic wrap.

    python3 reference.py out.bin betaX betaY betaZ continuous mag
"""
import sys

import numpy as np

fn = sys.argv[1]
beta = np.array([float(a) for a in sys.argv[2:5]])
continuous, mag = int(sys.argv[5]), int(sys.argv[6])

with open(fn, 'rb') as fh:
    Nx, Ny, Nz, Np, Ng = np.fromfile(fh, dtype=np.int32, count=5)
    x = np.fromfile(fh, dtype=np.float64, count=Nx)
    y = np.fromfile(fh, dtype=np.float64, count=Ny)
    z = np.fromfile(fh, dtype=np.float64, count=Nz)
    P = np.fromfile(fh, dtype=np.float64).reshape(Np, Nx, Ny, Nz)

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
B0, p0, rho0, rho1, yl, al = 0.001, 1.0, 2.0, 1.0, 0.25, 0.01
v0x, v0y, v0z, dvx, dvy, dvz = 0.5, 0.0, 0.5, 0.01, 0.1, 0.01
nx, ny, nz, sy2, sz2 = 11, 7, 5, 0.01, 0.1
sgn, ay = np.sign(Y), np.abs(Y)
layer = np.tanh((ay - yl) / al)
envy, envz = np.exp(-(ay - yl) ** 2 / sy2), np.exp(-(ay - yl) ** 2 / sz2)
if continuous:
    rho = rho0 + rho1 * layer
    v = np.stack([v0x * layer + dvx * np.sin(2 * np.pi * nx * X), v0y * layer + dvy * np.sin(2 * np.pi * ny * ay) * envy,
                  v0z * np.sin(2 * np.pi * Y) * envz + dvz * np.sin(2 * np.pi * nz * Z)], axis=-1)
else:
    rho = rho0 + rho1 * sgn * layer
    v = np.stack([v0x * sgn * layer + dvx * np.sin(2 * np.pi * nx * X), v0y * sgn * layer + dvy * np.sin(2 * np.pi * ny * Y) * sgn * envy,
                  v0z * sgn * envz + dvz * np.sin(2 * np.pi * nz * Z)], axis=-1)
B = np.zeros_like(v)
B[..., 0] = B0 * mag
E = -np.cross(v, B)
W = 1 / np.sqrt(1 - np.sum(v ** 2, axis=-1))
u = np.concatenate([W[..., None], W[..., None] * v], axis=-1)
F = np.zeros(v.shape[:-1] + (4, 4))
F[..., 0, 1:] = E
F[..., 1:, 0] = -E
F[..., 1:, 1:] = np.einsum('ijk,...k->...ij', np.array([[[0, 0, 0], [0, 0, 1], [0, -1, 0]], [[0, 0, -1], [0, 0, 0], [1, 0, 0]], [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]], dtype=float), B)

b2 = beta @ beta
if b2 > 0:
    g = 1 / np.sqrt(1 - b2)
    L = np.eye(4)
    L[0, 0], L[0, 1:], L[1:, 0] = g, -g * beta, -g * beta
    L[1:, 1:] += (g - 1) * np.outer(beta, beta) / b2
    u = u @ L.T
    F = np.einsum('ma,...ab,nb->...mn', L, F, L)
vb = u[..., 1:] / u[..., :1]
Bb = 0.5 * np.einsum('kij,...ij->...k', np.array([[[0, 0, 0], [0, 0, 1], [0, -1, 0]], [[0, 0, -1], [0, 0, 0], [1, 0, 0]], [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]], dtype=float), F[..., 1:, 1:])
ref = np.concatenate([rho[None], np.moveaxis(vb, -1, 0), np.full((1,) + rho.shape, p0), np.moveaxis(Bb, -1, 0)], axis=0)

fail = False
def check(name, ok, msg):
    global fail
    print("%-62s %s   %s" % (name, "ok" if ok else "FAIL", msg))
    fail |= not ok

err = np.abs(P[:8] - ref).max(axis=(1, 2, 3))
check("primitives against the numpy boost (max abs diff per variable)", err.max() < 1e-13 * max(1.0, np.abs(ref).max()), np.array2string(err, precision=1))
check("rho and p unchanged by the boost", np.abs(P[0] - rho).max() < 1e-15 and np.abs(P[4] - p0).max() < 1e-15, "")
vC = np.moveaxis(P[1:4], 0, -1)
check("subluminal everywhere", np.sum(vC ** 2, axis=-1).max() < 1, "max |v| %.6f, W max %.3f" % (np.sqrt(np.sum(vC ** 2, axis=-1).max()), 1 / np.sqrt(1 - np.sum(vC ** 2, axis=-1).max())))
Wc = 1 / np.sqrt(1 - np.sum(vC ** 2, axis=-1))
uc = np.concatenate([Wc[..., None], Wc[..., None] * vC], axis=-1)
Bc = np.moveaxis(P[5:8], 0, -1)
Ec = -np.cross(vC, Bc)
Fc = np.zeros(vC.shape[:-1] + (4, 4))
Fc[..., 0, 1:] = Ec
Fc[..., 1:, 0] = -Ec
Fc[..., 1:, 1:] = np.einsum('ijk,...k->...ij', np.array([[[0, 0, 0], [0, 0, 1], [0, -1, 0]], [[0, 0, -1], [0, 0, 0], [1, 0, 0]], [[0, 1, 0], [-1, 0, 0], [0, 0, 0]]], dtype=float), Bc)
eta = np.array([-1.0, 1, 1, 1])
check("ideal condition F^{mu nu} u_nu = 0 in the boosted frame", np.abs(np.einsum('...mn,...n->...m', Fc, uc * eta)).max() < 1e-12 * max(1, np.abs(Fc).max() * Wc.max()), "")

# continuity of the data at y = 0 and across the periodic wrap (interior cells only): the jump over the interface against the largest step elsewhere
i0 = slice(Ng, Ny - Ng)
yi = y[i0]
mid = np.argmin(np.abs(yi))
def jump(f):
    g = f[:, i0, :]
    return max(np.abs(g[:, mid] - g[:, mid - 1]).max(), np.abs(g[:, 0] - g[:, -1]).max())
def step(f):
    g = f[:, i0, :]
    d = np.abs(np.diff(g, axis=1))
    d[:, mid - 1] = 0
    return d.max()
if beta @ beta == 0:
    for name, f in (("rho", P[0]), ("v^x", P[1]), ("v^z", P[3])):
        j, s = jump(f), step(f)
        msg = "jump over y = 0 / the wrap %.3e, largest step elsewhere %.3e" % (j, s)
        check("continuity of %s" % name, (j < 0.5 * s) if continuous else (j > s), msg + ("" if continuous else "  (the Vigano data are discontinuous there)"))
sys.exit(1 if fail else 0)
