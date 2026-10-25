"""numpy reference for SubfilterGradient: same analytic fields, the formulas of Lagrangian-filtering OhmsRegression.em_flux_deficits
(validated against filtered analytic fields there), centred differences; compares with the C++ source array."""
import sys

import numpy as np

n, Ng = 32, 4
Nn = n + 2 * Ng
h = 1.0 / n
densityTerm, limiter, c, quadFactor = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
sigma2 = c * h * h
sigma2q = quadFactor * sigma2      # sigma^2 of the quadratic gradient parts of the Poynting and stress deficits

x = (np.arange(Nn) - Ng + 0.5) * h
X = np.stack(np.meshgrid(x, x, x, indexing='ij'), axis=-1)
ph = (0.3, 1.1, 2.0)
s = [np.sin(2 * np.pi * X[..., m] + ph[m]) for m in range(3)]
co = [np.cos(2 * np.pi * X[..., m] + ph[m]) for m in range(3)]
rho = 1.5 + 0.4 * s[0] * s[1] * co[2]
v = 0.35 * np.stack([s[1] * co[2], s[2] * co[0], s[0] * co[1]], axis=-1)
B = np.stack([1.0 + 0.5 * co[0] * s[1], 0.8 * co[1] * s[2], 0.6 * co[2] * s[0]], axis=-1)

g = lambda f, l: np.gradient(f, h, axis=l)
E = -np.cross(v, B)
D = rho / np.sqrt(1.0 - np.sum(v ** 2, axis=-1))
dv = [g(v, l) for l in range(3)]
dB = [g(B, l) for l in range(3)]
dD = [g(D, l) for l in range(3)]
dE = [-np.cross(dv[l], B) - np.cross(v, dB[l]) for l in range(3)]
emf = -sigma2 * sum(np.cross(dv[l], dB[l]) for l in range(3))
if densityTerm:
    q = sum(dD[l][..., None] * dv[l] for l in range(3)) / D[..., None]
    emf = emf + sigma2 * np.cross(q, B)
poy = np.cross(emf, B) + sigma2q * sum(np.cross(dE[l], dB[l]) for l in range(3))
eye = np.eye(3)
Eem = np.sum(E * emf, axis=-1)
quad = sum(-dE[l][..., :, None] * dE[l][..., None, :] - dB[l][..., :, None] * dB[l][..., None, :]
           + 0.5 * eye * (np.sum(dE[l] ** 2, -1) + np.sum(dB[l] ** 2, -1))[..., None, None] for l in range(3))
stress = -(E[..., :, None] * emf[..., None, :] + emf[..., :, None] * E[..., None, :]) + eye * Eem[..., None, None] + sigma2q * quad

curl = lambda F: np.stack([g(F[..., 2], 1) - g(F[..., 1], 2), g(F[..., 0], 2) - g(F[..., 2], 0), g(F[..., 1], 0) - g(F[..., 0], 1)], axis=-1)
src = np.zeros((9, Nn, Nn, Nn))
sB = -curl(emf)
sS = -np.stack([sum(g(stress[..., a, l], l) for l in range(3)) for a in range(3)], axis=-1)
sT = -sum(g(poy[..., l], l) for l in range(3))
J = curl(B)
exch = np.sum(J * emf, axis=-1)
keep = np.ones_like(exch, dtype=bool) if not limiter else (exch >= 0)
for a in range(3):
    src[1 + a] = np.where(keep, sS[..., a], 0.0)
    src[5 + a] = np.where(keep, sB[..., a], 0.0)
src[4] = np.where(keep, sT, 0.0)

cpp = np.fromfile(sys.argv[5], dtype=np.float64).reshape(9, Nn, Nn, Nn)
sl = (slice(None),) + (slice(Ng, Nn - Ng),) * 3
scale = np.max(np.abs(src[sl]))
err = np.max(np.abs(cpp[sl] - src[sl])) / scale
print("max |C++ - numpy| / max|source| = %.2e   (max |source| = %.3e)" % (err, scale))
sys.exit(0 if err < 1e-10 else 1)
