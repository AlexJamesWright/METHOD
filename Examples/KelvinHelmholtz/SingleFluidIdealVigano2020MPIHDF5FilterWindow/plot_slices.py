#!/usr/bin/env python
"""
Slice plots and simple turbulence diagnostics for METHOD 3D SRMHD HDF5 output
(e.g. the data_window_*.hdf5 files written by this example's main.cc).

Non-interactive: uses the Agg backend and writes PNG files only, so it is safe to
run in a batch job or on a login node. Only the requested slices (and a strided
subsample for volume statistics) are read from disk, never whole 3D arrays.

For each selected snapshot it writes slices_<snapshot>.png, a 2x5 figure:
  row 1: z-mid plane (x-y)    row 2: y-mid plane (x-z)
  cols : density, Bx, By, Bz, vorticity component normal to the plane
and it writes one spectrum.png comparing the planar velocity power spectra
of all selected snapshots. A table of statistics is printed to stdout.

Usage (from the directory holding the hdf5 files):
    python plot_slices.py
    python plot_slices.py /path/to/run --snapshots 0 4 9 --outdir plots
"""
import argparse
import glob
import os
import sys

import h5py
import numpy as np

import matplotlib
matplotlib.use("Agg")  # no display needed; must precede the pyplot import
import matplotlib.pyplot as plt

# plane name -> (first in-plane axis, second in-plane axis, normal axis)
PLANES = {"xy": (0, 1, 2), "xz": (0, 2, 1), "yz": (1, 2, 0)}
AXIS_NAMES = "xyz"


def scalar_attr(group, name, default=None):
    """METHOD writes attributes as shape-(1,) arrays."""
    if name not in group.attrs:
        return default
    return float(np.ravel(group.attrs[name])[0])


def domain_info(f):
    """Extents and grid spacings (from the Domain attributes, with sane fallbacks)."""
    shape = f["Primitive/rho"].shape
    dom = f["Domain"]
    info = {"shape": shape, "min": [], "max": [], "d": []}
    for ax, n in zip(AXIS_NAMES, shape):
        lo = scalar_attr(dom, ax + "min", 0.0)
        hi = scalar_attr(dom, ax + "max", float(n))
        info["min"].append(lo)
        info["max"].append(hi)
        info["d"].append((hi - lo) / n)
    return info


def read_plane(f, name, plane, index):
    """Read one 2D slice of a primitive variable. Returns array indexed [a, b]."""
    ds = f["Primitive/" + name]
    normal = PLANES[plane][2]
    sel = [slice(None)] * 3
    sel[normal] = index
    return ds[tuple(sel)]


def vorticity_normal(va, vb, da, db):
    """d_a v_b - d_b v_a for a 2D slice (sign set by the (a, b) plane orientation)."""
    return np.gradient(vb, da, axis=0) - np.gradient(va, db, axis=1)


def plane_spectrum(components):
    """Shell-summed power spectrum of the velocity components on a periodic 2D plane."""
    n0, n1 = components[0].shape
    power = np.zeros((n0, n1))
    for c in components:
        fc = np.fft.fft2(c - c.mean()) / (n0 * n1)
        power += 0.5 * np.abs(fc) ** 2
    k0 = np.fft.fftfreq(n0, d=1.0 / n0)
    k1 = np.fft.fftfreq(n1, d=1.0 / n1)
    kmag = np.rint(np.sqrt(k0[:, None] ** 2 + k1[None, :] ** 2)).astype(int)
    kmax = min(n0, n1) // 2
    spec = np.bincount(kmag.ravel(), weights=power.ravel(), minlength=kmax + 1)[: kmax + 1]
    return np.arange(kmax + 1), spec


def show(ax, fig, data, extent, title, symmetric, cmap):
    """Draw one panel. Colour limits are robust percentiles of the slice."""
    if symmetric:
        lim = np.percentile(np.abs(data), 99.5) or 1.0
        vmin, vmax = -lim, lim
    else:
        vmin, vmax = np.percentile(data, [0.5, 99.5])
        if vmin == vmax:
            vmin, vmax = vmin - 1e-12, vmax + 1e-12
    im = ax.imshow(data.T, origin="lower", extent=extent, cmap=cmap,
                   vmin=vmin, vmax=vmax, aspect="equal", interpolation="nearest")
    ax.set_title(title, fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.8)


def volume_stats(f, stride):
    """Statistics from a strided subsample of the full volume."""
    sub = (slice(None, None, stride),) * 3
    g = {n: f["Primitive/" + n][sub] for n in ("rho", "vx", "vy", "vz", "Bx", "By", "Bz")}
    b2 = g["Bx"] ** 2 + g["By"] ** 2 + g["Bz"] ** 2
    v2 = g["vx"] ** 2 + g["vy"] ** 2 + g["vz"] ** 2
    return {
        "rho_min": g["rho"].min(), "rho_max": g["rho"].max(),
        "rho_std": g["rho"].std(),
        "B_rms": np.sqrt(b2.mean()),
        "Bz_std": g["Bz"].std(),
        "v_rms": np.sqrt(v2.mean()),
        "v_max": np.sqrt(v2.max()),
        "finite": all(np.all(np.isfinite(a)) for a in g.values()),
    }


def make_slice_figure(f, info, frac, out_path, label, dpi):
    fig, axes = plt.subplots(2, 5, figsize=(21, 8), constrained_layout=True)
    for row, plane in enumerate(("xy", "xz")):
        a, b, n = PLANES[plane]
        index = int(np.clip(round(frac * info["shape"][n]), 0, info["shape"][n] - 1))
        pos = info["min"][n] + (index + 0.5) * info["d"][n]
        extent = [info["min"][a], info["max"][a], info["min"][b], info["max"][b]]
        tag = "%s-plane at %s=%.3f" % (plane, AXIS_NAMES[n], pos)

        rho = read_plane(f, "rho", plane, index)
        show(axes[row, 0], fig, rho, extent, "density, " + tag, False, "viridis")
        for col, comp in enumerate(("Bx", "By", "Bz"), start=1):
            data = read_plane(f, comp, plane, index)
            show(axes[row, col], fig, data, extent, comp + ", " + tag, True, "RdBu_r")
        va = read_plane(f, "v" + AXIS_NAMES[a], plane, index)
        vb = read_plane(f, "v" + AXIS_NAMES[b], plane, index)
        w = vorticity_normal(va, vb, info["d"][a], info["d"][b])
        show(axes[row, 4], fig, w, extent,
             "vorticity (normal to plane), " + tag, True, "RdBu_r")
        for col in range(5):
            axes[row, col].set_xlabel(AXIS_NAMES[a])
            axes[row, col].set_ylabel(AXIS_NAMES[b])
    fig.suptitle(label, fontsize=13)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("directory", nargs="?", default=".", help="directory holding the hdf5 files")
    p.add_argument("--pattern", default="data_window_*.hdf5")
    p.add_argument("--snapshots", type=int, nargs="+", default=[0, -1],
                   help="indices into the sorted file list (negative allowed); default: first and last")
    p.add_argument("--outdir", default=None, help="default: <directory>/plots")
    p.add_argument("--frac", type=float, default=0.5,
                   help="position of the slice along the normal axis, as a fraction of the domain")
    p.add_argument("--stride", type=int, default=4,
                   help="subsampling stride for the volume statistics")
    p.add_argument("--dpi", type=int, default=120)
    args = p.parse_args()

    files = sorted(glob.glob(os.path.join(args.directory, args.pattern)))
    if not files:
        sys.exit("No files matching %s in %s" % (args.pattern, args.directory))
    outdir = args.outdir or os.path.join(args.directory, "plots")
    os.makedirs(outdir, exist_ok=True)

    chosen = []
    for s in args.snapshots:
        try:
            chosen.append(files[s])
        except IndexError:
            sys.exit("Snapshot index %d out of range (found %d files)" % (s, len(files)))
    chosen = list(dict.fromkeys(chosen))  # drop duplicates, keep order

    print("%d files found; plotting %d. Output directory: %s\n" % (len(files), len(chosen), outdir))
    header = "%-26s %10s %9s %9s %9s %9s %9s %9s %9s %s" % (
        "file", "t", "rho_min", "rho_max", "rho_std", "B_rms", "Bz_std", "v_rms", "v_max", "finite")
    print(header)

    spectra = []
    for path in chosen:
        base = os.path.splitext(os.path.basename(path))[0]
        with h5py.File(path, "r") as f:
            t = scalar_attr(f, "t", float("nan"))
            info = domain_info(f)
            st = volume_stats(f, max(1, args.stride))
            print("%-26s %10.5f %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %9.3g %s" % (
                base, t, st["rho_min"], st["rho_max"], st["rho_std"], st["B_rms"],
                st["Bz_std"], st["v_rms"], st["v_max"], st["finite"]))

            label = "%s   t = %.5f   grid %dx%dx%d" % ((base, t) + tuple(info["shape"]))
            make_slice_figure(f, info, args.frac, os.path.join(outdir, "slices_%s.png" % base),
                              label, args.dpi)

            nk = info["shape"][2]
            kz = int(np.clip(round(args.frac * nk), 0, nk - 1))
            comps = [read_plane(f, v, "xy", kz) for v in ("vx", "vy", "vz")]
            spectra.append((base, t) + plane_spectrum(comps))

    fig, ax = plt.subplots(figsize=(6.5, 5), constrained_layout=True)
    for base, t, k, spec in spectra:
        ax.loglog(k[1:], spec[1:], label="%s (t=%.4f)" % (base, t))
    k_ref = np.arange(2, spectra[0][2].max() + 1)
    anchor = spectra[0][3][min(4, len(spectra[0][3]) - 1)]
    ax.loglog(k_ref, anchor * (k_ref / 4.0) ** (-5.0 / 3.0), "k--", lw=1, label="$k^{-5/3}$ (guide)")
    ax.set_xlabel("wavenumber k (cycles per box)")
    ax.set_ylabel("planar velocity power (xy mid-plane)")
    ax.legend(fontsize=8)
    ax.set_title("Velocity power spectrum, 2D slice")
    fig.savefig(os.path.join(outdir, "spectrum.png"), dpi=args.dpi)
    plt.close(fig)

    print("\nWrote %d slice figure(s) and spectrum.png to %s" % (len(chosen), outdir))


if __name__ == "__main__":
    main()
