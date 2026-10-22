#!/usr/bin/env python
"""
Table of the ratio of each run's kinetic and magnetic spectrum to a reference at chosen wavenumbers, from the .npz written by compare_spectra.py.

    python spectra_ratio_table.py spectra_stages.npz [--k 4 8 16 32 48] [--ref SUBSTRING] [--csv table.csv]

The .npz holds, for each panel (an E_kin level in the `stages` layout, a group in the `groups` layout) and each curve, the raw shell spectra E(k) of the convention of compare_spectra.py (so a ratio
does not depend on the normalisation). The reference is the curve whose label contains --ref; default: the filtered fine run (`..., filtered`), else the native one, else the first curve. Ratios are
taken at the integer wavenumbers --k (k in units of 2 pi/L; only those present in both curves: the filtered reference lives on the grid of the runs, the native one extends to higher k).
The match file `<name>_match.txt` next to the .npz, if present, supplies the E_kin level and the frame time of each run in the table. A ratio below 1 means less power than the reference.
"""
import argparse
import os
import re

import numpy as np


def load_npz(path):
    z = np.load(path)
    panels = {}
    for key in z.files:
        if not key.endswith('_k'):
            continue
        stem = key[:-2]
        pi, label = stem.split('_', 1)
        panels.setdefault(int(pi), {})[label] = (z[stem + '_k'], z[stem + '_Ek'], z[stem + '_Em'])
    return panels


def read_match(path):
    """[(level, {label: time})] from the match file of compare_spectra.py."""
    out = []
    if not os.path.exists(path):
        return out
    for line in open(path):
        if line.startswith('E_kin ='):
            out.append((float(line.split('=')[1]), {}))
        else:
            m = re.match(r'\s+(.+?)\s+t =\s+([0-9.]+)', line)
            if m and out:
                out[-1][1][m.group(1).strip()] = float(m.group(2))
    return out


def pick_reference(labels, ref):
    if ref:
        hits = [l for l in labels if ref in l]
        if not hits:
            raise SystemExit("no curve label contains %r (labels: %s)" % (ref, ", ".join(labels)))
        return hits[0]
    for pattern in (',_filtered', '(native)'):
        hits = [l for l in labels if pattern in l]
        if hits:
            return hits[0]
    return labels[0]


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument('npz')
    p.add_argument('--k', nargs='+', type=int, default=[4, 8, 16, 32, 48])
    p.add_argument('--ref')
    p.add_argument('--csv')
    a = p.parse_args()
    panels = load_npz(a.npz)
    match = read_match(os.path.splitext(a.npz)[0] + '_match.txt')
    rows = []
    for pi in sorted(panels):
        curves = panels[pi]
        labels = list(curves)
        ref = pick_reference(labels, a.ref)
        kr, ekr, emr = curves[ref]
        title = "panel %d" % pi
        times = {}
        if pi < len(match):
            title = "E_kin = %g" % match[pi][0]
            times = match[pi][1]
        print("\n== %s, reference: %s" % (title, ref.replace('_', ' ')))
        ks = [k for k in a.k if k in kr]
        for which, idx in (('kinetic', 1), ('magnetic', 2)):
            print("  %s: ratio to the reference at k = %s" % (which, ", ".join(str(k) for k in ks)))
            for label in labels:
                if label == ref:
                    continue
                kk, E = curves[label][0], curves[label][idx]
                R = (ekr if idx == 1 else emr)
                vals = []
                for k in ks:
                    if k in kk:
                        i, j = int(np.where(kk == k)[0][0]), int(np.where(kr == k)[0][0])
                        vals.append(E[i] / R[j])
                    else:
                        vals.append(float('nan'))
                short = label.replace('_', ' ')
                t = next((tm for lab, tm in times.items() if short.startswith(lab)), None)
                print("    %-28s %s%s" % (short, "".join("%10.3g" % v for v in vals), "   (t = %.2f)" % t if t is not None else ""))
                rows.append((title, which, short, t, vals))
    if a.csv:
        with open(a.csv, 'w') as fh:
            fh.write("panel,spectrum,curve,t," + ",".join("k%d" % k for k in a.k) + "\n")
            for title, which, short, t, vals in rows:
                fh.write("%s,%s,%s,%s,%s\n" % (title, which, short.replace(',', ';'), "" if t is None else "%.3f" % t, ",".join("%.6g" % v for v in vals)))
        print("\nwrote", a.csv)


if __name__ == '__main__':
    main()
