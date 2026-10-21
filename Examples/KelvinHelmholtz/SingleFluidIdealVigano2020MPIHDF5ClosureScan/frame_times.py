#!/usr/bin/env python
"""
Times of the 384^3 reference frames, as a comma-separated list for the ClosureScan runs (`outputTimes` argument of ./main).

    python frame_times.py --dir /scratch/ih3/METHOD_simulations/vigano2020_384_run1 [--tmin 4] [--tmax 14] > times.txt   # from the files
    python frame_times.py --formula [--every 370] [--tmax 14] > times.txt                                                   # before they exist

The files give the exact times (attribute `t`). While the job runs, the newest frame is still open and unreadable; it is reported on stderr and skipped (the final frame of a finished job is closed properly). The formula is for
the 384^3 run as set up (cfl 0.9, dx = 1/384, three dimensions, wave speed 1) with its slow start (the first 15 steps are 0.1 dt0 and the
next 15 are 0.25 dt0): t(n) = dt0 (n - 24.75) for n >= 30, dt0 = 0.9/(384 sqrt 3) = 1.35316e-3. It reproduced the logged
`diag t=0.10182564 iter=100` to 8 digits; use the files when they exist. Frames are written every 370 iterations.
"""
import argparse
import glob
import math
import sys

p = argparse.ArgumentParser()
g = p.add_mutually_exclusive_group(required=True)
g.add_argument('--dir')
g.add_argument('--formula', action='store_true')
p.add_argument('--every', type=int, default=370, help='iterations between frames (formula)')
p.add_argument('--tmin', type=float, default=0.0)
p.add_argument('--tmax', type=float, default=14.0)
a = p.parse_args()

times = []
if a.formula:
    dt0 = 0.9 / 384.0 / math.sqrt(3.0)
    n = a.every
    while True:
        t = dt0 * (n - 24.75)
        if t > a.tmax + 1e-9:
            break
        times.append(t)
        n += a.every
else:
    import h5py
    import numpy as np
    for f in glob.glob(a.dir.rstrip('/') + '/frame.checkpoint.*.hdf5'):
        try:
            with h5py.File(f, 'r') as h:
                times.append(float(np.ravel(h.attrs['t'])[0]))
        except OSError:
            # the newest frame is left open by METHOD (saveAll(true)) until the next one is written, so it cannot be read yet
            print("skipped (not readable yet; the newest frame of a running job): %s" % f, file=sys.stderr)
    times.sort()
if not times:
    print("no readable frames found", file=sys.stderr)
print(",".join(repr(t) for t in times if a.tmin <= t <= a.tmax))
