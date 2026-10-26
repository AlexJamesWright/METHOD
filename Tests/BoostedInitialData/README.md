# Boosted, continuous Vigano KHI initial data

Standalone serial check of `KHInstabilitySingleFluidVigano2020Boosted` (`Project/CPU/Src/initFunc.cc`), design doc section 69.

    ./run.sh build     # compile only (g++, no HDF5/MPI)
    ./run.sh           # build if needed and run all cases (needs numpy: python3 or PY=/path/to/python)

`test.cc` builds the initial data on a 32^3 grid (with ghosts), checks inside C++ that `continuous = 0`, `beta = 0` reproduces `KHInstabilitySingleFluidVigano2020` bitwise, and
writes the primitives; `reference.py` recomputes them in numpy (rest-frame formulas, Lorentz matrix on u^mu and F^{mu nu}), compares to 1e-13, and checks rho and p unchanged,
|v| < 1, F u = 0 in the boosted frame and, at beta = 0, the continuity at y = 0 and across the periodic wrap (jump below half the largest step elsewhere; the old data fail this).
Cases: beta_z = 0, 0.5, 0.9, 0.99, beta_x = 0.9, a general beta, the old data (continuous = 0) and no seed field.
