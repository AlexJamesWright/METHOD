Model
-----
  - Ideal MHD (SRMHD), 3D


Initial data
------------
  - `KHInstabilitySingleFluidCookBernuzzi`: the 3D special-relativistic Kelvin-Helmholtz initial data of Cook & Bernuzzi,
  arXiv:2606.21659, section II.5.3 (which follows Vigano et al.), on the box [-1, 1]^3. Deterministic, no random seed:
  shear layers with velocity perturbations in all three directions (n_x, n_y, n_z = 11, 7, 5), weak seed field B_x = 0.001,
  gamma = 4/3 here (the paper does not state its equation of state).
  - The formulas are implemented as written. Because of the sgn(y) factors the density and velocity are
  discontinuous at y = 0 and across the periodic boundary y = +-1, as well as at the smooth layers y = +-0.25. Check this is
  what is wanted before a long run (see the comment in `Include/initFunc.h`). Max |v| at t = 0 is about 0.91 (W ~ 2.4).
  - Not yet run: neither the initial data nor this example has been built or run.


Purpose
-------
  - Produces output for Lagrangian filtering, exactly as `../SingleFluidIdealRandomMPIHDF5FilterWindow`: spin up to `tStart`
  with no output, then write `nSnap` snapshots of all variables (cons, prims, aux), `data_window_00.hdf5` ..., one every
  `snapStride` timesteps; snapshot times are in `Domain/endTime`; `main` refuses a window that would start inside the
  30-iteration slow-start. See that example's README for the details.
  - `tStart = 5` is carried over from the unit-box example. Here the box side is 2 and the shear velocity 0.5, so the
  turnover time is longer: check with `plot_slices.py` that the flow is turbulent at `tStart` before using it.


Running
-------
  - `make main`, then `mpirun -np 128 ./main` (nxRanks x nyRanks x nzRanks = 8 x 8 x 2 must equal the number of processes;
  the grid must divide evenly: 256 / 8, 256 / 8, 256 / 2). No seed argument.

  - Size: 256^3 is about 4.5 GB per snapshot, so 50 snapshots are about 220 GB (the file quota on iridis is 2 TB).
  384^3 would be about 14 GB per snapshot and needs 8 x 8 x 3 = 192 ranks. For comparison, the earlier unit-box
  runs had dx = 1/576 = 1.7e-3 whereas 256^3 on this box has dx = 7.8e-3, so the dissipation range is wider in units
  of the box. The window-length and reach rules of the filtering pipeline (`make_mhd_config.py`) depend on dx and dt,
  not on the box size.

  - `plot_slices.py` is copied unchanged from the other example; check its slice coordinates assume nothing about [0, 1].
