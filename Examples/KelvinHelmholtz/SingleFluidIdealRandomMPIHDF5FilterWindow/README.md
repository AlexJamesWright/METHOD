Model
-----
  - Ideal MHD (SRMHD), 3D


Initial data
------------
  - Kelvin-Helmholtz instability with random perturbation and seed field Bz = 0.1


Purpose
-------
  - Produces output for Lagrangian filtering. The simulation is spun up to `tStart`
  (default 5.0, once turbulence has developed) with no output. It then writes
  `nSnap` (default 10) snapshots of all variables (cons, prims, aux):
  `data_window_00.hdf5` ... `data_window_09.hdf5`, one every `snapStride` timesteps.

  - Snapshot 00 is at t = tStart exactly. Every step in the window has the nominal
  CFL dt, so the snapshots are uniformly spaced in time (spacing = snapStride * dt,
  recorded in the `Domain/dt` attribute).

  - The window is written by calling `Simulation::updateTime()` directly, so dt is
  never clipped to an end time. `main` refuses to run if the window would start
  inside METHOD's 30-iteration slow-start.


Running
-------
  - `make main`, then `mpirun -np 8 ./main seed` (nxRanks x nyRanks x nzRanks must equal
  the number of processes).

  - Each snapshot is ~900 MB at 200x200x100, gathered to rank 0. For a first run
  consider nx=ny=128, nz=64.

  - Parameters (`tStart`, `nSnap`, `snapStride`, grid, ranks) are at the top of `main.cc`.
  The timespan of the window is only (nSnap-1)*snapStride*dt (~0.026 at 200x200x100 with
  the defaults). Increase `snapStride` if time derivatives of the filtered fields are noisy.
