Model
-----
  - Ideal MHD (SRMHD), 3D


Initial data
------------
  - `KHInstabilitySingleFluidVigano2020`: the 3D special-relativistic Kelvin-Helmholtz initial data of Vigano et al.,
  arXiv:2004.00870 (Eqs. 34-38), as used by Cook & Bernuzzi, arXiv:2606.21659 (II.5.3), on the periodic box [-1/2, 1/2]^3
  (L = 1). Deterministic, no random seed: shear layers at y = +-1/4 with velocity perturbations in all three directions
  (n_x, n_y, n_z = 11, 7, 5) and a weak seed field B_x = 0.001. gamma = 4/3 here.
  - Deviations from Cook & Bernuzzi as printed, following Vigano et al.: v^y_0 = 0 (they list v^i_0 = 0.5 for every
  component; v^y_0 = 0.5 gives head-on flows at y = 0 and across the periodic boundary, and the conserved-to-primitive
  recovery fails within a few dozen steps, as found in a 64^3 test), and the v^z envelope uses sigma_z^2 = 0.1.
  v^z_0 = 0.5 is kept from Cook & Bernuzzi (not among Vigano's standard values).
  - Max |v| at t = 0 is about 0.72 (W ~ 1.4).
  - As in both papers, because of the sgn(y) factors the density and velocity are discontinuous at y = 0 and across the
  periodic boundary y = +-1/2, as well as at the smooth layers.
  - Local 64^3 test (8 ranks, 100 steps to t ~ 1.2): no failed recoveries, mass, momentum, energy and B_x conserved to
  rounding; see the commit message. Resolved turbulence has not been tested.


Purpose
-------
  - Produces output for Lagrangian filtering, exactly as `../SingleFluidIdealRandomMPIHDF5FilterWindow`: spin up to `tStart`
  with no output, then write `nSnap` snapshots of all variables (cons, prims, aux), `data_window_00.hdf5` ..., one every
  `snapStride` timesteps; snapshot times are in `Domain/endTime`; `main` refuses a window that would start inside the
  30-iteration slow-start. See that example's README for the details.
  - `tStart = 5` is carried over from the earlier example. Vigano et al. show the 3D flow at t = 6 with the eddies only just
  forming, so check with `plot_slices.py` that the flow is turbulent at `tStart` before using it.


Running
-------
  - `make main`, then `mpirun -np 128 ./main` (nxRanks x nyRanks x nzRanks = 8 x 8 x 2 must equal the number of processes;
  the grid must divide evenly: 256 / 8, 256 / 8, 256 / 2). No seed argument.

  - Size: 256^3 is about 4.5 GB per snapshot, so 50 snapshots are about 220 GB (the file quota on iridis is 2 TB).
  384^3 would be about 14 GB per snapshot and needs 8 x 8 x 3 = 192 ranks. The earlier runs had dx = 1/576 = 1.7e-3
  on a unit box; 256^3 here has dx = 1/256 = 3.9e-3, so the dissipation range is wider in units of the box. The window-length
  and reach rules of the filtering pipeline (`make_mhd_config.py`) depend on dx and dt, not on the box size.

  - `plot_slices.py` is copied unchanged from the other example; check its slice coordinates assume nothing about [0, 1].


Figures
-------
  - `plot_fig3.py`: recreation of Fig. 3 of Vigano et al. (density above, B^z below, x-y plane, one column per snapshot;
  viridis and seismic colour maps). `plot_quantities.py RUN_DIR`: maps in the x-y and x-z planes of rho, p, |v|, W, vorticity,
  B_x, B_z, |B|, plasma beta, |J|, |E| and the sub-filter magnetic energy (top-hat), the energy history (cf. their Fig. 4),
  spectra, profiles in y, and the |J| distribution. Both use `viz_common.py` (needs h5py, numpy, scipy, matplotlib).
  - Periodic checkpoint frames: calling `sim.evolve(output, safety)` writes a full-primitive/conserved checkpoint file
  `<name>.checkpoint.<t>.hdf5` every `safety` iterations (the main in this directory uses `evolve(output)`, no frames). The
  last frame of a run is never closed by `saveAll(true)` and may be unreadable; the scripts skip it.
  - Local run used for the first figures: 128^3 (the N = 128, C = 0 case of the paper), 8 ranks, `tStart = 6`, frames every 100
  iterations (see the commit message for how it was configured).

Comparing a coarse closure run with the fine run
------------------------------------------------
  - `compare_coarse.py --fine FINE_FILE --coarse LABEL:FILE ... [--fine-frames DIR ...] [--out DIR]` (test:
  `../../../Tests/CompareCoarse`). The fine run (prims or restart frame) is block-averaged by r = N_fine / N_coarse in each direction,
  which restricts the conserved variables exactly onto the coarse grid and is a top-hat filter of one coarse cell (c = 1/12 in
  sigma^2 = c dx^2). Each coarse `data_final.hdf5` (run to the SAME time as the fine file; the script refuses a mismatch) is then
  scored by: relative error of each conserved variable, fluctuation amplitude ratio, resolved magnetic energy, shell spectra of B and
  S (`spectra.png`), the high-k tail of the B spectrum, mid-plane |B| maps and differences (`slices.png`), and the magnetic energy
  history against the block-averaged fine frames (`timeseries.png`, from the `diag` lines of the `run.log` beside each file).
  The first label is the baseline (the run without the closure). Output: `comparison.txt`, `comparison.npz` and the figures.
  - Lab-frame and flat: the real-data (Lagrangian-filtered 384^3) version swaps in the filtered pickles as the reference.

A-priori test of the closure's EMF on fine frames
-------------------------------------------------
  - `apriori_emf.py --block R --dirs RUN [RUN2 ...] [--tmin T] [--every N] [--out FILE]`: for each frame, block-averages the fine primitives over R^3 cells, computes the
  exact sub-filter EMF of E = -v x B (lab-frame and Favre-weighted) and the closure's second-order prediction (c = 1/12, with and without the density term), and prints R^2, the best
  scale (and `scale_c`, which adds back the variance inside one fine cell), the size relative to the resolved E, and the energy exchange J.E (net/abs, backscatter fraction) against time.
  `--block 3` for the 384^3 run (about 4 GB of memory: use a compute node), 2 for a 128^3 run. Test: `Tests/CompareCoarse`.

Budget of the resolved magnetic energy
--------------------------------------
  - `energy_budget.py --block R --dirs RUN [...] [--tmin T] [--out FILE]`: from frames, dE/dt = P - X - D with E = int |<B>|^2/2, the production P = int J.(<v> x <B>) by the resolved
  flow, the sub-filter drain X = int J.E_sub (0 for R = 1) and the remainder D (numerical dissipation and discretisation error); all printed as rates per unit time (divided by E) against time.
  Use `--block 3` on the 384^3 run (the filtered fine run) and `--block 1` on a coarse run to compare the closure's X with the coarse run's own D. Frames should be closer in time than 1/(growth rate).

Delay of one run's history against another
------------------------------------------
  - `time_shift.py REFERENCE_LOG OTHER_LOG`: from the `diag` lines of two run logs, the times at which each run reaches the same E_kin (falling) and E_mag (rising) levels and their difference.
  A constant difference is the same evolution, delayed. Between the 384^3 run and the 128^3 baseline: 2.7-3.0 for E_kin, about 4.7 for E_mag (design doc section 11).

Convergence of the energy histories
-----------------------------------
  - `convergence_order.py N1:LOG1 N2:LOG2 N3:LOG3 ...`: for each E_kin / E_mag level, fits `t_N = t_inf + A N^-p` to the crossing times of runs at several grids (three or more); prints p, the continuum `t_inf` and how
  late each grid is. The runs: `../SingleFluidIdealVigano2020MPIHDF5ClosureScan/submit_convergence.sh`. Test: `Tests/CompareCoarse/test_convergence_order.py`.

Comparing two run logs (works on a run in progress)
---------------------------------------------------
  - `compare_logs.py BASELINE_LOG CLOSURE_LOG [--dt 0.5]`: at matched times, E_mag and E_kin of the closure run relative to the baseline, the same for the perturbation energy (E_mag minus the 5e-7 of the uniform B_x seed, the
  sensitive measure early on), W_max and the failed-recovery counts of both, and the closure's exchange (`net/abs`, and the rate int J.EMF dV / E_mag) from the `closure` lines.

Summary of a scan
-----------------
  - `scan_summary.py REFERENCE_LOG BASELINE_LOG LABEL:LOG ...`: one line per run: the times at which E_mag reaches 1e-6 and 2e-5 and their delays against the reference, the mean growth rate between them, the time at which
  E_kin falls to 0.2 and its delay, E_mag at the end relative to the baseline, the failed recoveries and the peak closure drain X/E_mag.

  - `compare_scan.py --match ekin --fine-log LOG`: compare at equal E_kin instead of equal time (see the ClosureScan README); `--ekin-tol` sets the largest accepted mismatch.

Ensembles of runs
-----------------
  - `ensemble_stats.py NAME=LOG1,LOG2,... [NAME=...] [--reference LOG]`: mean, standard deviation and range over a group of nominally identical runs (the noise runs, or the runs at one closure setting) of the time at which E_mag reaches 1e-6 and 2e-5, the growth
  rate between them, the time at which E_kin falls to 0.2, the largest E_mag and E_mag at the end; with two groups, the difference of the means in standard errors.

The closure's exchange against the exact one, at equal E_kin
------------------------------------------------------------
  - `exchange_at_ekin.py --fine-log FINE_LOG --exact BUDGET_FILE [...] LABEL=C:RUN_LOG ...`: the closure's own `exchange / E_mag` (from the `closure` lines of a coarse run) at the stage where the coarse E_kin equals the fine E_kin, against the exact
  sub-filter exchange X/E of the filtered fine run (the `X/E` column of `energy_budget.py --block R` tables); `C_eff = C x exact / closure` is the C that reproduces the exact exchange at that stage (the closure's exchange is proportional to C).
  An amplitude of a flux, not a level of E_mag, so not limited by the run-to-run timing noise. Test: `Tests/CompareCoarse/test_exchange_at_ekin.py`.
