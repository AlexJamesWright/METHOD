Model
-----
  - Ideal MHD (SRMHD), 3D, gamma = 4/3, WENO5 + RK3


Initial data
------------
  - `KHInstabilitySingleFluidVigano2020` (Vigano et al., arXiv:2004.00870, Eqs. 34-38; see the `FilterWindow` example's README for
  the deviations from Cook & Bernuzzi's printed listing) on [-1/2, 1/2]^3. Deterministic.


Purpose
-------
  - A long run to a chosen end time with what is needed to pick the filtering window afterwards:
    - energy diagnostics in the log, `diag t=... iter=... Emag=... Ekin=... mass=... Wmax=...`, every `diagIters` iterations (global
    magnetic energy, kinetic energy rho W (W-1), mass and the largest Lorentz factor), so the growth phase can be read off the
    log without opening any file;
    - primitive-variable frames `frame.checkpoint.<t>.hdf5` for plotting (`../SingleFluidIdealVigano2020MPIHDF5FilterWindow/plot_quantities.py`);
    - restartable frames `restart.checkpoint.<t>.hdf5` (conserved + primitive, no auxiliary variables);
    - `data_final.hdf5` (everything) at the end.
  - The filtering window is then written by restarting from the last restart frame before the chosen time, with
  `../SingleFluidIdealVigano2020MPIHDF5Restart` and `nSnap`, `snapStride` (see its README). The window needs full output
  (`data_window_NN.hdf5` with auxiliary variables), about 14 GB per snapshot at 384^3.


Running
-------
  - `make main`, then `mpirun -np 192 ./main [endTime] [diagIters] [primFrameIters] [restartFrameIters] [closureC] [densityTerm] [limiter]`;
  defaults 14.0, 100, 370, 1480 iterations, 0 (no closure), 0, 0. At 384^3 the step is 0.52 dx = 1.35e-3, so 370 iterations are 0.5 and 1480 are 2 time units (0 switches a kind off).
  `submit_run.sh` is a Slurm script (60 h, 192 tasks); it sources `$HOME/METHOD/modules.sh` as in the earlier jobs.
  - Cost (an extrapolation from the 600x600x300 timing, 4e-5 core-s per cell-step; not measured at this size): about 12 s per step
  on 192 cores, 10,300 steps to t = 14, about 30 h.
  - Disk at 384^3: primitive frame about 3.6 GB (28 of them, about 100 GB), restart frame about 7.7 GB (7 of them, about 54 GB), final output
  about 14 GB. A 30-snapshot window is about 410 GB. Check the quota (2 TB total) before submitting.
  The window needs far fewer snapshots than 30: `make_mhd_config.py` requires the half-window (in dx) to cover the filter reach `gamma (1 + sqrt3 v) sqrt(3/5) max(W, L) / 2` plus the meso time-slice offset;
  with `snapStride = 20` the snapshots are about 10.4 dx apart (dt = 0.52 dx), so 7 snapshots (about 100 GB) give a half-window of 31 dx (enough for W = L = 12 cells at v = 0.87) and 9 (about 130 GB) 42 dx (enough at v = 0.94, the largest
  Lorentz factors seen, W = 12). Run `make_mhd_config.py` on the first few snapshots to check before committing the disk.
  - Chosen from the 128^3 laptop runs (energy history in `action_filtering/figures`): at 128^3 the magnetic energy grows from t = 5.6 to
  about 12 and the best-developed turbulence is at t ~ 8-12; the onset is expected earlier at higher resolution, so pick the
  window time from the `diag` lines.
  - Tested at 32^3 on 8 ranks (frames, diagnostics, clean exit); not run at 384^3.


Sub-filter closure
------------------
  - `closureC > 0` switches on `SubfilterGradient` (`Project/CPU/Include/subfilterGradient.h`): the second-order electromotive force, Maxwell
  stress and Poynting-flux deficits of a coarse run as sources for B, S and tau, with filter variance sigma^2 = closureC dx^2 (see the header for the
  formulas); `densityTerm = 1` adds the density-weighting term of the EMF; `limiter = 1` zeroes the source where the local energy exchange J.EMF
  (J = curl B) is negative. The log then also has `closure t=... exchange=... abs=... net/abs=...` lines.
  - Use it at coarse resolution (the 128^3 case of the paper; this example's 384^3 is the reference).
  - Tested only on the laptop at 32^3 and 64^3 (see the commit messages and planning/closure_stage2_design.md): stable at small `closureC`,
  unstable above a threshold unless the limiter is on.
