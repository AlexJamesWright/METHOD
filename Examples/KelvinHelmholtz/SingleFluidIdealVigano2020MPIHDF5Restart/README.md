Continue a run of `../SingleFluidIdealVigano2020MPIHDF5FilterWindow` (or any SRMHD 3D run with the same Data layout) from a
full-output HDF5 file to a later time.

    make main
    mpirun -np 8 ./main restart_file.hdf5 endTime [frameIters] [finalName] [nSnap] [snapStride] [diagIters]

  - The restart file needs the `Conserved` and `Primitive` groups (a window snapshot, a `data_final` or an `OUTPUT_ALL` frame);
  primitive-only frames (`OUTPUT_MINIMAL`) cannot be restarted from.
  - The number of ranks (`nxRanks`, `nyRanks`, `nzRanks`, top of `main.cc`) must equal that of the run that wrote the file.
  - Writes `frame.checkpoint.<t>.hdf5` every `frameIters` iterations (default 250; conserved + primitive, themselves restartable),
  prints energy diagnostics (`diag ...`) every `diagIters` iterations, and writes `<finalName>.hdf5` (default `data_final`, full output)
  at the end. The last frame is closed explicitly before exit.
  - **Filtering window:** with `nSnap > 0`, after `endTime` it writes `nSnap` full-output snapshots `data_window_00 ...`, one every
  `snapStride` steps, with each file's time in `Domain/endTime` (what the Lagrangian-filtering reader uses). To put the window at a chosen
  time, restart from the last `restart.checkpoint.<t>.hdf5` before it and set `endTime` to the window start (at least 30 iterations
  after the restart, METHOD's slow start; the program refuses otherwise).
  - Window and diagnostics checked at 32^3 (8 ranks): three snapshots with uniform spacing, all readable. The restart itself:
  checked on a 32^3 run (8 ranks): primitives read back exactly, conserved variables to 9e-8 (the restart recomputes them from the
  primitives), and a continuation to a later time runs without failed recoveries. The iteration counter restarts at zero, so the
  30-iteration slow-start applies again.
