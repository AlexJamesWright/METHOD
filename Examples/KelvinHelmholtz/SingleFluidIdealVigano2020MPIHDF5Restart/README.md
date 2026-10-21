Continue a run of `../SingleFluidIdealVigano2020MPIHDF5FilterWindow` (or any SRMHD 3D run with the same Data layout) from a
full-output HDF5 file to a later time.

    make main
    mpirun -np 8 ./main restart_file.hdf5 endTime [frameIters] [finalName]

  - The restart file needs the `Conserved` and `Primitive` groups (a window snapshot, a `data_final` or an `OUTPUT_ALL` frame);
  primitive-only frames (`OUTPUT_MINIMAL`) cannot be restarted from.
  - The number of ranks (`nxRanks`, `nyRanks`, `nzRanks`, top of `main.cc`) must equal that of the run that wrote the file.
  - Writes `frame.checkpoint.<t>.hdf5` every `frameIters` iterations (default 250), full output, themselves restartable, and
  `<finalName>.hdf5` (default `data_final`) at the end. The last frame is closed explicitly before exit.
  - Checked on a 32^3 run (8 ranks): primitives read back exactly, conserved variables to 9e-8 (the restart recomputes them from the
  primitives), and a continuation to a later time runs without failed recoveries. The iteration counter restarts at zero, so the
  30-iteration slow-start applies again.
