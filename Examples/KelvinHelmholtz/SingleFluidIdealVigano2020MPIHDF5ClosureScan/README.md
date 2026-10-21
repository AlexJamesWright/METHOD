Model
-----
  - Ideal MHD (SRMHD), 3D, gamma = 4/3, WENO5 + RK3, with the optional sub-filter gradient closure (`SubfilterGradient`)


Purpose
-------
  - Experiment A of `action_filtering/planning/closure_stage2_design.md`: coarse 128^3 runs of the Vigano et al. KHI, with and without the
  closure, to be compared with the 384^3 reference of `../SingleFluidIdealVigano2020MPIHDF5Run` filtered down to the 128^3 grid (a factor 3).
  Same initial data, numerics and parameters as the reference run; only the grid (128^3 on 4 x 4 x 4 = 64 ranks) and the closure differ.
  - Output is at EXACTLY the times given (`outputTimes`; the step is cut to land on each one), so every coarse frame has a reference frame of
  the same time: the 384^3 frames are written every 370 iterations and their times are not round numbers.


Running
-------
  1. On iridis, after the 384^3 run has written its frames, get the times: `python frame_times.py --dir <384 run dir> --tmin 4 > times.txt`
  (exact, from the files). Before they exist, `python frame_times.py --formula --tmin 4 > times.txt` gives the same times from the time step
  (checked against a logged `diag` line to 8 digits). Use `--tmin` to skip the early times (the field is flat before t ~ 5).
  2. `make main`; then in a scratch directory with `main` and `times.txt`: `sbatch submit_scan.sh`. This is an 11-task Slurm array, 64 tasks each,
  one configuration per task (c = 0, 1/12, 1/6, 1/3, 2/3, 1, 2; c = 1/3 with the density term; c = 1 and 2 with the backscatter limiter; c = 1e-6, the baseline up to rounding, whose difference from c = 0 is the run-to-run noise of the chaotic flow: a closure effect smaller than that is not measurable),
  each in its own directory `c<c>_d<densityTerm>_l<limiter>/` with `run.log`, `frame.checkpoint.<t>.hdf5` and `data_final.hdf5`.
  `submit_scan.sh` sources `$HOME/METHOD/modules.sh` as in the other jobs and requests 6 h per task.
  3. Compare (on the laptop or iridis; needs h5py, numpy, scipy, matplotlib):
  `python ../SingleFluidIdealVigano2020MPIHDF5FilterWindow/compare_scan.py --fine-dir <384 run dir> --coarse c=0:<scan>/c0_d0_l0 c=1/12:<scan>/c0.0833333333_d0_l0 ... --out cmp --tmin 4`
  Reads each 384^3 frame (about 3.6 GB) once per time and writes `cmp/summary.txt` (magnetic energy, B-spectrum tail and B_z error against time and
  closure setting) and per-time figures; see `compare_coarse.py` for what each number means.
  - Usage: `mpirun -np 64 ./main endTime closureC densityTerm limiter [diagIters] [outputTimes]` (`outputTimes`: a comma list or the name of a
  file holding them). Log lines `diag t=... Emag=...`, and with the closure `closure t=... net/abs=...`.
  - Cost (not measured): estimate 1.3 s per step on 64 cores; measured 1.0 s per step on iridis (first runs), about 3500 steps to t = 14, so about 1 h per run; the 64^3 laptop
  runs took 15-18 min each on 8 ranks. Disk: about 134 MB per frame, 28 frames, 3.8 GB per run, about 45 GB for the eleven.
  - Judge the closure on E_mag(t), the B spectrum and its high-k tail, not on pointwise errors (in the 64^3 rehearsal, design doc section 9, those were
  dominated by where the eddies are). Expect the coarse run to be most useful past the onset of field growth (the window time chosen from the 384^3 `diag` lines).
  - Tested at 32^3 on 8 ranks (frames at requested times to full precision, closure lines, clean exit), and the comparison chain on 32^3 against 64^3;
  not run at 128^3.
