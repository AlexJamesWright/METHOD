# CompareCoarse

Synthetic check of `compare_coarse.py` (in `Examples/KelvinHelmholtz/SingleFluidIdealVigano2020MPIHDF5FilterWindow`), the script that
compares coarse sub-filter-closure runs with the block-averaged fine run. No METHOD build and no real data.

```
python Tests/CompareCoarse/test_compare_coarse.py      # needs numpy, h5py, scipy, matplotlib
```

It checks `conserved_from_prims` against hand values, writes a smooth random fine state, builds an "exact" coarse state (the block average of
the fine conserved variables, error must be zero, energy ratio 1) and a noisy one (error must be small and positive), runs the script and
checks its `comparison.npz` and figures, and checks that a coarse file at a different time is refused. On real data `conserved_from_prims`
agrees with the stored `Conserved` group to about 1e-7 (the recovery tolerance).

It also runs `compare_scan.py` (the experiment-A driver) on two frame times, one of them missing in one run, and checks that the common time is
compared and the other is skipped. The example's `main.cc` (`../../Examples/KelvinHelmholtz/SingleFluidIdealVigano2020MPIHDF5ClosureScan`) was checked
separately on 32^3: frames written at the requested times to full precision, readable, closure lines present.

The equal-E_kin mode of `compare_scan.py` is checked with synthetic logs (coarse runs lagging by 0.2 in E_kin(t)): frames at 0.7 and 1.2 must be paired with the fine frames at 0.5 and 1.0, frames 6 % off in E_kin are refused at a 1 % tolerance and accepted at 10 %.

`test_exchange_at_ekin.py` checks `exchange_at_ekin.py` on synthetic logs: a coarse run lagging the fine one by one time unit whose closure exchange is half the exact one at the matched stage must give the same `C_eff` for two different `C`.

`compare_slices.py` is smoke-tested on synthetic frames and logs: one figure per level, the match file, and an empty panel (`no frame within`) for a level that no frame of a run is close to.

`compare_spectra.py` is smoke-tested in both layouts (stages and groups) with a reference: the figure, the `.npz` and the match file are written.
