Standalone check of `SubfilterGradient` (mock `Data`, no METHOD build): the real `sourceExtension` on analytic periodic fields (32^3, 4 ghost layers) against an independent numpy implementation of the closure formulas
(`reference.py`, mirroring `OhmsRegression.em_flux_deficits` in Lagrangian-filtering) for the plain closure, with the density term, with the backscatter limiter and with `quadFactor` 1.5 and 2 (agreement to 2e-15), and the conservation sums of the sources
(zero to rounding without the limiter). Run `./run.sh` (needs `g++` and `numpy`).
