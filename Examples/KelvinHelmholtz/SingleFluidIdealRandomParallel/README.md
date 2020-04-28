Model
-----
  - Ideal MHD



Initial data
------------
  - Kelvin-Helmholtz instability with random perturbation



InteractivePlot help
--------------------
  - This is a 2D problem, so we want to plot heatmaps. After running
  interactivePlot.py, use

      `Plot.plotHeatMaps()`

  to plot the primitive variables as heatmaps.

  - Use

      `help(Plot.plotHeatMaps)`

    for more information about plotting heatmaps.



Animation help
--------------
  - This is a great simulations to generate animations with. Running animation.py
  as main will save a gif of the density evolution in the Examples/ directory.



Notes
-----
  - As this is a 2D simulation, it will take significantly longer to run. Expect
  this to take five minutes or so to generate the timeseries data.

  - If you have the time, I recommend running a higher resolution simulation of
  this. It is aesthetically pleasing!

  - This simulation should be run using `make main` and then `./main seed` where
  `seed` is some integar seed for the random perturbation on the layer.


Parallel version -- current limitations
---------------------------------------
  - Save data does not include boundary cells
  - Save data gathers all state vectors to proc0 so while this version will be faster, it doesn't currently allow a larger problem to be simulated than will fit on one node
  - nx, ny, nz must currently be even multiple of nxRanks, nyRanks, nzRanks, and these must be specified as args to PlatformEnv manually rather than being calculated

