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
