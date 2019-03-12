Model
-----
  - Ideal MHD + REGIME



Initial data
------------
  - Brio-Wu shock tube



InteractivePlot help
--------------------
  - This is a 1D problem, so we want to plot slices. After running
  interactivePlot.py, use
      `Plot.plotSlice()`
  to plot the primitive variables.

  - Use
      `help(Plot.plotSlice)`
    for more information about plotting heatmaps.



Notes
-----
  - This simulation demonstrates the effectiveness of REGIME at capturing
  resistive effects. Feel free to change the conductivity (sigma) to see how
  this changes the result. It can be instructive to run this as an ideal
  simulation also, by using a conductivity of 9999999999, for example.
