Model
-----
  - Ideal MHD + REGIME



Initial data
------------
  - Self-similar current sheet



InteractivePlot help
--------------------
  - The current sheet problem has an exact solution, which we can plot against
  using
    `Plot.plotSingleFluidCurrentSheetAgainstExact()`

  - Use
      `help(Plot.plotSingleFluidCurrentSheetAgainstExact)`
    for more information about this plotting tool.



Notes
-----
  - This simulation demonstrates the effectiveness of REGIME at capturing
  resistive effects. Feel free to change the conductivity (sigma) to see how
  this changes the result. It can be instructive to run this as an ideal
  simulation also, by using a conductivity of 10000, for example (beware
  that in this limit, the resolution becomes the dominating factor for the
  diffusion of the magnetic field).
