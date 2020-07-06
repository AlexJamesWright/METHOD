Model
-----
  - Two-fluid EMHD



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
    for more information about plotting slices.



Notes
-----
  - The parameters we've chosen here mean that the resolution is sufficient to
  observe the two-fluid effect (i.e. charge separation). This is what the wiggles
  are in the data. If you increase the charge mass ratio to p/m 2000 for both
  species, the skin depth decreases and you would need a higher resolution to
  see the separation. At this point, the wiggles disappear and the result looks
  like a single fluid solution.
