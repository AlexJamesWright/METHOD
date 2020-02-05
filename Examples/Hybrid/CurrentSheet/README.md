Model
-----
  - Ideal MHD (+ REGIME) + Resistive MHD hybrid model
      This model uses a crossover between resistive MHD and ideal MHD+REGIME.
      The crossover switches between the models smoothly dependant upon the
      value of the conductivity using a tanh function. The speed/steepness of
      the crossover can be changed, but I have played a lot and found the values
      used here and set as the defaults seem to give good results.


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
  - This simulation demonstrates the effectiveness of the hybrid model at capturing
  resistive effects. Notice how including REGIME in the hybrid model (dubbed
  "hybrid+") improves the solution with respect to the exact solution for
  higher conductivities---see this by switching `useREGIME = false`  
